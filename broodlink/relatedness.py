"""Step 1: pairwise relatedness (theta, k0, k1, k2).

Source, in priority order:
  1. --relatedness, used as is;
  2. --bisnp genotype likelihoods through the external tool NGSrelate when it is installed
     (the original analysis used NGSrelate; the project has no own likelihood estimator of
     k0/k1/k2, so there is no Python fallback of that kind);
  3. a VanRaden genomic relationship matrix (GRM) from the densest genotype input,
     kinship phi = GRM / 2, as in compute_sv_kinship_v1.py (grm_acc) and grm_verify.py.
     No k2 exists in this mode.
"""
from __future__ import annotations
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

from .io import allele_freq, log, read_relatedness_file


def vanraden_grm(G, af_min=0.05, af_max=0.95, min_called_fraction=0.0, chunk=50000):
    """GRM_ij = sum_m Z_im Z_jm / (markers called in both), Z = (g - 2p) / sqrt(2p(1-p)), missing -> 0.

    Returns (GRM, number of markers used). Filters are strict inequalities, as in the original.
    """
    n = G.shape[1]
    K = np.zeros((n, n)); D = np.zeros((n, n)); used = 0
    nmin = min_called_fraction * n
    for a in range(0, G.shape[0], chunk):
        Gc = G[a:a + chunk]
        called = Gc >= 0
        nc = called.sum(1)
        p = np.where(called, Gc, 0).sum(1) / (2.0 * np.maximum(nc, 1))
        ok = (nc >= nmin) & (p > af_min) & (p < af_max) & (nc > 0)
        if not ok.any():
            continue
        Gc, called, p = Gc[ok], called[ok], p[ok]
        Z = np.where(called, (Gc - 2 * p[:, None]) / np.sqrt(2 * p * (1 - p))[:, None], 0.0).astype(np.float32)
        K += Z.T @ Z
        c = called.astype(np.float32)
        D += c.T @ c
        used += int(ok.sum())
    return K / np.maximum(D, 1), used


def grm_to_pairs(K, samples, source):
    iu, ju = np.triu_indices(len(samples), 1)
    s = np.asarray(samples)
    R = pd.DataFrame({"id1": s[iu], "id2": s[ju], "theta": K[iu, ju] / 2.0})
    a = np.where(R.id1 < R.id2, R.id1, R.id2); b = np.where(R.id1 < R.id2, R.id2, R.id1)
    R["id1"], R["id2"] = a, b
    R["k0"] = np.nan; R["k1"] = np.nan; R["k2"] = np.nan; R["source"] = source
    return R


def relatedness_from_markers(ms, cfg, label=None):
    g = cfg["relatedness"]["grm"]
    frac = g["min_called_fraction_snp"] if ms.klass in ("SNP", "INDEL", "SNP_INDEL") else g["min_called_fraction_sv"]
    K, used = vanraden_grm(ms.G, g["af_min"], g["af_max"], frac)
    R = grm_to_pairs(K, ms.samples, label or f"GRM:{ms.klass}")
    return R, used


def try_ngsrelate(beagle_files, n_ind, sample_names, cfg):
    """Run NGSrelate on one beagle file when the binary is installed. Returns (pairs or None, message)."""
    exe = shutil.which(cfg["relatedness"]["ngsrelate_binary"])
    if exe is None:
        return None, "NGSrelate is not installed (binary '%s' not on PATH)" % cfg["relatedness"]["ngsrelate_binary"]
    if len(beagle_files) != 1:
        return None, "NGSrelate needs one beagle file; several were given (concatenate them to use NGSrelate)"
    out = os.path.join(tempfile.mkdtemp(prefix="broodlink_ngsrelate_"), "relatedness.res")
    cmd = [exe, "-G", beagle_files[0], "-n", str(n_ind), "-O", out]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=24 * 3600)
        R, _ = read_relatedness_file(out, sample_names)
        return R, "NGSrelate: " + " ".join(os.path.basename(c) if os.path.sep in c else c for c in cmd)
    except Exception as exc:                                         # pragma: no cover - needs the external tool
        return None, f"NGSrelate failed ({type(exc).__name__}); falling back to the GRM"


def pair_matrices(R, ids=None):
    """Symmetric theta / k0 / k1 / k2 matrices (NaN where a pair is absent)."""
    ids = sorted(set(R.id1) | set(R.id2)) if ids is None else list(ids)
    ix = {s: i for i, s in enumerate(ids)}
    R = R[R.id1.isin(ix) & R.id2.isin(ix)]
    ii = R.id1.map(ix).to_numpy(); jj = R.id2.map(ix).to_numpy()
    out = {}
    for c in ("theta", "k0", "k1", "k2"):
        M = np.full((len(ids), len(ids)), np.nan)
        v = R[c].to_numpy(float)
        M[ii, jj] = v; M[jj, ii] = v
        out[c] = M
    return ids, out


def classify_pairs(R, cfg):
    """Relationship class of every pair, as common.classify (k2 used only where it exists)."""
    r = cfg["relatedness"]; f = cfg["families"]
    th = R.theta.to_numpy(float); k2 = R.k2.to_numpy(float)
    has = ~np.isnan(k2)
    cls = np.full(len(R), "unrelated", dtype=object)
    cls[th >= r["theta_third"]] = "third degree"
    cls[th >= r["theta_second"]] = "second degree"
    fd = th >= r["theta_first"]
    cls[fd & has & (k2 >= f["k2_min"])] = "full-sib"
    cls[fd & has & (k2 < f["k2_min"])] = "first degree, not full-sib (k2 below gate)"
    cls[fd & ~has] = "first degree (no k2)"
    cls[th >= r["theta_duplicate"]] = "identical genotypes (theta at or above the identity threshold)"
    return cls


def unrelated_subset(R, ids, cfg):
    """Greedy maximal set free of first-degree pairs (lowest degree first; deterministic).

    The original used an exact maximum independent set; the greedy set can be smaller.
    """
    t = cfg["relatedness"]["theta_first"]
    E = R[(R.theta >= t) & R.id1.isin(ids) & R.id2.isin(ids)]
    adj = {s: set() for s in ids}
    for a, b in zip(E.id1, E.id2):
        adj[a].add(b); adj[b].add(a)
    left = set(ids); keep = []
    while left:
        s = min(left, key=lambda x: (len(adj[x] & left), x))
        keep.append(s); left -= adj[s] | {s}
    return sorted(keep)
