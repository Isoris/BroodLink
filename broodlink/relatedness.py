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


def find_ngsrelate(cfg, override=None):
    """Path of an NGSrelate binary, or None. `override` (--ngsrelate) may be a path or a name;
    otherwise the names of relatedness.ngsrelate_binaries are looked up on the PATH, in order
    (ngsRelate-fast first: a fork of NGSrelate with per-chromosome downsampling of sites built
    in; with a standard NGSrelate binary BroodLink does that downsampling itself)."""
    names = [override] if override else list(cfg["relatedness"]["ngsrelate_binaries"])
    for n in names:
        if n and os.path.isfile(n):
            return os.path.abspath(n)
        exe = shutil.which(n) if n else None
        if exe:
            return exe
    return None


def write_beagle_frequencies(beagle_files, out_beagle, out_freq, sites_per_gb=0):
    """One Beagle file and one allele-frequency file for NGSrelate (-G, -f).

    Several Beagle files (one per chromosome) are concatenated under the header of the first.
    The frequency file has one line per site: the frequency of allele 2 of the Beagle file (the
    minor allele in ANGSD output, i.e. the knownEM column of the ANGSD .mafs file), estimated
    here from the genotype likelihoods by EM under Hardy-Weinberg.

    sites_per_gb > 0 thins the sites to a balanced density along every chromosome before
    NGSrelate sees them: the first site of each window of 1e9 / sites_per_gb bp is kept (marker
    names must be <chromosome>_<position>, as ANGSD writes them). This is the downsampling idea
    of the ngsRelate-fast fork, done here so that a standard NGSrelate binary can be used.
    Returns (sites written, sites read)."""
    import gzip
    from .io import _open
    n_sites = n_read = 0
    width = int(1e9 // sites_per_gb) if sites_per_gb and sites_per_gb > 0 else 0
    last = (None, -1)
    with gzip.open(out_beagle, "wt", encoding="utf-8", newline="\n", compresslevel=3) as ob, open(out_freq, "w", encoding="utf-8", newline="\n") as of:
        for k, path in enumerate(beagle_files):
            with _open(path) as fh:
                header = fh.readline()
                if k == 0:
                    ob.write(header if header.endswith("\n") else header + "\n")
                while True:
                    lines = fh.readlines(64 * 1024 * 1024)
                    if not lines:
                        break
                    lines = [l for l in lines if l.strip()]
                    n_read += len(lines)
                    if width:
                        kept = []
                        for l in lines:
                            mk = l[:l.index("\t")]
                            ch, _sep, ps = mk.rpartition("_")
                            try:
                                key = (ch, int(ps) // width)
                            except ValueError:
                                key = (mk, 0)                 # not <chromosome>_<position>: keep the site
                            if key != last:
                                kept.append(l); last = key
                        lines = kept
                    if not lines:
                        continue
                    rows = [l.rstrip("\n").split("\t") for l in lines]
                    L = np.array([r[3:] for r in rows], dtype=np.float32)
                    L = L.reshape(len(rows), -1, 3)
                    sm = L.sum(2, keepdims=True)
                    p = gl_allele_freq(L / np.where(sm > 0, sm, 1))
                    ob.write("".join(l if l.endswith("\n") else l + "\n" for l in lines))
                    of.write("".join(f"{x:.6f}\n" for x in p))
                    n_sites += len(rows)
    return n_sites, n_read


def try_ngsrelate(beagle_files, n_ind, sample_names, cfg, override=None, freqs=None, threads=None, workdir=None):
    """Relatedness (theta, k0, k1, k2) from Beagle genotype likelihoods with NGSrelate, as in the
    original study, when a binary is installed. The call is the one of the original cohort run:
    -G beagle -f frequencies -n individuals -z sample names -p threads -O output, plus -D
    (sites per Gb) when the binary is ngsRelate-fast. Returns (pairs or None, message)."""
    c = cfg["relatedness"]
    exe = find_ngsrelate(cfg, override)
    if exe is None:
        return None, "NGSrelate is not installed (looked for %s)" % (override or ", ".join(c["ngsrelate_binaries"]))
    fast_bin = "fast" in os.path.basename(exe).lower()
    wd = workdir or tempfile.mkdtemp(prefix="broodlink_ngsrelate_")
    os.makedirs(wd, exist_ok=True)
    out = os.path.join(wd, "relatedness.res")
    try:
        if freqs is None or len(beagle_files) != 1:
            beagle = os.path.join(wd, "all_sites.beagle.gz"); fq = os.path.join(wd, "allele2_frequencies.txt")
            thin = 0 if (fast_bin or freqs is not None) else int(c["ngsrelate_sites_per_gb"] or 0)
            n_sites, n_read = write_beagle_frequencies(beagle_files, beagle, fq, thin)
            if freqs is not None:
                fq = freqs
            how = "frequencies estimated from the likelihoods by BroodLink" if freqs is None else "frequencies supplied"
            if thin:
                how += f"; sites thinned by BroodLink to {thin:,} per Gb along each chromosome, {n_sites:,} of {n_read:,} kept"
        else:
            beagle, fq, n_sites, how = beagle_files[0], freqs, None, "frequencies supplied"
        zf = os.path.join(wd, "samples.txt")
        with open(zf, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(map(str, sample_names)) + "\n")
        th = int(threads or c["ngsrelate_threads"] or os.cpu_count() or 1)
        cmd = [exe, "-G", beagle, "-f", fq, "-n", str(n_ind), "-z", zf, "-p", str(th), "-O", out]
        fast = fast_bin
        if fast:
            cmd += ["-D", str(c["ngsrelate_fast_sites_per_gb"])]
        log(f"relatedness: running {os.path.basename(exe)} on {n_ind} individuals" + (f", {n_sites:,} sites" if n_sites else "") + f", {th} threads (this can take minutes to hours)")
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7 * 24 * 3600)
        if r.returncode != 0 or not os.path.exists(out) or os.path.getsize(out) == 0:
            tail = " | ".join((r.stderr or "").strip().splitlines()[-3:])
            return None, f"{os.path.basename(exe)} failed (exit {r.returncode}: {tail}); falling back to the GRM"
        R, _ = read_relatedness_file(out, list(sample_names))
        shown = " ".join(os.path.basename(x) if (os.path.sep in x or "/" in x) else x for x in cmd)
        return R, f"{os.path.basename(exe)} ({how}" + (f"; -D {c['ngsrelate_fast_sites_per_gb']} sites per Gb" if fast else "") + f"): {shown}"
    except Exception as exc:
        return None, f"{os.path.basename(exe)} could not be run ({type(exc).__name__}: {exc}); falling back to the GRM"


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


# --------------------------------------------------------------------------- likelihood estimator
def gl_from_genotypes(G, error=0.01):
    """Genotype likelihoods from hard genotypes: 1 - error on the call, error / 2 elsewhere; flat when missing."""
    m, n = G.shape
    L = np.full((m, n, 3), error / 2, np.float32)
    for g in (0, 1, 2):
        L[:, :, g][G == g] = 1 - error
    L[G < 0] = 1 / 3
    return L


def gl_allele_freq(L, iters=30):
    """Frequency of the counted allele per site, by EM on the genotype likelihoods under Hardy-Weinberg."""
    L = L.astype(np.float64)
    p = np.clip((L[:, :, 1] + 2 * L[:, :, 2]).sum(1) / np.maximum(2 * L.sum(2).sum(1), 1e-12), 0.01, 0.99)
    for _ in range(iters):
        pr = np.stack([(1 - p) ** 2, 2 * p * (1 - p), p ** 2], 1)[:, None, :]
        w = L * pr
        w /= np.maximum(w.sum(2, keepdims=True), 1e-300)
        p = np.clip((w[:, :, 1] + 2 * w[:, :, 2]).mean(1) / 2, 1e-4, 1 - 1e-4)
    return p


def _ibd_tables(p):
    """P(g1, g2 | z shared alleles identical by descent) for z = 0, 1, 2: three (sites, 3, 3) arrays."""
    q = 1 - p
    h = np.stack([q ** 2, 2 * p * q, p ** 2], 1)
    T0 = h[:, :, None] * h[:, None, :]
    T1 = np.zeros((len(p), 3, 3))
    T1[:, 0, 0] = q ** 3; T1[:, 0, 1] = p * q ** 2; T1[:, 1, 0] = p * q ** 2
    T1[:, 1, 1] = p * q; T1[:, 1, 2] = p ** 2 * q; T1[:, 2, 1] = p ** 2 * q; T1[:, 2, 2] = p ** 3
    T2 = np.zeros((len(p), 3, 3))
    for g in range(3):
        T2[:, g, g] = h[:, g]
    return T0, T1, T2


def ml_relatedness(L, samples, cfg, source, progress=None):
    """Maximum-likelihood k0, k1, k2 for every pair from genotype likelihoods (sites x samples x 3).

    The model is the one of NGSrelate for non-inbred individuals (Korneliussen and Moltke 2015):
    at each site the pair shares 0, 1 or 2 alleles identical by descent with probabilities
    k0, k1, k2; genotypes are summed over with their likelihoods; allele frequencies come from
    the same likelihoods; sites are treated as independent. k is fitted by EM from several
    starting points. theta = k1 / 4 + k2 / 2. Reimplemented here; not the NGSrelate code.
    Returns (pair table, number of sites used).
    """
    c = cfg["relatedness"]["ml"]
    L = np.asarray(L, np.float32)
    p = gl_allele_freq(L)
    info = (np.abs(L.max(2) - L.min(2)) > 1e-6).mean(1)
    ok = (p >= c["maf_min"]) & (p <= 1 - c["maf_min"]) & (info >= c["min_informative_fraction"])
    L, p = L[ok], p[ok]
    m, n = L.shape[0], L.shape[1]
    if m == 0:
        return None, 0
    T = [t.astype(np.float32) for t in _ibd_tables(p)]
    starts = np.array([[1 / 3, 1 / 3, 1 / 3], [0.98, 0.01, 0.01], [0.25, 0.5, 0.25], [0.5, 0.5, 0.0], [0.02, 0.02, 0.96]])
    rows = []
    for i in range(n - 1):
        J = np.arange(i + 1, n)
        # likelihood of each site under z = 0, 1, 2 for i against every later individual: (pairs, sites, 3)
        S = np.stack([np.einsum("sjb,sb->js", L[:, J, :], np.einsum("sa,sab->sb", L[:, i, :], Tz)) for Tz in T], 2).astype(np.float64)
        S = np.maximum(S, 1e-300)
        best_ll = np.full(len(J), -np.inf); best_k = np.zeros((len(J), 3))
        for k0 in starts:
            k = np.tile(k0, (len(J), 1))
            for _ in range(c["em_iterations"]):
                w = S * k[:, None, :]
                tot = w.sum(2)
                kn = (w / tot[:, :, None]).mean(1)
                done = np.abs(kn - k).max() < c["em_tolerance"]
                k = kn
                if done:
                    break
            ll = np.log((S * k[:, None, :]).sum(2)).sum(1)
            better = ll > best_ll
            best_ll[better] = ll[better]; best_k[better] = k[better]
        for jj, j in enumerate(J):
            rows.append((samples[i], samples[j], best_k[jj, 0], best_k[jj, 1], best_k[jj, 2]))
        if progress is not None:
            progress.update(i)
    R = pd.DataFrame(rows, columns=["id1", "id2", "k0", "k1", "k2"])
    a = np.where(R.id1 < R.id2, R.id1, R.id2); b = np.where(R.id1 < R.id2, R.id2, R.id1)
    R["id1"], R["id2"] = a, b
    R.insert(2, "theta", R.k1 / 4 + R.k2 / 2)
    R["source"] = source
    return R, int(m)
