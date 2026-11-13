"""Step 7: family effect, the variance explained by family (eta squared) per marker class.

Extracted from compute_association_eta2_v1.py (eta_vec, assoc). Only individuals assigned to a
family enter. Null: family labels permuted, genotype matrix fixed (keeps allele frequency,
missingness and family sizes); empirical P with +1; Benjamini-Hochberg within class.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .io import Progress, allele_freq, detail


def eta_vec(G, lab, min_called):
    ids, idx = np.unique(lab, return_inverse=True)
    M = G >= 0; V = np.where(M, G, 0).astype(float); n = M.sum(1).astype(float)
    Z = np.zeros((G.shape[1], len(ids))); Z[np.arange(G.shape[1]), idx] = 1
    S = V @ Z; N = M.astype(float) @ Z; tot = V.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ssb = np.where(N > 0, S ** 2 / np.where(N > 0, N, 1), 0).sum(1) - tot ** 2 / n
        sst = (V ** 2).sum(1) - tot ** 2 / n
        e = ssb / sst
    e[(n < min_called) | (sst <= 0)] = np.nan
    return e


def family_effect(ms, fam_of, cfg, seed, n_perm=None):
    """Returns (per-marker table, summary row) or (None, row with the reason)."""
    c = cfg["family_effect"]; n_perm = c["n_permutations"] if n_perm is None else n_perm
    cols = [i for i, s in enumerate(ms.samples) if fam_of.get(s, "unassigned") != "unassigned"]
    base = {"marker_class": ms.klass, "marker_set": ms.name, "n_individuals_in_families": len(cols)}
    if len(cols) < c["min_called"]:
        return None, {**base, "n_tested": 0, "note": f"fewer than {c['min_called']} individuals in families carry this input"}
    G = ms.G[:, cols]; lab = np.array([fam_of[ms.samples[i]] for i in cols])
    af = allele_freq(G)
    keep = np.flatnonzero((af > c["af_min"]) & (af < c["af_max"]) & ((G >= 0).sum(1) >= c["min_called"]))
    note = ""
    if len(keep) > c["max_markers_per_class"]:
        keep = np.sort(np.random.default_rng(seed + 101).choice(keep, c["max_markers_per_class"], replace=False))
        note = f"seeded subsample of {c['max_markers_per_class']} markers"
    G = G[keep]
    obs = eta_vec(G, lab, c["min_called"]); ok = ~np.isnan(obs)
    if not ok.any():
        return None, {**base, "n_tested": 0, "note": "no marker passes the filters"}
    Gk = G[ok]; o = obs[ok]
    rng = np.random.default_rng(seed)
    null = np.empty((n_perm, len(o)))
    detail(f"{ms.name}: {len(o):,} markers tested in {len(set(lab))} families")
    pr = Progress(f"{ms.name}: family-label permutations", n_perm, min_total=100)
    for i in range(n_perm):
        null[i] = eta_vec(Gk, rng.permutation(lab), c["min_called"])
        pr.update(i)
    pe = (1 + (null >= o[None, :]).sum(0)) / (n_perm + 1)
    d = pd.DataFrame({"marker_class": ms.klass, "marker_set": ms.name, "id": ms.ids[keep][ok], "chrom": ms.chrom[keep][ok],
                      "pos": ms.pos[keep][ok], "af": allele_freq(Gk), "eta2_obs": o, "null_median": np.nanmedian(null, 0),
                      "null_p95": np.nanpercentile(null, 95, axis=0), "null_max": np.nanmax(null, 0), "p_emp": pe})
    srt = np.sort(pe); m = len(srt); pas = srt <= np.arange(1, m + 1) / m * 0.05
    kmax = int(np.flatnonzero(pas).max() + 1) if pas.any() else 0
    row = {**base, "n_tested": m, "eta2_obs_median": float(np.median(o)),
           "null_median_of_medians": float(np.median(d.null_median)),
           "pct_above_null_median": 100 * float((d.eta2_obs > d.null_median).mean()),
           "pct_above_null_p95": 100 * float((d.eta2_obs > d.null_p95).mean()),
           "pct_above_all_permutations": 100 * float((d.eta2_obs > d.null_max).mean()),
           "n_BH_0.05": kmax, "n_permutations": n_perm, "note": note}
    return d, row
