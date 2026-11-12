"""Step 3: nomination of candidate half-sib links between units (screening only).

Extracted from sp_lib.pair_matrix: for every pair of units, the fraction of cross-unit pairs of
individuals at second-degree kinship or closer (theta >= theta_min; optional upper bound) is
compared with a size-matched null made by permuting the unit labels among individuals;
P = (exceedances + 1) / (permutations + 1); Benjamini-Hochberg over all unit pairs.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for i, ix in enumerate(o[::-1]):
        prev = min(prev, p[ix] * n / (n - i)); q[ix] = prev
    return q


def nominate(R, unit_of, cfg, seed, n_perm=None, marker_class="", n_markers=np.nan):
    """R: pair table; unit_of: individual -> unit name. Returns the unit-pair table."""
    c = cfg["nomination"]
    n_perm = c["n_permutations"] if n_perm is None else n_perm
    ids = sorted(unit_of); ix = {s: i for i, s in enumerate(ids)}; N = len(ids)
    R = R[R.id1.isin(ix) & R.id2.isin(ix)]
    H = np.zeros((N, N), np.float32)
    ii = R.id1.map(ix).to_numpy(); jj = R.id2.map(ix).to_numpy()
    th = R.theta.to_numpy(float)
    h = th >= c["theta_min"]
    if c["theta_max"] is not None:
        h &= th < c["theta_max"]
    h = h.astype(np.float32)
    H[ii, jj] = h; H[jj, ii] = h
    units = sorted(set(unit_of.values())); uix = {u: i for i, u in enumerate(units)}
    lab = np.array([uix[unit_of[s]] for s in ids]); U = len(units)

    def frac(l):
        L = np.zeros((N, U), np.float32); L[np.arange(N), l] = 1
        S = L.T @ H @ L
        n = L.sum(0)
        return S / np.maximum(np.outer(n, n), 1), S
    obs, S = frac(lab)
    n = np.bincount(lab, minlength=U)
    iu, ju = np.triu_indices(U, 1)
    rng = np.random.default_rng(seed)
    ge = np.zeros(len(iu), int); o = obs[iu, ju]
    for _ in range(n_perm):
        f, _s = frac(rng.permutation(lab))
        ge += (f[iu, ju] >= o - 1e-12)
    p = (ge + 1) / (n_perm + 1)
    un = np.array(units)
    g = pd.DataFrame({"unit_a": un[iu], "unit_b": un[ju], "n_a": n[iu], "n_b": n[ju], "n_cross_pairs": n[iu] * n[ju],
                      "n_pairs_second_degree_or_closer": S[iu, ju].astype(int), "fraction": o.astype(float), "p_perm": p})
    g["unassigned_sides"] = g.unit_a.str.startswith("S_").astype(int) + g.unit_b.str.startswith("S_").astype(int)
    g["q_value"] = bh(g.p_perm)
    g["min_n"] = g[["n_a", "n_b"]].min(axis=1)
    g["nominated"] = g.q_value < c["q_max"]
    g["marker_class"] = marker_class; g["n_markers"] = n_markers; g["n_permutations"] = n_perm
    return g


def split_markers(ms, cfg):
    """Single-input mode: nomination set and audit set by chromosome parity or fixed blocks
    (never at random inside linkage blocks). Returns two boolean masks."""
    c = cfg["nomination"]
    chroms = list(dict.fromkeys(ms.chrom.tolist()))
    if c["single_input_split"] == "chromosome_parity" and len(chroms) >= 2:
        rank = {ch: i for i, ch in enumerate(chroms)}
        par = np.array([rank[x] % 2 for x in ms.chrom])
        return par == 0, par == 1, f"chromosome parity ({(len(chroms) + 1) // 2} chromosomes for nomination, {len(chroms) // 2} for the audit)"
    blk = (ms.pos // c["block_size_bp"]).astype(np.int64)
    off = {ch: i for i, ch in enumerate(chroms)}
    par = (blk + np.array([off[x] for x in ms.chrom])) % 2
    return par == 0, par == 1, f"alternating fixed blocks of {c['block_size_bp']:,} bp"
