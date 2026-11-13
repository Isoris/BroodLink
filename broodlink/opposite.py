"""Opposite homozygotes by relationship class (panel "Opposite-homozygote gate" of the shared-parent figure).

Two quantities per pair of individuals, each summarised by relationship class and shown as a
ratio to the class of unrelated / unlinked pairs:

* R0 of the biSNPs. The original analysis read the R0 column of the NGSrelate table
  (`opposite_homozygote_test_v01.tsv`: columns ida, idb, theta, k0, k1, k2, R0, cls). BroodLink
  uses the R0 column of `--relatedness` when the table has one. Without it, and with SNP
  genotypes, R0 is counted from hard calls as (sites where the two individuals are opposite
  homozygotes) / (sites where both are heterozygous): the definition of R0, on hard calls
  instead of genotype likelihoods (this count is new in BroodLink).
* Opposite-homozygote rate of the deletions (`opposite_homozygote_DEL_v01.tsv`: columns
  oh_rate, n_shared, cls): sites where one individual is 0/0 and the other 1/1, over the sites
  genotyped in both.

The script that wrote the two original tables was not found. The class definitions are read
from the `cls` column of the tables and from the analysis note that describes them
(reimplemented from that note; no marker filter is recorded, so every marker of the class is
used):

  R0 classes:  within-family (full-sib) | cross-family, same shared-parent group (two families
               of one network component) | second degree or closer without such a link |
               third degree | unrelated
  DEL classes: within-family (full-sib) | directly linked families | same group, not directly
               linked | no link
"""
from __future__ import annotations
import numpy as np
import pandas as pd

R0_CLASSES = ["within-family (full-sib)", "cross-family, same shared-parent group", "second degree or closer, no link", "third degree", "unrelated"]
DEL_CLASSES = ["within-family (full-sib)", "directly linked families", "same group, not directly linked", "no link"]


def _pair_counts(G, chunk=20000):
    """Opposite homozygotes, both heterozygous and both called, for every pair of columns."""
    n = G.shape[1]
    oh = np.zeros((n, n)); hh = np.zeros((n, n)); cc = np.zeros((n, n))
    for a in range(0, G.shape[0], chunk):
        g = G[a:a + chunk]
        z = (g == 0).astype(np.float32); t = (g == 2).astype(np.float32); h = (g == 1).astype(np.float32); c = (g >= 0).astype(np.float32)
        x = z.T @ t
        oh += x + x.T; hh += h.T @ h; cc += c.T @ c
    return oh, hh, cc


def opposite_homozygotes(R, fam_of, links, comp_of, snp, dels, cfg):
    """R: relatedness pairs (id1, id2, theta, optional R0). fam_of: individual -> family or
    'unassigned'. links: kept (non-rejected) links between families (unit_a, unit_b).
    comp_of: family -> network component. snp / dels: MarkerSet or None.
    Returns (per-pair table, per-class summary, info dict); the tables are empty when neither
    quantity can be computed."""
    r = cfg["relatedness"]
    P = R[["id1", "id2", "theta"]].copy()
    fa = P.id1.map(lambda s: fam_of.get(s, "unassigned")).to_numpy(); fb = P.id2.map(lambda s: fam_of.get(s, "unassigned")).to_numpy()
    asg = (fa != "unassigned") & (fb != "unassigned")
    within = asg & (fa == fb)
    linked = {frozenset(x) for x in zip(links.unit_a, links.unit_b)} if links is not None and len(links) else set()
    ca = np.array([comp_of.get(f, 0) for f in fa]); cb = np.array([comp_of.get(f, 0) for f in fb])
    same_group = asg & ~within & (ca == cb) & (ca > 0)
    direct = np.array([bool(s) and frozenset((x, y)) in linked for s, x, y in zip(same_group, fa, fb)])
    th = P.theta.to_numpy(float)
    P["class_R0"] = np.select([within, same_group, th >= r["theta_second"], th >= r["theta_third"]], R0_CLASSES[:4], R0_CLASSES[4])
    P["class_DEL"] = np.select([within, direct, same_group], DEL_CLASSES[:3], DEL_CLASSES[3])
    info = {"R0_source": None, "n_snps": None, "n_deletions": None}
    P["R0"] = np.nan; P["DEL_opposite_homozygote_rate"] = np.nan; P["n_deletions_called_in_both"] = np.nan

    def per_pair(ms, M):
        ix = {s: i for i, s in enumerate(ms.samples)}
        i = P.id1.map(ix); j = P.id2.map(ix); ok = (i.notna() & j.notna()).to_numpy()
        out = np.full(len(P), np.nan)
        out[ok] = M[i[ok].astype(int).to_numpy(), j[ok].astype(int).to_numpy()]
        return out
    if "R0" in R.columns and R["R0"].notna().any():
        P["R0"] = R["R0"].to_numpy(float); info["R0_source"] = "R0 column of the relatedness table"
    elif snp is not None and snp.n:
        oh, hh, _cc = _pair_counts(snp.G)
        with np.errstate(divide="ignore", invalid="ignore"):
            P["R0"] = per_pair(snp, np.where(hh > 0, oh / hh, np.nan))
        info["R0_source"] = f"hard calls of {snp.name} (opposite homozygotes / both heterozygous)"; info["n_snps"] = int(snp.n)
    if dels is not None and dels.n:
        oh, _hh, cc = _pair_counts(dels.G)
        with np.errstate(divide="ignore", invalid="ignore"):
            P["DEL_opposite_homozygote_rate"] = per_pair(dels, np.where(cc > 0, oh / cc, np.nan))
        P["n_deletions_called_in_both"] = per_pair(dels, cc)
        info["n_deletions"] = int(dels.n)
    rows = []
    for col, ccol, classes, q in (("R0", "class_R0", R0_CLASSES, "biSNP R0"), ("DEL_opposite_homozygote_rate", "class_DEL", DEL_CLASSES, "DEL opposite-homozygote rate")):
        if not P[col].notna().any():
            continue
        base = P.loc[P[ccol] == classes[-1], col].median()
        for c in classes:
            v = P.loc[P[ccol] == c, col].dropna()
            rows.append({"quantity": q, "class": c, "n_pairs": int(len(v)), "median": float(v.median()) if len(v) else np.nan,
                         "ratio_to_reference_class": float(v.median() / base) if len(v) and base and np.isfinite(base) else np.nan, "reference_class": classes[-1]})
    S = pd.DataFrame(rows, columns=["quantity", "class", "n_pairs", "median", "ratio_to_reference_class", "reference_class"])
    if not len(S):
        P = P.iloc[0:0]
    return P.drop(columns=["id1", "id2"]).reset_index(drop=True), S, info
