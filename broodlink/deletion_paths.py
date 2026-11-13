"""Deletion genotypes along linked families (part of step 5, network).

Two things are computed for the figures "deletion genotypes along a path of linked families"
and "deletion genotypes, the other links":

1. Blocks. For a list of families (a chain = one family with its link partners, or the two
   families of one link) and one chromosome: the deletions that fit a Mendelian cross in every
   family of the list, the cross expected among the siblings at each of them (all alike, split
   1:1, split 1:2:1), the even subsample that is drawn, the order of the individuals, and the
   deletions homozygous for the deletion allele (1/1) in at least one individual of both
   families of a link.
   Extracted from HAT/code/build_supp_chromosome_DEL_heatmaps_v2.py (`fit_family`, `block2`,
   `strips2`, `summarise`) and HAT/code/chromDEL_common.py (`af_mask`, `order_fish`).
2. Sharing per chromosome. The Gate-2 sharing score of a link, evaluated on each chromosome,
   against the same score of the k non-nominated auditable pairs nearest in the genome-wide
   Gate-1 score (median, 5th and 95th percentile, maximum).
   Extracted from HAT/code/build_supp_chromosome_DEL_heatmaps_v1.py (`chrom_share_table`) and
   HAT/code/chromDEL_common.py (`share`, `matched_null`, `edges`).

Rule that differs from `audit.sharing`: as in the original per-chromosome scripts, the allele
frequency filter is computed on the individuals assigned to a family (the original genotype
matrix held the family individuals only), not on the whole analysis set.

New in BroodLink (the original typed these in): which families form a chain (`find_chains`),
the order of the families of a chain, and the chromosome shown (`choose_chromosome`: the
chromosome with the most deletions, the rule stated in the original script).
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .network import _fkey

# thresholds of the original scripts; cfg["deletion_paths"] overrides them
DEFAULTS = {
    "low_call": 0.30,            # individuals with more missing genotypes on the chromosome go below the dotted line and are not used for the fit
    "genotype_error": 0.02,      # error mixed into the expected genotype proportions of a cross
    "fit_p_min": 0.05,           # a deletion fits when the simulated P of its best cross is not below this
    "fit_draws": 300, "fit_seed": 1, "fit_min_called": 3,
    "min_call": 0.90,            # call rate among the well-called individuals of every family
    "thin_kb": 100, "thin_concordance": 0.90, "thin_min_called": 3,
    "max_columns": 160,          # deletions drawn (even subsample of the eligible ones)
    "max_chain_families": 6,     # new in BroodLink: a chain is cut to its centre and this many families less one
    "chain_links": "kept",       # new in BroodLink: "kept" = junction families, then families on two or more kept links; "supported" = junction families only
}
STATUS_RANK = {"supported": 0, "Gate-1 only": 1, "provisional": 2, "nominated, not audited": 3}
CROSS_NAMES = ["all alike", "split 1:1", "split 1:2:1"]
_PAIRS = [(a, b) for a in range(3) for b in range(a, 3)]           # parental genotypes (alt allele counts)
_CROSS = np.array([0, 1, 0, 2, 1, 0])                              # index of _PAIRS -> 0 all alike, 1 = 1:1, 2 = 1:2:1


def params(cfg):
    p = dict(DEFAULTS); p.update((cfg or {}).get("deletion_paths") or {})
    return p


def _cross_probs(eps):
    def pr(a, b):
        x, y = a / 2, b / 2
        return np.array([(1 - x) * (1 - y), x * (1 - y) + y * (1 - x), x * y])
    return np.array([pr(a, b) * (1 - eps) + eps / 3 for a, b in _PAIRS])


def fit_family(sub, p):
    """Per deletion: the best parental pair (as cross type) and the Mendelian-fit P of the
    genotype counts of one family (`sub` = deletions x individuals), by parametric simulation."""
    PM = _cross_probs(p["genotype_error"]); lp = np.log(PM)
    n = (sub >= 0).sum(1); C = np.stack([(sub == g).sum(1) for g in range(3)], 1)
    LL = C @ lp.T; best = LL.argmax(1); ll = LL.max(1)
    pv = np.full(len(sub), np.nan); rng = np.random.default_rng(p["fit_seed"]); nd = int(p["fit_draws"])
    for i in range(len(sub)):
        if n[i] < p["fit_min_called"]:
            continue
        sim = rng.multinomial(n[i], PM[best[i]], size=nd) @ lp[best[i]]
        pv[i] = (1 + (sim <= ll[i] + 1e-9).sum()) / (nd + 1)
    return _CROSS[best], pv


def order_individuals(sub):
    """Hierarchical order of the columns of a deletions x individuals matrix (average linkage on
    the mean absolute genotype difference over deletions called in both, optimal leaf order)."""
    n = sub.shape[1]
    if n <= 2:
        return np.arange(n)
    from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            ok = (sub[:, i] >= 0) & (sub[:, j] >= 0)
            D[i, j] = D[j, i] = (np.abs(sub[ok, i].astype(int) - sub[ok, j]).mean() / 2) if ok.any() else 1.0
    d = squareform(D, checks=False)
    return leaves_list(optimal_leaf_ordering(linkage(d, "average"), d))


def af_mask(G, lo, hi):
    called = G >= 0; n = called.sum(1)
    af = np.where(n > 0, np.where(called, G, 0).sum(1) / (2.0 * np.maximum(n, 1)), np.nan)
    return (af > lo) & (af < hi)


def chrom_key(c):
    """Natural sort key of a chromosome name."""
    import re
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", str(c))]


def choose_chromosome(ms):
    """The chromosome with the most deletions (first in natural order on a tie)."""
    if not ms.n:
        return None
    u, n = np.unique(ms.chrom.astype(str), return_counts=True)
    return sorted(zip(u, n), key=lambda x: (-x[1], chrom_key(x[0])))[0][0]


class Context:
    """Genotypes of one deletion set, the columns of every family and the caches of the fit."""

    def __init__(self, ms, members, family_units, cfg):
        self.ms, self.cfg, self.p = ms, cfg, params(cfg)
        ix = {s: i for i, s in enumerate(ms.samples)}
        self.cols = {str(f): np.array([ix[s] for s in members.get(f, []) if s in ix], int) for f in family_units}
        fam_cols = np.concatenate([v for v in self.cols.values()]) if self.cols else np.array([], int)
        self.fam_cols = np.sort(fam_cols)
        a = cfg["audit"]
        # allele frequency among the individuals assigned to a family (original: genotype matrix of the family individuals)
        self.af_ok = af_mask(ms.G[:, self.fam_cols], a["af_min"], a["af_max"]) if len(self.fam_cols) else np.zeros(ms.n, bool)
        self.chrom = ms.chrom.astype(str)
        self._fit, self._chrom = {}, {}

    def _c(self, chrom):
        if chrom not in self._chrom:
            m = np.where(self.chrom == chrom)[0]
            m = m[np.argsort(self.ms.pos[m], kind="stable")]
            Gc = self.ms.G[m]
            self._chrom[chrom] = (m, Gc, self.ms.pos[m], (Gc < 0).mean(0) if len(m) else np.zeros(Gc.shape[1]))
        return self._chrom[chrom]

    def block(self, chrom, fam_list):
        """Eligible and drawn deletions of `fam_list` on `chrom` (see the module docstring)."""
        p = self.p; m, Gc, pc, miss = self._c(chrom)
        groups = []
        for f in fam_list:
            ii = self.cols.get(str(f), np.array([], int))
            groups.append((str(f), ii[miss[ii] <= p["low_call"]], ii[miss[ii] > p["low_call"]]))
        ok = self.af_ok[m].copy(); anyalt = np.zeros(len(Gc), bool); cross = {}
        for f, good, low in groups:
            if not len(good):                                   # no well-called individual: nothing can be fitted
                ok[:] = False; cross[f] = np.zeros(len(Gc), int); continue
            sub = Gc[:, good]; ok &= ((sub >= 0).mean(1) >= p["min_call"]); anyalt |= (sub > 0).any(1)
            k = (chrom, f)
            if k not in self._fit:
                self._fit[k] = fit_family(sub, p)
            cr, pv = self._fit[k]
            with np.errstate(invalid="ignore"):
                ok &= ~(pv < p["fit_p_min"])
            cross[f] = cr
        rows = np.where(ok & anyalt)[0]
        ref = np.concatenate([g for _f, g, _l in groups]) if groups else np.array([], int)
        kept = []
        for r in rows:                                          # thinning: drop a deletion close to and concordant with the last kept one
            if kept and pc[r] - pc[kept[-1]] <= p["thin_kb"] * 1000:
                q = kept[-1]; o = (Gc[r, ref] >= 0) & (Gc[q, ref] >= 0)
                if o.sum() >= p["thin_min_called"] and (Gc[r, ref][o] == Gc[q, ref][o]).mean() >= p["thin_concordance"]:
                    continue
            kept.append(r)
        elig = np.array(kept, int); rows = elig
        if len(rows) > p["max_columns"]:                        # even subsample along the chromosome, for legibility
            rows = rows[np.unique(np.round(np.linspace(0, len(rows) - 1, int(p["max_columns"]))).astype(int))]
        out = []
        for f, good, low in groups:
            o = good[order_individuals(Gc[rows][:, good])] if len(good) > 2 and len(rows) else good
            ff = np.concatenate([o, low]).astype(int)
            out.append({"fam": f, "cols": ff, "n_good": len(o), "n_low": len(low), "cross": cross[f][rows], "cross_eligible": cross[f][elig],
                        "M": Gc[rows][:, ff].T, "M_eligible": Gc[elig][:, ff].T, "miss": miss[ff]})
        return {"chrom": chrom, "fams": [g["fam"] for g in out], "n_eligible": len(elig), "eligible": elig, "rows": rows, "drawn": np.isin(elig, rows),
                "pos": pc[rows], "pos_eligible": pc[elig], "end_eligible": self.ms.end[m][elig], "groups": out, "n_chrom": int(len(m))}


def both_hom(A, B):
    """Deletions (columns) that are 1/1 in at least one individual of both families (rows = individuals)."""
    if not A.size or not B.size:
        return np.zeros(A.shape[1] if A.ndim == 2 else 0, bool)
    return (A == 2).any(0) & (B == 2).any(0) & (A >= 0).any(0) & (B >= 0).any(0)


def summarise(B):
    """Per family of a block: individuals, low-call individuals, expected crosses among the drawn deletions."""
    return {g["fam"]: {"n": int(g["M"].shape[0]), "n_low": int(g["n_low"]), "all_alike": int((g["cross"] == 0).sum()),
                       "split_1_1": int((g["cross"] == 1).sum()), "split_1_2_1": int((g["cross"] == 2).sum())} for g in B["groups"]}


# --------------------------------------------------------------------------- chains
def kept_links(links):
    """Links that are not rejected, best status first."""
    if links is None or not len(links):
        return pd.DataFrame(columns=["unit_a", "unit_b", "status", "min_n"])
    L = links[links.status.astype(str) != "rejected"].copy()
    L["_r"] = [STATUS_RANK.get(str(s), 9) for s in L.status]
    L["_ka"] = [_fkey(a) for a in L.unit_a]; L["_kb"] = [_fkey(b) for b in L.unit_b]
    return L.sort_values(["_r", "_ka", "_kb"], kind="stable").drop(columns=["_r", "_ka", "_kb"]).reset_index(drop=True)


def find_chains(links, junction, sizes, max_families, chain_links="kept"):
    """Chains of linked families. A chain = a centre family with its link partners.
    First the junction families of the network (two or more supported links; partners = the
    supported partners), then the other families on two or more kept links (partners = all kept
    partners, best status first). Ordered by support and size; a chain whose families are all
    part of an earlier chain is left out. Families are ordered with the centre in the middle and
    the partners by family size (smallest first)."""
    L = kept_links(links)
    st, adj = {}, {}
    for r in L.itertuples():
        a, b = str(r.unit_a), str(r.unit_b); st[(a, b)] = st[(b, a)] = str(r.status)
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    cand = []
    for f, ps in (junction or {}).items():
        ps = [str(x) for x in ps if (str(f), str(x)) in st]
        if len(ps) >= 2:
            cand.append((0, str(f), ps))
    done = {c[1] for c in cand}
    for f, ps in adj.items():
        if chain_links == "kept" and f not in done and len(ps) >= 2:
            cand.append((1, f, sorted(ps, key=lambda x: (STATUS_RANK.get(st[(f, x)], 9), _fkey(x)))))
    def key(c):
        tier, f, ps = c; r = [STATUS_RANK.get(st[(f, x)], 9) for x in ps]
        return (tier, -sum(x == 0 for x in r), -sum(x == 1 for x in r), -len(ps), -sizes.get(f, 0), _fkey(f))
    out, seen = [], []
    for tier, f, ps in sorted(cand, key=key):
        n_all = len(ps); ps = ps[:max(int(max_families) - 1, 2)]
        fs = {f, *ps}
        if any(fs <= s for s in seen):
            continue
        seen.append(fs)
        po = sorted(ps, key=lambda x: (sizes.get(x, 0), _fkey(x))); h = len(po) // 2
        out.append({"centre": f, "partners": ps, "order": po[:h] + [f] + po[h:], "links": [(f, x, st[(f, x)]) for x in po],
                    "kind": "junction family, supported links" if tier == 0 else "family on two or more kept links", "n_partners_all": n_all})
    return out


# --------------------------------------------------------------------------- sharing per chromosome
def chrom_sharing(ctx, E, links, cfg):
    """Sharing score per link and chromosome against its matched non-nominated pairs."""
    cols = ["unit_a", "unit_b", "status", "chrom", "observed", "null_median", "null_p05", "null_p95", "null_max", "n_null_pairs", "n_deletions"]
    a = cfg["audit"]; ms = ctx.ms
    need = {"auditable", "LLR_one", "nominated"}
    if E is None or not need <= set(E.columns) or not len(links):
        return pd.DataFrame(columns=cols)
    pool = E[E.auditable.astype(bool) & E.LLR_one.notna() & ~E.nominated.astype(bool)]
    key = {(str(r.unit_a), str(r.unit_b)): r for r in E.itertuples()}
    todo = []
    for r in links.itertuples():
        e = key.get((str(r.unit_a), str(r.unit_b))) or key.get((str(r.unit_b), str(r.unit_a)))
        if e is None or not bool(e.auditable) or not np.isfinite(e.LLR_one) or not len(pool):
            continue
        d = (pool.LLR_one - e.LLR_one).abs().sort_values(kind="stable").index[:a["k_match"]]
        todo.append((str(r.unit_a), str(r.unit_b), str(r.status), [(str(x), str(y)) for x, y in zip(pool.loc[d, "unit_a"], pool.loc[d, "unit_b"])]))
    if not todo:
        return pd.DataFrame(columns=cols)
    G = ms.G[ctx.af_ok]; ch = ctx.chrom[ctx.af_ok]
    fams = {f for x, y, _s, nul in todo for f in (x, y, *[u for p_ in nul for u in p_])}
    called = {f: (G[:, ctx.cols[f]] >= 0).any(1) for f in fams if f in ctx.cols}
    carrier = {f: (G[:, ctx.cols[f]] > 0).any(1) for f in fams if f in ctx.cols}
    rows = []
    for c in sorted(set(ch), key=chrom_key):
        m = ch == c; cache = {}

        def sh(x, y):
            if (x, y) not in cache:
                if x not in called or y not in called:
                    cache[(x, y)] = np.nan
                else:
                    k = called[x][m] & called[y][m]
                    cache[(x, y)] = float((carrier[x][m] & carrier[y][m] & k).sum() / max(((carrier[x][m] | carrier[y][m]) & k).sum(), 1))
            return cache[(x, y)]
        for x, y, s, nul in todo:
            nv = np.array([sh(u, v) for u, v in nul], float); nv = nv[np.isfinite(nv)]
            if not len(nv):
                continue
            rows.append({"unit_a": x, "unit_b": y, "status": s, "chrom": c, "observed": sh(x, y), "null_median": float(np.median(nv)),
                         "null_p05": float(np.quantile(nv, .05)), "null_p95": float(np.quantile(nv, .95)), "null_max": float(nv.max()),
                         "n_null_pairs": int(len(nv)), "n_deletions": int(m.sum())})
    return pd.DataFrame(rows, columns=cols)


# --------------------------------------------------------------------------- step
def run(ms, members, family_units, sizes, E, links, junction, cfg, max_chains=None):
    """Tables of the step: blocks (per family), their links, their deletions, and the sharing
    per chromosome. Blocks are computed for the first `max_chains` chains and for every kept
    link. Returns (dict of DataFrames, summary dict)."""
    ctx = Context(ms, members, family_units, cfg); p = ctx.p
    chrom = choose_chromosome(ms)
    L = kept_links(links)
    chains = find_chains(links, junction, sizes, p["max_chain_families"], p["chain_links"]); n_chains = len(chains)
    chains = chains[:max_chains] if max_chains is not None else chains
    defs = [(f"chain{i + 1}", "chain", c["kind"], c["order"], c["links"], c["centre"]) for i, c in enumerate(chains)]
    defs += [(f"link{i + 1}", "link", "kept link", [str(r.unit_a), str(r.unit_b)], [(str(r.unit_a), str(r.unit_b), str(r.status))], "") for i, r in enumerate(L.itertuples())]
    fam_rows, link_rows, mk = [], [], []
    for name, kind, desc, fams, lks, centre in defs:
        if chrom is None or any(f not in ctx.cols for f in fams):
            continue
        B = ctx.block(chrom, fams); sm = summarise(B); g = {x["fam"]: x for x in B["groups"]}
        for f in fams:
            fam_rows.append({"block": name, "kind": kind, "definition": desc, "chrom": chrom, "family": f, "role": "centre" if f == centre else ("partner" if kind == "chain" else "linked family"),
                             "n_individuals": sm[f]["n"], "n_low_call": sm[f]["n_low"], "n_deletions_chromosome": B["n_chrom"], "n_eligible": B["n_eligible"],
                             "n_drawn": int(len(B["rows"])), "n_all_alike": sm[f]["all_alike"], "n_split_1_1": sm[f]["split_1_1"], "n_split_1_2_1": sm[f]["split_1_2_1"]})
            mk.append(pd.DataFrame({"block": name, "family": f, "chrom": chrom, "pos": B["pos_eligible"], "end": B["end_eligible"],
                                    "expected_cross": [CROSS_NAMES[k] for k in g[f]["cross_eligible"]], "drawn": B["drawn"]}))
        for x, y, s in lks:
            link_rows.append({"block": name, "kind": kind, "chrom": chrom, "unit_a": x, "unit_b": y, "status": s, "n_eligible": B["n_eligible"], "n_drawn": int(len(B["rows"])),
                              "n_hom_alt_in_both_drawn": int(both_hom(g[x]["M"], g[y]["M"]).sum()), "n_hom_alt_in_both_eligible": int(both_hom(g[x]["M_eligible"], g[y]["M_eligible"]).sum())})
    T = {"blocks": pd.DataFrame(fam_rows, columns=["block", "kind", "definition", "chrom", "family", "role", "n_individuals", "n_low_call", "n_deletions_chromosome", "n_eligible", "n_drawn",
                                                   "n_all_alike", "n_split_1_1", "n_split_1_2_1"]),
         "links": pd.DataFrame(link_rows, columns=["block", "kind", "chrom", "unit_a", "unit_b", "status", "n_eligible", "n_drawn", "n_hom_alt_in_both_drawn", "n_hom_alt_in_both_eligible"]),
         "markers": pd.concat(mk, ignore_index=True) if mk else pd.DataFrame(columns=["block", "family", "chrom", "pos", "end", "expected_cross", "drawn"]),
         "sharing": chrom_sharing(ctx, E, L, cfg)}
    S = T["sharing"]
    summ = {"chromosome_shown": None if chrom is None else str(chrom), "chromosome_rule": "the chromosome with the most deletions", "n_chains": n_chains, "n_chain_blocks": len(chains), "n_link_blocks": int(len(L)),
            "n_links_with_sharing_per_chromosome": int(S.drop_duplicates(["unit_a", "unit_b"]).shape[0]) if len(S) else 0,
            "n_cells_above_null_median": int((S.observed > S.null_median).sum()) if len(S) else 0, "n_cells_above_null_p95": int((S.observed > S.null_p95).sum()) if len(S) else 0,
            "n_cells": int(len(S)), "parameters": {k: p[k] for k in sorted(p)}, "deletion_set": ms.name}
    return T, summ
