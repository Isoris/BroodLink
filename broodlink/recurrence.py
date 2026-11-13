"""Step 7, third part: does one parental pair of a junction family recur on two of its links?

Extracted from `FAM/_project_docs/code/refine_LP.py` (function `main`; `load_gt_matrix` of
`common.py` for the allele-frequency filter). A family on two links has two parental slots.
Rare, well-called deletions (allele frequency 0.01-0.10 among the individuals in families,
called in >= 0.80 of them) are clustered by co-transmission inside the junction family: the
signs of the two leading eigenvectors of their correlation matrix split the informative
deletions into two pairs of founder groups ("slots"). The score of a slot in another family is
the mean, over that family's individuals, of the larger of the two group carrier rates (an
offspring inherits one homolog of the pair). Statistic = the smaller, over the two partner
families, of the best slot score. Null = pseudo-chains F_x - junction - F_y built from pairs
of families that are not linked to the junction (at least the minimum family size); P with +1.

Differences from the original: the junction family and its two partners come from the network
of the run (the first two partners in the order of the network summary; the original took the
first two rows of the edge table); the families of a null pair are taken in sorted order, and
when there are more pairs than `n_null_max` a seeded random subset is drawn (the original
stopped after the first `nperm` pairs of an unordered set).
"""
from __future__ import annotations
import itertools

import numpy as np
import pandas as pd

DEFAULT = {"af_min": 0.01, "af_max": 0.10, "min_call_rate": 0.80, "n_null_max": 200, "p_max": 0.05}
COLS = ["junction", "partner_1", "partner_2", "n_partners_of_junction", "n_rare_deletions", "n_individuals_junction", "n_informative_deletions",
        "slot1_partner_1", "slot2_partner_1", "slot1_partner_2", "slot2_partner_2", "same_slot_favoured", "statistic", "null_mean", "null_max",
        "n_null", "n_null_families", "P", "verdict", "note"]


def settings(cfg):
    c = dict(DEFAULT); c.update((cfg.get("barcode", {}) or {}).get("recurrence", {}) or {})
    return c


def _key(f):
    s = str(f); d = "".join(ch for ch in s if ch.isdigit())
    return (int(d) if d else 0, s)


def test_junction(ms, members, junction, partners, cfg, seed):
    c = settings(cfg); mn = int(cfg["audit"]["min_family_size"])
    fams = {f: v for f, v in members.items() if not str(f).startswith("S_")}
    si = {s: i for i, s in enumerate(ms.samples)}
    allc = [si[s] for f in fams for s in fams[f] if s in si]
    row = {k: np.nan for k in COLS}
    row.update({"junction": junction, "partner_1": partners[0], "partner_2": partners[1], "verdict": "not tested", "note": ""})
    if not allc:
        row["note"] = "no individual of a family in the deletion input"; return row
    G0 = ms.G[:, allc]; nc = (G0 >= 0).sum(1)
    af = np.where(nc > 0, np.where(G0 >= 0, G0, 0).sum(1) / (2.0 * np.maximum(nc, 1)), np.nan)
    keep = (af > c["af_min"]) & (af < c["af_max"]) & ((G0 >= 0).mean(1) >= c["min_call_rate"])
    G = ms.G[keep]
    J = [si[s] for s in fams[junction] if s in si]
    row["n_rare_deletions"] = int(keep.sum()); row["n_individuals_junction"] = len(J)
    sub = G[:, J] if J else np.zeros((G.shape[0], 0), np.int8)
    inf = np.where(((sub >= 0).sum(1) >= max(3, len(J) - 2)) & ((sub > 0).sum(1) >= 2) & ((sub > 0).sum(1) <= len(J) - 2))[0]
    row["n_informative_deletions"] = int(len(inf))
    if len(inf) < 4:
        row["note"] = "fewer than 4 informative rare deletions in the junction family"; return row
    C = (sub[inf] > 0).astype(float); C[sub[inf] < 0] = np.nan
    Z = np.nan_to_num(C - np.nanmean(C, 1, keepdims=True))
    Sm = Z @ Z.T; d = np.sqrt(np.diag(Sm))
    Rm = np.nan_to_num(Sm / np.outer(np.maximum(d, 1e-9), np.maximum(d, 1e-9)))
    _w, v = np.linalg.eigh(Rm)
    slot1, slot2 = v[:, -1] > 0, v[:, -2] > 0
    if min(slot1.sum(), (~slot1).sum(), slot2.sum(), (~slot2).sum()) == 0:
        row["note"] = "the informative deletions do not split into two groups per slot"; return row

    def pair_score(grpA, grpB, target):
        T = [si[s] for s in fams.get(target, []) if s in si]
        if not T:
            return np.nan
        x = G[np.ix_(inf, T)]; cc = (x > 0).astype(float); cc[x < 0] = np.nan
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return float(np.nanmean(np.maximum(np.nanmean(cc[grpA], 0), np.nanmean(cc[grpB], 0))))
    sc = [(pair_score(slot1, ~slot1, p), pair_score(slot2, ~slot2, p)) for p in partners[:2]]
    if any(np.isnan(x) for t in sc for x in t):
        row["note"] = "a partner family has no called genotype at the informative deletions"; return row
    fav = ["slot1" if a > b else "slot2" for a, b in sc]
    stat = float(min(max(t) for t in sc))
    linked = set(partners[:2]) | {junction}
    cand = sorted([f for f in fams if f not in linked and len(fams[f]) >= mn], key=_key)
    best = {}
    for f in cand:
        best[f] = max(pair_score(slot1, ~slot1, f), pair_score(slot2, ~slot2, f))
    pairs = [(x, y) for x, y in itertools.combinations(cand, 2) if not (np.isnan(best[x]) or np.isnan(best[y]))]
    if len(pairs) > c["n_null_max"]:
        rng = np.random.default_rng([int(seed), 11]); pick = sorted(rng.choice(len(pairs), int(c["n_null_max"]), replace=False))
        pairs = [pairs[i] for i in pick]
    null = np.array([min(best[x], best[y]) for x, y in pairs])
    row.update({"slot1_partner_1": sc[0][0], "slot2_partner_1": sc[0][1], "slot1_partner_2": sc[1][0], "slot2_partner_2": sc[1][1],
                "same_slot_favoured": bool(len(set(fav)) == 1), "statistic": stat, "n_null": int(len(null)), "n_null_families": len(cand)})
    if not len(null):
        row["note"] = f"no pair of families of at least {mn} individuals that are not linked to the junction: no null"; return row
    P = float((np.sum(null >= stat) + 1) / (len(null) + 1))
    row.update({"null_mean": float(null.mean()), "null_max": float(null.max()), "P": P, "verdict": "candidate" if P < c["p_max"] else "not supported"})
    return row


def run(ms, members, junctions, cfg, seed):
    """One row per junction family. ms: deletion marker set (None = no test)."""
    rows = []
    if ms is not None:
        for f, partners in (junctions or {}).items():
            if len(partners) < 2 or f not in members or any(p not in members for p in partners[:2]):
                continue
            r = test_junction(ms, members, f, list(partners), cfg, seed); r["n_partners_of_junction"] = len(partners)
            if len(partners) > 2:
                r["note"] = (r["note"] + "; " if r["note"] else "") + f"the junction family has {len(partners)} supported links; the first two partners are tested"
            rows.append(r)
    T = pd.DataFrame(rows, columns=COLS)
    ok = T[T.verdict != "not tested"] if len(T) else T
    return T, {"ran": bool(len(ok)), "n_junctions_tested": int(len(ok)), "n_candidate": int((ok.verdict == "candidate").sum()) if len(ok) else 0,
               "settings": settings(cfg), "why": "" if len(ok) else ("no deletion input" if ms is None else"no junction family with two supported links, or too few rare deletions")}
