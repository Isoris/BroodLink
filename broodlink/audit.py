"""Step 4: audit of nominated family pairs on markers not used for nomination.

Extracted from sp_lib.oh_llr_vec, run_audit.gate_stats / run and make_tables.build.
Gate 1: per-marker log-likelihood ratio of "one shared parent" and of "both parents shared",
        each against "no shared parent"; compatible when LLR_one > 0 and LLR_both < 0
        (an exclusion test: passing means one shared parent is not excluded).
Gate 2: sharing of inherited structural variants (carrier in both families / carrier in either),
        once per variant class, against the k non-nominated auditable pairs nearest in the
        Gate-1 score; a link is supported when its sharing exceeds the median of the
        Gate-1-compatible nominated pairs in at least one class.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .io import allele_freq, detail

SV_LIKE = ("DEL", "DUP", "INS", "INV", "SV")


def transmission_tensor(err=0.03):
    """T[gP, gQ, g_obs] = P(observed offspring genotype | parental genotypes), with genotype error."""
    alleles = {0: (0, 0), 1: (0, 1), 2: (1, 1)}
    T = np.zeros((3, 3, 3))
    for p in range(3):
        for q in range(3):
            for a in alleles[p]:
                for b in alleles[q]:
                    T[p, q, a + b] += 0.25
    E = np.array([[(1 - err) ** 2, 2 * err * (1 - err), err ** 2],
                  [err * (1 - err), (1 - err) ** 2 + err ** 2, err * (1 - err)],
                  [err ** 2, 2 * err * (1 - err), (1 - err) ** 2]])
    return np.einsum("pqg,gh->pqh", T, E)


def hwe(af):
    return np.stack([(1 - af) ** 2, 2 * af * (1 - af), af ** 2], 1)


def family_cells(G, idx, TE):
    """sites x 3 x 3: product over the offspring of T[gP, gQ, g_i] (missing genotype -> 1)."""
    TEp = np.concatenate([TE, np.ones((3, 3, 1))], axis=2)
    g = G[:, idx].astype(int); g = np.where(g < 0, 3, g)
    out = np.empty((G.shape[0], 3, 3))
    for gp in range(3):
        for gq in range(3):
            out[:, gp, gq] = np.prod(TEp[gp, gq][g], axis=1)
    return out


def _mean_llr(x, base):
    ok = (x > 0) & (base > 0)
    return float(np.where(ok, np.log(np.where(ok, x, 1) / np.where(ok, base, 1)), 0).sum() / len(x))


def gate1_scores(ms, fam_members, pairs, cfg, seed):
    """LLR_one and LLR_both per family pair on one marker set."""
    a = cfg["audit"]
    af = allele_freq(ms.G)
    keep = (af > a["af_min"]) & (af < a["af_max"])
    G = ms.G[keep]; af = af[keep]
    if G.shape[0] > a["gate1_max_markers"]:
        sel = np.sort(np.random.default_rng(seed).choice(G.shape[0], a["gate1_max_markers"], replace=False))
        G, af = G[sel], af[sel]
    if G.shape[0] == 0:
        return pd.DataFrame(columns=["unit_a", "unit_b", "LLR_one", "LLR_both"]), 0
    si = {s: i for i, s in enumerate(ms.samples)}
    TE = transmission_tensor(a["genotype_error"]); pri = hwe(af)
    cache = {}

    def fam(f):
        if f not in cache:
            C = family_cells(G, [si[s] for s in fam_members[f] if s in si], TE)
            cache[f] = (C, np.einsum("sq,spq->sp", pri, C))
        return cache[f]
    rows = []
    for x, y in pairs:
        CA, LA = fam(x); CB, LB = fam(y)
        one = (pri * LA * LB).sum(1)
        base = (pri * LA).sum(1) * (pri * LB).sum(1)
        both = np.einsum("sp,sq,spq->s", pri, pri, CA * CB)
        rows.append({"unit_a": x, "unit_b": y, "LLR_one": _mean_llr(one, base), "LLR_both": _mean_llr(both, base)})
    return pd.DataFrame(rows), int(G.shape[0])


def sharing(ms, fam_members, pairs, cfg):
    """Fraction of loci with a carrier in either family that have a carrier in both (loci called in both)."""
    a = cfg["audit"]
    af = allele_freq(ms.G)
    G = ms.G[(af > a["af_min"]) & (af < a["af_max"])]
    si = {s: i for i, s in enumerate(ms.samples)}
    called, carrier = {}, {}
    for f in {u for p in pairs for u in p}:
        idx = [si[s] for s in fam_members[f] if s in si]
        called[f] = (G[:, idx] >= 0).any(1); carrier[f] = (G[:, idx] > 0).any(1)
    out = []
    for x, y in pairs:
        k = called[x] & called[y]
        out.append(float((carrier[x] & carrier[y] & k).sum() / max(((carrier[x] | carrier[y]) & k).sum(), 1)))
    return np.array(out), int(G.shape[0])


def run_audit(nom, unit_of, audit_sets, cfg, seed, single_input=False):
    """nom: nomination table. audit_sets: MarkerSets never used for nomination.
    Returns (edge table, summary dict)."""
    a = cfg["audit"]; mn = a["min_family_size"]
    E = nom.copy()
    fam_members = {}
    for s, u in unit_of.items():
        fam_members.setdefault(u, []).append(s)
    isfam = (E.unassigned_sides == 0)
    E["auditable"] = isfam & (E.min_n >= mn)
    summ = {"min_family_size": mn, "audit_classes": [], "gate1_class": None, "gate2": {}}
    suffix = " (single-input, provisional)" if single_input else ""
    if not audit_sets:
        E["link_class"] = np.where(E.nominated, "nominated, not audited", "not nominated")
        E["gate1"] = np.nan; E["LLR_one"] = np.nan; E["LLR_both"] = np.nan; E["evidence"] = ""
        summ["ran"] = False
        return E, summ
    summ["ran"] = True
    pref = sorted(audit_sets, key=lambda m: (m.klass != "DEL", -m.n))
    g1set = pref[0]
    scored = E[isfam & (E.min_n >= 2)]
    pairs = list(zip(scored.unit_a, scored.unit_b))
    detail(f"Gate 1 (one shared parent vs none or both): scoring {len(pairs):,} family pairs on {g1set.name}")
    st, n1 = gate1_scores(g1set, fam_members, pairs, cfg, seed)
    summ["gate1_class"] = g1set.name; summ["gate1_markers"] = n1
    detail(f"Gate 1 used {n1:,} markers")
    E = E.merge(st, on=["unit_a", "unit_b"], how="left")
    E["gate1"] = np.where(E.LLR_one.notna(), (E.LLR_one > 0) & (E.LLR_both < 0), np.nan)
    g2sets = [m for m in audit_sets if m.klass in SV_LIKE]
    sup_cols = []
    for gi, m in enumerate(g2sets):
        k = m.klass if sum(x.klass == m.klass for x in g2sets) == 1 else f"{m.klass}_{gi + 1}"
        col = f"share_{k}"
        detail(f"Gate 2 (sharing of inherited {m.klass}): {m.name}, against matched non-nominated pairs")
        sh, n2 = sharing(m, fam_members, pairs, cfg)
        E = E.merge(pd.DataFrame({"unit_a": [p[0] for p in pairs], "unit_b": [p[1] for p in pairs], col: sh}), on=["unit_a", "unit_b"], how="left")
        pool = E[E.auditable & E.LLR_one.notna()]
        g1 = pool[pool.nominated & (pool.gate1 == 1)]; npool = pool[~pool.nominated]
        info = {"marker_set": m.name, "n_markers": n2, "n_gate1_pass": len(g1), "n_null_pool": len(npool)}
        if len(g1) and len(npool):
            matched, pairp = [], []
            for r in g1.itertuples():
                d = (npool.LLR_one - r.LLR_one).abs().sort_values(kind="stable").index[:a["k_match"]]
                v = npool.loc[d, col].to_numpy(); matched.extend(v.tolist())
                pairp.append((1 + (v >= getattr(r, col)).sum()) / (1 + len(v)))
            matched = np.array(matched); obs = g1[col].to_numpy(); omed = float(np.median(obs))
            info.update({"threshold_median_of_gate1_pass": omed, "null_median": float(np.median(matched)),
                         "n_null_values": len(matched), "p_one_sided": float((1 + (matched >= omed).sum()) / (1 + len(matched))),
                         "null_p99": float(np.quantile(matched, 0.99))})
            E.loc[g1.index, f"pairP_{k}"] = pairp
            E[f"above_median_{k}"] = (E[col] > omed) & E.nominated & (E.gate1 == 1) & isfam
            E[f"above_null_p99_{k}"] = (E[col] > info["null_p99"]) & E.nominated & (E.gate1 == 1) & isfam
            sup_cols.append((k, f"above_median_{k}"))
        elif len(g1):
            omed = float(np.median(g1[col]))
            info.update({"threshold_median_of_gate1_pass": omed, "note": "no non-nominated auditable pair: matched null not available"})
            E[f"above_median_{k}"] = (E[col] > omed) & E.nominated & (E.gate1 == 1) & isfam
            sup_cols.append((k, f"above_median_{k}"))
        summ["gate2"][k] = info
    summ["audit_classes"] = [m.name for m in audit_sets]

    def klass(r):
        if not r.nominated:
            return "not nominated"
        if not r.auditable:
            return "provisional (not testable: unit below %d individuals)" % mn
        if r.gate1 != 1:
            return "rejected (failed Gate 1)"
        ev = [k for k, c in sup_cols if bool(getattr(r, c))]
        return ("supported" if ev else "Gate-1 only") + suffix
    E["evidence"] = [" + ".join(k for k, c in sup_cols if bool(r[c]) and r["nominated"] and r["auditable"] and r["gate1"] == 1) for _, r in E.iterrows()] if sup_cols else ""
    E["link_class"] = [klass(r) for r in E.itertuples()]
    if not g2sets:
        summ["gate2_note"] = "no structural-variant class in the audit set: Gate 2 not run, links can reach 'Gate-1 only' at most"
    return E, summ


def matched_null(E, col, k_match):
    """Observed sharing of the Gate-1-compatible nominated auditable pairs and the pooled matched
    null (for each of them, the k_match non-nominated auditable pairs nearest in LLR_one), as
    computed inside run_audit. Used by the figures; returns (observed, null) arrays."""
    pool = E[E.auditable & E.LLR_one.notna() & E[col].notna()]
    g1 = pool[pool.nominated & (pool.gate1 == 1)]; npool = pool[~pool.nominated]
    matched = []
    for r in g1.itertuples():
        d = (npool.LLR_one - r.LLR_one).abs().sort_values(kind="stable").index[:k_match]
        matched.extend(npool.loc[d, col].tolist())
    return g1[col].to_numpy(float), np.array(matched, float)
