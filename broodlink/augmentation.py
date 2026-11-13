"""Step 7, second part: evidence on the shared parent from arrangement-tag deletions and
duplications, compared with the barcode and added to it; and the vote for a parent shared by
several links of one junction family.

Extracted from the original scripts (paths relative to the roots named in docs/EXTRACTED_METHOD.md):

- `FAM/_manuscript_figures_2026-10-28/code/compute_del_augmentation_v1.py` (`obligate`, `evidence`,
  `acc`) and its successor `HAT/per_inversion_first_pass/relabelled_augmentation/code/
  compute_augmentation_relabelled_v1.py` (`prep`, `ob`, `evidence(a, b, k, dmin, use)`, `acc`), which
  adds the duplications and the combined evidence with the same rule:
    * a deletion (duplication) is an OBLIGATE-CARRIER locus of a link when both families contain an
      individual homozygous for it and it is called in >= 2 individuals of each family; loci with an
      allele frequency outside 0.02-0.98 among the individuals in families are not used;
    * an arrangement-tag variant (carrier frequency differs between the two homokaryotype classes by
      at least delta) that is obligate for the link is evidence that the shared parent carries the
      tagged arrangement;
    * evidence class per link and interval: >= k obligate tag variants of class 0 -> "class 0",
      >= k of class 2 -> "class 2", both -> "class 0 + class 2", else "none"
      (k in 1, 2, 3; delta >= 0.6, 0.8, 0.9; setting shown: k = 2, delta >= 0.8);
    * comparison with the barcode state of the shared parent: supported (homokaryotype 0 with
      "class 0", homokaryotype 2 with "class 2", heterokaryotype with both), contradicted
      (a homokaryotype with evidence of the other arrangement), partial (heterokaryotype with one
      arrangement), and "fills" where the barcode is unresolved;
    * match rate on homozygous barcode states with evidence, against a null in which the evidence
      classes are permuted across the intervals inside each link (500 permutations, P with +1).
- `HAT/code/build_extdata_SP_DEL_augmentation_v5.py` and `HAT/figures/upgrades_2026-11-03/round_ed/
  ED_Fig_H5/ED_Fig_H5.py` (the vote, lines "G  vote logic unchanged from v5"): per interval, every
  link of the junction family votes with its barcode state and its evidence class; opposing votes
  are flagged, agreeing votes give a full state (0/0, 0/1, 1/1) or a partial one (at least one
  arrangement of class 0 / class 2).

The original labels STD / INV are the classes 0 / 2 of the structural-haplotype genotypes here;
they carry no ancestral meaning. The tag variants are those of step 6 (`inversions.tag_variants`),
so the smallest tag strength available is `inversions.tags.min_delta`.
"""
from __future__ import annotations
import re

import numpy as np
import pandas as pd

DEFAULT = {"k_sweep": [1, 2, 3], "delta_sweep": [0.6, 0.8, 0.9], "n_permutations": 500, "af_min": 0.02, "af_max": 0.98}
STATE = {"0/0": 0, "0/1": 1, "1/1": 2}
EVID = ["none", "class 0", "class 2", "class 0 + class 2"]
SUPER = {0: "0/0", 1: "0/1", 2: "1/1", 10: "at least one class 0", 12: "at least one class 2", -2: "opposing votes", -1: "no vote"}
EV_COLS = ["evidence_set", "k", "delta_min", "default_setting", "parent", "families", "link_status", "interval", "interval_number", "barcode_state",
           "n_obligate_tags_class_0", "n_obligate_tags_class_2", "evidence", "relation", "augmented_state"]
SW_COLS = ["evidence_set", "k", "delta_min", "default_setting", "links_pooled", "n_links", "supported", "contradicted", "partial", "n_informative", "unresolved",
           "filled", "n_homozygous_informative", "n_homozygous_match", "acc_homozygous", "null_mean", "p_perm", "n_permutations"]
SB_COLS = ["junction", "links", "interval", "interval_number", "junction_best_genotype", "junction_posterior", "junction_called", "votes_class_0",
           "votes_class_2", "votes_against_class_0", "votes_against_class_2", "super_state", "super_kind"]


def settings(cfg):
    """Setting shown (k, delta, calls per family) = config barcode.tag_evidence (shared with tag_evidence.py);
    sweep grid, permutations and allele-frequency window = DEFAULT, overridden by barcode.augmentation."""
    b = cfg.get("barcode", {}) or {}; e = b.get("tag_evidence", {}) or {}
    c = dict(DEFAULT); c.update(b.get("augmentation", {}) or {})
    c.update({"af_min": float(e.get("af_min", c["af_min"])), "af_max": float(e.get("af_max", c["af_max"])), "k": int(e.get("min_tags", 2)), "delta": float(e.get("min_delta", 0.8)), "min_called_per_family": int(e.get("min_called_per_family", 2))})
    c["k_sweep"] = sorted(set(list(c["k_sweep"]) + [c["k"]])); c["delta_sweep"] = sorted(set([float(x) for x in c["delta_sweep"]] + [c["delta"]]))
    return c


def interval_numbers(names):
    """Number of each interval in genome order (chromosome by its number, then start): {name: 1..n}."""
    def key(x):
        m = re.match(r"(?:scan_)?(.+):(\d+)-(\d+)$", str(x))
        if not m:
            return (2, 0, str(x), 0)
        d = re.findall(r"\d+", m.group(1))
        return (0 if d else 1, int(d[-1]) if d else 0, m.group(1), int(m.group(2)))
    return {n: i + 1 for i, n in enumerate(sorted(dict.fromkeys(names), key=key))}


def _family_tables(ms, members, c):
    """Per family: a homozygous carrier exists (per marker), number of called individuals; and the
    allele-frequency filter among the individuals in families (function `prep` of the original)."""
    si = {s: i for i, s in enumerate(ms.samples)}
    fams = [f for f in members if not str(f).startswith("S_")]
    cols = {f: [si[s] for s in members[f] if s in si] for f in fams}
    allc = [i for f in fams for i in cols[f]]
    if not allc:
        return None
    G = ms.G[:, allc]; nc = (G >= 0).sum(1)
    af = np.where(nc > 0, np.where(G >= 0, G, 0).sum(1) / (2.0 * np.maximum(nc, 1)), np.nan)
    sel = (af > c["af_min"]) & (af < c["af_max"])
    HOM = {f: (ms.G[:, cols[f]] == 2).any(1) if cols[f] else np.zeros(ms.n, bool) for f in fams}
    NC = {f: (ms.G[:, cols[f]] >= 0).sum(1) if cols[f] else np.zeros(ms.n, int) for f in fams}
    return sel, HOM, NC


def _tags_of_class(T, ms, numbers):
    """Tag variants of one SV class matched to the rows of its marker set: arrays (row, interval index, class, delta)."""
    key = {}
    for i, (ch, p, v) in enumerate(zip(ms.chrom, ms.pos, ms.ids)):
        key.setdefault((str(ch), int(p), str(v)), i)
    g = T[(T.svtype.astype(str) == ms.klass) & T.interval.isin(numbers)]
    row = np.array([key.get((str(ch), int(p), str(v)), -1) for ch, p, v in zip(g.chrom, g.start, g.variant_id)], int)
    j = np.array([numbers[x] - 1 for x in g.interval], int)
    cl = np.array([2 if str(x).strip().endswith("2") else 0 for x in g.tagged_arrangement], int)
    d = g.difference.to_numpy(float)
    ok = row >= 0
    return row[ok], j[ok], cl[ok], d[ok]


def _compare(p, cls):
    res = p >= 0; inf = res & (cls > 0)
    sup = ((p == 0) & (cls == 1)) | ((p == 2) & (cls == 2)) | ((p == 1) & (cls == 3))
    con = ((p == 0) & ((cls == 2) | (cls == 3))) | ((p == 2) & ((cls == 1) | (cls == 3)))
    return res, inf, sup, con


def _acc(store, rng=None):
    """Match rate on the homozygous barcode states with evidence, pooled over links (function `acc`)."""
    s_ = n_ = 0
    for p, cls in store:
        c2 = cls[rng.permutation(len(cls))] if rng is not None else cls
        r_ = (p >= 0) & (c2 > 0) & ((p == 0) | (p == 2))
        s_ += int((r_ & (((p == 0) & (c2 == 1)) | ((p == 2) & (c2 == 2)))).sum()); n_ += int(r_.sum())
    return s_, n_


def super_vote(src):
    """The vote of the original figure script, unchanged. src: list of (barcode states, evidence classes),
    one pair of arrays per link. Returns (state per interval, votes): state 0 / 1 / 2 full, 10 / 12 at least
    one arrangement of class 0 / 2, -2 opposing votes, -1 no vote."""
    n = len(src[0][0]) if src else 0
    SUP = np.full(n, -1); V = np.zeros((n, 4), int)
    for i in range(n):
        sv = iv = 0; sup_s = sup_i = 0
        for p, cls in src:
            if p[i] >= 0:
                sv += int(p[i] in (0, 1)); iv += int(p[i] in (1, 2)); sup_s += int(p[i] == 2); sup_i += int(p[i] == 0)
            sv += int(cls[i] in (1, 3)); iv += int(cls[i] in (2, 3))
        conflict = (sv > 0 and sup_s > 0) or (iv > 0 and sup_i > 0); hasS, hasI = sv > sup_s and sv > 0, iv > sup_i and iv > 0
        bc = [int(p[i]) for p, _ in src if p[i] >= 0]
        if conflict:
            SUP[i] = -2
        elif hasS and hasI:
            SUP[i] = 1
        elif hasS and 0 in bc:
            SUP[i] = 0
        elif hasI and 2 in bc:
            SUP[i] = 2
        elif hasS:
            SUP[i] = 10
        elif hasI:
            SUP[i] = 12
        V[i] = (sv, iv, sup_s, sup_i)
    return SUP, V


def augment(B, T, sv_sets, members, junctions, cfg, seed):
    """B: barcode table (step 7). T: tag variants (step 6). sv_sets: marker sets of class DEL / DUP.
    Returns (evidence table, sweep table, super-barcode table, summary)."""
    c = settings(cfg)
    EV = pd.DataFrame(columns=EV_COLS); SW = pd.DataFrame(columns=SW_COLS); SB = pd.DataFrame(columns=SB_COLS)
    summ = {"k": c["k"], "delta": c["delta"], "ran": False}
    if B is None or not len(B) or T is None or not len(T) or not sv_sets:
        summ["why"] = "no barcode, no tag variant or no deletion / duplication input"
        return EV, SW, SB, summ
    numbers = interval_numbers(B.interval.astype(str))
    names = sorted(numbers, key=numbers.get); n = len(names)
    L = B[B.kind.astype(str) == "link"]
    links = []
    for par, g in L.groupby("parent", sort=False):
        fams = str(g.families.iloc[0]).split(",")
        if len(fams) != 2 or fams[0] not in members or fams[1] not in members:
            continue
        st = dict(zip(g.interval.astype(str), g.call.astype(str)))
        p = np.array([STATE.get(st.get(x, "unresolved"), -1) for x in names], int)
        links.append({"parent": str(par), "a": fams[0], "b": fams[1], "status": str(g.link_status.iloc[0]), "p": p})
    if not links:
        summ["why"] = "no link in the barcode table"
        return EV, SW, SB, summ
    # ---- one marker set per SV class (the largest), its family tables and its tag variants
    per = {}
    for kl in ("DEL", "DUP"):
        ss = [m for m in sv_sets if m.klass == kl]
        if not ss:
            continue
        ms = max(ss, key=lambda m: m.n)
        ft = _family_tables(ms, members, c)
        if ft is None:
            continue
        row, j, cl, d = _tags_of_class(T, ms, numbers)
        if not len(row):
            continue
        per[kl] = {"sel": ft[0], "HOM": ft[1], "NC": ft[2], "row": row, "j": j, "cl": cl, "d": d}
    if not per:
        summ["why"] = "no tag variant matched to a deletion / duplication input"
        return EV, SW, SB, summ
    uses = [(k,) for k in per] + ([tuple(per)] if len(per) > 1 else [])
    mc = int(c["min_called_per_family"])

    def evidence(lk, k, dmin, use):
        n0 = np.zeros(n, int); n2 = np.zeros(n, int); nobl = 0
        for kl in use:
            q = per[kl]; a, b = lk["a"], lk["b"]
            ob = q["sel"] & (q["NC"][a] >= mc) & (q["NC"][b] >= mc) & q["HOM"][a] & q["HOM"][b]
            nobl += int(ob.sum())
            m = (q["d"] >= dmin) & ob[q["row"]]
            np.add.at(n0, q["j"][m & (q["cl"] == 0)], 1); np.add.at(n2, q["j"][m & (q["cl"] == 2)], 1)
        cls = np.where((n0 >= k) & (n2 >= k), 3, np.where(n2 >= k, 2, np.where(n0 >= k, 1, 0)))
        return cls, n0, n2, nobl
    sup_links = [lk for lk in links if lk["status"].startswith("supported")]
    pooled, pooled_label = (sup_links, "supported links") if sup_links else (links, "all links with a barcode (no supported link)")
    pooled_names = {lk["parent"] for lk in pooled}
    rng = np.random.default_rng([int(seed), 5]); nperm = int(c["n_permutations"])
    ev_rows, sw_rows, store0 = [], [], {}
    inv_state = {v: k for k, v in STATE.items()}
    for use in uses:
        uname = "+".join(use)
        for k in c["k_sweep"]:
            for dmin in c["delta_sweep"]:
                dflt = bool(k == c["k"] and abs(dmin - c["delta"]) < 1e-9)
                tot = {"supported": 0, "contradicted": 0, "partial": 0, "n_informative": 0, "unresolved": 0, "filled": 0}; st = []
                for lk in links:
                    cls, n0, n2, _nobl = evidence(lk, k, dmin, use); p = lk["p"]
                    res, inf, sup, con = _compare(p, cls)
                    if lk["parent"] in pooled_names:
                        st.append((p, cls))
                        tot["supported"] += int((inf & sup).sum()); tot["contradicted"] += int((inf & con).sum()); tot["partial"] += int((inf & ~sup & ~con).sum())
                        tot["n_informative"] += int(inf.sum()); tot["unresolved"] += int((~res).sum()); tot["filled"] += int((~res & (cls > 0)).sum())
                    if dflt:
                        store0[(uname, lk["parent"])] = cls
                    for i in range(n):
                        rel = "supported" if inf[i] and sup[i] else "contradicted" if inf[i] and con[i] else "partial" if inf[i] else "fills unresolved" if (not res[i] and cls[i] > 0) else ""
                        aug = inv_state[int(p[i])] if p[i] >= 0 else ("at least one " + EVID[cls[i]] if cls[i] in (1, 2) else "0/1 (class 0 + class 2)" if cls[i] == 3 else "no information")
                        ev_rows.append((uname, k, dmin, dflt, lk["parent"], lk["a"] + "," + lk["b"], lk["status"], names[i], i + 1, inv_state.get(int(p[i]), "unresolved"),
                                        int(n0[i]), int(n2[i]), EVID[cls[i]], rel, aug))
                s_, n_h = _acc(st)
                nl = np.array([(lambda x: x[0] / max(x[1], 1))(_acc(st, rng)) for _ in range(nperm)]) if nperm else np.array([])
                o = s_ / n_h if n_h else np.nan
                sw_rows.append({"evidence_set": uname, "k": k, "delta_min": dmin, "default_setting": dflt, "links_pooled": pooled_label, "n_links": len(st), **tot,
                                "n_homozygous_informative": n_h, "n_homozygous_match": s_, "acc_homozygous": o,
                                "null_mean": float(nl.mean()) if len(nl) else np.nan,
                                "p_perm": float((1 + (nl >= o).sum()) / (len(nl) + 1)) if n_h and len(nl) else np.nan, "n_permutations": nperm})
    EV = pd.DataFrame(ev_rows, columns=EV_COLS); SW = pd.DataFrame(sw_rows, columns=SW_COLS)
    # ---- super barcode: the links of a junction family vote (deletion evidence when there is a deletion input)
    vote_set = "DEL" if "DEL" in per else uses[0][0]
    by_fams = {frozenset((lk["a"], lk["b"])): lk for lk in sup_links}
    sb_rows, sb_sum = [], {}
    J = B[B.kind.astype(str) == "junction"]
    for f, partners in (junctions or {}).items():
        src_l = [by_fams[frozenset((f, q))] for q in partners if frozenset((f, q)) in by_fams and (vote_set, by_fams[frozenset((f, q))]["parent"]) in store0]
        if len(src_l) < 2:
            continue
        SUP, V = super_vote([(lk["p"], store0[(vote_set, lk["parent"])]) for lk in src_l])
        g = J[J.parent.astype(str) == "J_" + str(f)]
        bg = dict(zip(g.interval.astype(str), g.best_genotype.astype(str))); po = dict(zip(g.interval.astype(str), g.posterior.astype(float)))
        ca = dict(zip(g.interval.astype(str), g.call.astype(str)))
        code = np.array([STATE.get(bg.get(x, ""), -1) for x in names]); conf = np.array([ca.get(x, "unresolved") != "unresolved" for x in names])
        both = np.isin(SUP, (0, 1, 2)) & conf
        sb_sum[str(f)] = {"links": [lk["parent"] for lk in src_l], "evidence_set": vote_set, "n_full": int(np.isin(SUP, (0, 1, 2)).sum()), "n_partial": int(np.isin(SUP, (10, 12)).sum()),
                          "n_opposing_votes": int((SUP == -2).sum()), "n_no_vote": int((SUP == -1).sum()), "n_junction_called": int(conf.sum()),
                          "n_both": int(both.sum()), "n_agree": int((SUP[both] == code[both]).sum())}
        for i in range(n):
            sb_rows.append(("J_" + str(f), ",".join(lk["parent"] for lk in src_l), names[i], i + 1, bg.get(names[i], ""), po.get(names[i], np.nan), bool(conf[i]),
                            int(V[i, 0]), int(V[i, 1]), int(V[i, 2]), int(V[i, 3]), SUPER[int(SUP[i])],
                            "full" if SUP[i] in (0, 1, 2) else "partial" if SUP[i] in (10, 12) else "opposing votes" if SUP[i] == -2 else "no vote"))
    SB = pd.DataFrame(sb_rows, columns=SB_COLS)
    d0 = SW[SW.default_setting.astype(bool)]
    per_link = {}
    for lk in links:
        cls = store0.get((vote_set, lk["parent"]))
        if cls is not None:
            per_link[lk["parent"]] = {"barcode": int((lk["p"] >= 0).sum()), "added": int(((lk["p"] < 0) & (cls > 0)).sum()), "status": lk["status"]}
    summ.update({"ran": True, "evidence_sets": ["+".join(u) for u in uses], "n_intervals": n, "n_links": len(links), "links_pooled": pooled_label, "n_links_pooled": len(pooled),
                 "n_tag_variants_matched": {k: int(len(v["row"])) for k, v in per.items()},
                 "default": {r.evidence_set: {"n_homozygous_informative": int(r.n_homozygous_informative), "n_homozygous_match": int(r.n_homozygous_match),
                                              "acc_homozygous": None if pd.isna(r.acc_homozygous) else float(r.acc_homozygous),
                                              "null_mean": None if pd.isna(r.null_mean) else float(r.null_mean), "p_perm": None if pd.isna(r.p_perm) else float(r.p_perm),
                                              "filled": int(r.filled)} for r in d0.itertuples()},
                 "per_link": per_link, "super_barcode": sb_sum, "settings": c})
    return EV, SW, SB, summ


def add_to_run(run, B, members, junctions):
    """Called by step_barcode after the barcode: writes the evidence tables and the recurrence test of
    the junction families; returns the entries added to the summary of the step. Never stops a run."""
    from . import recurrence
    out = {}
    sv = [m for m in run.sets if m.klass in ("DEL", "DUP")]
    try:
        T = run.t.get("haplotypes", {}).get("tags")
        EV, SW, SB, s = augment(B, T, sv, members, junctions, run.cfg, run.seed)
        out["augmentation"] = s
        if s.get("ran"):
            run.put("barcode", "evidence_by_setting", EV, "tag_evidence_by_setting.tsv", "Evidence on the shared parent from arrangement-tag deletions / duplications, per evidence set (DEL, DUP, DEL+DUP), minimum number of tag variants k, minimum tag strength delta, link and interval: obligate tag variants of each class (both families hold a homozygous carrier), evidence class, relation to the barcode state (supported / contradicted / partial / fills unresolved) and the augmented state. default_setting marks the setting shown in the figure.")
            run.put("barcode", "evidence_sweep", SW, "tag_evidence_sweep.tsv", "Sweep over k and delta per evidence set, pooled over the links named in links_pooled: counts of supported / contradicted / partial / filled intervals, match rate on homozygous barcode states with evidence, mean of the permutation null (evidence classes permuted across intervals inside each link) and its P.")
            if len(SB):
                run.put("barcode", "super_barcode", SB, "junction_super_barcode.tsv", "Per junction family and interval: the vote of its supported links (barcode state and deletion evidence of each link) next to the joint barcode of the junction (J_ rows of the barcode table). super_state: full genotype, at least one arrangement of a class, opposing votes, or no vote.")
    except Exception as exc:                                    # an add-on of the barcode never stops the step
        run.note(f"barcode: tag-variant evidence not computed ({type(exc).__name__}: {exc})")
        out["augmentation"] = {"ran": False, "why": f"{type(exc).__name__}: {exc}"}
    try:
        dls = [m for m in run.sets if m.klass == "DEL"]
        RC, s = recurrence.run(max(dls, key=lambda m: m.n) if dls else None, members, junctions, run.cfg, run.seed)
        out["recurrence"] = s
        if len(RC):
            run.put("barcode", "recurrence", RC, "junction_recurrence.tsv", "Per junction family (its first two supported links): does one parental pair of the junction family recur in both partner families? Slot scores from rare deletions clustered by co-transmission inside the junction family, the statistic, and a null of pseudo-chains built from families not linked to the junction (mean, maximum, number, empirical P). A candidate, not a proof.")
    except Exception as exc:
        run.note(f"barcode: recurrence test of the junction families not computed ({type(exc).__name__}: {exc})")
        out["recurrence"] = {"ran": False, "why": f"{type(exc).__name__}: {exc}"}
    return out
