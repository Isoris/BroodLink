"""Step 7, second part: what the arrangement tag variants say about a shared parent.

Extracted from HAT/per_inversion_first_pass/relabelled_augmentation/code/
compute_augmentation_relabelled_v1.py (functions `prep`, `ob`, `evidence`; same rule as
compute_del_augmentation_v1.py) and from the per-link rows of
HAT/figures/upgrades_2026-11-03/round_ed/ED_Fig_H6/ED_Fig_H6.py:
  * obligate carrier: a tag variant is obligate for the shared parent of a link when BOTH
    families hold an individual homozygous for it (each parent of such an individual carries
    it, and the shared parent is a parent of both families); the variant must be called in at
    least `min_called_per_family` individuals of each family;
  * evidence class of a link at an interval, from the obligate tag variants whose carrier
    frequency differs between the homokaryotype classes by at least `min_delta`: at least
    `min_tags` tagging class 0 -> the parent carries a class-0 arrangement; at least `min_tags`
    tagging class 2 -> a class-2 arrangement; both -> both;
  * agreement with the barcode call: supported (0/0 with class 0, 1/1 with class 2, 0/1 with
    both), contradicted (0/0 with class 2 or both; 1/1 with class 0 or both), or the evidence
    fills a state the barcode left unresolved.
The original settings were chosen in a sweep (k = 1, 2, 3; difference >= 0.6, 0.8, 0.9) with
a permutation check; the sweep and the check are not part of the run. The original catalogue
required >= 20 and >= 10 called homokaryotypes per class and an allele-frequency window; here
the tag variants are those of step 6 (config inversions.tags).
"""
from __future__ import annotations
import numpy as np
import pandas as pd

CALL_P = {"0/0": 0, "0/1": 1, "1/1": 2}
CLASS_NAME = {0: "", 1: "class 0", 2: "class 2", 3: "class 0 + class 2"}
OB_COLS = ["link", "families", "interval", "variant_id", "svtype", "chrom", "start", "end", "length_bp", "tagged_arrangement", "difference", "counts_as_evidence"]
EV_COLS = ["link", "families", "interval", "svtype", "n_obligate_class_0", "n_obligate_class_2", "evidence_class", "barcode_call", "agreement", "augmented_state"]


def settings(cfg):
    e = cfg["barcode"].get("tag_evidence", {})
    return int(e.get("min_called_per_family", 2)), int(e.get("min_tags", 2)), float(e.get("min_delta", 0.8))


def obligate(tags, sv_sets, members, links, cfg):
    """One row per link and tag variant that is obligate for the shared parent of the link."""
    mc, _k, dmin = settings(cfg)
    if tags is None or not len(tags) or not links:
        return pd.DataFrame(columns=OB_COLS)
    rows = []
    for ms in sv_sets:
        tg = tags[tags.svtype == ms.klass]
        if not len(tg):
            continue
        key = {(str(c), int(p), str(i)): r for r, (c, p, i) in enumerate(zip(ms.chrom, ms.pos, ms.ids))}
        row = np.array([key.get((str(c), int(p), str(i)), -1) for c, p, i in zip(tg.chrom, tg.start, tg.variant_id)])
        tg = tg[row >= 0]; G = ms.G[row[row >= 0]]
        if not len(tg):
            continue
        six = {s: i for i, s in enumerate(ms.samples)}
        hom, nc = {}, {}
        for f in {x for a, b, _s in links for x in (a, b)}:
            c = [six[s] for s in members.get(f, []) if s in six]
            sub = G[:, c]; hom[f] = (sub == 2).any(1); nc[f] = (sub >= 0).sum(1)
        for a, b, _s in links:
            o = (nc[a] >= mc) & (nc[b] >= mc) & hom[a] & hom[b]
            for t in tg[o].itertuples():
                rows.append({"link": f"{a}x{b}", "families": f"{a},{b}", "interval": t.interval, "variant_id": t.variant_id, "svtype": t.svtype, "chrom": t.chrom, "start": int(t.start),
                             "end": int(t.end), "length_bp": int(t.length_bp), "tagged_arrangement": t.tagged_arrangement, "difference": float(t.difference),
                             "counts_as_evidence": bool(t.difference >= dmin)})
    O = pd.DataFrame(rows, columns=OB_COLS)
    return O.drop_duplicates(["link", "interval", "svtype", "variant_id", "start"]).reset_index(drop=True)


def evidence(O, B, links, intervals, cfg):
    """One row per link, interval and SV class (and 'all' when several classes tag): counts of
    obligate tag variants, evidence class, the barcode call and the augmented state."""
    _mc, k, _d = settings(cfg)
    if not links or not len(intervals):
        return pd.DataFrame(columns=EV_COLS)
    Bl = B[B.kind == "link"] if len(B) else B
    call = {(r.parent, r.interval): r.call for r in Bl.itertuples()}
    E = O[O.counts_as_evidence] if len(O) else O
    kinds = sorted(set(E.svtype)) if len(E) else []
    kinds = kinds + (["all"] if len(kinds) > 1 else [])
    rows = []
    for kind in kinds:
        e = E if kind == "all" else E[E.svtype == kind]
        n = e.groupby(["link", "interval", "tagged_arrangement"]).size().to_dict()
        for a, b, _s in links:
            nm = f"{a}x{b}"
            for iv in intervals:
                n0, n2 = int(n.get((nm, iv, "class 0"), 0)), int(n.get((nm, iv, "class 2"), 0))
                cls = 3 if (n0 >= k and n2 >= k) else 2 if n2 >= k else 1 if n0 >= k else 0
                c = call.get((nm, iv), "unresolved"); p = CALL_P.get(c, -1)
                if cls == 0:
                    agree = ""
                elif p < 0:
                    agree = "fills an unresolved state"
                elif (p == 0 and cls == 1) or (p == 2 and cls == 2) or (p == 1 and cls == 3):
                    agree = "supported"
                elif (p == 0 and cls in (2, 3)) or (p == 2 and cls in (1, 3)):
                    agree = "contradicted"
                else:
                    agree = "compatible (heterokaryotype call, one class seen)"
                aug = c if p >= 0 else (f"at least one {CLASS_NAME[cls]} arrangement" if cls in (1, 2) else "class 0 and class 2 arrangements" if cls == 3 else "no information")
                rows.append({"link": nm, "families": f"{a},{b}", "interval": iv, "svtype": kind, "n_obligate_class_0": n0, "n_obligate_class_2": n2, "evidence_class": CLASS_NAME[cls],
                             "barcode_call": c, "agreement": agree, "augmented_state": aug})
    return pd.DataFrame(rows, columns=EV_COLS)


def barcode_tables(run, links=None):
    """(key, table, file name, description) for step 7; empty when step 6 found no tag variant."""
    tags = run.t.get("haplotypes", {}).get("tags")
    B = run.t.get("barcode", {}).get("barcode")
    if tags is None or not len(tags) or B is None or not len(B):
        return []
    if links is None:
        Bl = B[B.kind == "link"].drop_duplicates("parent")
        links = [tuple(str(r.families).split(",")[:2]) + (r.link_status,) for r in Bl.itertuples()]
    _u, _f, members = run.units()
    sv = [m for m in run.sets if m.klass in set(tags.svtype)]
    ints = list(dict.fromkeys(B.interval))
    mc, k, dmin = settings(run.cfg)
    O = obligate(tags[tags.interval.isin(ints)], sv, members, links, run.cfg)
    E = evidence(O, B, links, ints, run.cfg)
    return [("obligate_tags", O, "obligate_tag_variants.tsv", f"Arrangement tag variants that the shared parent of a link must carry: both families hold an individual homozygous for the variant (called in at least {mc} individuals of each family). counts_as_evidence: carrier-frequency difference between the homokaryotype classes >= {dmin}."),
            ("tag_evidence", E, "barcode_tag_evidence.tsv", f"Per link, barcode interval and SV class: obligate tag variants tagging class 0 and class 2 (difference >= {dmin}), evidence class (at least {k} per class), the barcode call, their agreement, and the augmented state (the barcode call, or what the tag variants show where the barcode is unresolved).")]
