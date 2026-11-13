"""Comparison of results with a reference (the truth of the synthetic cohort, or a reference run)."""
from __future__ import annotations
import itertools

import pandas as pd


def partition(M):
    out = {}
    for s, f in zip(M.individual.astype(str), M.family.astype(str)):
        if f != "unassigned":
            out.setdefault(f, set()).add(s)
    return out


def family_scores(M, ref):
    """ref: dict individual -> reference family ('unassigned' or '' when none)."""
    A = {frozenset(v) for v in partition(M).values()}
    T = {}
    for s, f in ref.items():
        if f and f != "unassigned":
            T.setdefault(f, set()).add(s)
    T = {frozenset(v) for v in T.values() if len(v) > 1}
    pa = {frozenset(p) for v in A for p in itertools.combinations(sorted(v), 2)}
    pt = {frozenset(p) for v in T for p in itertools.combinations(sorted(v), 2)}
    return {"families_called": len(A), "families_reference": len(T), "families_identical": len(A & T),
            "pair_precision": round(len(pa & pt) / len(pa), 4) if pa else float("nan"), "pair_recall": round(len(pa & pt) / len(pt), 4) if pt else float("nan")}


def family_map(M, ref):
    """called family label -> reference family label (majority of its members)."""
    d = pd.DataFrame({"f": M.family.astype(str), "t": [ref.get(s, "unassigned") or "unassigned" for s in M.individual.astype(str)]})
    d = d[d.f != "unassigned"]
    return d.groupby("f").t.agg(lambda x: x.value_counts().index[0]).to_dict()


def link_set(table, fmap, mask=None):
    """Set of reference-family pairs behind the rows of a unit-pair table (families on both sides)."""
    t = table if mask is None else table[mask]
    out = set()
    for a, b in zip(t.unit_a, t.unit_b):
        if a in fmap and b in fmap and fmap[a] != "unassigned" and fmap[b] != "unassigned" and fmap[a] != fmap[b]:
            out.add(frozenset((fmap[a], fmap[b])))
    return out


def link_scores(called, ref_links):
    return {"links_called": len(called), "links_true": len(called & ref_links), "links_false": len(called - ref_links),
            "links_missed": len(ref_links - called), "links_reference": len(ref_links)}
