"""Step 2: full-sib families (relatedness-based groups, not pedigree records).

Extracted from get_FS.py: full-sib edge = theta >= 0.177 and k2 >= 0.15; a family is a block
(connected component) of full-sib edges; closure = share of within-family pairs that pass both
gates; core = maximal complete clique; attached = at least two full-sib edges into the core and
no first-degree edge outside the block. When no k2 exists (GRM mode) the edge rests on theta alone.
"""
from __future__ import annotations
import itertools

import networkx as nx
import numpy as np
import pandas as pd

HEADER = "Families are relatedness-based groups, not pedigree records."


def reconstruct(R, ids, cfg, given=None):
    """given: {family name: members} imported with --families. The families are then taken as
    given and only described (closure, status, membership class) with the same rules."""
    f = cfg["families"]
    has_k2 = bool(R.k2.notna().any())
    fd = R[R.theta >= cfg["relatedness"]["theta_first"]]
    fs = R[(R.theta >= f["theta_min"]) & (R.k2 >= f["k2_min"])] if has_k2 else R[R.theta >= f["theta_min"]]
    FS = {(a, b) for a, b in zip(fs.id1, fs.id2)}
    FD = {(a, b) for a, b in zip(fd.id1, fd.id2)}
    key = lambda a, b: (a, b) if a < b else (b, a)
    G = nx.Graph(); G.add_nodes_from(ids); G.add_edges_from(FS)
    if f["rule"] == "clique":                       # optional: greedy maximal cliques instead of blocks
        blocks, H = [], G.copy()
        while H.number_of_edges():
            c = max(nx.find_cliques(H), key=lambda c: (len(c), sorted(c)[::-1]))
            blocks.append(sorted(c)); H.remove_nodes_from(c)
    else:
        blocks = [sorted(c) for c in nx.connected_components(G) if len(c) > 1]
    blocks.sort(key=lambda b: (-len(b), b[0]))
    names = [f"F{bi}" for bi in range(1, len(blocks) + 1)]
    if given is not None:
        names = list(given); blocks = [sorted(given[n]) for n in names]
    rows, fam_rows = [], []
    theta = {key(a, b): t for a, b, t in zip(R.id1, R.id2, R.theta)}
    for name, blk in zip(names, blocks):
        sub = G.subgraph(blk)
        core = set(max(nx.find_cliques(sub), key=lambda c: (len(c), sum(sub.degree(x) for x in c), sorted(c)[::-1])))
        inblk = set(blk)
        npairs = len(blk) * (len(blk) - 1) // 2
        nfs = sum(1 for a, b in itertools.combinations(blk, 2) if key(a, b) in FS)
        th = [theta.get(key(a, b), np.nan) for a, b in itertools.combinations(blk, 2)]
        ext_fd = 0
        for x in blk:
            nin = sum(1 for c in core if c != x and key(x, c) in FS)
            ext = any((key(x, y) in FD) for y in ids if y not in inblk)
            ext_fd += int(ext)
            m = "core" if x in core else ("attached" if nin >= 2 and not ext else "unresolved member")
            rows.append({"individual": x, "unit": name, "family": name, "membership": m, "family_size": len(blk)})
        fam_rows.append({"family": name, "n": len(blk), "n_pairs": npairs, "n_pairs_full_sib": nfs, "closure": nfs / npairs,
                         "status": "closed" if nfs == npairs else "not closed", "n_core": len(core),
                         "n_with_first_degree_link_outside": ext_fd, "min_theta": float(np.nanmin(th)) if np.isfinite(th).any() else np.nan,
                         "median_theta": float(np.nanmedian(th)) if np.isfinite(th).any() else np.nan})
    infam = {r["individual"] for r in rows}
    for x in sorted(set(ids) - infam):
        rows.append({"individual": x, "unit": f"S_{x}", "family": "unassigned", "membership": "unassigned", "family_size": 1})
    M = pd.DataFrame(rows).sort_values("individual").reset_index(drop=True)
    F = pd.DataFrame(fam_rows)
    info = {"has_k2": has_k2, "n_full_sib_edges": len(FS), "n_first_degree_pairs": len(FD),
            "edge_rule": (f"theta >= {f['theta_min']} and k2 >= {f['k2_min']}" if has_k2 else
                          f"kinship >= {f['theta_min']} alone (no k2 can be estimated from a GRM)"),
            "group_rule": "connected blocks of full-sib edges" if f["rule"] != "clique" else "complete cliques of full-sib edges"}
    return M, F, info


def read_families(path, ids):
    """Import family labels from a previous run or another tool (columns: individual + family/unit)."""
    d = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    low = {c.lower(): c for c in d.columns}
    ic = next(low[c] for c in ("individual", "sample_id", "sample", "fish_id", "id") if c in low)
    fc = next(low[c] for c in ("family", "family_label", "unit") if c in low)
    fam = {}
    for s, v in zip(d[ic], d[fc]):
        v = "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)
        fam[s] = None if v in ("", "unassigned", "nan", "NA", "-1") or v.startswith("S_") else v
    rows = []
    sizes = pd.Series([v for s, v in fam.items() if v and s in set(ids)]).value_counts().to_dict()
    for x in ids:
        v = fam.get(x)
        if v and sizes.get(v, 0) > 1:
            rows.append({"individual": x, "unit": v, "family": v, "membership": "imported", "family_size": sizes[v]})
        else:
            rows.append({"individual": x, "unit": f"S_{x}", "family": "unassigned", "membership": "unassigned", "family_size": 1})
    return pd.DataFrame(rows).sort_values("individual").reset_index(drop=True)


def compare_partitions(M, truth):
    """Families recovered exactly, and pair-level precision / recall against a reference assignment."""
    a = {s: f for s, f in zip(M.individual, M.family) if f != "unassigned"}
    sets_a = {}
    for s, f in a.items():
        sets_a.setdefault(f, set()).add(s)
    sets_t = {}
    for s, f in truth.items():
        if f and f != "unassigned":
            sets_t.setdefault(f, set()).add(s)
    A = {frozenset(v) for v in sets_a.values()}; T = {frozenset(v) for v in sets_t.values() if len(v) > 1}
    pa = {frozenset(p) for v in A for p in itertools.combinations(sorted(v), 2)}
    pt = {frozenset(p) for v in T for p in itertools.combinations(sorted(v), 2)}
    return {"families_called": len(A), "families_reference": len(T), "families_identical": len(A & T),
            "individuals_called": sum(len(v) for v in A), "individuals_reference": sum(len(v) for v in T),
            "pair_precision": len(pa & pt) / max(len(pa), 1), "pair_recall": len(pa & pt) / max(len(pt), 1)}
