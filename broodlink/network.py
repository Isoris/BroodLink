"""Step 5: network of families and links, components, latent parents and the slot limit.

Extracted from make_tables.py (components, maximal cliques) and build_fig2_halfsib_v8.py (viol).
Latent parents are bookkeeping labels for maximal sets of mutually linked families; they are
not reconstructed parents. A family has two parental slots.
"""
from __future__ import annotations
import networkx as nx
import numpy as np
import pandas as pd


def _fkey(f):
    d = "".join(ch for ch in str(f) if ch.isdigit())
    return (int(d) if d else 0, str(f))


def slot_violations(edges, min_n):
    """Families whose link partners cannot be split into two groups of mutually linked families."""
    s = edges[edges.min_n >= min_n]
    adj = {}
    for r in s.itertuples():
        adj.setdefault(r.unit_a, set()).add(r.unit_b); adj.setdefault(r.unit_b, set()).add(r.unit_a)
    bad = []
    for f, partners in adj.items():
        P = sorted(partners, key=_fkey); col = {}; ok = True
        for s0 in P:
            if s0 in col:
                continue
            col[s0] = 0; st = [s0]
            while st:
                u = st.pop()
                for v in P:
                    if v != u and v not in adj[u]:
                        if v not in col:
                            col[v] = 1 - col[u]; st.append(v)
                        elif col[v] == col[u]:
                            ok = False
        if not ok:
            bad.append(f)
    return sorted(bad, key=_fkey)


def status_of(link_class):
    for k in ("rejected", "supported", "Gate-1 only", "provisional"):
        if str(link_class).startswith(k):
            return k
    return "nominated, not audited"


def build(E, cfg, sizes):
    """E: audited edge table. Returns dict of tables and summary numbers."""
    fam = E[(E.unassigned_sides == 0) & E.nominated].copy()
    fam["status"] = [status_of(x) for x in fam.link_class]
    kept = fam[fam.status != "rejected"]
    G = nx.Graph(); G.add_edges_from(zip(kept.unit_a, kept.unit_b))
    comps = sorted((sorted(c, key=_fkey) for c in nx.connected_components(G)),
                   key=lambda c: (-sum(sizes.get(f, 0) for f in c), _fkey(c[0])))
    comp_of = {f: i + 1 for i, c in enumerate(comps) for f in c}
    fam["component"] = [comp_of.get(a, 0) for a in fam.unit_a]
    # latent parents = maximal cliques of the kept links, once per minimum-family-size tier
    lp_rows, tier_summary = [], {}
    for t in cfg["network"]["slot_tiers"]:
        kt = kept[kept.min_n >= t]
        Gt = nx.Graph(); Gt.add_edges_from(zip(kt.unit_a, kt.unit_b))
        clt = sorted((tuple(sorted(c, key=_fkey)) for c in nx.find_cliques(Gt)),
                     key=lambda c: (comp_of[c[0]], [_fkey(x) for x in c])) if Gt.number_of_nodes() else []
        for i, c in enumerate(clt):
            lp_rows.append({"min_family_size": t, "latent_parent": f"H{i + 1}", "component": comp_of[c[0]],
                            "n_families": len(c), "families": ",".join(c)})
        tier_summary[t] = {"n_links": int(len(kt)), "n_families": Gt.number_of_nodes(), "n_latent_parents": len(clt),
                           "n_components": nx.number_connected_components(Gt) if Gt.number_of_nodes() else 0}
    LP = pd.DataFrame(lp_rows, columns=["min_family_size", "latent_parent", "component", "n_families", "families"])
    ft = cfg["network"]["figure_min_family_size"]
    cl = [tuple(x.split(",")) for x in LP[LP.min_family_size == ft].families] if ft in cfg["network"]["slot_tiers"] else []
    inc = {}
    for c in cl:
        for f in c:
            inc[f] = inc.get(f, 0) + 1
    sup = fam[fam.status == "supported"]
    Gs = nx.Graph(); Gs.add_edges_from(zip(sup.unit_a, sup.unit_b))
    scomps = sorted((sorted(c, key=_fkey) for c in nx.connected_components(Gs)), key=lambda c: _fkey(c[0]))
    nodes = pd.DataFrame([{"family": f, "n": sizes.get(f, 0), "component": comp_of[f], "n_links_kept": G.degree(f),
                           "n_links_supported": Gs.degree(f) if f in Gs else 0, "n_latent_parents": inc.get(f, 0)}
                          for c in comps for f in c],
                         columns=["family", "n", "component", "n_links_kept", "n_links_supported", "n_latent_parents"])
    viol = {t: slot_violations(kept, t) for t in cfg["network"]["slot_tiers"]}
    junction = {f: sorted(Gs.neighbors(f), key=_fkey) for f in sorted(Gs.nodes, key=_fkey) if Gs.degree(f) >= 2}
    summary = {"n_links_nominated_family_pairs": len(fam), "n_links_kept": len(kept), "n_links_supported": len(sup),
               "n_families_in_network": G.number_of_nodes(), "n_components": len(comps),
               "n_families_on_supported_links": Gs.number_of_nodes(), "n_components_supported": len(scomps),
               "n_latent_parents": len(cl), "latent_parent_tier": ft, "tiers": tier_summary, "slot_violations": viol, "junction_families": junction}
    return {"edges": fam, "nodes": nodes, "latent_parents": LP, "summary": summary, "graph": G,
            "supported_graph": Gs, "components": comps}
