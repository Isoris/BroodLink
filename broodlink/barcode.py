"""Step 6: barcode, the posterior genotype of a shared parent at each large structural haplotype.

Extracted from get_INV_barcode_SP_get_LP.py (per link) and lp01_profile.py (one parent shared by
several families). Pr(gP | D) is proportional to prior(gP) x product over families of L_F(gP),
with L_F(gP) = sum over the private parent gQ of that family of prior(gQ) x product over
offspring of T[gP, gQ, g]. Priors: Hardy-Weinberg from the arrangement frequency among the
individuals in families. Genotypes 0/1/2 are arrangement classes; the labels carry no
ancestral meaning.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .audit import transmission_tensor

CALLS = ["0/0", "0/1", "1/1"]


def _fam_joint(g, pri, TE):
    J = np.empty((3, 3))
    for gp in range(3):
        for gq in range(3):
            pr = pri[gq]
            for x in g:
                if x >= 0:
                    pr *= TE[gp, gq, x]
            J[gp, gq] = pr
    return J


def barcode(hap, fam_members, links, junctions, cfg, interval_info=None):
    """hap: MarkerSet (intervals x individuals). links: list of (family_a, family_b, status).
    junctions: dict family -> partner families on supported links."""
    b = cfg["barcode"]; TE = transmission_tensor(b["genotype_error"])
    si = {s: i for i, s in enumerate(hap.samples)}
    infam = [si[s] for f in fam_members if not f.startswith("S_") for s in fam_members[f] if s in si]
    Gf = hap.G[:, infam]
    nc = (Gf >= 0).sum(1)
    af = np.where(nc > 0, np.where(Gf >= 0, Gf, 0).sum(1) / (2.0 * np.maximum(nc, 1)), np.nan)
    info = interval_info or {}
    rows = []

    def add(kind, label, fams, status):
        for j in range(hap.n):
            if np.isnan(af[j]):
                continue
            pri = np.array([(1 - af[j]) ** 2, 2 * af[j] * (1 - af[j]), af[j] ** 2])
            post = pri.copy(); n = 0; Js = {}
            for f in fams:
                g = [int(hap.G[j, si[s]]) for s in fam_members[f] if s in si]
                n += sum(x >= 0 for x in g)
                Js[f] = _fam_joint(g, pri, TE)
                post = post * Js[f].sum(1)
            if not (post.sum() > 0):
                continue
            post = post / post.sum()
            gp = int(np.argmax(post)); pmax = float(post[gp])
            ii = info.get(hap.ids[j], {})
            r = {"parent": label, "kind": kind, "families": ",".join(fams), "link_status": status, "interval": hap.ids[j],
                 "p_00": post[0], "p_01": post[1], "p_11": post[2], "best_genotype": CALLS[gp], "posterior": pmax,
                 "call": CALLS[gp] if pmax >= b["posterior_call"] else "unresolved",
                 "called_strict": bool(pmax >= b["posterior_strict"]), "n_descendant_genotypes": n,
                 "arrangement_frequency": af[j], "interval_source": ii.get("source", hap.name),
                 "interval_tier": ii.get("tier", ""), "interval_label": ii.get("label", ""),
                 "genotype_markers": ii.get("markers", "")}
            if kind == "link":
                for f, tag in zip(fams, ("a", "b")):
                    v = Js[f][gp]
                    v = v / v.sum() if v.sum() > 0 else np.full(3, np.nan)
                    ok = not np.all(np.isnan(v))
                    r[f"private_parent_{tag}_call"] = CALLS[int(np.nanargmax(v))] if ok and np.nanmax(v) >= b["posterior_call"] else "unresolved"
                    r[f"private_parent_{tag}_posterior"] = float(np.nanmax(v)) if ok else np.nan
            rows.append(r)
    for a, c, status in links:
        add("link", f"{a}x{c}", [a, c], status)
    for f, partners in junctions.items():
        add("junction", "J_" + f, [f] + list(partners), "same parent assumed on every supported link of " + f + " (candidate, not tested)")
    return pd.DataFrame(rows)
