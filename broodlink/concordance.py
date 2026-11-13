"""Step 9: cross-marker concordance of kinship.

Extracted from two scripts of the original study:

* `compute_sv_kinship_v1.py`: VanRaden GRM per marker class (biSNPs; the same biSNPs with every
  SV region removed; deletions; duplications; structural haplotypes). "SV regions" = the
  intervals of the large structural haplotypes with a flank of 100 kb on each side, and the
  span of every deletion and duplication call (functions `merged`, `in_mask`).
* `build_extdata_SV_validation_v2.py`: rank AUC of within-family against between-family pairs of
  assigned individuals; recovery = share of within-family pairs above the 99.9th percentile of
  the between-family pairs of the same class; family-level bootstrap (`boot_ci`: the families
  are resampled with replacement, a within-family pair weighs m_f and a between-family pair
  m_f * m_g, 1,000 resamples, seed 7, percentile interval); Pearson correlation of each SV GRM
  with the GRM of the masked biSNPs over all pairs.

The families were built from one of the inputs, so this is concordance between marker
classes, not pedigree truth.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .io import MarkerSet, detail
from .relatedness import vanraden_grm

POINT = ("SNP", "INDEL", "SNP_INDEL", "SNP_MASKED")


def _merged(iv):
    iv = sorted(iv); out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return np.array(out)


def sv_region_mask(sets, hap, intervals, flank):
    """{chromosome: merged [start, end] array} of structural-haplotype intervals (+ flank) and DEL / DUP spans."""
    mask = {}
    if intervals is not None and len(intervals):
        for c, s, e in zip(intervals.chrom, intervals.start, intervals.end):
            mask.setdefault(str(c), []).append((max(0, int(s) - flank), int(e) + flank))
    elif hap is not None:
        for c, s, e in zip(hap.chrom, hap.pos, hap.end):
            if int(e) > 0:
                mask.setdefault(str(c), []).append((max(0, int(s) - flank), int(e) + flank))
    n_sv = 0
    for m in sets:
        if m.klass in ("DEL", "DUP"):
            n_sv += m.n
            for c, s, e in zip(m.chrom, m.pos, m.end):
                mask.setdefault(str(c), []).append((int(s), int(max(e, s))))
    return {c: _merged(v) for c, v in mask.items()}, n_sv


def in_mask(mask, chrom, pos):
    out = np.zeros(len(pos), bool)
    chrom = np.asarray(chrom).astype(str); pos = np.asarray(pos)
    for c in np.unique(chrom):
        if c not in mask:
            continue
        sel = chrom == c; m = mask[c]; p = pos[sel]
        i = np.searchsorted(m[:, 0], p, side="right") - 1
        out[sel] = (i >= 0) & (p <= m[np.clip(i, 0, None), 1])
    return out


def marker_sets(sets, hap, intervals, cfg):
    """The marker classes of the concordance, in table order: each input, the SNPs outside SV
    regions right after the SNPs (when SNPs and SV calls or intervals are both present), then the
    structural haplotypes. Returns ([(grm key, MarkerSet)], note)."""
    out, note, used = [], "", {}

    def key(k):
        used[k] = used.get(k, 0) + 1
        return k if used[k] == 1 else f"{k}_{used[k]}"
    snps = [m for m in sets if m.klass == "SNP"]
    main = max(snps, key=lambda m: (m.flag == "bisnp", m.n)) if snps else None
    for m in sets:
        out.append((key(m.klass), m))
        if m is main:
            mask, n_sv = sv_region_mask(sets, hap, intervals, int(cfg["concordance"]["sv_mask_flank_bp"]))
            if mask:
                keep = ~in_mask(mask, m.chrom, m.pos)
                mm = m.subset_markers(keep)
                mm.name = "SNP, SV regions masked" + (m.name[m.name.index(" ("):] if " (" in m.name else ""); mm.klass = "SNP_MASKED"
                out.append(("SNP_masked", mm))
                note = (f"{int((~keep).sum()):,} of {m.n:,} SNPs ({100 * float((~keep).mean()):.1f}%) lie in SV regions (structural-haplotype intervals with a flank of "
                        f"{cfg['concordance']['sv_mask_flank_bp']:,} bp, and the spans of {n_sv:,} deletion / duplication calls) and are left out of the masked SNP class")
    if hap is not None:
        out.append((key("HAPLOTYPE"), hap))
    return out, note


def auc(a, b):
    from scipy.stats import rankdata
    r = rankdata(np.concatenate([a, b]))
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))


def boot_ci(vw, fw, vb, fa_b, fb_b, n_fam, nb=1000, seed=7):
    """Family-level bootstrap of AUC and recovery (boot_ci of the original script).
    vw / fw: value and family index of the within-family pairs; vb / fa_b / fb_b: value and the
    two family indices of the between-family pairs. Returns ((auc_lo, auc_hi), (rec_lo, rec_hi))."""
    rg = np.random.default_rng(seed)
    o_ = np.argsort(vb); vb, fa_b, fb_b = vb[o_], fa_b[o_], fb_b[o_]; A, R_ = [], []
    lo = np.searchsorted(vb, vw, side="left"); hi = np.searchsorted(vb, vw, side="right")
    for _ in range(nb):
        m = np.bincount(rg.integers(0, n_fam, n_fam), minlength=n_fam).astype(float); ww = m[fw]; wb = m[fa_b] * m[fb_b]
        if ww.sum() == 0 or wb.sum() == 0:
            continue
        cw = np.cumsum(wb); tot = cw[-1]
        below = np.where(lo > 0, cw[np.maximum(lo - 1, 0)], 0.0); tied = np.where(hi > 0, cw[np.maximum(hi - 1, 0)], 0.0) - below
        A.append(float((ww * (below + .5 * tied)).sum() / (ww.sum() * tot)))
        k_ = np.searchsorted(cw, .999 * tot); thr = vb[min(k_, len(vb) - 1)]; R_.append(100 * float((ww * (vw > thr)).sum() / ww.sum()))
    if not A:
        return (np.nan, np.nan), (np.nan, np.nan)
    return tuple(np.percentile(A, [2.5, 97.5])), tuple(np.percentile(R_, [2.5, 97.5]))


def family_order(fam_of, analysis):
    """Individuals in family order: largest family first (ties by name), unassigned individuals last."""
    from .figures.style import fam_key
    fam = {s: fam_of.get(s, "unassigned") for s in analysis}
    size = pd.Series([f for f in fam.values() if f != "unassigned"]).value_counts()
    fl = sorted(size.index, key=lambda f: (-size[f], fam_key(f)))
    order = [s for f in fl for s in analysis if fam[s] == f] + [s for s in analysis if fam[s] == "unassigned"]
    return order, fl, fam


def run(sets, fam_of, analysis, cfg, primary_name=None):
    """sets: [(grm key, MarkerSet)]. Returns (summary table, pair table, {key: GRM in family order}, order table)."""
    g = cfg["relatedness"]["grm"]; cc = cfg["concordance"]
    order, fl, fam = family_order(fam_of, analysis)
    n = len(order); oix = {s: i for i, s in enumerate(order)}
    f = np.array([fam[s] for s in order]); asg = f != "unassigned"
    fidx = {x: i for i, x in enumerate(fl)}; fi = np.array([fidx.get(x, -1) for x in f])
    iu, ju = np.triu_indices(n, 1)
    both = asg[iu] & asg[ju]; same = (f[iu] == f[ju]) & both
    rows, pairs, grm = [], [], {}
    for key, ms in sets:
        frac = g["min_called_fraction_snp"] if ms.klass in POINT else g["min_called_fraction_sv"]
        K0, used = vanraden_grm(ms.G, g["af_min"], g["af_max"], frac)
        row = {"marker_set": ms.name, "marker_class": ms.klass, "grm_key": key, "n_markers": used}
        if used == 0:
            rows.append({**row, "note": "no marker passes the GRM filters"}); continue
        K = np.full((n, n), np.nan); ix = np.array([oix[s] for s in ms.samples if s in oix]); src = np.array([i for i, s in enumerate(ms.samples) if s in oix])
        K[np.ix_(ix, ix)] = K0[np.ix_(src, src)]
        grm[key] = K
        v = K[iu, ju]; fin = np.isfinite(v)
        w, b = v[same & fin], v[both & ~same & fin]
        row.update({"n_pairs_within_family": int(len(w)), "n_pairs_between_families": int(len(b)),
                    "is_source_of_families": bool(primary_name is not None and ms.name == primary_name), "note": ""})
        if len(w) and len(b):
            thr = np.percentile(b, 99.9)
            sw = same & fin; sb = both & ~same & fin
            detail(f"{ms.name}: {used:,} markers; bootstrap over {len(fl)} families ({cc['n_bootstrap']:,} resamples)")
            ci_a, ci_r = boot_ci(w, fi[iu][sw], b, fi[iu][sb], fi[ju][sb], len(fl), cc["n_bootstrap"], cc["bootstrap_seed"])
            row.update({"auc_within_vs_between": float(auc(w, b)), "recovery_pct": 100 * float((w > thr).mean()),
                        "median_kinship_within": float(np.median(w)) / 2.0, "median_kinship_between": float(np.median(b)) / 2.0,
                        "auc_ci95_low": float(ci_a[0]), "auc_ci95_high": float(ci_a[1]), "recovery_ci95_low": float(ci_r[0]), "recovery_ci95_high": float(ci_r[1])})
            if row["is_source_of_families"]:
                row["note"] = "the families were built from this input: this row is circular"
            pairs.append(pd.DataFrame({"marker_set": ms.name, "marker_class": ms.klass, "same_family": same[both & fin], "kinship": np.round(v[both & fin] / 2.0, 5)}))
        rows.append(row)
    cols = ["marker_set", "marker_class", "grm_key", "n_markers", "n_pairs_within_family", "n_pairs_between_families", "is_source_of_families", "note", "auc_within_vs_between",
            "auc_ci95_low", "auc_ci95_high", "recovery_pct", "recovery_ci95_low", "recovery_ci95_high", "median_kinship_within", "median_kinship_between", "pearson_r_vs_snp_outside_sv"]
    T = pd.DataFrame(rows)
    for c in cols:
        if c not in T:
            T[c] = np.nan if c not in ("note",) else ""
    ref = grm.get("SNP_masked")
    if ref is not None and len(T):
        xr = ref[iu, ju]
        for i, r in T.iterrows():
            if r.grm_key in grm and r.marker_class not in POINT:
                y = grm[r.grm_key][iu, ju]; ok = np.isfinite(xr) & np.isfinite(y)
                if ok.sum() > 2 and np.ptp(y[ok]) > 0:
                    T.loc[i, "pearson_r_vs_snp_outside_sv"] = float(np.corrcoef(xr[ok], y[ok])[0, 1])
    T = T[cols] if len(T) else pd.DataFrame(columns=cols)
    P = pd.concat(pairs, ignore_index=True) if pairs else pd.DataFrame(columns=["marker_set", "marker_class", "same_family", "kinship"])
    O = pd.DataFrame({"order": np.arange(1, n + 1), "individual": order, "family": f})
    return T, P, grm, O
