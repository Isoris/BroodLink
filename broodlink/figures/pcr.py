"""Figures of step 10 (PCR markers), ported from the figure scripts of the original study:

  fig_assay      the assay figure (schema of the flanking-pair assay, expected bands of four
                 markers as lane strips, and one row per primer pair: deletion frequency by
                 karyotype, designed products, in-silico search);
  fig_gels       the simulated gel sheets, paged (pcr.gels_per_page gels per page);
  fig_schematic  the generic schematic of the band-shift assay.

No gel is real: every band is computed from the designed product sizes, and the figures say so.
Colours follow the barcode figures: class 0 amber, heterokaryotype taupe, class 2 purple.
"""
from __future__ import annotations
import math

import numpy as np

from . import style
from .style import LED, NAVY

GREYT = "#6B7280"
SLATE = "#8A949C"
FWD, REV = style.ROH_BASE["teal"], style.ROH_BASE["brick"]
C0, CH, C2 = style.STATE0, style.STATE1, style.STATE2
FS, TITLE = 5.5, 6.2
SRC_P, SRC_I = "10_pcr_markers/pcr_primers.tsv", "10_pcr_markers/pcr_markers_per_interval.tsv"
YL = (50, 5200)


def _tint(c, f=0.85):
    return style.mixw(c, f)


BG = {0: _tint(C0), 1: _tint(CH), 2: _tint(C2)}
LANES = (("class 0", 0, C0), ("het", 1, CH), ("class 2", 2, C2))


def _primers(run):
    P = run.t.get("pcr", {}).get("primers")
    if P is None or not len(P):
        return None
    return P.reset_index(drop=True)


def _lanes(r):
    """Expected bands per karyotype class: {0: [...], 1: [...], 2: [...]}, short band, long band."""
    d, u = int(r.product_with_deletion_bp), int(r.product_without_deletion_bp)
    if r.arrangement_with_deletion == "class 2":
        return {0: [u], 1: [d, u], 2: [d]}, d, u
    return {0: [d], 1: [d, u], 2: [u]}, d, u


def _freqs(r):
    return {0: float(r.deletion_frequency_class_0), 1: float(r.deletion_frequency_heterokaryotypes), 2: float(r.deletion_frequency_class_2)}


def _f2(v):
    if v != v:
        return "–"
    return "1.0" if v >= 0.995 else f"{v:.2f}".replace("0.", ".")


def _ins(P):
    return bool(len(P)) and P.in_silico_products_at_locus.notna().any() and str(P.in_silico_products_at_locus.iloc[0]) != ""


def _num(v, default=0):
    try:
        return default if v != v or v == "" else float(v)
    except (TypeError, ValueError):
        return default


def _mark(r, ins):
    if not ins:
        return ""
    return ("†" if _num(r.in_silico_products_at_locus) == 0 else "") + ("‡" if _num(r.in_silico_products_off_target) > 0 else "")


def _chrom(c):
    c = str(c)                                                   # long sequence names: the part after the last underscore
    if len(c) > 7 and "_" in c and c.rsplit("_", 1)[1]:
        c = c.rsplit("_", 1)[1]
    return c if len(c) <= 7 else c[-7:]


def _led_marker(panel, r):
    LED(panel, f"{r.marker} {r.variant_id}: deletion bp; on arrangement; products with / without the deletion",
        f"{int(r.length_bp)}; {r.arrangement_with_deletion}; {int(r.product_with_deletion_bp)}/{int(r.product_without_deletion_bp)}",
        f"{SRC_P}: length_bp, arrangement_with_deletion, product_with_deletion_bp, product_without_deletion_bp")
    f = _freqs(r)
    LED(panel, f"{r.marker} deletion frequency class 0, heterokaryotypes, class 2; difference", f"{f[0]:.3g}, {f[1]:.3g}, {f[2]:.3g}; {float(r.frequency_difference):.3g}",
        f"{SRC_P}: deletion_frequency_class_0, deletion_frequency_heterokaryotypes, deletion_frequency_class_2, frequency_difference")


def _examples(P):
    """Markers opened as examples (the original opened four by hand): the first whose deletion is
    on class 2, the first on class 0, the longest deletion, the shortest; duplicates removed."""
    idx = []
    for sel in (P.index[P.arrangement_with_deletion == "class 2"][:1], P.index[P.arrangement_with_deletion == "class 0"][:1],
                [P.length_bp.astype(int).idxmax()], [P.length_bp.astype(int).idxmin()]):
        for i in sel:
            if i not in idx:
                idx.append(i)
    for i in P.index:
        if len(idx) >= 4:
            break
        if i not in idx:
            idx.append(i)
    return [P.loc[i] for i in idx[:4]]


def _noun(run):
    return run.cfg["figures"].get("individual_noun", "individuals")


# =============================================================================== assay figure
def assay_pages(n_pairs, rp=2.75):
    """Rows of the per-pair panels on each page: the first page also holds panels A-E, the
    following pages only the rows; every page stays within A4."""
    first = int((style.A4_H - 18.0 - 53.6) // rp); later = int((style.A4_H - 18.0 - 12.9) // rp)
    out = [min(n_pairs, first)]
    while sum(out) < n_pairs:
        out.append(min(n_pairs - sum(out), later))
    return out


def fig_assay(run, plt):
    """Port of the assay figure of the original study (its Extended Data figure of the deletion
    markers): A schema for one marker, B-E expected lane strips, F-H one row per primer pair.
    One page when the pairs fit; otherwise the rows continue on following pages."""
    P = _primers(run)
    if P is None:
        return "no primer pair was designed"
    pages = assay_pages(len(P)); out = []; k = 0
    for pg, n in enumerate(pages):
        fig, title, desc, note = _assay_page(run, plt, P, P.iloc[k:k + n], pg, len(pages)); k += n
        if len(pages) == 1:
            return fig, title, desc, note
        out.append((f"page{pg + 1}of{len(pages)}", fig, title + f", page {pg + 1} of {len(pages)}", desc, (f"page {pg + 1} of {len(pages)}" + ("; " + note if note else ""))))
    return out


def _assay_page(run, plt, P, D, pg, npage):
    from matplotlib.patches import Rectangle
    ins = _ins(P); noun = _noun(run); first = pg == 0
    W = 183.0; rp = 2.75; y2 = 41.0 if first else 0.3; ttl = y2 + 6.0; t0 = ttl + 6.6
    n_all = len(P); n = len(D)
    bot = t0 + n * rp; H = bot + 18.0
    cv = style.canvas(plt, W, H, style.SZ_ED, thin=True); C = cv.overlay

    def tx(x, y, s, **kw):
        kw.setdefault("fontsize", FS); kw.setdefault("va", "center"); kw.setdefault("color", NAVY)
        return C.text(x, y, s, **kw)

    def arr(x0, x1, y, c, lw=0.6, st="-|>", ms=4):
        C.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle=st, color=c, lw=lw, shrinkA=0, shrinkB=0, mutation_scale=ms))

    ex = _examples(P); r8 = ex[0]
    if first:
        _assay_top(run, cv, C, tx, arr, ex, r8, noun)
    return _assay_rows(run, cv, C, tx, P, D, ins, pg, npage, y2, ttl, t0, rp, n, n_all, bot)


def _assay_top(run, cv, C, tx, arr, ex, r8, noun):
    from matplotlib.patches import Rectangle
    cv.banner(0.3, "1", f"Expected, not run: one primer pair and one reaction per {noun.rstrip('s') if noun.endswith('s') else noun} should read the arrangement as long band, short band, or both", fs=7.5)
    # ---- A schema
    on2 = r8.arrangement_with_deletion == "class 2"
    c_abs, c_del = (C0, C2) if on2 else (C2, C0); n_abs, n_del = ("class 0", "class 2") if on2 else ("class 2", "class 0")
    C.text(0.5, 7.2, "A", fontsize=8, fontweight="bold", va="center", color=NAVY)
    C.text(4.7, 7.2, f"Flanking-pair assay, example {r8.marker} on {_chrom(r8.chrom)}", fontsize=TITLE, fontweight="bold", va="center", color=NAVY)
    WB = dict(facecolor="white", edgecolor="none", pad=0.6)
    tx(4.7, 10.3, f"{n_abs} arrangement: deletion absent", fontweight="bold")
    arr(12, 22, 12.6, FWD); arr(76, 66, 12.6, REV); tx(10.6, 12.6, "F", ha="right"); tx(77.4, 12.6, "R")
    C.plot([5, 83], [14.7, 14.7], color=c_abs, lw=0.8); C.add_patch(Rectangle((33, 13.5), 22, 2.4, fc=_tint(c_abs, 0.6), ec=c_abs, lw=0.3, zorder=3))
    tx(44, 14.7, f"{int(r8.length_bp):,} bp present", ha="center", zorder=4)
    arr(12, 76, 17.6, NAVY, lw=0.35, st="<->", ms=3); tx(44, 17.6, f"product without the deletion {int(r8.product_without_deletion_bp):,} bp", ha="center", fontweight="bold", bbox=WB, zorder=4)
    tx(4.7, 20.7, f"{n_del} arrangement: deletion present", fontweight="bold")
    arr(26, 36, 23.0, FWD); arr(62, 52, 23.0, REV); tx(24.6, 23.0, "F", ha="right"); tx(63.4, 23.0, "R"); tx(44, 23.0, "junction", ha="center")
    C.plot([19, 69], [25.1, 25.1], color=c_del, lw=0.8); C.plot([44, 44], [24.1, 26.1], color=c_del, lw=0.6)
    arr(26, 62, 28.0, NAVY, lw=0.35, st="<->", ms=3); tx(44, 28.0, f"product with the deletion {int(r8.product_with_deletion_bp):,} bp", ha="center", fontweight="bold", bbox=WB, zorder=4)
    C.plot([4.7, 7.7], [31.3, 31.3], color=FWD, lw=0.7); tx(8.5, 31.3, f"F primer, + strand, 5′-{r8.forward_primer}-3′")
    C.plot([4.7, 7.7], [33.5, 33.5], color=REV, lw=0.7); tx(8.5, 33.5, f"R primer, − strand, 5′-{r8.reverse_primer}-3′")
    tx(4.7, 35.7, f"{str(r8.variant_id)[:40]}; not to scale; both primers lie in sequence kept on both alleles.")
    LED("A", f"{r8.marker} forward, reverse primer", f"{r8.forward_primer}, {r8.reverse_primer}", f"{SRC_P}: forward_primer, reverse_primer")
    _led_marker("A", r8)
    # ---- B-E expected lane strips (a size axis with one tick per expected band; not a gel image)
    FX0, FX1, FY0, FY1 = 90.2, 182.3, 6.0, 39.2
    C.add_patch(Rectangle((FX0, FY0), FX1 - FX0, FY1 - FY0, fill=False, ls=(0, (3, 2)), lw=0.45, ec=NAVY))
    tx(FX0 + 1.3, 8.1, "Expected, not run: bands computed from the designed products", fontweight="bold")

    def strip(lab, r, x0, y0):
        L, d, u = _lanes(r); tw = 25.0; tx0 = x0 + 9.8; f = _freqs(r)
        xm = lambda b: tx0 + tw * np.log10(max(b, 50) / 50.0) / 2.0     # noqa: E731
        C.text(x0, y0, lab, fontsize=8, fontweight="bold", va="center", color=NAVY); tx(x0 + 3.5, y0, f"{r.marker} · {_chrom(r.chrom)}", fontweight="bold")
        tx(x0 + 18.5, y0, f"{int(r.length_bp):,} bp on {r.arrangement_with_deletion}")
        pitch = 2.3; yl = y0 + 1.7
        for i, (nm, k, _c) in enumerate(LANES):
            yc = yl + (i + 0.5) * pitch; C.add_patch(Rectangle((tx0, yc - pitch / 2 + 0.12), tw, pitch - 0.24, fc=BG[k], ec="none", zorder=1)); tx(x0, yc, nm)
            for b in L[k]:
                xb = min(xm(b), tx0 + tw)
                C.plot([xb, xb], [yc - 0.85, yc + 0.85], color=NAVY, lw=0.7, zorder=3); tx(xb + 0.45, yc, f"{b:,}", zorder=4)
            tx(tx0 + tw + 5.0, yc, _f2(f[k]), ha="right")
        ya = yl + 3 * pitch + 0.15; C.plot([tx0, tx0 + tw], [ya, ya], color=NAVY, lw=0.35)
        for v in (100, 1000):
            C.plot([xm(v)] * 2, [ya, ya + 0.8], color=NAVY, lw=0.35); tx(xm(v), ya + 1.9, f"{v:,}", ha="center")
        tx(tx0 + tw, ya + 1.9, "bp", ha="right")
        _led_marker(lab, r)

    for (lab, x0, y0), r in zip((("B", 91.5, 11.0), ("C", 137.5, 11.0), ("D", 91.5, 24.3), ("E", 137.5, 24.3)), ex):
        strip(lab, r, x0, y0)
    tx(FX0 + 1.3, 37.6, "Tick: expected band, size beside it. Right column: deletion frequency of the karyotype class in the cohort.")


def _assay_rows(run, cv, C, tx, P, D, ins, pg, npage, y2, ttl, t0, rp, n, n_all, bot):
    from matplotlib.patches import Rectangle
    # ---- row 2
    dmin, dmax = float(P.frequency_difference.min()), float(P.frequency_difference.max())
    S = run.s.get("pcr", {})
    if ins:
        n_prod = int((P.in_silico_products_at_locus.astype(float) > 0).sum()); n_off = int((P.in_silico_products_off_target.astype(float) > 0).sum())
        B2 = f"In silico, {LED('banner 2', 'pairs with a product at the locus', n_prod, SRC_P + ': in_silico_products_at_locus')} of {n_all} pairs give a product at the locus in the reference, {n_all - n_prod} none; {LED('banner 2', 'pairs with off-target products', n_off, SRC_P + ': in_silico_products_off_target')} with off-target products"
    else:
        B2 = (f"{LED('banner 2', 'primer pairs', n_all, SRC_P + ': rows')} primer pairs for {LED('banner 2', 'intervals with a marker', S.get('n_intervals_with_marker', n_all), SRC_I + ': n_primer_pairs > 0')} of "
              f"{LED('banner 2', 'intervals', S.get('n_intervals', n_all), SRC_I + ': rows')} {run.cfg['figures'].get('haplotype_label', 'structural haplotypes')}; no reference given, so no in-silico search")
    cv.banner(y2, "2", B2 + (f" (rows {D.marker.iloc[0]} to {D.marker.iloc[-1]})" if npage > 1 else ""), fs=7.5)
    C.text(0.5, ttl, "F", fontsize=8, fontweight="bold", va="center", color=NAVY)
    C.text(4.7, ttl, f"Deletion frequency by karyotype, difference {dmin:.2f}–{dmax:.2f}", fontsize=TITLE, fontweight="bold", va="center", color=NAVY)
    LED("F title", "frequency difference range", f"{dmin:.2f}-{dmax:.2f}", f"{SRC_P}: frequency_difference")
    C.text(70.0, ttl, "G", fontsize=8, fontweight="bold", va="center", color=NAVY)
    C.text(74.2, ttl, "Designed products" + (" and predicted product" if ins else ""), fontsize=TITLE, fontweight="bold", va="center", color=NAVY)
    C.text(127.5, ttl, "H", fontsize=8, fontweight="bold", va="center", color=NAVY)
    C.text(131.7, ttl, "Primers" + (" and specificity against the reference" if ins else ""), fontsize=TITLE, fontweight="bold", va="center", color=NAVY)
    for i, r in enumerate(D.itertuples()):
        y = t0 + i * rp + rp / 2
        if i % 2 == 0:
            C.add_patch(Rectangle((1.0, y - rp / 2), 181.5, rp, fc="#F6F7F9", ec="none", zorder=0))
        tx(1.5, y, f"{r.marker}{_mark(r, ins)}", fontweight="bold"); tx(8.2, y, _chrom(r.chrom)); tx(19.5, y, f"{int(r.length_bp):,}", ha="right")
        C.add_patch(Rectangle((21.0, y - 0.7), 2.0, 1.4, fc=C0 if r.arrangement_with_deletion == "class 0" else C2, ec="none"))
        _led_marker("F-H", r)
    for x_, s_, ha in ((1.5, "pair", "left"), (8.2, "seq.", "left"), (19.5, "bp", "right"), (20.6, "on", "left")):
        tx(x_, t0 - 1.6, s_, fontweight="bold", ha=ha)
    fx = cv.ax(25.0, t0, 42.0, n * rp); fx.patch.set_alpha(0); fx.set_xlim(-0.02, 1.02); fx.set_ylim(n, 0); fx.set_yticks([])
    fx.set_xticks([0, 0.25, 0.5, 0.75, 1]); fx.set_xticklabels(["0", ".25", ".5", ".75", "1"], fontsize=FS); fx.set_xlabel("deletion frequency (cohort)", fontsize=FS, labelpad=1)
    fx.spines["left"].set_visible(False); fx.axvline(0.5, color=GREYT, lw=0.25, ls=(0, (3, 2)))
    for i, r in enumerate(D.itertuples()):
        f = _freqs(r)
        fx.plot([min(f[0], f[2]), max(f[0], f[2])], [i + 0.5] * 2, color="#C9D1D8", lw=0.25, zorder=1)
        for v, c in ((f[0], C0), (f[1], CH), (f[2], C2)):
            if v == v:
                fx.scatter([v], [i + 0.5], s=4.5, c=c, linewidths=0, zorder=3)
    import matplotlib
    xmax = max(2400, int(D.product_without_deletion_bp.astype(int).max() * 1.15))
    gx = cv.ax(70.0, t0, 52.0, n * rp); gx.patch.set_alpha(0); gx.set_xscale("log"); gx.set_xlim(50, xmax); gx.set_ylim(n, 0); gx.set_yticks([])
    gx.set_xticks([100, 200, 500, 1000, 2000]); gx.set_xticklabels(["100", "200", "500", "1,000", "2,000"], fontsize=FS); gx.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    gx.set_xlabel("product size (bp)", fontsize=FS, labelpad=1); gx.spines["left"].set_visible(False)
    for i, r in enumerate(D.itertuples()):
        gx.plot([r.product_with_deletion_bp, r.product_without_deletion_bp], [i + 0.5] * 2, color="#C9D1D8", lw=0.25, zorder=1)
        gx.scatter([r.product_with_deletion_bp], [i + 0.5], s=4.5, c=NAVY, linewidths=0, zorder=3); gx.scatter([r.product_without_deletion_bp], [i + 0.5], s=4.5, c=SLATE, linewidths=0, zorder=3)
        if ins and _num(r.in_silico_products_at_locus) > 0:
            gx.scatter([_num(r.in_silico_product_at_locus_bp)], [i + 0.5], s=11, facecolors="none", edgecolors=NAVY, linewidths=0.35, zorder=4)
    hc = dict(tf=127.5, tr=135.0, bd=146.0, bu=155.5, ins=157.2, dif=176.0, off=181.8)
    heads = [("tf", "Tm F", "left"), ("tr", "Tm R", "left"), ("bd", "with", "right"), ("bu", "without", "right")] + ([("ins", "at locus", "left"), ("dif", "− design", "right"), ("off", "off", "right")] if ins else [])
    for k_, s_, ha in heads:
        tx(hc[k_], t0 - 1.6, s_, ha=ha, fontweight="bold")
    for i, r in enumerate(D.itertuples()):
        y = t0 + i * rp + rp / 2
        tx(hc["tf"], y, f"{float(r.tm_forward):.1f}"); tx(hc["tr"], y, f"{float(r.tm_reverse):.1f}")
        tx(hc["bd"], y, f"{int(r.product_with_deletion_bp):,}", ha="right"); tx(hc["bu"], y, f"{int(r.product_without_deletion_bp):,}", ha="right")
        LED("H", f"{r.marker} Tm forward, reverse (degrees C)", f"{float(r.tm_forward):.1f}, {float(r.tm_reverse):.1f}", f"{SRC_P}: tm_forward, tm_reverse")
        if ins:
            if _num(r.in_silico_products_at_locus) == 0:
                tx(hc["ins"], y, "none"); tx(hc["dif"], y, "–", ha="right", color=GREYT)
            else:
                dd = int(_num(r.in_silico_product_minus_design_bp)); tx(hc["ins"], y, f"{int(_num(r.in_silico_product_at_locus_bp)):,}")
                tx(hc["dif"], y, f"{dd:+d}".replace("-", "−") if dd else "0", ha="right")
            no = int(_num(r.in_silico_products_off_target)); tx(hc["off"], y, f"{no}", ha="right", color=NAVY if no else GREYT)
            LED("H,G", f"{r.marker} in silico: products at the locus; product bp; off-target products", f"{int(_num(r.in_silico_products_at_locus))}; {r.in_silico_product_at_locus_bp}; {no}",
                f"{SRC_P}: in_silico_products_at_locus, in_silico_product_at_locus_bp, in_silico_products_off_target")
    yk = bot + 7.8
    xk = 1.5; tx(xk, yk, "F:", fontweight="bold"); xk += 3.0
    for nm, c in (("class 0", C0), ("heterokaryotypes", CH), ("class 2", C2)):
        C.scatter([xk + 0.8], [yk], s=5, c=c, linewidths=0); tx(xk + 2.0, yk, nm); xk += 4.5 + cv.text_width(nm, FS)
    tx(xk, yk, "dashed, 0.5")
    C.add_patch(Rectangle((1.5, yk + 1.5), 2.0, 1.4, fc=C0, ec="none")); C.add_patch(Rectangle((4.0, yk + 1.5), 2.0, 1.4, fc=C2, ec="none")); tx(6.7, yk + 2.2, "column on: deletion carried by class 0, class 2")
    if ins:
        tx(1.5, yk + 4.4, "† no product predicted at the locus in the reference."); tx(1.5, yk + 6.6, "‡ off-target products predicted.")
    x2 = 70.0; tx(x2, yk, "G:", fontweight="bold")
    C.scatter([x2 + 3.8], [yk], s=5, c=NAVY, linewidths=0); tx(x2 + 5, yk, "designed, with deletion"); C.scatter([x2 + 28.8], [yk], s=5, c=SLATE, linewidths=0); tx(x2 + 30, yk, "without")
    if ins:
        C.scatter([x2 + 3.8], [yk + 2.2], s=11, facecolors="none", edgecolors=NAVY, linewidths=0.35); tx(x2 + 5, yk + 2.2, "predicted in the reference, at the locus")
    tx(x2, yk + 4.4, "All products are predictions from sequence; none run.")
    tx(127.5, yk, "H: Tm of the forward and reverse primer (°C); designed product")
    tx(127.5, yk + 2.2, "with and without the deletion (bp)" + ("; product at the locus in the" if ins else "."))
    if ins:
        tx(127.5, yk + 4.4, "reference, its difference from the designed long product,")
        tx(127.5, yk + 6.6, "off-target products.")
    note = "" if ins else "panels G and H without the in-silico product and specificity columns (no --reference)"
    desc = ("A the flanking-pair assay for one marker; B-E expected bands of up to four markers per karyotype class (first marker with the deletion on class 2, first on class 0, longest and shortest deletion), "
            "with the deletion frequency of each class; F-H every primer pair: deletion frequency by karyotype, designed products" + (", product predicted in the reference and off-target products" if ins else "") +
            ". Bands are computed from the designed product sizes: these are not gel images and no laboratory test was made.")
    return cv.fig, "PCR markers of the structural haplotypes: assay, expected bands, all primer pairs", desc, note


# =============================================================================== simulated gels
def _gel(cv, l, t, w, h, r, ladder, ylab=True, heads=True, dashed=False):
    """Simulated gel: ladder and three lanes, band sizes in the right gutter, deletion frequency under each lane."""
    from matplotlib.patches import Rectangle
    ax = cv.ax(l, t, w, h); L, d, u = _lanes(r); f = _freqs(r)
    ax.set_yscale("log"); ax.set_ylim(*YL); ax.set_xlim(-0.55, 3.5)
    for i, (_nm, k, _c) in enumerate(LANES):
        x = i + 1; ax.add_patch(Rectangle((x - 0.46, YL[0]), 0.92, YL[1] - YL[0], facecolor=BG[k], edgecolor="none", zorder=0))
        for b in L[k]:
            ax.hlines(min(max(b, YL[0]), YL[1]), x - 0.4, x + 0.4, color=NAVY, lw=0.9, zorder=3)
    for b in ladder:
        ax.hlines(b, -0.4, 0.4, color=GREYT, lw=0.55 if b in (500, 1000) else 0.3, zorder=3)
    tk = [100, 500, 1000, 2000, 4000]; ax.set_yticks(tk); ax.set_yticks([], minor=True); ax.set_yticklabels([f"{v:,}" if ylab else "" for v in tk], fontsize=FS); ax.set_xticks([])
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    if ylab:
        ax.set_ylabel("size (bp)", fontsize=FS, labelpad=1)
    if heads:
        for x, s_ in ((0, "ladder"), (1, "class 0"), (2, "het"), (3, "class 2")):
            ax.text(x, YL[1] * 1.06, s_, rotation=90, ha="center", va="bottom", fontsize=FS, color=NAVY)
    for b in sorted({d, u}):
        ax.text(3.56, b, f"{b:,}", ha="left", va="center", fontsize=FS, clip_on=False, color=NAVY)
    for x, k in ((1, 0), (2, 1), (3, 2)):
        ax.text(x, YL[0] * 0.80, _f2(f[k]), ha="center", va="center", fontsize=FS, clip_on=False, color=NAVY)
    if dashed:
        ax.add_patch(Rectangle((0.5, YL[0]), 3.0, YL[1] - YL[0], fill=False, ls=(0, (3, 2)), lw=0.45, ec=NAVY, zorder=5))
    return ax


def gel_layout(per_page):
    """(gels per page actually used, columns, rows, gel width, gel height, column pitch, row pitch) in mm.
    Up to 9 gels per page: three columns of large gels; more: six columns as in the published sheets.
    The number of rows is capped so that a page stays within A4."""
    per_page = max(1, int(per_page))
    if per_page <= 9:
        ncol, gw, gh, cp, rpitch = 3, 34.0, 40.0, 57.0, 62.0
    else:
        ncol, gw, gh, cp, rpitch = 6, 17.5, 28.0, 28.3, 45.0
    nrow = int(math.ceil(per_page / ncol))
    max_rows = int((style.A4_H - 19.5 - 14.0 + (rpitch - gh - 8.0)) // rpitch)
    nrow = max(1, min(nrow, max_rows))
    return min(per_page, ncol * nrow), ncol, nrow, gw, gh, cp, rpitch


def fig_gels(run, plt):
    P = _primers(run)
    if P is None:
        return "no primer pair was designed"
    from matplotlib.patches import Rectangle
    ins = _ins(P); ladder = [int(x) for x in run.cfg["pcr"]["ladder_bp"]]
    per, ncol, nrow, gw, gh, cp, rpitch = gel_layout(run.cfg["pcr"]["gels_per_page"])
    npage = int(math.ceil(len(P) / per)); out = []
    S = run.s.get("pcr", {}); without = S.get("intervals_without_marker") or {}
    for pg in range(npage):
        sub = P.iloc[pg * per:(pg + 1) * per]
        rows_here = int(math.ceil(len(sub) / ncol))
        W = 183.0; t0 = 19.5; H = t0 + rows_here * rpitch - (rpitch - gh - 8.0) + 14.0
        cv = style.canvas(plt, W, H, style.SZ_ED, thin=True); C = cv.overlay

        def tx(x, y, s, **kw):
            kw.setdefault("fontsize", FS); kw.setdefault("va", "center"); kw.setdefault("color", NAVY)
            return C.text(x, y, s, **kw)
        B = f"Simulated gels, page {pg + 1} of {npage}: {LED('banner', f'page {pg + 1} primer pairs', len(sub), SRC_P + ': rows')} primer pairs, {sub.marker.iloc[0]} to {sub.marker.iloc[-1]}"
        if ins:
            nop = sub[sub.in_silico_products_at_locus.astype(float) == 0]; off = sub[sub.in_silico_products_off_target.astype(float) > 0]
            B += f"; {LED('banner', f'page {pg + 1} pairs with a product at the locus', len(sub) - len(nop), SRC_P + ': in_silico_products_at_locus')} predict a product at the locus, {len(nop)} none"
            if len(off):
                B += "; " + ", ".join(f"{r.marker} {int(float(r.in_silico_products_off_target))} off-target" for r in off.itertuples())
        cv.banner(0.3, "1", B, fs=7.5)
        lx0 = 12.0 if ncol == 3 else 11.0
        for i, r in enumerate(sub.itertuples()):
            cx, cy = i % ncol, i // ncol; x0 = lx0 + cx * cp; yt = t0 + cy * rpitch
            noprod = ins and _num(r.in_silico_products_at_locus) == 0
            _gel(cv, x0, yt, gw, gh, r, ladder, ylab=(cx == 0), heads=(cy == 0), dashed=bool(noprod))
            _led_marker(f"gel p{pg + 1}", r)
            yl = yt - 4.9 if cy > 0 else yt - 13.2
            tx(x0, yl, f"{r.marker} · {_chrom(r.chrom)}{_mark(r, ins)}", fontweight="bold"); tx(x0, yl + 2.2, f"{int(r.length_bp):,} bp on {r.arrangement_with_deletion}")
            yb = yt + gh + 4.2
            if ins:
                on = "none" if noprod else f"{int(_num(r.in_silico_product_at_locus_bp)):,}"
                no = int(_num(r.in_silico_products_off_target))
                tx(x0, yb, f"in silico {on}" + (f"; {no} off" if no else ""))
                LED(f"gel p{pg + 1}", f"{r.marker} in silico: product at the locus bp; off-target products", f"{r.in_silico_product_at_locus_bp}; {no}", f"{SRC_P}: in_silico_product_at_locus_bp, in_silico_products_off_target")
        yk = t0 + rows_here * rpitch - (rpitch - gh - 8.0) + 1.5
        tx(1.5, yk, "Lanes:")
        xk = 10.0
        for nm, k, c in (("class 0", 0, C0), ("heterokaryotype", 1, CH), ("class 2", 2, C2)):
            C.add_patch(Rectangle((xk, yk - 0.8), 3.0, 1.6, fc=BG[k], ec=c, lw=0.35)); tx(xk + 3.8, yk, nm); xk += 7.0 + cv.text_width(nm, FS)
        C.plot([xk, xk + 3], [yk, yk], color=NAVY, lw=0.9); tx(xk + 3.8, yk, "expected band (bp at right)"); xk += 34.0
        C.plot([xk, xk + 3], [yk, yk], color=GREYT, lw=0.35)
        tx(xk + 3.8, yk, f"ladder, {min(ladder):,}–{max(ladder):,} bp")
        tx(1.5, yk + 2.4, "Under the lanes: deletion frequency in the karyotype class. The heterokaryotype shows both bands; the short band is the allele with the deletion." + (" In silico: product at the locus in the reference (bp), off-target products." if ins else ""))
        if ins:
            C.plot([1.5, 4.5], [yk + 4.8, yk + 4.8], color=NAVY, lw=0.45, ls=(0, (3, 2))); tx(5.3, yk + 4.8, "† no product predicted at the locus in the reference.  ‡ off-target products predicted.")
        tx(1.5, yk + (7.2 if ins else 4.8), "Gels are simulated from the designed product sizes; no laboratory test." + (f"  No designed pair: {len(without)} interval(s)." if without and pg == npage - 1 else ""), fontweight="bold")
        if without and pg == npage - 1:
            LED("key", "intervals without a designed pair", len(without), f"{SRC_I}: n_primer_pairs == 0")
        note = f"page {pg + 1} of {npage}" + ("" if ins else "; in-silico line under each gel left out (no --reference)")
        out.append((f"page{pg + 1}of{npage}", cv.fig, f"Simulated gels of the PCR markers, page {pg + 1} of {npage}",
                    f"{len(sub)} primer pairs ({sub.marker.iloc[0]} to {sub.marker.iloc[-1]}), {per} gels per page (pcr.gels_per_page). Each gel: ladder and the expected bands of class 0, heterokaryotype and class 2 individuals, "
                    "band sizes at right, deletion frequency of the class under each lane. Simulated from the designed product sizes; no laboratory test.", note))
    return out


# =============================================================================== schematic
def fig_schematic(run, plt):
    """Generic schematic of the band-shift assay (no number of the run)."""
    if _primers(run) is None:
        return "no primer pair was designed"
    from matplotlib.patches import Rectangle
    W, H = 120.0, 52.0
    cv = style.canvas(plt, W + 10.0, H, style.SZ_ED, thin=True); C = cv.overlay; C.set_xlim(-9.0, W + 1.0)

    def tx(x, y, s, **kw):
        kw.setdefault("fontsize", 6.0); kw.setdefault("va", "center"); kw.setdefault("color", NAVY)
        return C.text(x, y, s, **kw)

    def arr(x0, x1, y, c):
        C.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=0.9, shrinkA=0, shrinkB=0, mutation_scale=5))
    tx(W / 2, 4.0, "Flanking-pair deletion assay (band shift)", ha="center", fontsize=7.5, fontweight="bold")
    seg, dl = "#DFE7EA", _tint(style.DEL_C, 0.55)
    tx(10.5, 15.0, "allele without\nthe deletion", ha="right", fontsize=5.5)
    C.add_patch(Rectangle((12, 13.5), 96, 3.0, fc=seg, ec="none")); C.add_patch(Rectangle((50.4, 13.5), 26.4, 3.0, fc=dl, ec=style.DEL_C, lw=0.3))
    tx(63.6, 15.0, "deleted segment", ha="center", fontsize=5.5)
    arr(24, 36, 11.2, FWD); arr(96, 84, 11.2, REV); tx(22.8, 11.2, "F", ha="right"); tx(97.2, 11.2, "R")
    tx(W / 2, 20.5, "product = flanks + deleted segment: long band", ha="center")
    tx(10.5, 32.0, "allele with\nthe deletion", ha="right", fontsize=5.5)
    C.add_patch(Rectangle((12, 30.5), 43.2, 3.0, fc=seg, ec="none")); C.add_patch(Rectangle((55.2, 30.5), 52.8, 3.0, fc=seg, ec="none")); C.plot([55.2, 55.2], [29.8, 34.2], color=NAVY, lw=0.6)
    arr(24, 36, 28.2, FWD); arr(79.2, 67.2, 28.2, REV); tx(22.8, 28.2, "F", ha="right"); tx(80.4, 28.2, "R"); tx(55.2, 27.0, "junction", ha="center", fontsize=5.5)
    tx(W / 2, 37.5, "product = flanks only: short band", ha="center")
    C.plot([12, 15], [43.5, 43.5], color=FWD, lw=0.9); tx(16, 43.5, "forward primer", fontsize=5.5); C.plot([36, 39], [43.5, 43.5], color=REV, lw=0.9); tx(40, 43.5, "reverse primer", fontsize=5.5)
    C.add_patch(Rectangle((60, 42.6), 3.0, 1.8, fc=seg, ec="none")); tx(64, 43.5, "sequence kept on both alleles", fontsize=5.5)
    C.add_patch(Rectangle((98, 42.6), 3.0, 1.8, fc=dl, ec=style.DEL_C, lw=0.3)); tx(102, 43.5, "deleted", fontsize=5.5)
    tx(W / 2, 48.5, "Both primers lie in sequence kept on both alleles (the split-read consensus across the junction, or reference flanks). Not to scale.", ha="center", fontsize=5.5)
    return cv.fig, "Schematic of the flanking-pair deletion assay", "One primer pair flanks the deletion: the allele without it gives the long product, the allele with it the short product, a heterozygote both. Schematic, not to scale; no number of the run."
