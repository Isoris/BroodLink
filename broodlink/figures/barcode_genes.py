"""Barcode at gene level (needs the optional input --genes).

Port of HAT/figures/upgrades_2026-11-03/round_ed/ED_Fig_H6/ED_Fig_H6.py (layout in mm, label
placement of `prepare`, cell colours, keys):
  A  barcode of the shared parents, one row per link and one column per barcode interval in
     genome order; plain cells = barcode call, hatched cells = states the barcode left
     unresolved and the obligate tag variants fill (tables of `broodlink/tag_evidence.py`);
     the intervals opened in B are boxed.
  B  for each opened interval: scale, genes (outlined = hit by an arrangement tag variant),
     the tag variants that hit genes ("calls", with their identifiers on leaders / brackets,
     adjacent calls merged), the tagging variants that count as evidence, and per link the
     obligate ones with three cells (barcode, tag-variant evidence, augmented) and the counts.

Which intervals are opened: the original shows four inversions chosen by the authors. Here:
the ids of config figures.gene_level_intervals when given; otherwise the
figures.gene_level_max_intervals barcode intervals with the most called parent states over
the links (ties: more genes, then genome order). They are drawn in genome order.

Left out, because a run has no generic input for them: gene families and their colours and
labels (NLR, TRIM, IFN / IL-1 of the original come from cohort inventories), the candidate
centromeric satellite array on the scale, and the names given to the links (SP1 ...).
"Hit" means coding sequence when the gene file has CDS features, else the gene span.
"""
from __future__ import annotations
import re

import numpy as np

from .. import genes as GN
from .. import tag_evidence as TE
from . import style as S
from .kinship import hap_label
from .shared_parent import noun  # noqa: F401  (same vocabulary as the other ported figures)

L = S.LED
SRC_B = "07_barcode/shared_parent_barcode.tsv"
SRC_E = "07_barcode/barcode_tag_evidence.tsv"
SRC_O = "07_barcode/obligate_tag_variants.tsv"
SRC_G = "07_barcode/barcode_interval_genes.tsv"
SRC_H = "07_barcode/barcode_interval_tag_gene_hits.tsv"
SRC_T = "06_structural_haplotypes/arrangement_tag_variants.tsv"
# colours of the original script (sky / lavender / navy family)
GREY, GENE, C0, C2, HET, BOX, GREYT = "#EEF1F4", "#C9D1D8", "#5E9BD1", S.ROH_BASE["navy"], "#9B84BF", "#d9482f", "#3D4650"
NAVY = S.NAVY
HATCH = "/////"
FS, CWM, ROWS, LANE = 5.5, 1.0, 2.25, 2.45
GX0, GW = 25.0, 124.0
CLS = {"": 0, "class 0": 1, "class 2": 2, "class 0 + class 2": 3}


def size_label(b):
    """del_gene_loci_common.size_label of the original."""
    return f"{b:,.0f} bp" if b < 1000 else (f"{b / 1e3:,.0f} kb" if b < 1e6 else f"{b / 1e6:.2f} Mb")


def cell_colour(p, cls):
    if p >= 0:
        return {0: C0, 1: HET, 2: C2}[int(p)], None
    if cls > 0:
        return {1: C0, 2: C2, 3: HET}[int(cls)], HATCH
    return GREY, None


def _tables(run):
    """Tables of step 7 (computed here when the run was made before --genes was given)."""
    t = dict(run.t.get("barcode", {}))
    for mod, keys in ((TE, ("obligate_tags", "tag_evidence")), (GN, ("genes", "gene_hits"))):
        if any(k not in t for k in keys):
            for key, df, _f, _d in mod.barcode_tables(run):
                t.setdefault(key, df)
    return t


def _symbol(gene, name):
    """A gene symbol worth printing: a name that is not the identifier and not a placeholder (LOC<number>)."""
    n = str(name or "")
    return n if n and n != str(gene) and not re.fullmatch(r"LOC\d+", n) else ""


def _tick_step(length_bp):
    """Largest step of 1 or 5 x 10^k Mb that leaves more than five intervals on the scale
    (0.5 Mb up to 5 Mb and 1 Mb above, as in the original)."""
    mb = length_bp / 1e6; best = None
    for k in range(-4, 4):
        for m in (1, 5):
            s = m * 10.0 ** k
            if mb / s > 5:
                best = s
    return best or mb


def _prepare(iv, gg, calls):
    """Positions of genes and call labels of one interval (function `prepare` of the original):
    labels are placed first so that each interval takes only the height it needs."""
    lo, hi = int(iv.start), int(iv.end); L_ = max(hi - lo, 1)
    X = lambda pos: GX0 + (min(max(pos, lo), hi) - lo) / L_ * GW  # noqa: E731
    sym = {g: _symbol(g, n) for g, n in zip(gg.gene, gg.name)}
    grp = []
    for cid, c in calls.groupby("call_id", sort=False):
        hg = gg[gg.gene.isin(set(c.gene))]
        if not len(hg):
            continue
        grp.append(dict(x1=X(hg.start.min()), x2=X(hg.end.max()), ids=[cid], genes=list(dict.fromkeys(hg.gene)), ticks=[X((r.start + r.end) / 2) for r in hg.itertuples()],
                        size=int(c.length_bp.iloc[0])))
    grp.sort(key=lambda g: g["x1"]); merged = []
    for g in grp:
        if merged and g["x1"] - merged[-1]["x2"] <= 8.0:
            m = merged[-1]; m["x2"] = max(m["x2"], g["x2"]); m["ids"] += g["ids"]; m["genes"] = list(dict.fromkeys(m["genes"] + g["genes"])); m["ticks"] += g["ticks"]
        else:
            merged.append(dict(g))
    items = []
    for m in merged:
        gs = m["genes"]; named = [sym[g] for g in gs if sym.get(g)]; idtxt = m["ids"][0] if len(m["ids"]) == 1 else f"{m['ids'][0]}–{m['ids'][-1]}"
        if len(m["ids"]) > 1:
            txt = f"{idtxt}: {len(gs)} gene{'s' if len(gs) != 1 else ''}"
        elif len(gs) > 2:
            txt = f"{idtxt}: {len(gs)} genes ({size_label(m['size'])})"
        elif len(named) == len(gs):
            txt = f"{idtxt} {', '.join(named)}"
        else:
            txt = idtxt
        items.append(((m["x1"] + m["x2"]) / 2 if len(gs) > 2 else m["x1"], m["x1"], m["x2"], txt, len(gs) > 1, m["ticks"]))
    items.sort(key=lambda t: t[1]); placed = []; out = []
    for anc, xa, xb, txt, br, ticks in items:
        lx = max(GX0, min(anc if not br else xa, GX0 + GW - len(txt) * CWM)); ext = (lx - .6, lx + len(txt) * CWM + .6); lane = 0; tries = 0
        while tries < 80:
            tries += 1
            cover = [e2 for l2, e2, _ in placed if l2 < lane and e2[0] - .3 <= lx <= e2[1] + .3]
            if cover:
                lx = max(e[1] for e in cover) + .6; ext = (lx - .6, lx + len(txt) * CWM + .6); continue
            if any(l2 == lane and not (ext[1] < e2[0] or ext[0] > e2[1]) for l2, e2, _ in placed):
                lane += 1; continue
            break
        if lx + len(txt) * CWM > GX0 + GW + 14:
            lx = GX0 + GW + 14 - len(txt) * CWM; ext = (lx - .6, lx + len(txt) * CWM + .6); lane = max([l2 for l2, _, _ in placed] + [0]) + 1
        placed.append((lane, ext, txt)); out.append((xa, xb, lx, lane, txt, br, ticks))
    nl = (max(l for l, _, _ in placed) + 1) if placed else 0
    return dict(iv=iv, lo=lo, hi=hi, L=L_, X=X, gg=gg, labels=out, nl=nl)


def fig_barcode_gene_level(run, plt):
    """Barcode at gene level: A barcode of the links with the tag-variant evidence, B genes and
    tag variants of the opened intervals. Intervals opened: config figures.gene_level_intervals,
    else the figures.gene_level_max_intervals intervals with the most called parent states
    (ties: more genes, then genome order)."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    if getattr(run, "genes", None) is None:
        return "no --genes input"
    B = run.t.get("barcode", {}).get("barcode")
    hap, _info = run.hap_set()
    if B is None or not len(B) or hap is None:
        return "no barcode (no structural-haplotype genotypes or no link between families)"
    Bl = B[B.kind == "link"]
    if not len(Bl):
        return "the barcode holds no link between families"
    T = _tables(run)
    I = GN.barcode_intervals(hap)
    I = I[I.interval.isin(set(Bl.interval))].reset_index(drop=True)
    if not len(I) or not I.has_coordinates.any():
        return "the barcode intervals have no coordinates (interval ids are not chrom:start-end)"
    I["number"] = np.arange(1, len(I) + 1)                           # column of panel A (the tables number the same way when every interval is in the barcode)
    num_tab = dict(zip(I.interval, I.number)); N = len(I)
    GT_, HT = T.get("genes"), T.get("gene_hits")
    if GT_ is None or HT is None:
        return "the genes of the barcode intervals could not be tabulated"
    notes = []
    fig_cfg = run.cfg["figures"]; pc = run.cfg["barcode"]["posterior_call"]
    _mc, kmin, dmin = TE.settings(run.cfg)
    tags = run.t.get("haplotypes", {}).get("tags")
    O, E = T.get("obligate_tags"), T.get("tag_evidence")
    kinds = [k for k in (sorted(set(tags.svtype)) if tags is not None and len(tags) else [])]
    kind = "DEL" if "DEL" in kinds else (max(kinds, key=lambda k: int((tags.svtype == k).sum())) if kinds else None)
    has_ev = kind is not None and E is not None and O is not None
    kw = S.class_label(kind, plural=False) if kind else ""
    if has_ev:
        tg = tags[(tags.svtype == kind) & tags.interval.isin(set(I.interval))]
        E = E[E.svtype == kind]; O = O[(O.svtype == kind) & O.counts_as_evidence.astype(bool)]; HT = HT[HT.svtype == kind]
        if len(kinds) > 1:
            notes.append(f"tag variants of class {kw} only (as the original); the other classes ({', '.join(S.class_label(k, plural=False) for k in kinds if k != kind)}) are in the tables")
    else:
        HT = HT.iloc[0:0]
        notes.append("no arrangement tag variant in the run (no DEL / DUP input, or none passes inversions.tags): hatched cells, tag-variant tracks and obligate rows are left out")
    hit_word = "coding sequence" if len(HT) and (HT.hit == "coding sequence").all() else ("coding sequence" if getattr(run, "genes_cds", None) is not None and len(run.genes_cds) else "a gene")
    # ---- links (rows): supported first, then by family number
    links = Bl.drop_duplicates("parent")
    order = sorted(range(len(links)), key=lambda i: (not str(links.link_status.iloc[i]).startswith("supported"), [S.fam_key(f) for f in str(links.families.iloc[i]).split(",")]))
    links = links.iloc[order]
    rows_all = [(r.parent, str(r.families).split(",")) for r in links.itertuples()]
    call = {(r.parent, r.interval): {"0/0": 0, "0/1": 1, "1/1": 2}.get(r.call, -1) for r in Bl.itertuples()}
    ev = {(r.link, r.interval): (CLS.get(str(r.evidence_class), 0), int(r.n_obligate_class_0), int(r.n_obligate_class_2)) for r in E.itertuples()} if has_ev else {}
    # ---- intervals opened in B
    ncalled = {iv: sum(call.get((nm, iv), -1) >= 0 for nm, _f in rows_all) for iv in I.interval}
    ngenes = GT_.groupby("interval").gene.nunique().to_dict() if len(GT_) else {}
    cand = I[I.has_coordinates]
    asked = [str(x) for x in (fig_cfg.get("gene_level_intervals") or [])]
    nmax = max(1, int(fig_cfg.get("gene_level_max_intervals", 4)))
    if asked:
        sel = [x for x in asked if x in set(cand.interval)]
        if len(sel) < len(asked):
            notes.append(f"{len(asked) - len(sel)} of the intervals named in figures.gene_level_intervals are not barcode intervals with coordinates")
        rule = "the intervals named in figures.gene_level_intervals"
    else:
        sel = []
    if not sel:
        sel = sorted(cand.interval, key=lambda iv: (-ncalled[iv], -ngenes.get(iv, 0), num_tab[iv]))[:nmax]
        rule = f"the {len(sel)} barcode intervals with the most called parent states over the links (ties: more genes, then genome order)"
    sel = sorted(sel, key=lambda iv: num_tab[iv])
    # ---- prepare the zooms, then fit the page (A4): fewer link rows first, then fewer intervals
    def prep(iv_name):
        iv = I[I.interval == iv_name].iloc[0]
        gg = GT_[GT_.interval == iv_name].drop_duplicates(["gene", "start"]).sort_values("start")
        return _prepare(iv, gg, HT[HT.interval == iv_name])
    PR = [prep(x) for x in sel]
    R = len(rows_all); rh = 2.7
    def zoom_h(p, r):
        return 2.8 + 2.6 + p["nl"] * LANE + 2.3 + (2.6 + 2.6 if has_ev else 0) + 2.5 + r * ROWS + 2.2 + 2.6
    def height(r, prs):
        return 10.2 + r * (rh + .4) + 5.3 + 3.2 + 7.4 + sum(zoom_h(p, r) for p in prs) + (8.0 if has_ev else 5.5)
    while height(R, PR) > S.A4_H and R > 4:
        R -= 1
    while height(R, PR) > S.A4_H and len(PR) > 1:
        PR = PR[:-1]
    while height(R, PR) > S.A4_H and R > 1:
        R -= 1
    if R < len(rows_all):
        notes.append(f"{R} of {len(rows_all)} links drawn (page height; supported links first)")
    if len(PR) < len(sel):
        notes.append(f"{len(PR)} of {len(sel)} chosen intervals opened (page height)")
        sel = sel[:len(PR)]
    rows = rows_all[:R]
    yA = 10.2; yEnd = yA + R * (rh + .4); yKey = yEnd + 5.3; yB = yKey + 3.2; yy = yB + 7.4
    W = S.WIDTH_FULL; H = height(R, PR)
    cv = S.canvas(plt, W, H, S.SZ_ED)
    import matplotlib
    matplotlib.rcParams.update({"hatch.linewidth": .35})
    ax_ = cv.overlay

    def T_(x, y, s, fs=FS, b=False, c=NAVY, ha="left", va="center", **k):
        return ax_.text(x, y, s, fontsize=fs, fontweight="bold" if b else "normal", color=c, ha=ha, va=va, **k)

    def BAN(y, n, text):
        fs = 6.8
        while fs > FS and cv.text_width(text, fs, bold=True) > W - 6.0:
            fs -= .2
        T_(1, y, n, 8, True); T_(5, y, text, fs, True); ax_.add_line(Line2D([.5, W - .5], [y + 2.4] * 2, color=NAVY, lw=.35))

    hl = hap_label(run)
    fills = {nm: sum((call.get((nm, iv), -1) >= 0) or (ev.get((nm, iv), (0, 0, 0))[0] > 0) for iv in I.interval) for nm, _f in rows}
    lo_f, hi_f = min(fills.values()), max(fills.values())
    L("banner1", "intervals informed per link (min, max)", f"{lo_f}, {hi_f}", SRC_B + (" + " + SRC_E if has_ev else ""))
    rng_txt = f"{lo_f}–{hi_f}" if lo_f != hi_f else f"{lo_f}"
    n_open = len(sel)
    BAN(2.2, "1", (f"With {kw} evidence each shared-parent link informs " if has_ev else "The barcode calls ") + f"{rng_txt} of the {L('banner1', 'barcode intervals', N, SRC_B)} {hl}"
        + (" per link" if not has_ev else "") + f"; the {n_open} boxed {'are' if n_open != 1 else 'is'} opened below")
    T_(.5, 7.0, "A", 8, True); T_(4.1, 7.0, ("Augmented barcode" if has_ev else "Barcode") + " of the shared parents", 6.5, True)
    # ---- panel A
    x0 = 25.0; cw = min(140.0 / N, 12.0); wA = cw * N
    for r, (nm, fams) in enumerate(rows):
        y = yA + r * (rh + .4)
        T_(x0 - 1.0, y + rh / 2, " × ".join(fams), ha="right")
        for i, iv in enumerate(I.interval):
            c, ht = cell_colour(call.get((nm, iv), -1), ev.get((nm, iv), (0, 0, 0))[0])
            ax_.add_patch(Rectangle((x0 + i * cw, y), cw, rh, fc=c, ec="white", lw=.15, hatch=ht, zorder=5))
        ax_.add_patch(Rectangle((x0, y), wA, rh, fill=False, ec=NAVY, lw=.3, zorder=6))
        T_(x0 + wA + 1.2, y + rh / 2, f"{L('A', f'intervals informed, {nm}', fills[nm], SRC_B + (' + ' + SRC_E if has_ev else ''))} of {N}", c=GREYT)
    seln = {int(num_tab[iv]) for iv in sel}
    step = max(1, int(np.ceil(6.0 / cw)))
    for i in range(N):
        if i % step == 0 or (i + 1) in seln:
            near = any(abs(i + 1 - s) < step and (i + 1) != s for s in seln)
            if near and (i + 1) not in seln and cw * 1 < 3.0:
                continue
            T_(x0 + (i + .5) * cw, yEnd + .1, f"{i + 1}", ha="center", va="top", b=(i + 1) in seln, c=BOX if (i + 1) in seln else GREYT)
    for s in seln:
        ax_.add_patch(Rectangle((x0 + (s - 1) * cw - .1, yA - .5), cw + .2, yEnd - yA + .3, fill=False, ec=BOX, lw=.7, zorder=8))
    key = [(C0, None, "0/0"), (HET, None, "0/1"), (C2, None, "1/1")]
    if has_ev:
        key += [(C0, HATCH, f"≥ 1 class 0 ({kw})"), (C2, HATCH, f"≥ 1 class 2 ({kw})"), (HET, HATCH, f"class 0 + class 2 ({kw})")]
    key += [(GREY, None, "no information")]
    x = x0
    for c, ht, lab in key:
        ax_.add_patch(Rectangle((x, yKey - .85), 2.8, 1.7, fc=c, ec=NAVY if c == GREY else ("white" if ht else c), lw=.25, hatch=ht)); T_(x + 3.5, yKey, lab)
        x += 3.5 + cv.text_width(lab, FS) + 4.0
    T_(x0 - 1, yEnd + 1.1, "interval", ha="right", c=GREYT)
    # ---- banner 2 and the zooms
    if has_ev:
        cs = HT[HT.interval.isin(set(sel))].drop_duplicates("call_id")
        n3 = int(len(cs)); nw = int((cs.difference < dmin).sum())
        L("banner2", "calls hitting genes in the opened intervals", n3, SRC_H); L("banner2", f"of them with difference < {dmin:g}", nw, SRC_H)
        BAN(yB, "2", f"The {kw}s that fill the cells are small; " + (f"{nw} of {n3} gene-hitting calls in the opened intervals are weakly tagged and stay out of the evidence" if n3 else
                                                                 f"no tag {kw} hits {hit_word} in the opened intervals"))
        T_(.5, yB + 4.6, "B", 8, True); T_(4.1, yB + 4.6, f"Inside the cell: genes and the {kw}s that fill it", 6.5, True)
    else:
        BAN(yB, "2", f"Genes of the opened {hl} and the barcode call of every link")
        T_(.5, yB + 4.6, "B", 8, True); T_(4.1, yB + 4.6, "Inside the cell: genes", 6.5, True)
    top = yy
    for pz in PR:
        iv, lo, hi, X, gg, Lz = pz["iv"], pz["lo"], pz["hi"], pz["X"], pz["gg"], pz["L"]; name = iv.interval; num = int(num_tab[name])
        ng = int(gg.gene.nunique()); hitg = set(HT[HT.interval == name].gene) if has_ev else set()
        head = f"{iv.chrom}:{lo / 1e6:.2f}–{hi / 1e6:.2f} Mb:  {L('B', f'genes, interval {num}', ng, SRC_G)} gene{'s' if ng != 1 else ''}"
        if has_ev:
            cd = tg[(tg.interval == name) & (tg.difference >= dmin)]
            ncall = int(HT[HT.interval == name].call_id.nunique()); n0, n2 = int((cd.tagged_arrangement == "class 0").sum()), int((cd.tagged_arrangement == "class 2").sum())
            head += (f",  {L('B', f'calls hitting genes, interval {num}', ncall, SRC_H)} {kw} call{'s' if ncall != 1 else ''} hit{'s' if ncall == 1 else ''} {hit_word},  "
                     f"{L('B', f'class-0-tagged, interval {num}', n0, SRC_T)} class-0- and {L('B', f'class-2-tagged, interval {num}', n2, SRC_T)} class-2-tagged {kw}s")
        ax_.add_patch(Rectangle((.8, top - .1), 2.0, 2.0, fc="none", ec=BOX, lw=.7)); T_(3.8, top + .9, f"{num}", 6.5, True); T_(10.5, top + .9, head)
        ya = top + 2.8 + 2.6; ax_.plot([GX0, GX0 + GW], [ya, ya], color=NAVY, lw=.35); stp = _tick_step(Lz)
        k0 = int(np.ceil(lo / 1e6 / stp - 1e-9))
        for kk in range(k0, int(np.floor(hi / 1e6 / stp + 1e-9)) + 1):
            mb = kk * stp
            if mb * 1e6 < lo or mb * 1e6 > hi:
                continue
            ax_.plot([X(mb * 1e6)] * 2, [ya - .7, ya], color=NAVY, lw=.3); T_(X(mb * 1e6), ya - .8, f"{round(mb, 6):g}", ha="center", va="bottom", c=GREYT)
        T_(GX0 - 1, ya, "Mb", ha="right", c=GREYT)
        yg = ya + 2.6 + pz["nl"] * LANE + 1.2
        for r in gg.itertuples():
            hit = r.gene in hitg
            ax_.add_patch(Rectangle((X(r.start), yg - .9), max(X(r.end) - X(r.start), .15), 1.8, fc=GENE, ec=NAVY if hit else "none", lw=.35, zorder=3 if hit else 2))
        for xa, xb, lx, lane, txt, br, ticks in pz["labels"]:
            yb = yg - 3.0 - lane * LANE
            if br:
                ax_.plot([xa, xb], [yg - 1.9, yg - 1.9], color=NAVY, lw=.4, zorder=6)
                for tx in ticks:
                    ax_.plot([tx, tx], [yg - 1.9, yg - 1.4], color=NAVY, lw=.3, zorder=6)
                if xa - .01 <= lx <= xb + .01:
                    ax_.plot([lx, lx], [yg - 1.9, yb + .2], color=NAVY, lw=.3, zorder=6)
                else:
                    ax_.plot([xa, lx], [yg - 1.9, yb + .2], color=NAVY, lw=.3, zorder=6)
            elif abs(ticks[0] - lx) > .05:
                ax_.plot([ticks[0], lx, lx], [yg - 1.4, yg - 2.1, yb + .2], color=NAVY, lw=.3, zorder=6)
            else:
                ax_.plot([lx, lx], [yg - 1.4, yb + .2], color=NAVY, lw=.3, zorder=6)
            T_(lx + .3, yb, txt, b=True, c=NAVY, va="bottom", zorder=7); L("B label", f"interval {num}", txt, SRC_H + " / " + SRC_G)
        T_(GX0 - 1, yg, "genes", ha="right")
        ye = yg + 2.5
        if has_ev:
            ycl = yg + 2.6; ax_.add_patch(Rectangle((GX0, ycl - .95), GW, 1.9, fc="white", ec="#C9D1D8", lw=.2, zorder=1))
            for cr in HT[HT.interval == name].drop_duplicates("call_id").itertuples():
                x1, x2 = X(cr.start), X(cr.end); ax_.add_patch(Rectangle((x1, ycl - .95), max(x2 - x1, .25), 1.9, fc=C0 if cr.tagged_arrangement == "class 0" else C2, ec="none", zorder=4))
            T_(GX0 - 1, ycl, f"{kw} calls hitting genes", ha="right")
            yd = ycl + 2.6; ax_.add_patch(Rectangle((GX0, yd - .95), GW, 1.9, fc="white", ec="#C9D1D8", lw=.2, zorder=1))
            for r in cd.itertuples():
                x1, x2 = X(r.start), X(r.end); ax_.add_patch(Rectangle((x1, yd - .95), max(x2 - x1, .2), 1.9, fc=C0 if r.tagged_arrangement == "class 0" else C2, ec="none", zorder=4))
            T_(GX0 - 1, yd, f"tagging {kw}s", ha="right")
            if len(cd):
                T_(GX0 + GW + 3, yd, "median " + size_label(L("B", f"median length of the tagging {kw}s, interval {num}", float(np.median(np.abs(cd.end - cd.start))), SRC_T)), c=GREYT)
            ye = yd + 2.5
        for r_, (nm, fams) in enumerate(rows):
            y = ye + r_ * ROWS + .2
            pv = call.get((nm, name), -1); cvv, n0o, n2o = ev.get((nm, name), (0, 0, 0)); xr = GX0 + GW + 3.0
            lab = "×".join(fams)
            if has_ev:
                ax_.add_patch(Rectangle((GX0, y - .9), GW, 1.8, fc="white", ec="#C9D1D8", lw=.2, zorder=1))
                for rr in O[(O.link == nm) & (O.interval == name)].itertuples():
                    x1, x2 = X(rr.start), X(rr.end); ax_.add_patch(Rectangle((x1, y - .9), max(x2 - x1, .2), 1.8, fc=C0 if rr.tagged_arrangement == "class 0" else C2, ec="none", zorder=4))
                T_(GX0 - 1, y, f"obligate, {lab}", ha="right")
                c1 = {0: C0, 1: HET, 2: C2}.get(pv, GREY); c2 = {1: C0, 2: C2, 3: HET}.get(cvv, GREY); c3, h3 = cell_colour(pv, cvv)
                for off, (cf, hh) in enumerate(((c1, None), (c2, HATCH if cvv > 0 else None), (c3, h3))):
                    ax_.add_patch(Rectangle((xr + off * 4.3, y - .9), 3.9, 1.8, fc=cf, ec="white", lw=.2, hatch=hh, zorder=5))
                T_(xr + 13.3, y, f"{L('B', f'obligate class-0 tags, interval {num}, {nm}', n0o, SRC_E)} c0, {L('B', f'obligate class-2 tags, interval {num}, {nm}', n2o, SRC_E)} c2", c=GREYT)
            else:
                T_(GX0 + 16, y, lab, ha="right")
                ax_.add_patch(Rectangle((GX0 + 17.5, y - .9), 3.9, 1.8, fc={0: C0, 1: HET, 2: C2}.get(pv, GREY), ec="white", lw=.2, zorder=5))
                L("B", f"barcode call, interval {num}, {nm}", {0: "0/0", 1: "0/1", 2: "1/1"}.get(pv, "unresolved"), SRC_B)
        top += zoom_h(pz, R)
    # ---- key
    yl = top + .9; xk = 25.0
    items = [(GENE, "gene")] + ([(C0, f"class-0-tagged {kw}"), (C2, f"class-2-tagged {kw}")] if has_ev else [(C0, "0/0"), (HET, "0/1"), (C2, "1/1"), (GREY, "unresolved")])
    for c, lab in items:
        ax_.add_patch(Rectangle((xk, yl - .85), 2.4, 1.7, fc=c, ec="none")); T_(xk + 3.0, yl, lab); xk += 3.0 + cv.text_width(lab, FS) + 4.0
    if has_ev:
        T_(xk + 2, yl, "0/0, 1/1: homokaryotypes of class 0, class 2 (labels without ancestral meaning)", c=GREYT)
        ax_.add_patch(Rectangle((25, yl + 1.65), 2.4, 1.7, fc=GENE, ec=NAVY, lw=.4)); T_(28, yl + 2.5, f"outlined: gene hit in {hit_word.replace('a gene', 'its span')} by an arrangement tag {kw}; brackets and leaders: call ID")
        T_(25, yl + 5.0, f"Three cells per link, left to right: barcode, {kw} evidence, augmented (colours as in A); c0, c2: class-0- and class-2-tagged obligate {kw}s of that link")
    else:
        T_(25, yl + 2.5, "One cell per link: barcode call of the shared parent (0/0, 1/1: homokaryotypes of class 0, class 2)")
    what = GN.drawn_genes(run.genes)[1]
    notes = [f"intervals opened: {rule}", f"genes drawn: {what}; a gene belongs to an interval when it overlaps it"] + \
            ([f"hit = a tag variant overlaps {hit_word.replace('a gene', 'the gene span (no CDS in the gene file)')}; calls = tag variants of step 6 (difference >= {run.cfg['inversions']['tags']['min_delta']:g}); "
              f"evidence = obligate tag variants with difference >= {dmin:g}, at least {kmin} per class"] if has_ev else []) + notes + \
            ["left out (no generic input): gene-family colours and labels, candidate centromeric satellite array, names of the links"]
    title = "Barcode at gene level"
    desc = (f"A barcode of the shared parents per link (barcode call at posterior >= {pc:g}" + (f"; hatched: state filled by obligate tag {kw}s" if has_ev else "") + f"), {N} {hl} in genome order, "
            f"{len(sel)} boxed and opened in B. B genes of each opened interval" + (f", the tag {kw}s that hit {hit_word}, the tagging {kw}s, and per link the obligate ones with the cells barcode / {kw} evidence / augmented." if has_ev else " and the barcode call of every link."))
    return cv.fig, title, desc, "; ".join(notes)
