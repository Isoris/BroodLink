"""Relatedness and the full-sib families: port of panels A-E of the published relatedness figure
(`build_hatchery_fig_A_v6.py`): A IBD coefficient simplex, B relatedness-class matrix,
C family sizes with the pairs per kinship band and the closure sensitivity to the k2 cut,
D within- against between-family kinship, E full-sib closure per family.

Drawn from `01_relatedness/relatedness_pairs.tsv`, `02_families/family_membership.tsv` and
`02_families/families.tsv`. Panels that need k0 / k1 / k2 (A, the k2 inset of C) are left out
when relatedness comes from a GRM. The ideogram and the SV panels of the published figure are
other figures (or out of scope: a genome display).
"""
from __future__ import annotations
import itertools

import numpy as np

from . import style as S
from .shared_parent import noun

L = S.LED
SRC_R = "01_relatedness/relatedness_pairs.tsv"
SRC_M = "02_families/family_membership.tsv"
SRC_F = "02_families/families.tsv"


def _stat_note(ax, text, sz, loc="upper left"):
    x, ha = (0.97, "right") if "right" in loc else (0.03, "left")
    ax.text(x, 0.98, text, transform=ax.transAxes, fontsize=sz["note"], ha=ha, va="top", color=S.NOTE_GREY)


def fig_relatedness(run, plt):
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Rectangle
    P = run.t["relatedness"]["pairs"].copy()
    if "membership" not in run.t.get("families", {}) or not len(P):
        return None
    M = run.t["families"]["membership"]; F = run.t["families"]["families"]
    cfg = run.cfg; r = cfg["relatedness"]; fcfg = cfg["families"]
    BRK = [r["theta_third"], r["theta_second"], r["theta_first"], r["theta_duplicate"]]
    fish = M.individual.astype(str).tolist(); fam = dict(zip(fish, M.family.astype(str)))
    P = P[P.id1.isin(fam) & P.id2.isin(fam)].copy()
    P["fa"], P["fb"] = P.id1.map(fam), P.id2.map(fam)
    P["within"] = (P.fa == P.fb) & (P.fa != "unassigned"); P["band"] = np.digitize(P.theta, BRK)
    has_k = bool(P.k2.notna().any()) and bool(P.k1.notna().any())
    size = M[M.family != "unassigned"].groupby("family").size()
    fams = sorted(size.index, key=S.fam_key); nfam = len(fams); n = len(fish); n_un = int((M.family == "unassigned").sum())
    MIN_FAM = int(cfg["audit"]["min_family_size"])
    fcol = S.family_colours(fams); notes = []
    SAND, BRICK = S.ROH_BASE["sand"], S.ROH_BASE["brick"]
    sz = dict(S.SZ_MAIN); W = S.WIDTH_FULL
    # closure of every within-family pair
    Wc = P[P.within].copy()
    Wc["ok"] = (Wc.theta >= fcfg["theta_min"]) & ((Wc.k2 >= fcfg["k2_min"]) if has_k else True)
    okset = {frozenset((a_, b_)): bool(o) for a_, b_, o in zip(Wc.id1, Wc.id2, Wc.ok)}
    mem_of = {f: M.individual[M.family == f].astype(str).tolist() for f in fams}
    ncol = 13; rows_e = int(np.ceil(max(nfam, 1) / ncol)); max_rows = 8
    if rows_e > max_rows:
        notes.append(f"closure graphs drawn for the {max_rows * ncol} largest of {nfam} families"); rows_e = max_rows
    rh = 11.5 if rows_e <= 4 else 10.0
    eh = max(42.0, rows_e * rh - 2.0)
    H = 78 + eh + 12.0
    cv = S.canvas(plt, W, H, sz); fig = cv.fig; letters = iter("ABCDEFG")
    dx = 0.0 if has_k else -64.0
    # ---------------------------------------------------------------- A simplex
    if has_k:
        cv.title(next(letters), "IBD coefficient simplex", 4.0, 1.0)
        SQ3 = np.sqrt(3) / 2; UNIT = 46.0 / SQ3
        ax = cv.ax(3.7, 12 - (1.045 - SQ3) * UNIT, 1.2 * UNIT, 1.2 * UNIT)
        tern = lambda k1, k2: (np.asarray(k1) + .5 * np.asarray(k2), SQ3 * np.asarray(k2))
        ks = (P.k0 + P.k1 + P.k2).to_numpy(); ks = np.where(ks > 0, ks, 1.0); P["n1"], P["n2"] = P.k1 / ks, P.k2 / ks
        ax.plot(*tern([0, 1, 0, 0], [0, 0, 1, 0]), color=S.NAVY, lw=S.LW["axes"], zorder=1)
        for q in (.25, .5, .75):
            for a_, b_ in ((tern(q, 0), tern(q, 1 - q)), (tern(0, q), tern(1 - q, q)), (tern(0, 1 - q), tern(q, 1 - q))):
                ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color="#E3E8EC", lw=.3, zorder=0)
        cls = [("below 3rd", P.band == 0, "#C9D1D8"), ("3rd degree", P.band == 1, S.PAL[3]), ("2nd degree", P.band == 2, SAND),
               ("1st, between fam.", (P.band >= 3) & ~P.within, S.PAL[4]), ("within full-sib", P.within, BRICK)]
        for lab, m, c in cls:
            g = P[m]
            x_, y_ = tern(g.n1, g.n2)
            ax.scatter(x_, y_, s=1.8 if lab != "within full-sib" else 2.6, c=c, lw=0, alpha=.7 if lab != "below 3rd" else .45,
                       label=f"{lab} ({L('A', 'pairs ' + lab, len(g), SRC_R):,})", rasterized=True, zorder=2 + (lab == "within full-sib"))
        po = P[(P.theta >= r["theta_first"]) & (P.k2 < .02)]
        x_, y_ = tern(po.n1, po.n2); ax.scatter(x_, y_, s=22, marker="+", c=S.RED, lw=0.9, zorder=7)
        k1m, k2m = cfg["qc"]["po_k1_min"], cfg["qc"]["po_k2_max"]
        ax.fill(*tern([k1m, 1, 1 - k2m, k1m], [0, 0, k2m, k2m]), fc=S.RED, alpha=.12, ec=S.RED, lw=.4, ls="--", zorder=1)
        for k1, k2 in ((0, 0), (1, 0), (.5, .25), (.5, 0), (.25, 0)):
            xx, yy = tern(k1, k2); ax.scatter([xx], [yy], marker="*", s=34, fc="white", ec=S.NAVY, lw=.6, zorder=8)
        for m_, c in ((P[P.within], BRICK), (P[P.band == 2], "#C9A227"), (P[P.band == 1], S.PAL[3]), (P[P.band == 0], "#8A949C")):
            if len(m_):
                xx, yy = tern(m_.n1.median(), m_.n2.median()); ax.scatter([xx], [yy], marker="o", s=16, fc=c, ec=S.NAVY, lw=.6, zorder=9)
        fs_x, fs_y = tern(.5, .25); ax.text(fs_x - .14, fs_y + .17, "full-sib", ha="right", fontsize=sz["small"], fontweight="bold", color=S.dark_text(BRICK))
        hx, hy = tern(.5, 0); ax.text(hx, hy + .045, "half-sib", fontsize=sz["small"], ha="center")
        tx, ty = tern(.25, 0); ax.text(tx, ty + .045, "3rd", fontsize=sz["small"], ha="center")
        ax.text(*tern(0, 1.0), "k$_2$ = 1", fontsize=sz["small"], ha="center", va="bottom")
        ax.text(-.02, -.05, "k$_0$ = 1\nunrelated", fontsize=sz["small"], ha="center", va="top")
        ax.text(1.05, -.05, "k$_1$ = 1\nparent–offspring", fontsize=sz["small"], ha="right", va="top", color=S.RED, fontweight="bold")
        ax.text(1.085, .62, "open star: expected\nfilled circle:\nobserved median", fontsize=sz["small"] - .3, ha="right", va="center", color=S.NOTE_GREY)
        ax.set_xlim(-.10, 1.10); ax.set_ylim(-.155, 1.045); ax.set_aspect("equal"); ax.axis("off")
        hs, ls = ax.get_legend_handles_labels()
        fig.legend(hs, ls, loc="upper left", bbox_to_anchor=(6.0 / W, 1 - 64.0 / H), bbox_transform=fig.transFigure, fontsize=sz["small"] - .5, markerscale=2.0, handletextpad=.1, labelspacing=.15, columnspacing=.7, ncol=3, borderaxespad=0, frameon=False)
    else:
        notes.append("IBD simplex and the k2 sensitivity inset omitted: relatedness comes from a GRM (no k0 / k1 / k2)")
    # ---------------------------------------------------------------- B class matrix
    ix = {s: i for i, s in enumerate(fish)}
    T = np.zeros((n, n)); a = P.id1.map(ix).to_numpy(); b = P.id2.map(ix).to_numpy()
    T[a, b] = P.theta.to_numpy(); T[b, a] = P.theta.to_numpy(); np.fill_diagonal(T, np.nan)
    fi = {f: [ix[s] for s in mem_of[f]] for f in fams}
    if nfam >= 3:
        with np.errstate(invalid="ignore"):
            FM = np.array([[np.nanmean(T[np.ix_(fi[f], fi[g])]) for g in fams] for f in fams]); FM = np.nan_to_num(FM); np.fill_diagonal(FM, 0); FM = np.clip(FM, 0, None)
        w_, v_ = np.linalg.eigh(np.diag(FM.sum(1)) - FM); fam_order = [fams[i] for i in np.argsort(v_[:, 1], kind="stable")]
    else:
        fam_order = list(fams)

    def leaf(idx):
        if len(idx) < 3:
            return list(idx)
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform
        sub = np.nan_to_num(T[np.ix_(idx, idx)], nan=1.0); mx = sub.max() or 1.0; d = squareform(1 - sub / mx, checks=False)
        return [idx[i] for i in leaves_list(linkage(np.clip(d, 0, None), "average"))]
    order = [i for f in fam_order for i in leaf(fi[f])]
    un = [ix[s] for s in fish if fam[s] == "unassigned"]
    un.sort(key=lambda i: -np.nanmax(np.delete(T[i], i)) if n > 1 else 0)
    order += un
    cv.title(next(letters), "Relatedness-class matrix", 4.0, 64.0 + dx)
    axm = cv.ax(70 + dx, 12, 50, 46)
    R = T[np.ix_(order, order)]
    cm = ListedColormap(S.BAND_COL); nm = BoundaryNorm([-9, *BRK, 9], cm.N)
    axm.imshow(np.ma.masked_invalid(R), cmap=cm, norm=nm, interpolation="nearest", aspect="auto", rasterized=True)
    axm.set_xticks([]); axm.set_yticks([]); S.frame(axm, S.LW["axes"])
    start = 0; cb_x = axm.inset_axes([0, 1.012, 1, .035]); cb_y = axm.inset_axes([-.04, 0, .03, 1])
    cb_x.set_xlim(0, n); cb_x.set_ylim(0, 1); cb_y.set_ylim(n, 0); cb_y.set_xlim(0, 1)
    for ax_ in (cb_x, cb_y):
        ax_.axis("off")
    for f in fam_order:
        k = len(fi[f]); c = fcol[f]
        cb_x.add_patch(Rectangle((start, 0), k, 1, fc=c, ec="none")); cb_y.add_patch(Rectangle((0, start), 1, k, fc=c, ec="none"))
        axm.add_patch(Rectangle((start - .5, start - .5), k, k, fill=False, ec=c, lw=.7))
        if k >= 5 and start > 0 and k / n * 50 >= 1.6:
            axm.text(start + k / 2 - .5, start - .006 * n, f, fontsize=3.9, ha="center", va="bottom", color=S.NAVY)
        start += k
    if un:
        cb_x.add_patch(Rectangle((start, 0), len(un), 1, fc=S.GREY_FAM, ec="none")); cb_y.add_patch(Rectangle((0, start), 1, len(un), fc=S.GREY_FAM, ec="none"))
        axm.axvline(start - .5, color=S.NOTE_GREY, ls=":", lw=S.LW["thin"]); axm.axhline(start - .5, color=S.NOTE_GREY, ls=":", lw=S.LW["thin"])
    axm.set_xlabel(f"{L('B', 'individuals', n, SRC_M)} {noun(run)}: {L('B', 'families', nfam, SRC_F)} families, then {L('B', 'unassigned', len(un), SRC_M)} unassigned", fontsize=sz["label"], labelpad=1.5)
    kb = axm.inset_axes([.565, .66, .43, .32]); kb.set_xlim(0, 1); kb.set_ylim(0, 1); kb.axis("off")
    kb.add_patch(Rectangle((0, 0), 1, 1, fc="white", ec=S.NAVY, lw=.3, zorder=0))
    g3 = lambda v: f"{v:.3f}".rstrip("0")
    labs = [f"θ ≥ {g3(BRK[3])}", f"1st degree {g3(BRK[2])}–{g3(BRK[3])}", f"2nd degree {g3(BRK[1])}–{g3(BRK[2])}", f"3rd degree {g3(BRK[0])}–{g3(BRK[1])}", f"below 3rd < {g3(BRK[0])}"]
    for k, (lab, col) in enumerate(zip(labs, S.BAND_COL[::-1])):
        yy = .85 - k * .17
        kb.add_patch(Rectangle((.05, yy - .057), .11, .115, fc=col, ec=S.NAVY, lw=.3)); kb.text(.22, yy, lab, fontsize=sz["small"] - 1, va="center")
    # ---------------------------------------------------------------- C family sizes
    cv.title(next(letters), "Family sizes", 4.0, 127.0 + dx)
    axc = cv.ax(135 + dx, 12, 46, 46)
    szc = size.value_counts().sort_index(); xmax = max(int(size.max()) if nfam else 2, 11); xs = np.arange(2, xmax + 1); vals = [int(szc.get(k, 0)) for k in xs]
    axc.bar(xs, vals, color=S.PAL[0], width=.78, lw=0)
    for x, v in zip(xs, vals):
        if v:
            axc.text(x, v + .014 * max(max(vals) * 1.14, max([v2 for x2, v2 in zip(xs, vals) if x2 > (2 + xmax) / 2 - .5] + [0]) / .17), str(L("C", f"families of {x}", v, SRC_F + " n")), ha="center", fontsize=sz["small"], fontweight="bold")
    axc.set_xticks(xs if len(xs) <= 14 else xs[::2]); axc.set_xlabel(f"{noun(run)} per full-sib family"); axc.set_ylabel("families"); tall_right = max([v for x, v in zip(xs, vals) if x > (2 + xmax) / 2 - .5] + [0])          # bars under the insets must stay below them
    axc.set_ylim(0, max(max(vals + [1]) * 1.14, tall_right / .17)); axc.set_xlim(1.4, xmax + .6)
    axc.axvline(MIN_FAM - .5, color=S.RED, ls="--", lw=S.LW["ref"]); axc.text(MIN_FAM - .5, axc.get_ylim()[1], f"≥ {MIN_FAM}", fontsize=sz["small"] + .3, color=S.RED, fontweight="bold", ha="left", va="top")
    axc.text(.98, .575, f"{nfam} families, {L('C', 'individuals in families', int(size.sum()), SRC_M)} {noun(run)}\n{n_un} {noun(run)} unassigned", transform=axc.transAxes, ha="right", va="center", fontsize=sz["small"], color=S.NOTE_GREY)
    ins = axc.inset_axes([.5, .71, .48, .26]); cnt = [int((P.band == k).sum()) for k in range(3)] + [int((P.band >= 3).sum())]
    ins.bar(range(4), cnt, color=S.BAND_COL[:4], ec=S.NAVY, lw=.4, width=.78); ins.set_yscale("log"); ins.set_ylim(1, max(1e5, max(cnt) * 5))
    for k, v in enumerate(cnt):
        ins.text(k, max(v, 1) * 1.3, f"{L('C inset', 'pairs in band ' + ['<3rd', '3rd', '2nd', '1st'][k], v, SRC_R + ' theta'):,}", ha="center", fontsize=sz["small"] - .4, fontweight="bold")
    ins.set_xticks(range(4)); ins.set_xticklabels(["<3rd", "3rd", "2nd", "1st"], fontsize=sz["small"]); ins.tick_params(labelsize=sz["small"], length=1.5)
    ins.set_ylabel("pairs", fontsize=sz["small"], labelpad=1); ins.set_title(f"pairs per θ band ({len(P):,})", fontsize=sz["small"] + .3, fontweight="bold", pad=2)
    if has_k and len(Wc):
        cuts = sorted({.05, .08, .10, .12, float(fcfg["k2_min"]), .18, .20})
        fclosed = [float(np.mean([bool(((g.theta >= fcfg["theta_min"]) & (g.k2 >= c_)).all()) for _, g in Wc.groupby("fa")])) for c_ in cuts]
        ik = axc.inset_axes([.62, .205, .36, .20]); ik.plot(cuts, fclosed, color="#555555", lw=.9, marker="o", ms=2.4); ik.axvline(fcfg["k2_min"], color=S.RED, ls="--", lw=.8)
        ik.set_ylim(min(.3, min(fclosed) - .05), 1.0); ik.set_xticks([.05, .10, .15, .20])
        ik.tick_params(labelsize=sz["small"] - .3, length=1.2, pad=.8); ik.set_xlabel("k$_2$ cut", fontsize=sz["small"] - .3, labelpad=.5); ik.set_ylabel("frac. closed", fontsize=sz["small"] - .3, labelpad=.8)
        ik.set_title("QC: closure vs classifier cut\n(sensitivity, not the decision)", fontsize=sz["small"] - .6, fontweight="bold", pad=3.4, loc="right")
    # ---------------------------------------------------------------- D within vs between
    cv.title(next(letters), "Within- vs between-family θ", 71.0, 1.0)
    axd = cv.ax(11, 80, 60, 38)
    Wn = P[P.within]; B = P[~P.within]; B_lo = B[B.theta < BRK[1]]; B_hi = B[B.theta >= BRK[1]]
    hi_x = max(.335, float(np.ceil(P.theta.quantile(.999) * 20) / 20) if len(P) else .335); lo_x = min(-.005, float(P.theta.min()) if len(P) else -.005)
    bins = np.linspace(lo_x, hi_x, 70)
    axd.hist(B_lo.theta, bins=bins, density=True, color="#C9D1D8", label=f"between < {g3(BRK[1])} ({L('D', 'between-family pairs below second degree', len(B_lo), SRC_R):,})", lw=0)
    if len(B_hi):
        axd.hist(B_hi.theta, bins=bins, density=True, color=SAND, label=f"between ≥ {g3(BRK[1])} ({L('D', 'between-family pairs at second degree or closer', len(B_hi), SRC_R):,})", lw=0, alpha=.95)
    if len(Wn):
        axd.hist(Wn.theta, bins=bins, density=True, color=BRICK, alpha=.92, label=f"within family ({L('D', 'within-family pairs', len(Wn), SRC_R):,})", lw=0)
    axd.set_yscale("log"); axd.set_ylim(1e-1, 8e2); axd.set_xlim(lo_x, hi_x)
    axd.axvline(BRK[2], color=S.NOTE_GREY, ls=":", lw=S.LW["thin"])
    axd.set_xlabel("kinship θ"); axd.set_ylabel("density (log)"); axd.legend(loc="upper right", fontsize=sz["small"] - .4, handlelength=.9, handletextpad=.4, labelspacing=.25, bbox_to_anchor=(1.0, 1.02))
    if len(Wn):
        axd.axvline(Wn.theta.median(), ymax=.70, color=S.RED, ls="--", lw=S.LW["ref"])
        nb = int((B.theta >= BRK[2]).sum()); nw = int((Wn.theta < BRK[2]).sum())
        _stat_note(axd, f"within median {L('D', 'within-family median theta', float(Wn.theta.median()), SRC_R):.3f} (red dashed)\nwithin < {g3(BRK[2])}: {L('D', 'within-family pairs below first degree', nw, SRC_R)} of {len(Wn)}\nbetween ≥ {g3(BRK[2])}: {L('D', 'between-family pairs at first degree', nb, SRC_R)}", sz)
        fm_ = Wn.groupby("fa").theta.agg(["median", lambda s: s.quantile(.25), lambda s: s.quantile(.75)])
        fm_.columns = ["med", "q1", "q3"]; fm_ = fm_[size.reindex(fm_.index) >= MIN_FAM].sort_values("med")
        if len(fm_) >= 2:
            id_ = axd.inset_axes([.2, .42, .3, .26]); xx = np.arange(len(fm_))
            id_.vlines(xx, fm_.q1, fm_.q3, color=S.NOTE_GREY, lw=.5); id_.scatter(xx, fm_.med, s=5, c=BRICK, zorder=3, lw=0)
            id_.axhline(.25, color=S.NAVY, ls="--", lw=.5); id_.set_xticks([]); id_.tick_params(labelsize=sz["small"], length=1.5)
            id_.set_ylabel("family median θ", fontsize=sz["small"], labelpad=1); id_.set_xlabel(f"{len(fm_)} families ≥ {MIN_FAM}", fontsize=sz["small"] - .4, labelpad=.5)
            id_.set_facecolor("white"); id_.set_ylim(min(.15, float(fm_.q1.min()) - .01), max(.34, float(fm_.q3.max()) + .01))
    # ---------------------------------------------------------------- E closure graphs
    status = dict(zip(F.family.astype(str), F.status.astype(str))) if len(F) else {}
    n_closed_ = int(sum(v == "closed" for v in status.values()))
    cv.title(next(letters), f"Full-sib closure per family ({L('E', 'closed families', n_closed_, SRC_F + ' status')} of {nfam} closed)", 71.0, 72.0)
    axg = cv.ax(72, 78, 110, eh); axg.set_xlim(0, 110); axg.set_ylim(eh, 0); axg.axis("off")
    fl = sorted(fams, key=lambda f: (-len(mem_of[f]), S.fam_key(f)))[:rows_e * ncol]; cw = 110 / ncol

    def bridged(u, v, mem):
        return any(okset.get(frozenset((u, w)), False) and okset.get(frozenset((v, w)), False) for w in mem if w not in (u, v))
    tri = {"closed": 0, "intransitive": 0, "one edge": 0}; n_br = n_unbr = 0; trio = {}
    for f_ in fams:
        mem_ = mem_of[f_]; c_ = {"closed": 0, "intransitive": 0, "one edge": 0}
        if len(mem_) <= 40:
            for a_, b_, d_ in itertools.combinations(mem_, 3):
                ne = okset.get(frozenset((a_, b_)), False) + okset.get(frozenset((a_, d_)), False) + okset.get(frozenset((b_, d_)), False)
                if ne:
                    c_["closed" if ne == 3 else "intransitive" if ne == 2 else "one edge"] += 1
        for a_, b_ in itertools.combinations(mem_, 2):
            if not okset.get(frozenset((a_, b_)), False):
                if bridged(a_, b_, mem_):
                    n_br += 1
                else:
                    n_unbr += 1
        for k_ in tri:
            tri[k_] += c_[k_]
        if len(mem_) == 3:
            trio[f_] = "closed" if c_["closed"] else "intransitive"
    for q, f in enumerate(fl):
        r_, c_ = divmod(q, ncol); cx, cy = (c_ + .5) * cw, r_ * rh + 4.0
        mem = mem_of[f]; k = len(mem); rad = 1.15 if k == 2 else min(3.2, .75 + .24 * k)
        ang = np.linspace(0, 2 * np.pi, k, endpoint=False) + (np.pi / 2 if k > 2 else 0); px = cx + rad * np.cos(ang); py = cy + rad * np.sin(ang)
        for u in range(k):
            for v in range(u + 1, k):
                o = okset.get(frozenset((mem[u], mem[v])), False); br = (not o) and bridged(mem[u], mem[v], mem)
                axg.plot([px[u], px[v]], [py[u], py[v]], color=S.NAVY if o else (S.PAL[4] if br else S.RED), lw=.22 if o else .55, ls="-" if o else (0, (1.6, 1.0)), alpha=.55 if o else 1, zorder=1 if o else 2)
        axg.scatter(px, py, s=3.2 if k <= 14 else 1.6, c=fcol[f], ec=S.NAVY, lw=.2, zorder=3)
        st = status.get(f, "")
        axg.text(cx, cy + rad + 1.0, f + ("" if st in ("closed", "imported") else "*"), fontsize=3.6, ha="center", va="top", color=S.NAVY if st in ("closed", "imported") else S.RED, fontweight="bold" if st not in ("closed", "imported") else "normal")
    rule = f"θ ≥ {g3(fcfg['theta_min'])}" + (f" and k$_2$ ≥ {fcfg['k2_min']:g}" if has_k else " (no k$_2$: GRM)")
    ntrio_c = sum(v == "closed" for v in trio.values()); intr = [k_ for k_, v in trio.items() if v == "intransitive"]
    n_ok = int(Wc.ok.sum()) if len(Wc) else 0
    txt = (f"Solid edges: {rule} ({L('E', 'within-family pairs passing the rule', n_ok, SRC_R)} of {L('E', 'within-family pairs', len(Wc), SRC_R)} within-family pairs). Dashed: below the rule; terracotta = bridged by a third full-sib confirmed to both "
           f"({L('E', 'bridged pairs', n_br, SRC_R)}), red = unbridged ({L('E', 'unbridged pairs', n_unbr, SRC_R)}). * = not closed ({nfam - n_closed_}; {n_closed_} closed). "
           f"Triads: {tri['closed']} closed, {tri['intransitive']} intransitive (2 of 3 edges), {tri['one edge']} with one edge; 3-{noun(run, 1)} families: {ntrio_c} of {len(trio)} closed"
           + (f", {', '.join(intr[:4])}{' ...' if len(intr) > 4 else ''} intransitive" if intr else "") + ". Node colour = family colour of the matrix; families ordered by size.")
    import textwrap
    cv.text(72, 78 + eh + 2.3, "\n".join(textwrap.wrap(txt, 150)), fs=sz["note"] - .3, color=S.NOTE_GREY, va="top", linespacing=1.25)
    desc = ("Port of panels A-E of the published relatedness figure. A k0 / k1 / k2 of every pair on the IBD simplex, coloured by kinship band, with expected positions (stars), observed medians (circles) and the parent-offspring "
            "region. B kinship band of every pair, individuals ordered by family (families ordered by similarity), then unassigned individuals. C families by size, with the number of pairs per kinship band and the share of "
            "families that stay closed as the k2 cut moves. D kinship of within- and between-family pairs. E one graph per family: solid = pair passes the full-sib rule, dashed = it does not. "
            "Families are relatedness-based groups, not pedigree records.")
    return fig, "Relatedness and the full-sib families", desc, "; ".join(notes)
