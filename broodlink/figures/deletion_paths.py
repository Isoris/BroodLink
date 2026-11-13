"""Deletion genotypes along linked families: ports of the published supplementary figures
"Deletion genotypes along a path of linked families" (`build_Supp_Fig_S-H6.py`, panels A-B,
drawing functions `draw_block`, `qc_cols`, `strips2` of `build_supp_chromosome_DEL_heatmaps_v2.py`)
and "Deletion genotypes, the other links" (`Supp_Fig_S-H7.py`, one block per link and the
per-chromosome summary).

The blocks (which families, which chromosome, how many deletions fit) and the sharing per
chromosome are tables of step 5 (`broodlink/deletion_paths.py`); the genotype matrices are
re-read from the deletion input with the same function, so the figure shows what the tables
count. The quality track of the originals had two bars per individual (% blank genotypes on
the chromosome, and the read support reported by the SV caller); BroodLink has no input for
the second one (sequencing depth is a different quantity), so only "% blank" is drawn.
"""
from __future__ import annotations
import os

import numpy as np

from .. import deletion_paths as DP
from . import style as S
from .shared_parent import family_sizes, noun

L = S.LED
SRC_B = "05_network/deletion_path_blocks.tsv"
SRC_L = "05_network/deletion_path_links.tsv"
SRC_S = "05_network/deletion_sharing_per_chromosome.tsv"
NA_C = "#AEB6BD"                                                   # missing genotype (grey, as in the original heat maps)
CROSS_COL = [S.ROH_BASE["sand"], S.ROH_BASE["teal"], S.ROH_BASE["navy"]]
CROSS_LAB = ["sibs all alike", "sibs split 1:1", "sibs split 1:2:1"]
BAD_C, OK_C = "#1F2A33", "#B7BEC6"
BOTH_C = "#7A4E9F"                                                 # 1/1 in both families (link figure); navy in the chain figure, as published
TIER_BLUE = ["#E3E8EC", "#A9CBE8", "#2F7FC1"]
TIER_BRICK = ["#DCE3EC", "#E3B193", S.ROH_BASE["brick"]]
MAX_SHARE_ROWS, MAX_SHARE_COLS, VALUE_COLS = 12, 60, 32


# ------------------------------------------------------------------ data
def _inputs(run):
    """(context, tables) or a reason string."""
    net = run.t.get("network", {})
    dls = [m for m in run.sets if m.klass == "DEL"]
    if not dls:
        return "no deletion genotypes in the inputs"
    if "links" not in net or not len(DP.kept_links(net["links"])):
        return "no kept link between families"
    if "deletion_blocks" not in net or not len(net["deletion_blocks"]):
        return "the tables of the deletion blocks are absent (re-run the network step)"
    _u, _f, members = run.units()
    sizes = family_sizes(run)
    ctx = DP.Context(max(dls, key=lambda m: m.n), members, list(sizes), run.cfg)
    return ctx, net


def chrom_labels(names):
    """Short labels of chromosome names: (label in titles, label on a crowded axis). The prefix
    common to all names is cut at its last separator, then the common letters are dropped."""
    names = [str(x) for x in names]
    p = os.path.commonprefix(names) if len(names) > 1 else ""
    i = max((k for k, ch in enumerate(p) if not ch.isalnum()), default=-1)
    mid = [n[i + 1:] or n for n in names]
    q = (os.path.commonprefix(mid) if len(mid) > 1 else "").rstrip("0123456789")
    return dict(zip(names, mid)), dict(zip(names, [m[len(q):] or m for m in mid]))


def _all_chroms(ctx):
    return sorted(set(ctx.chrom), key=DP.chrom_key)


def mb_ticks(pos):
    """(column index, label) of round megabase positions along the drawn deletions."""
    if not len(pos):
        return []
    mx = float(pos.max())
    step = 10e6 if mx > 40e6 else 5e6 if mx > 10e6 else next((s for s in (2e6, 1e6, .5e6, .2e6, .1e6) if mx / s >= 3), .1e6)
    out = []
    t = step
    while t < mx:
        k = int(np.searchsorted(pos, t))
        if k < len(pos):
            out.append((k, f"{t / 1e6:g}"))
        t += step
    return out


def _blocks(net, kind):
    """Blocks of one kind in table order: list of dicts (name, chrom, fams, rows of the family table, links)."""
    T = net["deletion_blocks"]; K = net.get("deletion_block_links")
    out = []
    for name in dict.fromkeys(T.block[T.kind == kind]):
        t = T[T.block == name]
        lk = K[K.block == name] if K is not None else None
        out.append({"name": str(name), "chrom": str(t.chrom.iloc[0]), "fams": t.family.astype(str).tolist(), "t": t.set_index(t.family.astype(str)),
                    "centre": next(iter(t.family[t.role == "centre"].astype(str)), ""), "definition": str(t.definition.iloc[0]),
                    "links": [] if lk is None else [(str(r.unit_a), str(r.unit_b), str(r.status), int(r.n_hom_alt_in_both_drawn)) for r in lk.itertuples()]})
    return out


def _sharing(net, pairs):
    """Per-chromosome sharing of the given links: (chromosomes, list of (a, b, status, observed, tier))."""
    T = net.get("deletion_sharing")
    if T is None or not len(T):
        return [], []
    chs = sorted(set(T.chrom.astype(str)), key=DP.chrom_key); rows = []
    for a, b in pairs:
        t = T[((T.unit_a.astype(str) == a) & (T.unit_b.astype(str) == b)) | ((T.unit_a.astype(str) == b) & (T.unit_b.astype(str) == a))]
        if not len(t):
            continue
        t = t.set_index(t.chrom.astype(str)).reindex(chs)
        tier = np.where(t.observed > t.null_p95, 2, np.where(t.observed > t.null_median, 1, 0)).astype(float)
        tier[t.observed.isna().to_numpy()] = np.nan
        rows.append((str(t.unit_a.dropna().iloc[0]), str(t.unit_b.dropna().iloc[0]), str(t.status.dropna().iloc[0]), t.observed.to_numpy(float), tier,
                     int(np.nanmax(t.n_null_pairs.to_numpy(float)))))
    return chs, rows


def _share_columns(net, chs, notes):
    """Chromosomes drawn in the sharing panel: all, or the MAX_SHARE_COLS with the most deletions."""
    if len(chs) <= MAX_SHARE_COLS:
        return chs
    T = net["deletion_sharing"]; n = T.groupby(T.chrom.astype(str)).n_deletions.max()
    keep = set(n.sort_values(ascending=False, kind="stable").index[:MAX_SHARE_COLS])
    notes.append(f"sharing panel: the {MAX_SHARE_COLS} chromosomes with the most deletions are drawn, {len(chs) - MAX_SHARE_COLS} are not")
    return [c for c in chs if c in keep]


def _share_links(net):
    """Links of the sharing panel: the supported links; when there is none, every link with a sharing score."""
    K = DP.kept_links(net["links"])
    sup = [(str(r.unit_a), str(r.unit_b)) for r in K.itertuples() if str(r.status) == "supported"]
    if sup:
        return sup, "supported"
    return [(str(r.unit_a), str(r.unit_b)) for r in K.itertuples()], "kept"


# ------------------------------------------------------------------ drawing helpers
def _clean(ax):
    for s in ax.spines.values():
        s.set_visible(False)


def _heat(ax, M):
    from matplotlib.colors import ListedColormap
    _clean(ax); ax.set_xticks([])
    if not M.size:
        return
    ax.imshow(np.where(M < 0, 3, M), cmap=ListedColormap([S.REF_C, S.HET_C, S.HOM_C, NA_C]), vmin=-.5, vmax=3.5, aspect="auto", interpolation="nearest")
    _clean(ax); ax.set_xticks([])


def _strip(ax, v, cols, vmax):
    from matplotlib.colors import ListedColormap
    _clean(ax); ax.set_xticks([])
    if not len(v):
        return
    ax.imshow(np.asarray(v, float)[None, :], cmap=ListedColormap(cols), vmin=-.5 if vmax > 1 else 0, vmax=vmax + .5 if vmax > 1 else 1, aspect="auto", interpolation="nearest")
    _clean(ax); ax.set_xticks([])


def _blank_bar(cv, x, y, w, h, miss, low, lw=.4):
    a = cv.ax(x, y, w, h); n = len(miss)
    a.barh(np.arange(n), miss * 100, height=.75, color=[BAD_C if b else OK_C for b in miss > low], lw=0)
    a.set_xlim(0, 100); a.set_ylim(n - .5, -.5); a.set_yticks([]); a.set_xticks([]); a.axvline(low * 100, color=S.NAVY, lw=lw, ls=(0, (2, 2)))
    for sp in ("top", "right", "bottom"):
        a.spines[sp].set_visible(False)
    a.spines["left"].set_linewidth(lw)
    return a


def _names(fams):
    return fams[0] if len(fams) == 1 else ", ".join(fams[:-1]) + " and " + fams[-1]


def _check(B, blk, notes):
    """The block redrawn from the genotypes must be the block of the table."""
    t = blk["t"]
    ok = int(t.n_eligible.iloc[0]) == B["n_eligible"] and int(t.n_drawn.iloc[0]) == len(B["rows"]) and all(int(t.loc[g["fam"], "n_individuals"]) == g["M"].shape[0] for g in B["groups"])
    if not ok:
        notes.append(f"{blk['name']}: the genotypes read for the figure give other counts than {SRC_B} (inputs differ from those of the network step); the figure shows the genotypes read now")
    return ok


def _share_panel(cv, top, left, width, rh, chs, rows, short, cols, panel, sz, values, label_fmt, grid_lw):
    """Links x chromosomes squares in three tiers; returns (cells above the median, above the 95th percentile, cells)."""
    from matplotlib.colors import ListedColormap
    st = np.array([r[4] for r in rows]); ob = np.array([r[3] for r in rows]); nr, nc = st.shape
    ax = cv.ax(left, top, width, nr * rh)
    cm = ListedColormap(cols).with_extremes(bad="white")
    ax.imshow(np.ma.masked_invalid(st), cmap=cm, vmin=-.5, vmax=2.5, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(-.5, nc, 1), minor=True); ax.set_yticks(np.arange(-.5, nr, 1), minor=True); ax.grid(which="minor", color="white", lw=grid_lw); ax.tick_params(which="minor", length=0)
    cnt = []
    for i, r in enumerate(rows):
        cnt.append((int(np.nansum(st[i] >= 1)), int(np.nansum(st[i] == 2))))
        L(panel, f"{r[0]} x {r[1]} chromosomes above null median", cnt[-1][0], SRC_S + " observed > null_median")
        L(panel, f"{r[0]} x {r[1]} chromosomes above null p95", cnt[-1][1], SRC_S + " observed > null_p95")
        for j, c in enumerate(chs):
            if values and np.isfinite(ob[i, j]):
                t = f"{ob[i, j]:.2f}"; t = t[1:] if t.startswith("0") else t
                ax.text(j, i, t, fontsize=sz["small"], ha="center", va="center", color="white" if st[i, j] == 2 else S.NAVY)
                L(panel, f"{r[0]} x {r[1]} {c} observed share", f"{ob[i, j]:.2f}", SRC_S + " observed")
    ax.set_xticks(range(nc)); ax.set_xticklabels([short[c] for c in chs], fontsize=sz["tick"] if nc <= 40 else sz["small"], rotation=0 if max(len(short[c]) for c in chs) <= 3 else 90)
    ax.set_xlabel("chromosome", fontsize=sz["tick"], labelpad=1.5); ax.tick_params(axis="x", length=1.2, pad=1)
    ax.set_yticks(range(nr)); ax.set_yticklabels([label_fmt(r, c) for r, c in zip(rows, cnt)], fontsize=sz["small"]); ax.tick_params(axis="y", length=0, pad=1.5)
    _clean(ax)
    return sum(c[0] for c in cnt), sum(c[1] for c in cnt), int(np.isfinite(st).sum())


# ------------------------------------------------------------------ figure: chains
def _chain_figure(run, plt, ctx, net, blk, share, notes):
    p = ctx.p; W = S.WIDTH_FULL; nn = noun(run); n1 = noun(run, 1)
    B = ctx.block(blk["chrom"], blk["fams"])
    if not len(B["rows"]):
        return None
    _check(B, blk, notes)
    mid, short = chrom_labels(_all_chroms(ctx)); cl = mid[blk["chrom"]]
    fams = B["fams"]; nfish = sum(g["M"].shape[0] for g in B["groups"]); nr, ne = len(B["rows"]), B["n_eligible"]
    links = blk["links"]; chs, srows, scope = share
    STRIP_H, LEFT, HW = 3.2, 24.0, 145.0                            # heat map 9 mm wider than published: one quality bar instead of two
    nb = len(srows)
    fixA = 13.5 + len(fams) * (STRIP_H + .5 + 1.8) + 3.4 * len(links) + 9
    fixB = (13 + nb * 4.6 + 9 + 4) if nb else 0
    rh = min(2.5, (S.A4_H - fixA - fixB) / max(nfish, 1))
    if rh < 2.5:
        notes.append(f"{blk['name']}: rows of {rh:.1f} mm instead of 2.5 mm to keep the page within A4")
    H = fixA + rh * nfish + fixB + (2 if not nb else 0)
    cv = S.canvas(plt, W, H, S.SZ_ED); sz = cv.sz; fig = cv.fig
    who = _names(fams) if len(fams) <= 3 else f"{len(fams)} linked families"
    cv.banner(0, "1", f"Across {who} ({L('1', 'individuals drawn', nfish, SRC_B + ' n_individuals')} {nn}), {L('1', 'DELs drawn', nr, SRC_B + ' n_drawn')} DELs on {cl} that fit a Mendelian cross show the split expected between siblings")
    L("A", "DELs eligible", ne, SRC_B + " n_eligible")
    for g in B["groups"]:
        c = g["cross"]
        L("A", f"{g['fam']} expected cross 1:1 / 1:2:1 / all alike", f"{int((c == 1).sum())}/{int((c == 2).sum())}/{int((c == 0).sum())}", SRC_B + " n_split_1_1, n_split_1_2_1, n_all_alike")
        L("A", f"{g['fam']} {nn} drawn", int(g["M"].shape[0]), SRC_B + " n_individuals"); L("A", f"{g['fam']} {nn} > {p['low_call']:.0%} blank on {cl}", int(g["n_low"]), SRC_B + " n_low_call")
    ttl = f"Linked families {' | '.join(fams)}, DELs along {cl}"
    cv.title("A", ttl, 9.5, .5, dx=4)
    x = max(82.0, 4.5 + cv.text_width(ttl, sz["title"], True) + 4)
    key = [(S.REF_C, "0/0"), (S.HET_C, "0/1"), (S.HOM_C, "1/1"), (NA_C, "missing")] + list(zip(CROSS_COL, CROSS_LAB))
    need = sum(2.6 + .7 + cv.text_width(lab, sz["small"]) + 2.5 for _c, lab in key)
    ky = 9.5
    if x + need > W - 1:                                            # long title: key on its own line is not possible in the fixed layout, so it is packed tighter from the right
        x = max(W - 1 - need, 30.0)
    for c, lab in key:
        cv.swatch(x, ky, c, lab, ec=S.BACKBONE); x += 2.6 + .7 + cv.text_width(lab, sz["small"]) + 2.5
    y = 13.5; last_axes = None
    for g in B["groups"]:
        sa = cv.ax(LEFT, y, HW, STRIP_H); _strip(sa, g["cross"], CROSS_COL, 2)
        sa.set_yticks([0]); sa.set_yticklabels(["expected\nin sibs"], fontsize=sz["small"]); sa.tick_params(axis="y", length=0, pad=2)
        y += STRIP_H + .5; n = g["M"].shape[0]; h = rh * n
        ax = cv.ax(LEFT, y, HW, h); _heat(ax, g["M"])
        ax.set_yticks([(n - 1) / 2]); ax.set_yticklabels([f"{g['fam']} (n = {n})"], fontsize=sz["small"]); ax.tick_params(axis="y", length=0, pad=2)
        qa = _blank_bar(cv, LEFT + HW + 1.5, y, 7.5, h, g["miss"], p["low_call"]); qa.set_title("% blank", fontsize=sz["note"], pad=1, loc="left")
        if g["n_low"] and g["n_good"]:
            ax.axhline(g["n_good"] - .5, color=S.NAVY, lw=.5, ls=(0, (2, 2)))
        if g["n_low"] and rh >= 2.0 and (g["n_good"] >= 2 or n == g["n_low"] > 1):
            ax.text(-.004, 1 - (g["n_good"] + (g["n_low"] - 1) / 2 + .5) / n, f"> {p['low_call']:.0%} missing", fontsize=sz["note"], va="center", ha="right", color=S.NOTE_GREY, transform=ax.transAxes)
        y += h + 1.8
    gm = {g["fam"]: g for g in B["groups"]}
    for k, (a, b, status, n_tab) in enumerate(links):
        st = DP.both_hom(gm[a]["M"], gm[b]["M"])
        ax = cv.ax(LEFT, y + k * 3.4, HW, 2.6); _strip(ax, st, ["white", S.ROH_BASE["navy"]], 1)
        ax.set_yticks([0]); ax.set_yticklabels([f"{a} x {b}"], fontsize=sz["small"]); ax.tick_params(axis="y", length=0, pad=2); last_axes = ax
        L("A", f"{a} x {b} DELs 1/1 in both families, of those drawn", int(st.sum()), SRC_L + " n_hom_alt_in_both_drawn")
    if last_axes is not None:
        tk = mb_ticks(B["pos"]); last_axes.set_xticks([q for q, _ in tk]); last_axes.set_xticklabels([t for _, t in tk], fontsize=sz["tick"]); last_axes.tick_params(axis="x", length=1.5, pad=1)
        last_axes.set_xlabel(f"{nr} DELs in position order (every {L('A', 'step of the subsample', round(ne / nr, 1), SRC_B + ' n_eligible / n_drawn'):.1f}th of {ne} that fit); ticks = Mb; "
                             f"navy strip = DEL 1/1 in at least one {n1} of both families; right: % blank genotypes per {n1}", fontsize=sz["small"], labelpad=1.5)
    t0 = y + 3.4 * len(links) + 9
    if nb:
        ttop = t0 + 13
        k_null = max(r[5] for r in srows)
        tot1, tot2, cells = _share_panel(cv, ttop, 39, 143, 4.6, chs, srows, short, TIER_BLUE, "B", sz, len(chs) <= VALUE_COLS,
                                         lambda r, c: f"{r[0]} x {r[1]}  {c[0]} | {c[1]}", .8)
        cv.banner(t0, "2", f"Linked families share more DELs than matched unlinked pairs: {L('2', 'squares above the null median', tot1, SRC_S)} of {L('2', 'squares', cells, SRC_S)} squares above the median, "
                           f"{L('2', 'squares above the null p95', tot2, SRC_S)} above the 95th percentile")
        cv.title("B", "Sharing score per link and chromosome (all DELs of the chromosome)" + ("" if scope == "supported" else "; no supported link: all kept links with a score"), t0 + 9.5, .5, dx=4)
        ky = ttop + nb * 4.6 + 9
        x = 39.0
        for c_, lab in ((TIER_BLUE[0], f"below the median of {L('B', 'matched unlinked pairs per link', k_null, SRC_S + ' n_null_pairs')} matched unlinked pairs"), (TIER_BLUE[1], "above that median"), (TIER_BLUE[2], "above their 95th percentile")):
            cv.swatch(x, ky, c_, lab, ec=S.BACKBONE); x += 2.6 + .7 + cv.text_width(lab, sz["small"]) + 8
        cv.text(1, ky, "row label: above median | above p95", fs=sz["small"], color=S.INK)
    return fig


def fig_chains(run, plt):
    """One figure per chain of linked families (port of the published supplementary figure):
    A genotypes of the deletions that fit a Mendelian cross in every family of the chain,
    B sharing score per link and chromosome against the matched non-nominated pairs."""
    inp = _inputs(run)
    if isinstance(inp, str):
        return inp
    ctx, net = inp
    chains = _blocks(net, "chain")
    dsum = run.s.get("network", {}).get("deletion_paths", {}) or {}
    if not chains:
        return "no chain: no family is on two or more kept links"
    cap = int(run.cfg["figures"].get("max_chain_figures", 4))
    n_all = int(dsum.get("n_chains", len(chains)))
    base = []
    if n_all > min(cap, len(chains)):
        base.append(f"{min(cap, len(chains))} of {n_all} chains drawn (figures.max_chain_figures); the links of the others are in the figure of the other links")
    pairs, scope = _share_links(net)
    chs, srows = _sharing(net, pairs)
    if len(srows) > MAX_SHARE_ROWS:
        base.append(f"sharing panel: the first {MAX_SHARE_ROWS} of {len(srows)} links are drawn"); srows = srows[:MAX_SHARE_ROWS]
    if srows:
        keep = _share_columns(net, chs, base); ix = [chs.index(c) for c in keep]
        srows = [(r[0], r[1], r[2], r[3][ix], r[4][ix], r[5]) for r in srows]; chs = keep
        if len(chs) > VALUE_COLS:
            base.append(f"sharing panel: the score is not printed in the squares (more than {VALUE_COLS} chromosomes)")
    else:
        base.append("sharing panel omitted: no link with a Gate-1 score and matched non-nominated pairs (no audit, or no auditable pair)")
    base.append("quality track: % blank genotypes only; the read-support bar of the published figure has no input in BroodLink (sequencing depth is a different quantity)")
    rule = dsum.get("chromosome_rule", "the chromosome with the most deletions")
    out = []
    for blk in chains[:cap]:
        notes = [f"chromosome shown: {blk['chrom']} ({rule})", f"chain: {blk['definition']}; families ordered with the centre ({blk['centre']}) in the middle and the partners by family size",
                 "links: " + ", ".join(f"{a} x {b} ({st})" for a, b, st, _n in blk["links"])]
        fig = _chain_figure(run, plt, ctx, net, blk, (chs, srows, scope), notes)
        if fig is None:
            base.append(f"{blk['name']} ({' | '.join(blk['fams'])}) not drawn: no deletion fits a Mendelian cross in all its families"); continue
        desc = ("Port of the published supplementary figure. A: deletions on one chromosome that fit a Mendelian cross in every family of the chain (even subsample in position order); strip above each family = "
                "cross expected among the siblings (all alike, split 1:1, split 1:2:1); rows = individuals in hierarchical order, those with too many blank genotypes below the dotted line; "
                "strips below = deletions 1/1 in at least one individual of both families of a link. B: sharing score per link and chromosome against the median and the 95th percentile of the matched non-nominated pairs.")
        out.append((f"{blk['name']}_{blk['centre']}", fig, f"Deletion genotypes along linked families: {' | '.join(blk['fams'])}", desc, notes))
    if not out:
        return "; ".join(base) or None
    return [(sfx, fig, title, desc, "; ".join(notes + base)) for sfx, fig, title, desc, notes in out]


# ------------------------------------------------------------------ figure: the other links
RH, STRIP, LEFT7, HW7, QX7 = 1.25, 1.5, 20.0, 151.3, 173.3     # heat map 9.3 mm wider than published: one quality bar instead of two


def _link_height(n_a, n_b):
    return 3.3 + (STRIP + .3 + RH * n_a + 1.0) + (STRIP + .3 + RH * n_b + 1.0) - .2 + 1.7 + 6.0


def _links_page(run, plt, ctx, net, blks, page_notes):
    p = ctx.p; W = S.WIDTH_FULL; n1 = noun(run, 1); nn = noun(run)
    mn = int(run.cfg["audit"]["min_family_size"])
    Bs = [(b, ctx.block(b["chrom"], b["fams"])) for b in blks]
    mid, short = chrom_labels(_all_chroms(ctx)); cl = mid[blks[0]["chrom"]]
    chs, srows = _sharing(net, [(b["links"][0][0], b["links"][0][1]) for b in blks if b["links"]])
    if srows:
        keep = _share_columns(net, chs, page_notes); ix = [chs.index(c) for c in keep]
        srows = [(r[0], r[1], r[2], r[3][ix], r[4][ix], r[5]) for r in srows]; chs = keep
    ns = len(srows)
    SUMH = (7 + 3.4 + ns * 2.5 + 6.2 + 4.5) if ns else 1.0
    H = 11.5 + sum(_link_height(*(g["M"].shape[0] for g in B["groups"])) for _b, B in Bs) + SUMH
    cv = S.canvas(plt, W, H, S.SZ_ED, thin=True); sz = cv.sz; fig = cv.fig
    from matplotlib.patches import Rectangle

    def key(x, col, lab, adv):
        fig.add_artist(Rectangle((x / W, 1 - (8.0 + 1.9) / H), 2.4 / W, 1.9 / H, transform=fig.transFigure, fc=col, ec=S.BACKBONE, lw=.2)); cv.text(x + 3, 8.95, lab, fs=5.5)
        return x + (adv if adv else 3 + cv.text_width(lab, 5.5) + 3)
    x = LEFT7
    for c, lab in [(S.REF_C, "0/0"), (S.HET_C, "0/1"), (S.HOM_C, "1/1"), (NA_C, "missing")]:
        x = key(x, c, lab, 4.2 + 1.15 * len(lab) + 1.6)
    x += 2
    for c, lab in zip(CROSS_COL, CROSS_LAB):
        x = key(x, c, lab, 4.2 + 1.0 * len(lab) + 2)
    x = key(x, BOTH_C, "DEL 1/1 in both families", 28)
    key(x, BAD_C, f"> {p['low_call']:.0%} blank", 0)
    y = 11.5; both, drawn = [], []
    letters = [chr(65 + k) if k < 26 else "A" + chr(65 + k - 26) for k in range(len(Bs) + 1)]
    for k, (blk, B) in enumerate(Bs):
        _check(B, blk, page_notes)
        a, b, status, _n = blk["links"][0]; let = letters[k]
        small = [g for g in B["groups"] if g["M"].shape[0] < mn]
        head = f"{a} x {b}: {status}" + (" (" + ", ".join(f"{g['fam']} has {g['M'].shape[0]} {nn}" for g in small) + ")" if small else "")
        cv.text(.8, y + .05, let, fs=8, bold=True, va="top"); cv.text(4.4, y + .25, head, fs=6.5, bold=True, va="top")
        y += 3.3
        for j, g in enumerate(B["groups"]):
            sa = cv.ax(LEFT7, y, HW7, STRIP); _strip(sa, g["cross"], CROSS_COL, 2)
            sa.set_yticks([0]); sa.set_yticklabels(["expected"], fontsize=5.5); sa.tick_params(axis="y", length=0)
            y += STRIP + .3; nf = g["M"].shape[0]; h = RH * nf
            ax = cv.ax(LEFT7, y, HW7, h); _heat(ax, g["M"])
            ax.set_yticks([(nf - 1) / 2]); ax.set_yticklabels([f"{g['fam']} (n = {nf})"], fontsize=5.5); ax.tick_params(axis="y", length=0)
            if g["n_low"] and g["n_good"]:
                ax.axhline(g["n_good"] - .5, color=S.NAVY, lw=.35, ls=(0, (2, 2)))
            _blank_bar(cv, QX7, y, 8, h, g["miss"], p["low_call"], lw=.3)
            if j == 0 and k == 0:
                cv.text(QX7, y - .75, "% blank", fs=5.5)
            L(let, f"{g['fam']} {nn} drawn", nf, SRC_B + " n_individuals"); L(let, f"{g['fam']} {nn} > {p['low_call']:.0%} blank on {cl}", int(g["n_low"]), SRC_B + " n_low_call")
            y += h + 1.0
        y -= .2
        st = DP.both_hom(B["groups"][0]["M"], B["groups"][1]["M"])
        ob = cv.ax(LEFT7, y, HW7, 1.7); _strip(ob, st, ["white", BOTH_C], 1)
        ob.set_yticks([0]); ob.set_yticklabels(["1/1 in both"], fontsize=5.5); ob.tick_params(axis="y", length=0)
        tk = mb_ticks(B["pos"]); ob.set_xticks([q for q, _ in tk]); ob.set_xticklabels([t for _, t in tk], fontsize=5.5); ob.tick_params(axis="x", pad=.6, length=1.2)
        nr, ne = len(B["rows"]), B["n_eligible"]
        L(let, "DELs drawn", nr, SRC_B + " n_drawn"); L(let, "DELs eligible", ne, SRC_B + " n_eligible"); L(let, "step of subsample", round(ne / max(nr, 1), 1), SRC_B + " n_eligible / n_drawn")
        L(let, "DELs 1/1 in both families, of those drawn", int(st.sum()), SRC_L + " n_hom_alt_in_both_drawn")
        ob.set_xlabel((f"{cl} (Mb), {nr} DELs, every {ne / nr:.1f}th of the {ne} that fit; {int(st.sum())} are 1/1 in both families" if nr else f"{cl}: no DEL fits a Mendelian cross in both families"), fontsize=5.5, labelpad=.5)
        y += 1.7 + 6.0
        both.append(int(st.sum())); drawn.append(nr)
    rng = f"{min(both)}-{max(both)}" if min(both) != max(both) else f"{both[0]}"
    per = f"the {drawn[0]} per link" if len(set(drawn)) == 1 else "those drawn per link"
    cv.banner(.3, "1", f"On {cl} every DEL drawn fits a Mendelian cross in both families of its link; {rng} of {per} are 1/1 in both" if len(Bs) > 1 else
              f"On {cl} every DEL drawn fits a Mendelian cross in both families of the link; {rng} of the {drawn[0]} are 1/1 in both")
    if ns:
        t0 = y - .3; t1 = t0 + 7.0; let = letters[len(Bs)]
        tot_med, tot_95, cells = _share_panel(cv, t1 + 3.4, 52, 130, 2.5, chs, srows, short, TIER_BRICK, let, sz, False,
                                              lambda r, c: f"{r[0]} x {r[1]}: {c[0]} and {c[1]} of {int(np.isfinite(r[4]).sum())}", .5)
        L(let, "cells above null median", tot_med, SRC_S); L(let, "cells above null p95", tot_95, SRC_S); L(let, "cells", cells, SRC_S)
        lk = "This link exceeds" if ns == 1 else f"These {ns} links exceed"
        cv.banner(t0, "2", f"{lk} the median of matched pairs in {tot_med} of {cells} chromosome cells, their 95th percentile in {tot_95}")
        cv.text(.8, t1, let, fs=8, bold=True, va="top"); cv.text(4.4, t1 + .2, "DEL sharing against matched non-nominated pairs, per chromosome", fs=6.5, bold=True, va="top")
        x = 52.0; yk = t1 + 3.4 + ns * 2.5 + 6.2
        for c, lab in zip(TIER_BRICK, ["below the null median", "above the null median", "above the null 95th percentile"]):
            fig.add_artist(Rectangle((x / W, 1 - (yk + 1.9) / H), 2.4 / W, 1.9 / H, transform=fig.transFigure, fc=c, ec=S.BACKBONE, lw=.2)); cv.text(x + 3, yk + .95, lab, fs=5.5); x += 40
    else:
        page_notes.append("per-chromosome summary panel omitted: none of these links has a sharing score with matched non-nominated pairs")
    return fig


def fig_other_links(run, plt):
    """The kept links that are not part of a drawn chain, one block per link (port of the
    published supplementary figure), and the per-chromosome summary of their sharing."""
    inp = _inputs(run)
    if isinstance(inp, str):
        return inp
    ctx, net = inp
    cap = int(run.cfg["figures"].get("max_chain_figures", 4))
    in_chain = set()
    for blk in _blocks(net, "chain")[:cap]:
        if len(ctx.block(blk["chrom"], blk["fams"])["rows"]):
            in_chain |= {frozenset(l[:2]) for l in blk["links"]}
    blks = [b for b in _blocks(net, "link") if b["links"] and frozenset(b["links"][0][:2]) not in in_chain]
    if not blks:
        return "every kept link is part of a drawn chain"
    sizes = family_sizes(run)
    fixed = 11.5 + 7 + 3.4 + 6.2 + 4.5                              # top, and the summary panel without its rows
    pages, cur, h = [], [], fixed
    for b in blks:
        hb = _link_height(*(sizes.get(f, 0) for f in b["fams"])) + 2.5
        if cur and h + hb > S.A4_H:
            pages.append(cur); cur, h = [], fixed
        cur.append(b); h += hb
    if cur:
        pages.append(cur)
    base = []
    if len(pages) > cap:
        n_left = sum(len(pg) for pg in pages[cap:])
        base.append(f"{sum(len(pg) for pg in pages[:cap])} of {len(blks)} links drawn on {cap} page{"s" if cap != 1 else ""} (figures.max_chain_figures), best status first; {n_left} links not drawn")
        pages = pages[:cap]
    base.append("quality track: % blank genotypes only; the read-support bar of the published figure has no input in BroodLink (sequencing depth is a different quantity)")
    dsum = run.s.get("network", {}).get("deletion_paths", {}) or {}
    base.append(f"chromosome shown: {blks[0]['chrom']} ({dsum.get('chromosome_rule', 'the chromosome with the most deletions')})")
    desc = ("Port of the published supplementary figure. One block per kept link that is not part of a drawn chain: deletions on one chromosome that fit a Mendelian cross in both families (even subsample), "
            "expected cross among the siblings above each family, individuals with too many blank genotypes below the dotted line, deletions 1/1 in at least one individual of both families below. "
            "Last panel: chromosomes on which the sharing score of each link exceeds the median and the 95th percentile of its matched non-nominated pairs.")
    out = []
    for i, pg in enumerate(pages):
        notes = []
        fig = _links_page(run, plt, ctx, net, pg, notes)
        sfx = "" if len(pages) == 1 else f"page{i + 1}"
        out.append((sfx, fig, "Deletion genotypes, the other links" + ("" if len(pages) == 1 else f" ({i + 1} of {len(pages)})"), desc, "; ".join(notes + base)))
    return out
