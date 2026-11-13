"""Figures of a run (PNG + PDF), drawn from the tables of the step folders only.

Colours: Okabe-Ito palette (colour-blind safe); every colour used is keyed in a legend.
Font: Arial when installed, DejaVu Sans otherwise.
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd

from .io import log, write_tsv

BLUE, ORANGE, GREEN, RED, PURPLE, SKY, YELLOW, GREY, INK = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999", "#222222"
STATUS_COL = {"supported": BLUE, "Gate-1 only": SKY, "provisional": ORANGE, "rejected": RED, "nominated, not audited": GREY}
GENO_COL = {0: BLUE, 1: GREY, 2: ORANGE, -1: "#DDDDDD"}
GENO_LAB = {0: "class 0 (more frequent homokaryotype)", 1: "class 1 (middle cluster)", 2: "class 2", -1: "ambiguous (missing)"}


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams.update({"font.family": "Arial" if "Arial" in have else "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 7,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    return plt


def _fkey(f):
    d = "".join(ch for ch in str(f) if ch.isdigit())
    return (0 if not str(f).startswith("S_") else 1, int(d) if d else 0, str(f))


def _status(s):
    for k in STATUS_COL:
        if str(s).startswith(k):
            return k
    return "nominated, not audited"


def _handles(plt, items, marker="s"):
    from matplotlib.lines import Line2D
    return [Line2D([], [], marker=marker, ls="", color=c, markersize=5, label=l) for l, c in items]


def _panel(ax, letter):
    ax.text(-0.09, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")


# ------------------------------------------------------------------ panels (each draws on one axis)
def p_callrate(run, plt, ax):
    PS = run.t["qc"]["per_sample"]; X = run.t["qc"]["flags"]
    low = set(X[X.reason == "low coverage"]["sample"]); idn = set(X[X.reason.str.startswith("pair of samples")]["sample"])
    inputs = list(dict.fromkeys(PS.input))
    for k, name in enumerate(inputs):
        g = PS[PS.input == name]
        x = k + (np.random.default_rng(k).random(len(g)) - 0.5) * 0.5
        col = np.where(g["sample"].isin(low), RED, np.where(g["sample"].isin(idn), ORANGE, GREY))
        ax.scatter(x, g.call_rate, s=8, c=col, lw=0, alpha=0.9)
        thr = run.cfg["qc"]["low_coverage_factor"] * g.call_rate.median()
        ax.plot([k - 0.35, k + 0.35], [thr, thr], color=INK, lw=0.8, ls="--")
    ax.set_xticks(range(len(inputs))); ax.set_xticklabels([n.split(" (")[0] for n in inputs])
    ax.set_ylabel("call rate"); ax.set_title("Call rate per sample and input")
    ax.legend(handles=_handles(plt, [("kept, no flag", GREY), ("low coverage (flag)", RED), ("pair with identical genotypes (flag)", ORANGE)], "o") +
              [plt.Line2D([], [], color=INK, ls="--", lw=0.8, label="low-coverage rule")], loc="lower right")


def p_het(run, plt, ax):
    PS = run.t["qc"]["per_sample"]
    cols = [BLUE, ORANGE, GREEN, PURPLE, SKY, RED]
    for k, (name, g) in enumerate(PS.groupby("input", sort=False)):
        ax.scatter(g.call_rate, g.heterozygosity, s=8, color=cols[k % len(cols)], lw=0, label=name.split(" (")[0], alpha=0.8)
    ax.set_xlabel("call rate"); ax.set_ylabel("heterozygosity (share of called genotypes)"); ax.set_title("Heterozygosity against call rate"); ax.legend()


def p_depth(run, plt, ax):
    PS = run.t["qc"]["per_sample"]
    if "depth_mask" not in PS or PS.depth_mask.notna().sum() == 0:
        ax.text(0.5, 0.5, "no --depth input", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    d = PS.drop_duplicates("sample")
    ax.scatter(d.depth_assembly, d.depth_mask, s=8, color=BLUE, lw=0, label="one sample")
    ax.set_xlabel("depth over the assembly (x)"); ax.set_ylabel("depth over the callable mask (x)"); ax.set_title("Depth: two denominators"); ax.legend()


def p_perchrom(run, plt, ax):
    PC = run.t["qc"]["per_chrom"]
    cols = [BLUE, ORANGE, GREEN, PURPLE, SKY, RED]
    chroms = list(dict.fromkeys(PC.chrom)); inputs = list(dict.fromkeys(PC.input)); w = 0.8 / max(len(inputs), 1)
    y = "markers_per_Mb" if PC.markers_per_Mb.notna().any() else "n_markers"
    for k, name in enumerate(inputs):
        g = PC[PC.input == name].set_index("chrom").reindex(chroms)
        ax.bar(np.arange(len(chroms)) + k * w, g[y] / g[y].max(), w * 0.9, color=cols[k % len(cols)], label=name.split(" (")[0])
    ax.set_xticks(np.arange(len(chroms)) + 0.4 - w / 2); ax.set_xticklabels(chroms, rotation=90 if len(chroms) > 10 else 0, fontsize=6)
    ax.set_ylabel(("markers per Mb" if y == "markers_per_Mb" else "markers") + "\n(relative to the densest chromosome)"); ax.set_title("Marker density per chromosome"); ax.legend(ncols=3)


def p_theta_k2(run, plt, ax):
    R = run.t["relatedness"]["pairs"]; cfg = run.cfg; r = cfg["relatedness"]
    order = [("unrelated", "#CCCCCC"), ("third degree", YELLOW), ("second degree", ORANGE), ("first degree, not full-sib (k2 below gate)", PURPLE),
             ("first degree (no k2)", BLUE), ("full-sib", BLUE), ("identical genotypes (theta at or above the identity threshold)", RED)]
    if R.k2.notna().any():
        for name, col in order:
            g = R[R["class"] == name]
            if len(g):
                ax.scatter(g.theta, g.k2, s=5, color=col, lw=0, label=f"{name} ({len(g):,})", alpha=0.8, rasterized=True)
        ax.axhline(cfg["families"]["k2_min"], color=INK, lw=0.7, ls="--"); ax.set_ylabel("k2")
        ax.text(ax.get_xlim()[1], cfg["families"]["k2_min"], f" k2 = {cfg['families']['k2_min']}", va="bottom", ha="right", fontsize=6)
    else:
        bins = np.linspace(min(R.theta.min(), -0.05), max(R.theta.max(), 0.3), 80)
        for name, col in order:
            g = R[R["class"] == name]
            if len(g):
                ax.hist(g.theta, bins=bins, color=col, label=f"{name} ({len(g):,})")
        ax.set_yscale("log"); ax.set_ylabel("pairs (log scale)")
    for v, lab in ((r["theta_second"], "second degree"), (r["theta_first"], "first degree"), (r["theta_duplicate"], "identical")):
        ax.axvline(v, color=INK, lw=0.7, ls=":")
        ax.text(v, ax.get_ylim()[1], f"{lab}\n{v}", fontsize=6, ha="center", va="bottom")
    ax.set_xlabel("theta (kinship)"); ax.legend(loc="center right" if R.k2.notna().any() else "upper right", fontsize=6)


def p_family_sizes(run, plt, ax):
    F = run.t["families"]["families"].sort_values("n", ascending=False)
    if not len(F):
        ax.text(0.5, 0.5, "no family", ha="center", va="center", transform=ax.transAxes); return
    closed = F.status == "closed"
    ax.bar(np.arange(len(F)), F.n, color=np.where(closed, BLUE, ORANGE), width=0.8)
    mn = run.cfg["audit"]["min_family_size"]
    ax.axhline(mn - 0.5, color=INK, ls="--", lw=0.7)
    ax.set_xlabel("families, by size"); ax.set_ylabel("individuals"); ax.set_xticks([])
    nun = int((run.t["families"]["membership"].family == "unassigned").sum())
    ax.set_title(f"{len(F)} families holding {int(F.n.sum())} individuals; {nun} unassigned")
    ax.legend(handles=_handles(plt, [("closed (every pair passes)", BLUE), ("not closed (a pair fails a threshold)", ORANGE)]) +
              [plt.Line2D([], [], color=INK, ls="--", lw=0.7, label=f"audit needs at least {mn}")])


def p_closure(run, plt, ax):
    F = run.t["families"]["families"]
    if "closure" not in F or F.closure.notna().sum() == 0:
        ax.text(0.5, 0.5, "closure not available (imported families)", ha="center", va="center", transform=ax.transAxes); return
    ax.scatter(F.n, F.closure, s=14, color=np.where(F.status == "closed", BLUE, ORANGE), lw=0)
    ax.set_xlabel("family size"); ax.set_ylabel("closure (share of pairs passing)"); ax.set_ylim(0, 1.05)
    ax.set_title("Families are connected blocks, not complete cliques")
    ax.legend(handles=_handles(plt, [("closed", BLUE), ("not closed", ORANGE)], "o"), loc="lower right")


def p_matrix(run, plt, ax, families_only=None):
    from matplotlib.colors import LinearSegmentedColormap
    E = run.t["audit"]["edges"] if "edges" in run.t.get("audit", {}) else run.t["nomination"]["pairs"]
    units = sorted(set(E.unit_a) | set(E.unit_b), key=_fkey)
    if families_only is None:
        families_only = len(units) > 120
    if families_only:
        units = [u for u in units if not u.startswith("S_")]
    ix = {u: i for i, u in enumerate(units)}; n = len(units)
    M = np.full((n, n), np.nan)
    cmap = LinearSegmentedColormap.from_list("seq", ["#F2F2F2", "#9CC3DC", "#0072B2", "#08306B"])
    up = []
    for r in E.itertuples():
        if r.unit_a in ix and r.unit_b in ix:
            i, j = sorted((ix[r.unit_a], ix[r.unit_b]))
            M[j, i] = r.fraction
            if r.nominated:
                up.append((i, j, _status(getattr(r, "link_class", "nominated, not audited"))))
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    for i, j, st in up:
        ax.scatter(j, i, s=max(4, 900 / max(n, 1) ** 1.3), marker="s", color=STATUS_COL[st], lw=0)
    ax.plot([-0.5, n - 0.5], [-0.5, n - 0.5], color=INK, lw=0.5)
    nf = sum(1 for u in units if not u.startswith("S_"))
    if nf < n:
        ax.axhline(nf - 0.5, color=INK, lw=0.4, ls=":"); ax.axvline(nf - 0.5, color=INK, lw=0.4, ls=":")
    if n <= 45:
        ax.set_xticks(range(n)); ax.set_yticks(range(n)); ax.set_xticklabels(units, rotation=90, fontsize=5); ax.set_yticklabels(units, fontsize=5)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"{n} units{' (families only)' if families_only else ''}: nomination (lower) and audit (upper)")
    used = sorted({s for _, _, s in up}, key=list(STATUS_COL).index)
    ax.legend(handles=_handles(plt, [(k, STATUS_COL[k]) for k in used]), loc="upper left", bbox_to_anchor=(1.01, 1.0), title="upper triangle: nominated link", title_fontsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02, location="bottom")
    cb.set_label("lower triangle: fraction of cross-unit pairs at second-degree kinship or closer", fontsize=7)


def p_gate1(run, plt, ax):
    E = run.t["audit"]["edges"]; a = E[E.auditable & E.LLR_one.notna()]
    non = a[~a.nominated]
    ax.scatter(non.LLR_one, non.LLR_both, s=7, color="#CCCCCC", lw=0, label=f"not nominated ({len(non)})", rasterized=True)
    for st in STATUS_COL:
        g = a[a.nominated & (a.link_class.map(_status) == st)]
        if len(g):
            ax.scatter(g.LLR_one, g.LLR_both, s=22, color=STATUS_COL[st], edgecolor="white", lw=0.4, label=f"nominated, {st} ({len(g)})")
    ax.axvline(0, color=INK, lw=0.6, ls="--"); ax.axhline(0, color=INK, lw=0.6, ls="--")
    ax.set_xlabel("LLR, one shared parent against none"); ax.set_ylabel("LLR, both parents shared against none")
    ax.set_title("Gate 1: compatible in the lower right quadrant"); ax.legend(fontsize=6)


def p_gate2(run, plt, ax):
    E = run.t["audit"]["edges"]; G2 = run.t["audit"]["gate2"]
    cols = [c for c in E.columns if c.startswith("share_")]
    if not cols:
        ax.text(0.5, 0.5, "Gate 2 not run\n(no structural-variant class in the audit set)", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    a = E[E.auditable & E.LLR_one.notna()]
    for k, c in enumerate(cols):
        non = a[~a.nominated][c].dropna(); g1 = a[a.nominated & (a.gate1 == 1)]
        if len(non) > 2:
            vp = ax.violinplot([non], positions=[k - 0.2], widths=0.35, showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor("#CCCCCC"); b.set_alpha(1)
        for st in STATUS_COL:
            g = g1[g1.link_class.map(_status) == st]
            if len(g):
                ax.scatter(np.full(len(g), k + 0.2) + (np.random.default_rng(3).random(len(g)) - 0.5) * 0.15, g[c], s=18, color=STATUS_COL[st], edgecolor="white", lw=0.4)
        row = G2[G2["class"] == c.replace("share_", "")]
        if len(row) and "threshold_median_of_gate1_pass" in row:
            ax.plot([k, k + 0.4], [row.threshold_median_of_gate1_pass.iloc[0]] * 2, color=INK, lw=0.8, ls="--")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.replace("share_", "") for c in cols])
    ax.set_ylabel("sharing (carrier in both / carrier in either)"); ax.set_title("Gate 2: sharing of inherited structural variants")
    used = [k for k in STATUS_COL if (a[a.nominated & (a.gate1 == 1)].link_class.map(_status) == k).any()]
    ax.legend(handles=_handles(plt, [("non-nominated auditable pairs", "#CCCCCC")]) + _handles(plt, [(f"nominated, {k}", STATUS_COL[k]) for k in used], "o") +
              [plt.Line2D([], [], color=INK, ls="--", lw=0.8, label="support bar (median of Gate-1-compatible pairs)")], fontsize=6, loc="lower right")


def p_network(run, plt, ax, min_n=None):
    import networkx as nx
    L = run.t["network"]["links"]; N = run.t["network"]["nodes"]
    min_n = run.cfg["network"]["figure_min_family_size"] if min_n is None else min_n
    L = L[L.min_n >= min_n]
    if not len(L):
        ax.text(0.5, 0.5, f"no link between families of at least {min_n} individuals", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    G = nx.Graph()
    for r in L.itertuples():
        G.add_edge(r.unit_a, r.unit_b, status=r.status)
    size = dict(zip(N.family, N.n))
    comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c, key=_fkey)[0]))
    pos = {}; ncol = int(np.ceil(np.sqrt(len(comps))))
    for k, c in enumerate(comps):
        sub = G.subgraph(sorted(c, key=_fkey))
        p = nx.spring_layout(sub, seed=run.seed, k=1.2 / np.sqrt(max(len(c), 1)), iterations=200) if len(c) > 2 else {n_: (i - 0.5, 0) for i, n_ in enumerate(sorted(c, key=_fkey))}
        arr = np.array(list(p.values()), float); arr = arr - arr.mean(0); sc = np.abs(arr).max() or 1
        for n_, v in zip(p, arr / sc * 0.4):
            pos[n_] = (k % ncol + v[0], -(k // ncol) + v[1])
    for st, col in STATUS_COL.items():
        ed = [(u, v) for u, v, d in G.edges(data=True) if d["status"] == st]
        if ed:
            nx.draw_networkx_edges(G, pos, edgelist=ed, edge_color=col, width=2.2 if st == "supported" else 1.2, style="dashed" if st in ("provisional", "rejected", "nominated, not audited") else "solid", ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=[60 + 22 * size.get(n_, 1) for n_ in G.nodes], node_color="white", edgecolors=INK, linewidths=0.8, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=6, ax=ax)
    used = [k for k in STATUS_COL if any(d["status"] == k for _, _, d in G.edges(data=True))]
    ax.legend(handles=[plt.Line2D([], [], color=STATUS_COL[k], lw=2, ls="--" if k in ("provisional", "rejected", "nominated, not audited") else "-", label=k) for k in used], loc="lower left", bbox_to_anchor=(0, -0.08), ncols=3, fontsize=6)
    ax.set_title(f"Network of links between families of at least {min_n} individuals (node area ~ family size)"); ax.set_axis_off(); ax.margins(0.12)


def p_incidence(run, plt, ax):
    LP = run.t["network"]["latent_parents"]; t = run.cfg["network"]["figure_min_family_size"]
    LP = LP[LP.min_family_size.astype(int) == t]
    if not len(LP):
        ax.text(0.5, 0.5, "no latent parent", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    fams = sorted({f for x in LP.families for f in x.split(",")}, key=_fkey)
    viol = set(str(run.s["network"]["slot_violations"].get(str(t), run.s["network"]["slot_violations"].get(t, []))).strip("[]").replace("'", "").replace(" ", "").split(","))
    for j, r in enumerate(LP.itertuples()):
        ys = [fams.index(f) for f in r.families.split(",")]
        ax.plot([j, j], [min(ys), max(ys)], color="#BBBBBB", lw=1, zorder=1)
        ax.scatter([j] * len(ys), ys, s=26, color=[RED if fams[y] in viol else BLUE for y in ys], zorder=2, lw=0)
    ax.set_xticks(range(len(LP))); ax.set_xticklabels(LP.latent_parent, rotation=90, fontsize=6)
    ax.set_yticks(range(len(fams))); ax.set_yticklabels(fams, fontsize=6); ax.invert_yaxis()
    ax.set_xlabel("latent parent (bookkeeping label, not a reconstructed parent)"); ax.set_title(f"Latent-parent incidence (families of at least {t})")
    ax.legend(handles=_handles(plt, [("family within the two-parent slot limit", BLUE)] + ([("family violating the slot limit", RED)] if viol & set(fams) else []), "o"), loc="upper left", bbox_to_anchor=(0, -0.22), ncols=2, fontsize=6)


def p_barcode(run, plt, ax):
    B = run.t["barcode"]["barcode"]
    parents = list(dict.fromkeys(B.parent)); ints = list(dict.fromkeys(B.interval))
    col = {"0/0": BLUE, "0/1": GREY, "1/1": ORANGE}
    for r in B.itertuples():
        x, y = ints.index(r.interval), parents.index(r.parent)
        called = r.call != "unresolved"
        ax.scatter(x, y, s=70 if called else 26, marker="s", color=col[r.best_genotype] if called else "white", edgecolor=col[r.best_genotype], lw=1.0)
    ax.set_xticks(range(len(ints))); ax.set_xticklabels(ints, rotation=90, fontsize=5 if len(ints) > 12 else 6)
    ax.set_yticks(range(len(parents))); ax.set_yticklabels(parents, fontsize=6); ax.invert_yaxis()
    ax.set_xlim(-0.6, len(ints) - 0.4); ax.set_ylim(len(parents) - 0.4, -0.6)
    pc = run.cfg["barcode"]["posterior_call"]
    ax.set_title(f"Barcode: genotype of the shared parent per structural haplotype (filled = posterior >= {pc})")
    ax.legend(handles=_handles(plt, [("0/0", BLUE), ("0/1", GREY), ("1/1", ORANGE)]) + [plt.Line2D([], [], marker="s", ls="", markerfacecolor="white", markeredgecolor=INK, markersize=5, label="unresolved (best genotype outlined)")],
              loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=6)


def p_eta(run, plt, ax):
    P = run.t["family_effect"].get("per_marker"); T = run.t["family_effect"]["summary"]
    if P is None or not len(P):
        ax.text(0.5, 0.5, "no marker tested", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    sets = list(dict.fromkeys(P.marker_set))
    for k, name in enumerate(sets):
        g = P[P.marker_set == name]
        for off, colname, c in ((-0.2, "null_median", "#CCCCCC"), (0.2, "eta2_obs", BLUE)):
            v = g[colname].dropna().to_numpy()
            if len(v) > 3 and np.ptp(v) > 0:
                vp = ax.violinplot([v], positions=[k + off], widths=0.36, showextrema=False)
                for b in vp["bodies"]:
                    b.set_facecolor(c); b.set_alpha(1)
            else:
                ax.scatter(np.full(len(v), k + off), v, s=8, color=c)
            ax.plot([k + off - 0.12, k + off + 0.12], [np.median(v)] * 2, color=INK, lw=1)
        ax.text(k, 1.02, f"n = {len(g):,}", ha="center", fontsize=6)
    ax.set_xticks(range(len(sets))); ax.set_xticklabels([s.split(" (")[0] for s in sets]); ax.set_ylim(0, 1.08)
    ax.set_ylabel("eta squared (variance explained by family)"); ax.set_title("Family effect per marker class")
    ax.legend(handles=_handles(plt, [("null (family labels permuted; per-marker median)", "#CCCCCC"), ("observed", BLUE)]) + [plt.Line2D([], [], color=INK, lw=1, label="median")], loc="lower right", fontsize=6)


def p_concordance(run, plt, ax):
    P = run.t["concordance"].get("pairs"); T = run.t["concordance"]["summary"]
    if P is None or not len(P):
        ax.text(0.5, 0.5, "no family or no genotype input", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off(); return
    sets = list(dict.fromkeys(P.marker_set))
    for k, name in enumerate(sets):
        g = P[P.marker_set == name]
        for off, flag, c in ((-0.2, False, "#CCCCCC"), (0.2, True, BLUE)):
            v = g[g.same_family == flag].kinship.to_numpy()
            if len(v) > 3 and np.ptp(v) > 0:
                vp = ax.violinplot([v], positions=[k + off], widths=0.36, showextrema=False)
                for b in vp["bodies"]:
                    b.set_facecolor(c); b.set_alpha(1)
        row = T[T.marker_set == name].iloc[0]
        ax.text(k, ax.get_ylim()[1], f"{int(row.n_markers):,} markers\nAUC {row.auc_within_vs_between:.3f}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(range(len(sets))); ax.set_xticklabels([s.split(" (")[0] for s in sets])
    ax.set_ylabel("kinship (GRM / 2)"); ax.set_title("Kinship per marker class, within and between families", pad=22)
    ax.legend(handles=_handles(plt, [("pairs in different families", "#CCCCCC"), ("pairs in the same family", BLUE)]), fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1))


# ------------------------------------------------------------------ figures
def f_qc(run, plt):
    fig, ax = plt.subplots(2, 2, figsize=(9.5, 7), constrained_layout=True)
    p_callrate(run, plt, ax[0, 0]); p_het(run, plt, ax[0, 1]); p_depth(run, plt, ax[1, 0]); p_perchrom(run, plt, ax[1, 1])
    for a, l in zip(ax.ravel(), "ABCD"):
        _panel(a, l)
    return fig, "QC panels", "A call rate per sample and input with the low-coverage rule; B heterozygosity against call rate; C depth over the mask against depth over the assembly; D marker density per chromosome."


def f_rel(run, plt):
    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    p_theta_k2(run, plt, ax)
    return fig, "Relatedness: theta and k2 with thresholds", "Every pair of individuals; dotted lines are the kinship thresholds, the dashed line the k2 gate of full sibs. Without k2 (GRM) a histogram of theta is drawn."


def f_fam(run, plt):
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    p_family_sizes(run, plt, ax[0]); p_closure(run, plt, ax[1]); _panel(ax[0], "A"); _panel(ax[1], "B")
    return fig, "Family sizes and closure", "A families by size, coloured by closure; B closure against family size. Families are relatedness-based groups, not pedigree records."


def f_matrix(run, plt):
    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    p_matrix(run, plt, ax)
    return fig, "Nomination matrix with audit triangle", "Lower triangle: fraction of cross-unit pairs at second-degree kinship or closer. Upper triangle: nominated pairs coloured by audit outcome. Dotted lines separate families from unassigned individuals."


def f_gates(run, plt):
    if not run.s.get("audit", {}).get("ran"):
        return None
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    p_gate1(run, plt, ax[0]); p_gate2(run, plt, ax[1]); _panel(ax[0], "A"); _panel(ax[1], "B")
    return fig, "Audit gates", "A Gate 1 scores of every auditable family pair; B Gate 2 sharing per SV class, nominated Gate-1-compatible pairs against the non-nominated auditable pairs."


def f_incidence(run, plt):
    if "latent_parents" not in run.t.get("network", {}) or not len(run.t["network"]["latent_parents"]):
        return None
    LP = run.t["network"]["latent_parents"]
    n = (LP.min_family_size.astype(int) == run.cfg["network"]["figure_min_family_size"]).sum()
    fig, ax = plt.subplots(figsize=(max(4.5, 0.35 * n + 2.5), 5), constrained_layout=True)
    p_incidence(run, plt, ax)
    return fig, "Latent-parent incidence", "Each column is a latent parent (a maximal set of mutually linked families); a family can sit on at most two."


def f_network(run, plt):
    if "links" not in run.t.get("network", {}) or not len(run.t["network"]["links"]):
        return None
    fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
    p_network(run, plt, ax)
    return fig, "Network of families and links", "Nodes are families, edges are nominated links coloured by audit outcome; layout seeded."


def f_hap(run, plt):
    S = run.t.get("haplotypes", {}).get("summary"); I = run.t.get("haplotypes", {}).get("per_individual")
    if S is None or I is None or not len(I):
        return None
    S = S[S.tier != "not genotyped"].head(24)
    n = len(S); nc = min(4, n); nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.7 * nr + 0.5), constrained_layout=True, squeeze=False)
    for ax in axs.ravel()[n:]:
        ax.set_axis_off()
    for ax, r in zip(axs.ravel(), S.itertuples()):
        g = I[I.interval == r.interval]
        g = g[g.markers == ("SNPs" if r.genotyped_from == "SNPs" else "deletions")]
        for k in (0, 1, 2, -1):
            x = g[g.genotype.astype(int) == k]
            ax.scatter(x.pc1, x.het_inside, s=9, color=GENO_COL[k], lw=0.3, edgecolor=INK if k == -1 else "none")
        ax.set_title(f"{r.interval}\n{r.source.split(' (')[0]}, {r.tier.replace('_', ' ')}, {r.genotyped_from}", fontsize=6.5)
        ax.set_xlabel("PC1"); ax.set_ylabel("heterozygosity inside")
    fig.legend(handles=_handles(plt, [(GENO_LAB[k], GENO_COL[k]) for k in (0, 1, 2, -1)], "o"), loc="outside lower center", ncols=4)
    return fig, "Structural haplotypes: PC1 against heterozygosity per interval", "One panel per interval (first 24). Three clusters on PC1 with excess heterozygosity in the middle cluster are the signature of a non-recombining structural haplotype. Labels 0 / 2 carry no ancestral meaning."


def f_scan(run, plt):
    h = run.t.get("haplotypes", {})
    W1, W2 = h.get("scan_windows_snp"), h.get("scan_windows_del")
    if W1 is None and W2 is None:
        return None
    S = h.get("summary")
    rows = [w for w in (W1, W2) if w is not None and len(w)]
    fig, axs = plt.subplots(len(rows), 1, figsize=(10, 2.8 * len(rows)), constrained_layout=True, squeeze=False)
    for ax, W in zip(axs[:, 0], rows):
        snp = "z1" in W
        chroms = list(dict.fromkeys(W.chrom)); off = {}; x0 = 0
        for c in chroms:
            off[c] = x0; x0 += W[W.chrom == c].end.max() * 1.02
        x = (W.start + W.end) / 2 + W.chrom.map(off)
        y = np.fmax(W.z1.abs(), W.z2.abs().fillna(0)) if snp else W.LD_mean_abs_r
        ax.scatter(x[~W.outlier], y[~W.outlier], s=5, color=GREY, lw=0, label="window")
        ax.scatter(x[W.outlier], y[W.outlier], s=9, color=RED, lw=0, label="outlier window")
        if S is not None and len(S):
            for r in S.itertuples():
                if r.chrom in off:
                    sup = str(r.source).startswith("supplied")
                    ax.axvspan(off[r.chrom] + r.start, off[r.chrom] + r.end, ymin=0.0 if sup else 0.94, ymax=0.06 if sup else 1.0, color=BLUE if sup else ORANGE, lw=0)
        for c in chroms:
            ax.axvline(off[c], color="#DDDDDD", lw=0.5, zorder=0)
        ax.set_xticks([off[c] + W[W.chrom == c].end.max() / 2 for c in chroms]); ax.set_xticklabels(chroms, fontsize=6, rotation=90 if len(chroms) > 12 else 0)
        ax.set_ylabel("max |z| of scaling coordinates 1-2" if snp else "LD among deletions (mean |r|)")
        ax.set_title("SNP scan: local PCA in windows (reimplemented from Methods)" if snp else "Deletion scan: linkage disequilibrium in fixed windows")
        ax.legend(handles=_handles(plt, [("window", GREY), ("outlier window", RED)], "o") + _handles(plt, [("supplied interval (bottom bar)", BLUE), ("scan candidate, not validated (top bar)", ORANGE)]), ncols=4, fontsize=6, loc="upper right")
    return fig, "Scan for large structural haplotypes", "Window statistics along the genome with outlier windows, supplied intervals and scan candidates. Candidates are detected automatically and not validated."


def f_tags(run, plt):
    T = run.t.get("haplotypes", {}).get("tags")
    if T is None or not len(T):
        return None
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for k, (cls, col, mk) in enumerate((("DEL", BLUE, "o"), ("DUP", ORANGE, "^"))):
        g = T[T.svtype == cls]
        ax[0].scatter(g.carrier_frequency_class_0, g.carrier_frequency_class_2, s=14, color=col, marker=mk, lw=0, alpha=0.8, label=f"{cls} ({len(g)})")
    ax[0].plot([0, 1], [0, 1], color="#DDDDDD", lw=0.6, zorder=0)
    ax[0].set_xlabel("carrier frequency in class 0"); ax[0].set_ylabel("carrier frequency in class 2"); ax[0].set_title("Arrangement-tag variants"); ax[0].legend()
    ints = list(dict.fromkeys(T.interval)); c = T.groupby(["interval", "fully_private"]).size().unstack(fill_value=0).reindex(ints).fillna(0)
    a = c.get(False, pd.Series(0, index=ints)); b = c.get(True, pd.Series(0, index=ints))
    ax[1].barh(range(len(ints)), a, color=SKY, label="tag variant"); ax[1].barh(range(len(ints)), b, left=a, color=BLUE, label="fully private to one class")
    ax[1].set_yticks(range(len(ints))); ax[1].set_yticklabels(ints, fontsize=6); ax[1].invert_yaxis(); ax[1].set_xlabel("variants"); ax[1].set_title("Tag variants per interval"); ax[1].legend()
    _panel(ax[0], "A"); _panel(ax[1], "B")
    return fig, "Arrangement-tag deletions and duplications", "A carrier frequency of each tag variant in the two homokaryotype classes; B tag variants per interval, with those fully private to one class."


def f_barcode(run, plt):
    if "barcode" not in run.t.get("barcode", {}) or not len(run.t["barcode"]["barcode"]):
        return None
    B = run.t["barcode"]["barcode"]
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * B.interval.nunique() + 4), max(2.6, 0.3 * B.parent.nunique() + 1.8)), constrained_layout=True)
    p_barcode(run, plt, ax)
    return fig, "Barcode of the shared parents", "Rows are links (family x family) and junction families (J_); columns are structural haplotypes."


def f_eta(run, plt):
    if "summary" not in run.t.get("family_effect", {}):
        return None
    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    p_eta(run, plt, ax)
    return fig, "Family effect (eta squared)", "Observed eta squared per marker against the per-marker null median, per marker class."


def f_conc(run, plt):
    if "pairs" not in run.t.get("concordance", {}):
        return None
    fig, ax = plt.subplots(figsize=(8.5, 4.2), constrained_layout=True)
    p_concordance(run, plt, ax)
    return fig, "Cross-marker kinship concordance", "Kinship from each marker class for pairs in the same family and in different families (families from the relatedness source; concordance, not pedigree truth)."


def f_overview(run, plt):
    if not run.s.get("audit", {}).get("ran") or "links" not in run.t.get("network", {}) or not len(run.t["network"]["links"]):
        return None
    fig = plt.figure(figsize=(13, 13.5), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1, 1])
    ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]
    p_matrix(run, plt, ax[0], families_only=True); p_network(run, plt, ax[1]); p_gate1(run, plt, ax[2]); p_gate2(run, plt, ax[3]); p_incidence(run, plt, ax[4])
    if "barcode" in run.t.get("barcode", {}) and len(run.t["barcode"]["barcode"]):
        p_barcode(run, plt, ax[5])
    else:
        ax[5].text(0.5, 0.5, "no barcode (no structural-haplotype genotypes)", ha="center", va="center", transform=ax[5].transAxes); ax[5].set_axis_off()
    for a, l in zip(ax, "ABCDEF"):
        _panel(a, l)
    return fig, "Shared-parent overview", "A nomination matrix (families only) with audit triangle; B network; C Gate 1; D Gate 2; E latent-parent incidence; F barcode. A compact composite in the spirit of the published shared-parent figure, not a panel-by-panel port."


FIGURES = [("fig01_qc_panels", f_qc, "qc"), ("fig02_relatedness_theta_k2", f_rel, "relatedness"), ("fig03_family_sizes_closure", f_fam, "families"),
           ("fig04_nomination_matrix_audit", f_matrix, "nomination"), ("fig05_audit_gates", f_gates, "audit"), ("fig06_latent_parent_incidence", f_incidence, "network"),
           ("fig07_network", f_network, "network"), ("fig08_structural_haplotypes_pc1", f_hap, "haplotypes"), ("fig09_structural_haplotype_scan", f_scan, "haplotypes"),
           ("fig10_arrangement_tag_variants", f_tags, "haplotypes"), ("fig11_barcode", f_barcode, "barcode"), ("fig12_family_effect", f_eta, "family_effect"),
           ("fig13_cross_marker_concordance", f_conc, "concordance"), ("fig14_shared_parent_overview", f_overview, "network")]


def make_figures(run):
    plt = _plt()
    d = os.path.join(run.out, "figures"); os.makedirs(d, exist_ok=True)
    dpi = run.cfg["figures"]["dpi"]; rows = []
    for stem, fn, step in FIGURES:
        try:
            if step not in run.t or (step == "nomination" and "pairs" not in run.t["nomination"]):
                res = None
            else:
                res = fn(run, plt)
        except Exception as exc:                                  # a figure never stops a run
            log(f"WARNING: figure {stem} failed ({type(exc).__name__}: {exc})")
            rows.append({"file": stem, "title": "", "status": f"failed: {type(exc).__name__}: {exc}", "description": ""}); plt.close("all"); continue
        if res is None:
            rows.append({"file": stem, "title": "", "status": "not drawn: the step did not run or produced nothing to draw", "description": ""}); continue
        fig, title, desc = res
        try:
            fig.canvas.draw(); fig.set_layout_engine("none")          # lay out once, then save twice
            fig.savefig(os.path.join(d, stem + ".png"), dpi=dpi); fig.savefig(os.path.join(d, stem + ".pdf"))
        except Exception as exc:                                  # e.g. a write error: report it, keep the run
            log(f"WARNING: figure {stem} could not be saved ({type(exc).__name__}: {exc})")
            rows.append({"file": stem, "title": title, "status": f"failed: {type(exc).__name__}: {exc}", "description": desc}); plt.close("all"); continue
        plt.close(fig)
        rows.append({"file": stem, "title": title, "status": "drawn (PNG + PDF)", "description": desc})
    M = pd.DataFrame(rows)
    write_tsv(M, os.path.join(d, "MANIFEST.tsv"), run.header("Figures of the run: file stem, title, status, description."))
    return M
