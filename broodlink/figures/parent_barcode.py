"""Two figures on the barcode of a shared parent.

`fig_augmentation`: port of the published Extended Data figure "Barcode of a shared parent with
deletion augmentation" (`HAT/figures/upgrades_2026-11-03/round_ed/ED_Fig_H5/ED_Fig_H5.py`,
panels A-H; panel I, the held-out prediction of withheld offspring, is a validation outside
the scope of the tool and is not drawn).
`fig_funnel`: port of the published Supplementary figure "Funnel from offspring relationships to
an unsampled parent" (`.../round_ed/Supp_Fig_S-H4/Supp_Fig_S-H4.py`, panels A-I).

Same panels, layout (millimetres), encodings and style as the originals; every name, count and
threshold comes from the tables of the run. The original arrangement labels STD / INV are the
classes 0 / 2 of the structural-haplotype genotypes (no ancestral meaning); links are named by
their two families and a parent shared by several links by its junction name (J_<family>).
"""
from __future__ import annotations
import numpy as np

from ..augmentation import interval_numbers
from ..network import status_of
from . import style as S
from .kinship import hap_label
from .shared_parent import audited, edges, noun

L = S.LED
FS = 5.5
GREYT = "#3D4650"
C0, C1, C2, NOI = S.STATE0, S.STATE1, S.STATE2, S.NO_INFO
CODE = {"0/0": 0, "0/1": 1, "1/1": 2}
ST_LAB = {0: "class 0 hom.", 1: "het.", 2: "class 2 hom."}
SRC_B, SRC_EV, SRC_SW, SRC_SB, SRC_T, SRC_RC = ("07_barcode/shared_parent_barcode.tsv", "07_barcode/tag_evidence_by_setting.tsv", "07_barcode/tag_evidence_sweep.tsv",
                                                "07_barcode/junction_super_barcode.tsv", "06_structural_haplotypes/arrangement_tag_variants.tsv", "07_barcode/junction_recurrence.tsv")
SRC_F, SRC_M, SRC_E, SRC_G2 = "02_families/families.tsv", "02_families/family_membership.tsv", "04_audit/audit_edges.tsv", "04_audit/audit_gate2_summary.tsv"


# ------------------------------------------------------------------ drawing helpers (same as the original scripts)
class _Kit:
    def __init__(self, cv):
        self.cv, self.W, self.H, self.fig, self.ov = cv, cv.W, cv.H, cv.fig, cv.overlay

    def AX(self, l, t, w, h):
        ax = self.cv.ax(l, t, w, h); ax.tick_params(labelsize=FS); return ax

    def TX(self, x, y, s, fs=FS, b=False, c=S.NAVY, ha="left", va="center", **k):
        return self.fig.text(x / self.W, 1 - y / self.H, s, fontsize=fs, fontweight="bold" if b else "normal", color=c, ha=ha, va=va, **k)

    def PT(self, letter, title, x, y, sub=None):
        self.TX(x, y, letter, 8, True); self.TX(x + 3.6, y, title, 6.5, True)
        if sub:
            self.TX(x + 3.6, y + 2.0, sub, c=GREYT, va="top", linespacing=1.05)

    def BAN(self, y, n, text):
        from matplotlib.lines import Line2D
        self.TX(1, y, n, 8, True); t = self.TX(5, y, text, 6.8, True)
        w = t.get_window_extent(self.fig.canvas.get_renderer()).width / self.fig.dpi * 25.4
        if w > self.W - 6:                                     # a long sentence is set smaller, never cut
            t.set_fontsize(max(FS, 6.8 * (self.W - 6) / w))
        self.fig.add_artist(Line2D([.5 / self.W, (self.W - .5) / self.W], [1 - (y + 2.4) / self.H] * 2, color=S.NAVY, lw=.35, transform=self.fig.transFigure))

    def SW(self, x, y, c, lab, ht=None, ec="none", w=2.2, h=1.7, alpha=None):
        from matplotlib.patches import Rectangle
        self.ov.add_patch(Rectangle((x, y - h / 2), w, h, fc=c, ec=(ec if ec != "none" else ("white" if ht else "none")), lw=.3, hatch=ht, alpha=alpha))
        self.TX(x + w + .8, y, lab)

    def strip(self, x0, y0, w, h, vals, colmap, hatchset=(), alpha=None, outline=None):
        from matplotlib.patches import Rectangle
        n = len(vals)
        for i, v in enumerate(vals):
            v = int(v)
            self.ov.add_patch(Rectangle((x0 + i * w / n, y0), w / n, h, fc=colmap.get(v, NOI), ec="white", lw=.12, hatch="//////" if v in hatchset else None,
                                        alpha=(1 if alpha is None or alpha[i] else .38), zorder=5))
            if outline is not None and v == outline:
                self.ov.add_patch(Rectangle((x0 + i * w / n, y0), w / n, h, fill=False, ec=S.CONFLICT, lw=.5, zorder=7))
        self.ov.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=S.NAVY, lw=.3, zorder=6))

    def numbers(self, x0, w, y, n):
        """Interval numbers under a strip; thinned when a cell is too narrow for a label."""
        cell = w / max(n, 1); step = 1 if cell >= 3.0 else 2 if cell >= 1.6 else int(np.ceil(4.0 / cell / 5.0) * 5)
        for i in range(0, n, step):
            self.TX(x0 + (i + .5) * w / n, y, f"{i + 1}", ha="center")


def _xname(parent):
    return str(parent).replace("x", " × ", 1) if "x" in str(parent) and not str(parent).startswith("J_") else str(parent)


def _singular(w):
    return w[:-1] if w.endswith("s") else w


def _vectors(run):
    """Barcode table as vectors in genome order: (names, {parent: {column: codes}}, link list, junction list)."""
    B = run.t["barcode"]["barcode"]
    num = interval_numbers(B.interval.astype(str)); names = sorted(num, key=num.get)
    out, links, juncs = {}, [], []
    for par, g in B.groupby("parent", sort=False):
        d = {}
        for col in ("call", "best_genotype", "private_parent_a_call", "private_parent_b_call"):
            if col in g:
                m = dict(zip(g.interval.astype(str), g[col].astype(str)))
                d[col] = np.array([CODE.get(m.get(x, ""), -1) for x in names], int)
        d["families"] = str(g.families.iloc[0]).split(","); d["status"] = str(g.link_status.iloc[0]); d["kind"] = str(g.kind.iloc[0])
        out[str(par)] = d
        (links if d["kind"] == "link" else juncs).append(str(par))
    return names, out, links, juncs


def choose_junction(run, V, links):
    """The junction family shown: most supported links, ties by the number of called states of its joint
    barcode, then by name. Returns (family, [its supported links in the order of the link table]) or (None, [])."""
    junc = run.s.get("network", {}).get("junction_families", {}) or {}
    best = None
    for f, partners in junc.items():
        ls = [p for p in links if f in V[p]["families"] and set(V[p]["families"]) - {f} <= set(partners) and V[p]["status"].startswith("supported")]
        if len(ls) < 2 or "J_" + str(f) not in V:
            continue
        key = (-len(ls), -int((V["J_" + str(f)]["call"] >= 0).sum()), S.fam_key(f))
        if best is None or key < best[0]:
            best = (key, str(f), ls)
    return (best[1], best[2]) if best else (None, [])


def _pct(x):
    return f"{x:.0%}"


# ====================================================================== figure: barcode with tag-variant evidence
def fig_augmentation(run, plt):
    """A tag variants per interval, B tag strength, C evidence class against barcode state, D sweep over k
    and delta, E barcode before and after, F information per link, G vote of a junction's links,
    H evidence per SV class."""
    from matplotlib.patches import Rectangle
    bt = run.t.get("barcode", {})
    if "barcode" not in bt or not len(bt["barcode"]):
        return None
    if "evidence_by_setting" not in bt or "evidence_sweep" not in bt or not len(bt["evidence_sweep"]):
        return "no arrangement-tag deletion / duplication matched to the barcode intervals (step 6 found no tag variant, or there is no deletion / duplication input)"
    EV, SWt = bt["evidence_by_setting"], bt["evidence_sweep"]
    SBt = bt.get("super_barcode")
    names, V, links, _j = _vectors(run); n = len(names); num = {x: i for i, x in enumerate(names)}
    sets = list(dict.fromkeys(SWt.evidence_set.astype(str)))
    main = "DEL" if "DEL" in sets else sets[0]
    kl = S.class_label(main); kl1 = S.class_label(main, plural=False)
    notes = []
    dflt = SWt.default_setting.astype(str) == "True"
    r0 = SWt[dflt & (SWt.evidence_set.astype(str) == main)].iloc[0]; K0, DM0 = int(r0.k), float(r0.delta_min)
    e0 = EV[(EV.default_setting.astype(str) == "True") & (EV.evidence_set.astype(str) == main)]
    ECODE = {"none": 0, "class 0": 1, "class 2": 2, "class 0 + class 2": 3}
    cls = {}
    for par, g in e0.groupby("parent", sort=False):
        m = dict(zip(g.interval.astype(str), g.evidence.astype(str)))
        cls[str(par)] = np.array([ECODE.get(m.get(x, "none"), 0) for x in names], int)
    links = [p for p in links if p in cls]
    sup = [p for p in links if V[p]["status"].startswith("supported")]
    pooled = sup or links
    hl = hap_label(run); hl1 = _singular(hl)
    # ---- tag variants of the main class inside the barcode intervals
    T = run.t.get("haplotypes", {}).get("tags")
    T = T[(T.svtype.astype(str) == main) & T.interval.astype(str).isin(num)] if T is not None and len(T) else None
    has_tags = T is not None and len(T) > 0
    # ---- which links / junction are shown
    jf, jl = choose_junction(run, V, links)
    has_G = jf is not None and SBt is not None and len(SBt) and ("J_" + jf) in set(SBt.junction.astype(str))
    if has_G:
        e_links = jl[:3]
        if len(jl) > 3:
            notes.append(f"panel E shows the first 3 of the {len(jl)} supported links of junction family {jf}")
        nj = len(run.s.get("network", {}).get("junction_families", {}) or {})
        notes.append(f"panels E and G show junction family {jf}" + (f" (of {nj} junction families: the one with the most supported links, ties by the number of called states of its joint barcode)" if nj > 1 else ""))
    else:
        e_links = sorted(pooled, key=lambda p: (-int((V[p]["call"] >= 0).sum()), links.index(p)))[:2]
        notes.append("panel G omitted: no junction family with two supported links; panel E shows the " + ("supported " if sup else "") + "links with the most called states")
    if not sup:
        notes.append("no supported link: panels C, D and F pool all links that have a barcode")
    f_links = pooled[:14]
    if len(pooled) > 14:
        notes.append(f"panel F shows the first 14 of {len(pooled)} links")
    multi = len(sets) > 1
    if not multi:
        notes.append(f"panel H omitted: one SV class ({main}) gives tag variants")
    notes.append("panel I of the published figure (held-out prediction of withheld offspring) is a validation outside the tool and is not drawn")
    # ---- geometry
    W = 183.0
    YA = 6.7 if has_tags else None
    Y_B = 30.0 if has_tags else 6.7                             # row B-D title line
    Y2 = Y_B + 35.5
    yyE = Y2 + 6.7; yk = yyE + 9.4 * len(e_links) + 3.0
    mlen = max([len(_xname(p).replace(" ", "")) for p in f_links] + [6])
    TF = yk + 8; yfk = TF + 28.5 + ((mlen - 6) * .6 + 1.0 if mlen > 6 else 0); Y3 = yfk + 5
    H = (Y3 + 34.0) if multi else (yfk + 4.0)
    cv = S.canvas(plt, W, H, S.SZ_ED, thin=True); K = _Kit(cv); ov = K.ov
    import matplotlib
    matplotlib.rcParams.update({"hatch.linewidth": .3})
    # ---- pooled table (rows: barcode state, columns: evidence class)
    tab = np.zeros((3, 4), int)
    for p in pooled:
        for pv, cvv in zip(V[p]["call"], cls[p]):
            if pv >= 0:
                tab[pv, cvv] += 1
    hom_ok, hom_n = int(tab[0, 1] + tab[2, 2]), int(tab[0, 1] + tab[0, 2] + tab[2, 1] + tab[2, 2])
    acc = hom_ok / hom_n if hom_n else float("nan")
    # ---- banner 1
    b1 = (f"{L('banner1', f'arrangement-tag {kl}', int(len(T)), SRC_T):,} arrangement-tag {kl}; " if has_tags else "")
    if hom_n:
        b1 += (f"{kl1} evidence matches the barcode at {_pct(L('banner1', 'match rate, homozygous calls', acc, SRC_EV))} ({L('banner1', 'matching homozygous calls', hom_ok, SRC_EV)} of "
               f"{L('banner1', 'homozygous calls with evidence', hom_n, SRC_EV)} homozygous calls) against a null of {_pct(L('banner1', 'null mean', float(r0.null_mean), SRC_SW))}")
    else:
        b1 += f"no homozygous barcode state carries {kl1} evidence at k = {K0}, δ ≥ {DM0:g}"
    K.BAN(2.2, "1", b1)
    # ---- A, B
    if has_tags:
        j = np.array([num[x] for x in T.interval.astype(str)]); t2 = T.tagged_arrangement.astype(str).str.endswith("2").to_numpy(); dl = T.difference.to_numpy(float)
        nI = np.bincount(j[t2], minlength=n); nS = np.bincount(j[~t2], minlength=n); n80j = np.bincount(j[dl >= DM0], minlength=n); n80 = int((dl >= DM0).sum())
        K.PT("A", f"Arrangement-tag {kl} per {hl1}", .5, YA)
        K.TX(72, YA, f"median {int(L('A', 'median tag variants per interval', float(np.median(nI + nS)), SRC_T))} per interval", c=GREYT)
        K.SW(101, YA, C2, f"class 2-tagged ({L('A', 'class 2-tagged', int(nI.sum()), SRC_T):,})"); K.SW(128, YA, C0, f"class 0-tagged ({L('A', 'class 0-tagged', int(nS.sum()), SRC_T):,})")
        ov.scatter([155.5], [YA], s=2, c=S.NAVY, lw=0); K.TX(157, YA, f"tag strength ≥ {DM0:g} ({L('A', f'tag strength >= {DM0:g}', n80, SRC_T):,})")
        ax = K.AX(9, 9, 172, 15); x = np.arange(n)
        ax.bar(x, nI, color=C2, width=.82); ax.bar(x, nS, bottom=nI, color=C0, width=.82); ax.scatter(x, n80j, s=1.6, c=S.NAVY, lw=0, zorder=4)
        step = 1 if n <= 45 else int(np.ceil(n / 45))
        ax.set_xticks(x[::step]); ax.set_xticklabels([f"{i + 1}" for i in x[::step]], fontsize=FS); ax.set_ylabel(kl, labelpad=.8, fontsize=FS); ax.set_xlim(-.7, n - .3); ax.tick_params(axis="x", length=0)
        dmin_cfg = float(run.cfg["inversions"]["tags"]["min_delta"])
        K.PT("B", f"Tag strength of the {kl}", .5, Y_B, f"cut δ ≥ {DM0:g}: {n80:,} of {len(T):,}"); ax = K.AX(9, Y_B + 6, 36, 21)
        lo = min(np.floor(dl.min() * 20 + 1e-9) / 20, DM0 - .05); bins = np.arange(lo, 1.0001, .05)
        ax.hist([dl[t2], dl[~t2]], bins=bins, stacked=True, color=[C2, C0], lw=0); ax.axvline(DM0, color=S.NAVY, ls=(0, (3, 2)), lw=.35)
        ax.set_xlabel("tag strength (δ)", labelpad=.8, fontsize=FS); ax.set_ylabel(kl, labelpad=.8, fontsize=FS)
        notes.append(f"panels A and B: tag variants of step 6 inside the barcode intervals (carrier-frequency difference of at least {dmin_cfg:g}, the configured minimum); interval numbers follow genome order")
        xC, xD = 50, 112
    else:
        notes.append("panels A and B omitted: no tag-variant table"); xC, xD = .5, 62
    # ---- C
    psub = f"P = {L('C', 'permutation P', float(r0.p_perm), SRC_SW):.3f}" if hom_n and r0.p_perm == r0.p_perm else "no test"
    K.PT("C", f"{kl1} evidence against barcode state of P", xC, Y_B, f"{L('C', 'links pooled', len(pooled), SRC_EV)} link{'s' if len(pooled) != 1 else ''} pooled; outlined, expected; shade, count; {psub}")
    ax = K.AX(xC + 17, Y_B + 6, 42, 21)
    ax.pcolormesh(np.arange(5) - .5, np.arange(4) - .5, np.ma.masked_where(tab == 0, tab), cmap=S.ros_seq(), vmin=0, vmax=max(tab.max(), 1), lw=0, antialiased=False); ax.set_xlim(-.5, 3.5); ax.set_ylim(2.5, -.5)
    exp = ((0, 1), (2, 2), (1, 3)); ev_lab = ["none", "class 0", "class 2", "class 0 + 2"]
    for i in range(3):
        for q in range(4):
            ax.text(q, i, str(L("C", f"{ST_LAB[i]} x {ev_lab[q]}", int(tab[i, q]), SRC_EV + f" pooled, k={K0}, delta>={DM0:g}")), ha="center", va="center", fontsize=FS,
                    color="white" if tab[i, q] > tab.max() * .55 else S.NAVY, fontweight="bold" if (i, q) in exp else "normal")
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set_xticks(range(4)); ax.set_xticklabels(ev_lab, fontsize=FS); ax.set_yticks(range(3)); ax.set_yticklabels([ST_LAB[i] for i in range(3)], fontsize=FS)
    ax.set_xlabel(f"{kl1} evidence in the interval", labelpad=.8, fontsize=FS); ax.set_ylabel("barcode state", labelpad=.8, fontsize=FS); ax.tick_params(length=0)
    for (i, q) in exp:
        ax.add_patch(Rectangle((q - .5, i - .5), 1, 1, fill=False, ec=S.NAVY, lw=.6))
    # ---- D
    Sm = SWt[SWt.evidence_set.astype(str) == main]; ok = Sm[Sm.n_homozygous_informative > 0]
    dsub = None
    if len(ok):
        amin, amax, nmin, nmax = float(ok.acc_homozygous.min()), float(ok.acc_homozygous.max()), float(ok.null_mean.min()), float(ok.null_mean.max())
        dsub = (f"{_pct(L('D', 'lowest match rate', amin, SRC_SW))}–{_pct(L('D', 'highest match rate', amax, SRC_SW))} match against nulls "
                f"{_pct(L('D', 'lowest null mean', nmin, SRC_SW))}–{_pct(L('D', 'highest null mean', nmax, SRC_SW))}")
    K.PT("D", "Sweep over k and δ", xD, Y_B, dsub); ax = K.AX(xD + 9, Y_B + 6, 36, 21)
    dms = sorted(set(Sm.delta_min.astype(float))); pal = ["#8A949C", S.NAVY, "#2B6A6A", "#D08C6E", "#B8A9CC"]
    DC = [(dm, pal[i % len(pal)]) for i, dm in enumerate(dms)]
    for dm, col in DC:
        g = Sm[np.isclose(Sm.delta_min.astype(float), dm)].sort_values("k"); gg = g[g.n_homozygous_informative > 0]
        ax.plot(gg.k, gg.acc_homozygous, marker="o", ms=1.6, mew=0, lw=.45, color=col); ax.plot(gg.k, gg.null_mean, ls=(0, (1, 1)), lw=.4, color=col)
    ks = sorted(set(Sm.k.astype(int)))
    ax.set_xticks(ks); ax.set_xlabel(f"minimum tagging {kl} (k)", labelpad=.8, fontsize=FS); ax.set_ylabel("match, homozygous calls", labelpad=.8, fontsize=FS); ax.set_ylim(0, 1.02)
    for k_, (dm, col) in enumerate(DC):
        ov.plot([xD + 47.5, xD + 50], [Y_B + 9 + k_ * 2.4] * 2, color=col, lw=.6); K.TX(xD + 51, Y_B + 9 + k_ * 2.4, f"δ ≥ {dm:g}")
    yn = Y_B + 9 + len(DC) * 2.4 + .2
    ov.plot([xD + 47.5, xD + 50], [yn] * 2, color=S.NAVY, lw=.45, ls=(0, (1, 1))); K.TX(xD + 51, yn, "null")
    # ---- E
    K.PT("E", f"Barcode of P before and after {kl1} augmentation ({hl}, intervals 1–{n})", .5, Y2 + 4.5); yy = yyE; XS, WS = 25, 145; rowsE = []
    B_ = {0: C0, 1: C1, 2: C2}; DE = {1: C0, 2: C2, 3: C1}; AU = {**B_, 11: C0, 12: C2, 13: C1}
    for p in e_links:
        p_a, c_a = V[p]["call"], cls[p]; aug = np.where(p_a >= 0, p_a, np.where(c_a > 0, 10 + c_a, -1))
        n1, n2, n3 = int((p_a >= 0).sum()), int((c_a > 0).sum()), int((aug >= 0).sum()); rowsE.append((n1, n3))
        K.TX(1, yy + 1, _xname(p), b=True)
        for k_, (l_, v_, cm, hs, nn) in enumerate((("barcode", p_a, B_, (), n1), (f"{kl1} evidence", c_a, DE, (1, 2, 3), n2), ("augmented", aug, AU, (11, 12, 13), n3))):
            y_ = yy + k_ * 2.9
            K.TX(XS - 1, y_ + 1.2, l_, ha="right", b=(k_ == 2), c=S.NAVY if k_ == 2 else GREYT); K.strip(XS, y_, WS, 2.4, v_, cm, hs)
            K.TX(XS + WS + 1.5, y_ + 1.2, f"{L('E', f'{p} {l_}', nn, SRC_EV if k_ else SRC_B)} / {n}", b=(k_ == 2))
        yy += 9.4
    K.numbers(XS, WS, yy + .3, n)
    for k_, (c, ht, lab) in enumerate(((C0, None, ST_LAB[0]), (C1, None, ST_LAB[1]), (C2, None, ST_LAB[2]), (C0, "//////", f"≥ 1 class 0 ({kl1})"), (C2, "//////", f"≥ 1 class 2 ({kl1})"),
                                       (C1, "//////", f"class 0 + 2 ({kl1})"), (NOI, None, "no information"))):
        K.SW(XS + k_ * 21, yk, c, lab, ht, ec=S.NAVY if ht is None else "none")
    # ---- F
    K.PT("F", "Information on P per link", .5, yk + 4); wF = 48.0 if len(f_links) > 3 else 16.0 * max(len(f_links), 1); ax = K.AX(9, TF, wF, 22); xs = np.arange(len(f_links))
    b_ = [int((V[p]["call"] >= 0).sum()) for p in f_links]; d_ = [int(((V[p]["call"] < 0) & (cls[p] > 0)).sum()) for p in f_links]; u_ = [n - b - d for b, d in zip(b_, d_)]
    ax.bar(xs, b_, color=S.NAVY, width=.74); ax.bar(xs, d_, bottom=b_, color="#95BADD", width=.74, hatch="//////", ec="white", lw=0); ax.bar(xs, u_, bottom=np.array(b_) + np.array(d_), color=NOI, width=.74, lw=0)
    for i, (bb, dd) in enumerate(zip(b_, d_)):
        L("F", f"{f_links[i]} barcode", bb, SRC_B); L("F", f"{f_links[i]} barcode + {kl1} evidence", bb + dd, SRC_EV)
        ax.text(i, n * 1.025, f"{bb}→{bb + dd}", ha="center", va="bottom", fontsize=FS if len(f_links) <= 7 else 4.6)
    from matplotlib.ticker import MaxNLocator
    ax.set_xticks(xs); ax.set_xticklabels([_xname(p).replace(" ", "") for p in f_links], fontsize=FS, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylabel(f"intervals (of {n})", labelpad=.8, fontsize=FS); ax.set_ylim(0, n * 47 / 41); ax.set_xlim(-.6, len(f_links) - .4); ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
    K.SW(1, yfk, S.NAVY, "barcode"); K.SW(14, yfk, "#95BADD", f"added by {kl1}", "//////"); K.SW(33, yfk, NOI, "none", ec=S.NAVY)
    # ---- G: the vote of the junction's links, next to the joint barcode of the junction
    b2_tail = ""
    if has_G:
        g = SBt[SBt.junction.astype(str) == "J_" + jf]
        st = dict(zip(g.interval.astype(str), g.super_state.astype(str)))
        SCODE = {"0/0": 0, "0/1": 1, "1/1": 2, "at least one class 0": 10, "at least one class 2": 12, "opposing votes": -2, "no vote": -1}
        SUP = np.array([SCODE.get(st.get(x, "no vote"), -1) for x in names]); lp_code = V["J_" + jf]["best_genotype"]; lp_conf = V["J_" + jf]["call"] >= 0
        full = int(np.isin(SUP, (0, 1, 2)).sum()); part = int(np.isin(SUP, (10, 12)).sum()); cflt = int((SUP == -2).sum()); none_ = int((SUP == -1).sum())
        both = np.isin(SUP, (0, 1, 2)) & lp_conf; agree = int((SUP[both] == lp_code[both]).sum())
        GX, GW = 86, 80
        K.PT("G", f"J_{jf} super barcode: barcode, {kl1} evidence and votes", 60, yk + 4,
             f"{L('G', 'agreeing intervals', agree, SRC_SB)} of {L('G', 'intervals with a full vote and a called joint state', int(both.sum()), SRC_SB)} intervals genotyped in both agree; "
             f"{L('G', 'opposing votes', cflt, SRC_SB)} opposing votes, {L('G', 'no vote', none_, SRC_SB)} without a vote")
        yg = TF + 2.5; SC = {0: C0, 1: C1, 2: C2, 10: C0, 12: C2, -2: "#FFFFFF", -1: NOI}
        K.TX(GX - 1, yg + 1.8, f"J_{jf} joint barcode", ha="right", c=GREYT); K.strip(GX, yg, GW, 3.6, lp_code, SC, alpha=lp_conf)
        K.TX(GX + GW + 1, yg + 1.8, f"{L('G', 'called states of the joint barcode', int(lp_conf.sum()), SRC_B)} called")
        K.TX(GX - 1, yg + 6.3, "super barcode", ha="right", b=True); K.strip(GX, yg + 4.5, GW, 3.6, SUP, SC, hatchset=(10, 12), outline=-2)
        K.TX(GX + GW + 1, yg + 6.3, f"{L('G', 'full votes', full, SRC_SB)} + {L('G', 'partial votes', part, SRC_SB)} partial", b=True)
        K.numbers(GX, GW, yg + 9.6, n)
        for k_, (c, ht, lab, ec) in enumerate(((C0, None, ST_LAB[0], S.NAVY), (C1, None, ST_LAB[1], S.NAVY), (C2, None, ST_LAB[2], S.NAVY), (C0, "//////", "≥ 1 class 0", "none"),
                                               (C2, "//////", "≥ 1 class 2", "none"), ("#FFFFFF", None, "opposing votes", S.CONFLICT), (NOI, None, "no vote", S.NAVY))):
            K.SW(63 + (k_ % 4) * 27, yg + 13 + (k_ // 4) * 2.6, c, lab, ht, ec=ec)
        K.SW(63 + 3 * 27, yg + 15.6, C0, "pale, not called (joint barcode)", alpha=.38)
        b2_tail = f" and agrees with the joint barcode of J_{jf}" if both.sum() and agree == both.sum() else f"; the vote agrees with the joint barcode of J_{jf} at {agree} of {int(both.sum())} intervals"
    # ---- banner 2
    if rowsE:
        ext = " and ".join(f"{a} → {b}" for a, b in rowsE)
        where = (f"the links of {jf}" if has_G else ("link " if len(e_links) == 1 else "links ") + " and ".join(_xname(p) for p in e_links))
        K.BAN(Y2, "2", f"{kl1} evidence extends the barcode of P ({ext} of {n} intervals in {where}){b2_tail}")
    # ---- banner 3, H: evidence per SV class
    if multi:
        ref = min(dms, key=lambda d: abs(d - DM0)); A = SWt[np.isclose(SWt.delta_min.astype(float), ref)]
        singles = [s for s in sets if "+" not in s]; comb = [s for s in sets if "+" in s]
        def nk(s):
            g = A[A.evidence_set.astype(str) == s].sort_values("k"); return [int(x) for x in g.n_homozygous_informative]
        K.BAN(Y3, "3", "Evidence per SV class at tag strength ≥ " + f"{ref:g}: " + "; ".join(f"{S.class_label(s)} inform {' / '.join(str(x) for x in nk(s))} homozygous calls" for s in singles) + f" at k = {' / '.join(str(k) for k in ks)}")
        K.PT("H", ", ".join(S.class_label(s, plural=False) for s in singles) + " and combined evidence", .5, Y3 + 4.4, f"tag strength ≥ {ref:g}; n = informative cells"); ax = K.AX(9, Y3 + 10, 48, 18)
        same = False
        if comb and main in singles:
            ga = A[A.evidence_set.astype(str) == comb[0]].sort_values("k"); gb = A[A.evidence_set.astype(str) == main].sort_values("k")
            same = len(ga) == len(gb) and np.allclose(ga.acc_homozygous.fillna(-1), gb.acc_homozygous.fillna(-1)) and list(ga.n_homozygous_informative) == list(gb.n_homozygous_informative)
        items = []
        for q, s in enumerate(singles + comb):
            col = S.NAVY if "+" in s else S.class_colour(s); ls_ = (0, (2, 1.5)) if "+" in s else "-"; lw_ = .4 if "+" in s else (.8 if s == main else .5)
            g = A[(A.evidence_set.astype(str) == s) & (A.n_homozygous_informative > 0)].sort_values("k")
            ax.plot(g.k, g.acc_homozygous, marker="o", ms=.9 if "+" in s else 1.8, mew=0, lw=lw_, ls=ls_, color=col)
            if "+" not in s:
                first = s == main
                for kk, aa, nn in zip(g.k, g.acc_homozygous, g.n_homozygous_informative):
                    ax.text(kk + (0 if first else .09), aa + (-.09 if first else .02), f"n = {L('H', f'{s} k={int(kk)} informative homozygous calls', int(nn), SRC_SW)}", fontsize=FS,
                            ha="center" if first else "left", va="top" if first else "bottom", color=S.NAVY)
            items.append((col, ls_, " + ".join(S.class_label(x, plural=False) for x in s.split("+")) + (f" (= {kl1})" if "+" in s and same else "")))
        gD = A[(A.evidence_set.astype(str) == main) & (A.n_homozygous_informative > 0)].sort_values("k")
        ax.plot(gD.k, gD.null_mean, ls=(0, (1, 1)), lw=.4, color="#6B7280"); items.append(("#6B7280", (0, (1, 1)), f"permutation null ({kl1})"))
        ax.set_xticks(ks); ax.set_xlabel("minimum tagging loci (k)", labelpad=.8, fontsize=FS); ax.set_ylabel("match rate", labelpad=.8, fontsize=FS); ax.set_ylim(0, 1.2); ax.set_yticks([0, .5, 1]); ax.set_xlim(min(ks) - .2, max(ks) + .2)
        for k_, (c, ls_, lab) in enumerate(items):
            ov.plot([60, 62.5], [Y3 + 13 + k_ * 2.4] * 2, color=c, lw=.6, ls=ls_); K.TX(63.2, Y3 + 13 + k_ * 2.4, lab)
    title = f"Barcode of a shared parent with {kl1} augmentation"
    desc = (f"A arrangement-tag {kl} per interval (class 0- and class 2-tagged; dot, tag strength at least {DM0:g}); B tag strength; C {kl1} evidence class against the barcode state of the shared parent P, pooled over links "
            f"(k = {K0}, tag strength at least {DM0:g}; permutation of the evidence across intervals inside each link); D match rate on homozygous calls over k and tag strength, with the permutation null; "
            "E barcode, evidence and augmented barcode of the links shown; F intervals with a barcode state and intervals added by the evidence, per link; G vote of the links of a junction family next to its joint barcode; "
            "H match rate per SV class and combined. Classes 0 / 2 are the two homokaryotype classes (no ancestral meaning).")
    return cv.fig, title, desc, "; ".join(notes)


# ====================================================================== figure: funnel to an unsampled parent
def _bands(mn, word):
    """Size bands of the smaller unit: (label, lowest, highest)."""
    w1 = word if word == "fish" else (word[:-1] if word.endswith("s") else word)
    out = [(f"≥ {mn} {word}, auditable", mn, 10 ** 9)]
    if mn - 1 >= 2:
        if mn - 1 <= 3:
            out += [(f"{s} {word}, provisional", s, s) for s in range(mn - 1, 1, -1)]
        else:
            out.append((f"2–{mn - 1} {word}, provisional", 2, mn - 1))
    out.append((f"single {w1}, not testable", 0, 1))
    return out


def fig_funnel(run, plt):
    """A family size by closure, B closure of each family, C nomination fraction of every unit pair,
    D funnel by size of the smaller unit, E sharing of two SV classes per compatible link, F audit
    outcome, G parental barcode of one link, H the barcodes of two links of a junction family and
    its joint barcode, I pseudo-chain null of the recurrence test."""
    if "families" not in run.t.get("families", {}) or "pairs" not in run.t.get("nomination", {}):
        return None
    F = run.t["families"]["families"]; M = run.t["families"]["membership"]
    if not len(F):
        return "no family"
    E = edges(run); aud = audited(run); mn = int(run.cfg["audit"]["min_family_size"]); word = noun(run); notes = []
    q_max = run.cfg["nomination"]["q_max"]
    nomi = E.nominated.astype(bool).to_numpy()
    fam_pair = (E.unassigned_sides.astype(int) == 0).to_numpy() if "unassigned_sides" in E else np.ones(len(E), bool)
    min_n = E.min_n.astype(int).to_numpy()
    st = np.array([status_of(x) for x in E.link_class]) if "link_class" in E else np.array(["nominated, not audited"] * len(E))
    g1 = (E.gate1.fillna(0).astype(float) == 1).to_numpy() if aud and "gate1" in E else np.zeros(len(E), bool)
    Tt = nomi & fam_pair & (min_n >= 2)                          # nominated pairs of families
    OH = Tt & g1; SUPP = Tt & (st == "supported")
    g2 = run.t.get("audit", {}).get("gate2") if aud else None
    g2c = [str(c) for c in g2["class"] if f"share_{c}" in E.columns] if g2 is not None and len(g2) else []
    has_EF = aud and OH.any()
    # ---- barcode part
    bt = run.t.get("barcode", {}); has_bar = "barcode" in bt and len(bt["barcode"]) > 0
    names, V, links, jf, jl = [], {}, [], None, []
    if has_bar:
        names, V, links, _j = _vectors(run); jf, jl = choose_junction(run, V, links)
    n = len(names)
    g_link = None
    if has_bar and links:
        sup = [p for p in links if V[p]["status"].startswith("supported")]
        g_link = jl[0] if jf else (max(sup or links, key=lambda p: (int((V[p]["call"] >= 0).sum()), -links.index(p))))
    RC = bt.get("recurrence") if has_bar else None
    rc = None
    if jf and RC is not None and len(RC):
        r_ = RC[(RC.junction.astype(str) == jf) & (RC.verdict.astype(str) != "not tested")]
        if len(r_) and {str(r_.iloc[0].partner_1), str(r_.iloc[0].partner_2)} == {x for p in jl[:2] for x in V[p]["families"]} - {jf}:
            rc = r_.iloc[0]
    # ---- geometry (sections close up when their input is absent)
    W = 183.0; Y1 = 2.2; Y2 = 33.5
    Y3 = Y2 + 36.0
    Y4 = Y3 + 45.0 if has_EF else Y3
    yG = Y4 + 6.6; ykG = yG + 3 * 3.1 + 3
    Y5 = (ykG + 4) if g_link else Y4
    H = (Y5 + 23.0) if jf else Y5 + (1.0 if (g_link or has_EF) else 0.5)
    cv = S.canvas(plt, W, H, S.SZ_ED, thin=True); K = _Kit(cv); ov = K.ov
    CL = {"closed": "#2A9D8F", "not closed": "#A7AEB3"}; VC = "#2F7FC1"
    # ---- 1: families
    asg = M.family.astype(str) != "unassigned"; n_ind, n_fam, n_asg, n_un = len(M), len(F), int(asg.sum()), int((~asg).sum())
    stF = F.status.astype(str); has_cl = stF.isin(list(CL)).any() and F.closure.notna().any()
    closed_f = set(F.family.astype(str)[stF == "closed"]); n_closed = len(closed_f); ind_closed = int(M.family.astype(str).isin(closed_f).sum())
    K.BAN(Y1, "1", f"{L('1', 'individuals in the analysis set', n_ind, SRC_M)} {word} in the analysis set: {L('1', 'full-sib families', n_fam, SRC_F)} full-sib families ({L('1', 'individuals in families', n_asg, SRC_M)} {word})"
          + (f", {L('1', 'closed families', n_closed, SRC_F)} closed ({L('1', 'individuals in closed families', ind_closed, SRC_M)} {word})" if has_cl else "") + f", {L('1', 'unassigned individuals', n_un, SRC_M)} {word} unassigned")
    ST = ["closed", "not closed"] if has_cl else sorted(set(stF))
    colF = {s: CL.get(s, "#A7AEB3") for s in ST}
    nF = F.n.astype(int).to_numpy(); sizes = np.arange(2, max(nF.max(), 2) + 1)
    K.PT("A", "Family size by closure outcome" if has_cl else "Family size", .5, 6.8); ax = K.AX(9, 10, 62, 17); bot = np.zeros(len(sizes))
    for s in ST:
        h_ = np.array([int(((nF == s_) & (stF == s)).sum()) for s_ in sizes]); ax.bar(sizes, h_, bottom=bot, color=colF[s], width=.8, lw=0); bot += h_
    for s_, b_ in zip(sizes, bot):
        if b_ > 0:
            ax.text(s_, b_ + .4, str(L("A", f"families of {s_}", int(b_), SRC_F)), ha="center", va="bottom", fontsize=FS)
    stp = 1 if len(sizes) <= 16 else int(np.ceil(len(sizes) / 12))
    ax.set_xticks(sizes[::stp]); ax.set_ylim(0, max(bot.max(), 1) * 1.2); ax.set_xlabel(f"{word} per family", labelpad=.8, fontsize=FS); ax.set_ylabel("families", labelpad=.8, fontsize=FS)
    if has_cl:
        K.PT("B", "Closure of each family", 76, 6.8); ax = K.AX(87, 10, 50, 17); rng = np.random.default_rng(1); clo = F.closure.astype(float).to_numpy()
        for s in ST:
            m = (stF == s).to_numpy() & ~np.isnan(clo); ax.scatter(nF[m] + rng.normal(0, .1, int(m.sum())), clo[m], s=2.6, c=colF[s], lw=0, zorder=3)
        ax.set_xlabel(f"{word} per family", labelpad=.8, fontsize=FS); ax.set_ylabel("closure", labelpad=.8, fontsize=FS); ax.set_xticks(sizes[::max(2, stp)]); ax.set_ylim(min(.0, np.nanmin(clo)) - .02, 1.05)
        notes.append("panels A and B: two closure classes (closed / not closed); the published figure split 'not closed' into unresolved and incompatible by a rule that is not part of the tool")
    else:
        notes.append("panel B omitted and panel A without closure classes: the families carry no closure (imported labels without a relatedness table)")
    yy = 12
    for s in ST:
        K.SW(141, yy, colF[s], f"{s}: {L('A,B', f'families {s}', int((stF == s).sum()), SRC_F)}"); yy += 2.6
    if has_cl:
        K.TX(141, yy + .6, f"the other {L('A,B', 'assigned individuals outside closed families', n_asg - ind_closed, SRC_M)} assigned {word} are", c=GREYT); K.TX(141, yy + 3, "not necessarily misassigned", c=GREYT)
    # ---- 2: nomination and funnel
    n_pairs, n_nom, n_units = len(E), int(nomi.sum()), int(M.unit.nunique()); n_oh, n_val = int(OH.sum()), int(SUPP.sum())
    b2 = f"{L('2', 'unit pairs', n_pairs, SRC_E):,} pairs of {L('2', 'units', n_units, SRC_M)} units: {L('2', 'nominated', n_nom, SRC_E)} nominated at BH q < {q_max:g}"
    if aud:
        b2 += f", {L('2', 'nominated family pairs compatible at Gate 1', n_oh, SRC_E)} family pairs compatible at Gate 1 (opposite homozygotes), {L('2', 'supported', n_val, SRC_E)} supported"
    else:
        b2 += "; no audit was run (links are nominated, not audited)"; notes.append("panels E and F omitted and the funnel shows nominations only: no audit was run")
    K.BAN(Y2, "2", b2)
    K.PT("C", "Half-sib fraction of every unit pair", .5, Y2 + 4.3); ax = K.AX(9, Y2 + 7, 62, 21)
    xc = E.n_cross_pairs.astype(float).to_numpy(); fr = E.fraction.astype(float).to_numpy()
    ax.scatter(xc[~nomi], fr[~nomi], s=.5, c="#A7AEB3", lw=0, zorder=2, rasterized=True); ax.scatter(xc[nomi], fr[nomi], s=2.2, c=VC, lw=0, zorder=3); ax.set_xscale("log")
    tk = [t for t in (1, 10, 100, 1000, 10000) if t <= max(xc.max(), 10) * 1.5]; ax.set_xticks(tk); ax.set_xticklabels([str(t) for t in tk], fontsize=FS); ax.minorticks_off()
    ax.set_xlabel(f"cross-unit {noun(run, 1)} pairs", labelpad=.8, fontsize=FS); ax.set_ylabel("second-degree fraction", labelpad=.8, fontsize=FS); ax.set_ylim(-.03, 1.05)
    K.SW(49, Y2 + 4.3, VC, "nominated"); K.SW(62, Y2 + 4.3, "#A7AEB3", "other")
    K.PT("D", "Funnel by size of the smaller unit", 76, Y2 + 4.3); ax = K.AX(107, Y2 + 7, 42, 21); band = _bands(mn, word); ys = np.arange(len(band))[::-1]; hh = .25
    FC = (("nominated", "#BFD6EA"), ("Gate-1 compatible", "#6FA6D6"), ("supported", S.NAVY)) if aud else (("nominated", "#BFD6EA"),)
    vmax = 1
    for (lab, lo, hi), y_ in zip(band, ys):
        inb = (min_n >= lo) & (min_n <= hi)
        vals = [int((nomi & inb).sum())] + ([int((OH & inb).sum()), int((SUPP & inb).sum())] if aud and lo >= 2 else [])
        for k_, val in enumerate(vals):
            ax.barh(y_ + (1 - k_) * hh, val, height=hh * .9, color=FC[k_][1], lw=0); vmax = max(vmax, val)
            ax.text(val, y_ + (1 - k_) * hh, " " + str(L("D", f"{lab}: {FC[k_][0]}", val, SRC_E)), va="center", fontsize=FS)
    ax.set_yticks(ys); ax.set_yticklabels([b_[0] for b_ in band], fontsize=FS); ax.set_xlabel("unit pairs", labelpad=.8, fontsize=FS); ax.tick_params(axis="y", length=0); ax.set_xlim(0, vmax * 1.17)
    yy = Y2 + 9
    for lab, c in FC:
        K.SW(155, yy, c, lab); yy += 2.6
    # ---- 3: audit
    if has_EF:
        aud_nom = nomi & E.auditable.astype(bool).to_numpy(); n_aud = int(aud_nom.sum()); n_aud_ok = int((aud_nom & g1).sum())
        K.BAN(Y3, "3", f"Structural variants not used to nominate support the links: {L('3', 'supported links', n_val, SRC_E)} supported; "
              f"{L('3', 'compatible links without support', n_oh - n_val, SRC_E)} of {n_oh} compatible links remain unresolved")
        if g2c:
            kx = g2c[0]; ky = g2c[1] if len(g2c) > 1 else None; row = {str(r["class"]): r for _i, r in g2.iterrows()}
            lx = S.class_label(kx.split("_")[0], plural=False); ly = S.class_label(ky.split("_")[0], plural=False) if ky else "second-degree fraction"
            K.PT("E", f"{lx} and {ly} sharing per compatible link" if ky else f"{lx} sharing per compatible link", .5, Y3 + 4.3, f"dashed, median of the n ≥ {mn} compatible pairs" + ("; × matched null" if ky else ""))
            ax = K.AX(9, Y3 + 10, 66, 24)
            X = E[f"share_{kx}"].astype(float).to_numpy(); Yv = E[f"share_{ky}"].astype(float).to_numpy() if ky else fr
            okp = OH & ~np.isnan(X) & ~np.isnan(Yv); o_ = okp & ~SUPP; s_ = okp & SUPP
            ax.scatter(X[o_], Yv[o_], s=1.8, c="#A7AEB3", lw=0, zorder=2); ax.scatter(X[s_], Yv[s_], s=4, c=S.NAVY, lw=0, zorder=3)
            tx = float(row[kx]["threshold_median_of_gate1_pass"]); ax.axvline(L("E", f"support threshold {kx}", tx, SRC_G2), color=S.NAVY, ls=(0, (3, 2)), lw=.35)
            if ky:
                ty = float(row[ky]["threshold_median_of_gate1_pass"]); ax.axhline(L("E", f"support threshold {ky}", ty, SRC_G2), color=S.NAVY, ls=(0, (3, 2)), lw=.35)
                nx_, ny_ = float(row[kx]["null_median"]), float(row[ky]["null_median"])
                if nx_ == nx_ and ny_ == ny_:
                    ax.scatter([L("E", f"matched-null median {kx}", nx_, SRC_G2)], [L("E", f"matched-null median {ky}", ny_, SRC_G2)], marker="x", s=7, c="#D9482F", lw=.5, zorder=5)
            x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); y1n = y1 + (y1 - y0) * .34; ax.set_ylim(y0, y1n)
            idx = np.flatnonzero(s_); idx = idx[np.argsort(X[idx])]
            if len(idx) > 8:
                notes.append(f"panel E: the 8 supported links with the highest {lx} sharing are named ({len(idx)} supported)"); idx = idx[-8:]
            if len(idx):
                ua, ub = E.unit_a.astype(str).to_numpy(), E.unit_b.astype(str).to_numpy()
                pitch = (x1 - x0) * .095; xr = X[idx].mean() + (np.arange(len(idx)) - (len(idx) - 1) / 2) * pitch
                xr += min(0, x1 - pitch * .55 - xr.max()) if len(idx) else 0; xr += max(0, x0 + pitch * .55 - xr.min())
                for q, i in enumerate(idx):
                    ax.annotate(f"{ua[i]}×{ub[i]}", (X[i], Yv[i]), xytext=(xr[q], y1 + (y1n - y1) * (.12 if q % 2 == 0 else .56)), fontsize=FS, ha="center", va="bottom", zorder=4,
                                arrowprops=dict(arrowstyle="-", lw=.15, color="#8A949C", shrinkA=0, shrinkB=.3))
            ax.set_xlabel(f"{lx} sharing", labelpad=.8, fontsize=FS); ax.set_ylabel(f"{ly} sharing" if ky else ly, labelpad=.8, fontsize=FS)
            K.SW(1, Y3 + 41, S.NAVY, "supported"); K.SW(15, Y3 + 41, "#A7AEB3", "compatible, unresolved")
            if ky:
                ov.scatter([44], [Y3 + 41], marker="x", s=7, c="#D9482F", lw=.5); K.TX(45.5, Y3 + 41, "matched null (medians)")
            if len(g2c) > 2:
                notes.append(f"panel E shows the first two of {len(g2c)} SV classes of Gate 2")
            if not ky:
                notes.append("panel E: one SV class in Gate 2; the vertical axis is the nomination fraction")
        else:
            notes.append("panel E omitted: no SV class in the audit (no Gate 2)")
        K.PT("F", f"Audit of the pairs of at least {mn} {word}", 80, Y3 + 4.3, f"{L('F', 'auditable nominated pairs', n_aud, SRC_E)} auditable, {L('F', 'auditable pairs compatible at Gate 1', n_aud_ok, SRC_E)} compatible; an exclusion test, not proof")
        ev = E.evidence.fillna("").astype(str).to_numpy() if "evidence" in E else np.array([""] * len(E))
        cats = []
        for e in sorted({x for x in ev[aud_nom & (st == "supported")]}, key=lambda x: (-x.count("+"), x)):
            cats.append((f"{e} support" if e else "supported", int((aud_nom & (st == "supported") & (ev == e)).sum()), S.NAVY))
        cats.append(("passed Gate 1, no SV support", int((aud_nom & (st == "Gate-1 only")).sum()), "#A7AEB3")); cats.append(("failed Gate 1", int((aud_nom & (st == "rejected")).sum()), "#D9482F"))
        other = n_aud - sum(c[1] for c in cats)
        if other:
            cats.append(("other", other, NOI))
        ax = K.AX(111, Y3 + 10, 30, 24); yv = np.arange(len(cats))[::-1]
        ax.barh(yv, [c[1] for c in cats], color=[c[2] for c in cats], height=.7, lw=0)
        for y_, c in zip(yv, cats):
            ax.text(c[1], y_, " " + str(L("F", c[0], c[1], SRC_E + " link_class, evidence")), va="center", fontsize=FS)
        ax.set_yticks(yv); ax.set_yticklabels([c[0] for c in cats], fontsize=FS); ax.tick_params(axis="y", length=0); ax.set_xlabel("pairs", labelpad=.8, fontsize=FS); ax.set_xlim(0, max(max(c[1] for c in cats), 1) * 1.3)
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
        if g2c:
            K.TX(146, Y3 + 11, "observed vs matched null", b=True); yy = Y3 + 14
            for c in g2c[:3]:
                r = row[c]; lab = S.class_label(c.split("_")[0], plural=False)
                o_, nl_, p_ = float(r["threshold_median_of_gate1_pass"]), float(r["null_median"]), float(r["p_one_sided"])
                K.TX(146, yy, f"{lab} {L('F', f'{c} observed median', o_, SRC_G2):.3f} vs " + (f"{L('F', f'{c} matched-null median', nl_, SRC_G2):.3f}" if nl_ == nl_ else "no null"))
                if p_ == p_:
                    K.TX(146, yy + 2.4, f"P = {L('F', f'{c} one-sided P', p_, SRC_G2):.4f}", c=GREYT)
                yy += 5.4
    elif aud:
        notes.append("panels E and F omitted: no nominated family pair passed Gate 1")
    # ---- 4: the three parental vectors of one link
    SC = {0: C0, 1: C1, 2: C2, -1: NOI}
    if g_link:
        fa, fb = V[g_link]["families"][:2]
        K.BAN(Y4, "4", f"One link constrains three parental vectors over the {L('4', 'barcode intervals', n, SRC_B)} {hap_label(run)} ({_xname(g_link)}): three vectors of one link, not three latent parents")
        K.PT("G", "Parental barcode", .5, Y4 + 4.3); XS, WS = 25, 140; y_ = yG
        for lab, col in (("P, shared", "call"), (f"QA, {fa} private", "private_parent_a_call"), (f"QB, {fb} private", "private_parent_b_call")):
            v_ = V[g_link].get(col, np.full(n, -1))
            K.TX(XS - 1, y_ + 1.3, lab, ha="right", b=lab.startswith("P")); K.strip(XS, y_, WS, 2.6, v_, SC); K.TX(XS + WS + 1.5, y_ + 1.3, f"{L('G', f'{g_link} {lab}: called states', int((v_ >= 0).sum()), SRC_B)} / {n}"); y_ += 3.1
        K.numbers(XS, WS, y_ + .6, n)
        for k_, (g, lab) in enumerate(((0, ST_LAB[0]), (1, ST_LAB[1]), (2, ST_LAB[2]), (-1, "unresolved"))):
            K.SW(XS + k_ * 22, ykG, SC[g], lab, ec=S.NAVY if g < 0 else "none")
        if jf:
            nj = len(run.s.get("network", {}).get("junction_families", {}) or {})
            notes.append(f"panels G-I show junction family {jf}" + (f" (of {nj}: the one with the most supported links, ties by the number of called states of its joint barcode)" if nj > 1 else "") + "; G is its first supported link in the link table, H its first two")
        else:
            notes.append(f"panel G shows the {'supported ' if V[g_link]['status'].startswith('supported') else ''}link with the most called states; panels H and I omitted: no junction family with two supported links")
    else:
        notes.append("panels G-I omitted: no barcode (no link between families of the minimum size, or no structural-haplotype genotypes)")
    # ---- 5: the same barcode on two links of a junction family
    if jf:
        pa, pb = jl[0], jl[1]; P1, P4 = V[pa]["call"], V[pb]["call"]; both = (P1 >= 0) & (P4 >= 0); nmatch, nboth = int((P1[both] == P4[both]).sum()), int(both.sum())
        b5 = f"{jf} bridges two links: {L('5', 'matching jointly resolved intervals', nmatch, SRC_B)} of {L('5', 'jointly resolved intervals', nboth, SRC_B)} jointly resolved intervals match"
        if rc is not None and rc.n_null == rc.n_null and float(rc.n_null) > 0:
            nn_ = int(rc.n_null); beat = float(rc.statistic) > float(rc.null_max)
            b5 += ("; slot agreement beats all " if beat else "; slot agreement against ") + f"{L('5', 'pseudo-chains', nn_, SRC_RC)} pseudo-chains (P = {L('5', 'empirical P', float(rc.P), SRC_RC):.3f}): " + \
                  (f"J_{jf} is a candidate" if str(rc.verdict) == "candidate" else "recurrence not supported")
        K.BAN(Y5, "5", b5)
        K.PT("H", "Does the same barcode recur?", .5, Y5 + 4.3); XS2, WS2 = 25, 105; y_ = Y5 + 6.6
        K.TX(XS2 - 1, y_ + 1.3, f"P, {_xname(pa)}", ha="right"); K.strip(XS2, y_, WS2, 2.6, P1, SC); K.TX(XS2 - 1, y_ + 5.3, f"P, {_xname(pb)}", ha="right"); K.strip(XS2, y_ + 4, WS2, 2.6, P4, SC)
        for i in np.where(both)[0]:
            ov.plot([XS2 + (i + .5) * WS2 / n] * 2, [y_ + 2.7, y_ + 3.9], color=S.NAVY if P1[i] == P4[i] else "#D9482F", lw=.45, zorder=7)
        lp_code, lp_conf = V["J_" + jf]["best_genotype"], V["J_" + jf]["call"] >= 0
        K.TX(XS2 - 1, y_ + 9.3, f"J_{jf}", ha="right", b=True); K.strip(XS2, y_ + 8, WS2, 2.6, lp_code, SC, alpha=lp_conf); K.numbers(XS2, WS2, y_ + 11.8, n)
        cnt = [int((lp_conf & (lp_code == g)).sum()) for g in range(3)]
        K.TX(XS2, y_ + 14.2, f"J_{jf}: {L('H', 'called states of the joint barcode', int(lp_conf.sum()), SRC_B)} called = {L('H', 'called class 0 hom.', cnt[0], SRC_B)} {ST_LAB[0]}, {L('H', 'called het.', cnt[1], SRC_B)} {ST_LAB[1]}, "
              f"{L('H', 'called class 2 hom.', cnt[2], SRC_B)} {ST_LAB[2]} (pale, not called); ticks, jointly resolved: navy match" + (", red mismatch" if nmatch < nboth else ""))
        if len(jl) > 2:
            notes.append(f"panel H shows the first two of the {len(jl)} supported links of {jf}; the joint barcode J_{jf} uses all of them")
        if rc is not None and rc.null_mean == rc.null_mean:
            stt, nmean, nmax = float(rc.statistic), float(rc.null_mean), float(rc.null_max)
            K.PT("I", "Pseudo-chain null", 134, Y5 + 4.3, "matched-slot agreement; grey line,\nnull mean to maximum"); ax = K.AX(137.5, Y5 + 11, 42, 4); ax.set_xlim(0, max(stt, nmax) * 1.125); ax.set_ylim(0, 1); ax.set_yticks([]); ax.spines["left"].set_visible(False)
            ax.plot([nmean, nmax], [.35, .35], color="#A7AEB3", lw=.8); ax.plot([nmean], [.35], "o", ms=2, mew=0, color="#6B7280"); ax.plot([nmax], [.35], "|", ms=4, mew=.5, color="#6B7280"); ax.plot([stt], [.35], "D", ms=2.4, mew=0, color=S.NAVY)
            from matplotlib.ticker import MaxNLocator
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ov.plot([138], [Y5 + 19], "D", ms=2.4, mew=0, color=S.NAVY); K.TX(139.5, Y5 + 19, f"observed {L('I', 'statistic', stt, SRC_RC):.3f}")
            ov.plot([159], [Y5 + 19], "o", ms=2, mew=0, color="#6B7280"); K.TX(160.5, Y5 + 19, f"null mean {L('I', 'null mean', nmean, SRC_RC):.3f}")
            ov.plot([138], [Y5 + 21.4], "|", ms=4, mew=.5, color="#6B7280"); K.TX(139.5, Y5 + 21.4, f"null max {L('I', 'null maximum', nmax, SRC_RC):.3f}; {int(rc.n_null)} pseudo-chains")
        else:
            notes.append("panel I omitted: the recurrence test of the junction family was not run (no deletion input, too few rare deletions, or no family pair for the null)")
    title = "Funnel from offspring relationships to an unsampled parent"
    desc = ("A family sizes by closure; B closure of each family; C fraction of cross-unit pairs at second degree or closer for every unit pair (nominated pairs in blue); D nominated, Gate-1-compatible and supported pairs by size of the smaller unit; "
            "E sharing of the SV classes of Gate 2 per Gate-1-compatible nominated family pair, with the support thresholds (dashed) and the matched-null medians; F audit outcome of the auditable nominated pairs; "
            "G barcode of the shared parent P and of the two private parents of one link; H barcodes of two links of a junction family, ticks at jointly resolved intervals, and the joint barcode of the junction; "
            "I recurrence statistic against pseudo-chains of families not linked to the junction. Classes 0 / 2 are the two homokaryotype classes (no ancestral meaning).")
    return cv.fig, title, desc, "; ".join(notes)
