"""Shared-parent figure (nomination, gates, audited network) and the matrix of all units.

Ports of the published figure scripts `build_fig3_sharedparent_v13.py` (panels A-J) and
`build_suppfig_matrix79_v1.py`. Same panels, encodings and style; every family name, count
and threshold comes from the tables of the run. Differences from the originals that follow
from being generic are listed in docs/DISPLAY_ITEMS.md.
"""
from __future__ import annotations
import itertools

import numpy as np

from ..audit import matched_null
from ..network import status_of
from . import style as S

SYM_BASE = {"Gate-1 only": ("white", "o", 11, "#1B1B2F"), "rejected": (S.RED, "X", 16, None), "provisional": ("#A7AEB3", "^", 6.5, "#555555"),
            "nominated, not audited": ("#A7AEB3", "D", 6, "#555555")}
PROV = ("#AEB8C0", .6, (0, (2, 2)))
CCOL = ["#283B7C", "#D08C6E", "#71A29C", "#B8A9CC", "#C9A227", "#83251D"]


# ------------------------------------------------------------------ data helpers
def edges(run):
    a = run.t.get("audit", {})
    return a["edges"] if "edges" in a else run.t["nomination"]["pairs"]


def audited(run):
    return bool(run.s.get("audit", {}).get("ran"))


def family_sizes(run):
    F = run.t["families"]["families"]
    return dict(zip(F.family.astype(str), F.n.astype(int)))


def noun(run, n=2):
    w = run.cfg["figures"].get("individual_noun", "individuals")
    return w if n != 1 else (w[:-1] if w.endswith("s") and w != "fish" else w)


def evidence_classes(r):
    ev = str(getattr(r, "evidence", "") or "")
    return [x.strip() for x in ev.split("+") if x.strip() and x.strip() != "nan"]


def link_key(r):
    """(status, evidence label) of a nominated pair; evidence label is 'DEL + DUP', 'DEL', ... for supported links."""
    st = status_of(getattr(r, "link_class", "nominated, not audited"))
    return st, " + ".join(evidence_classes(r)) if st == "supported" else ""


def evidence_style(ev):
    """Colour of a supported link by the SV classes that support it: two or more = navy."""
    k = [x.split("_")[0] for x in ev.split(" + ") if x]
    if len(k) >= 2:
        return S.ROH_BASE["navy"]
    return S.class_colour(k[0]) if k else S.ROH_BASE["navy"]


def symbol(st, ev, scale=1.0):
    if st == "supported":
        return evidence_style(ev), "s", (14 if " + " in ev else 13) * scale, None
    c, m, s, e = SYM_BASE[st]
    return c, m, s * scale, e


def order_linked_first(fams, sizes, pairs, optimise=True):
    """Families with a nominated link first (largest group first; linked families adjacent inside
    a group when the group is small enough to try every order), then the others by size."""
    par = {f: f for f in fams}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in pairs:
        par[find(a)] = find(b)
    groups = {}
    for f in fams:
        groups.setdefault(find(f), []).append(f)
    fkey = lambda f: (-sizes[f], S.fam_key(f))
    comps = sorted([g for g in groups.values() if len(g) > 1], key=lambda g: (-len(g), -max(sizes[x] for x in g), min(S.fam_key(x) for x in g)))
    iso = sorted([g[0] for g in groups.values() if len(g) == 1], key=fkey)

    def best(g):
        g0 = sorted(g, key=fkey)
        if not optimise or len(g0) > 8:
            return g0
        eg = [(a, b) for a, b in pairs if a in g and b in g]
        return list(min(itertools.permutations(g0), key=lambda p: (sum(abs(p.index(a) - p.index(b)) for a, b in eg), [fkey(x) for x in p])))
    fl = [f for g in comps for f in best(g)] + iso
    return fl, sum(len(g) for g in comps)


# ------------------------------------------------------------------ panel A / all-units matrix
def _matrix_core(ax, E, fl, sizes, n_comp, size_min_label, sym_scale, tick_fs, full_dotted):
    from matplotlib.patches import Rectangle
    pos = {f: i for i, f in enumerate(fl)}; nf = len(fl)
    G = np.full((nf, nf), np.nan)
    sub = E[E.unit_a.isin(pos) & E.unit_b.isin(pos)]
    for a, b, fr in zip(sub.unit_a, sub.unit_b, sub.fraction):
        i, j = pos[a], pos[b]; G[max(i, j), min(i, j)] = fr
    cm = S.nomination_cmap()
    ax.imshow(np.ma.masked_invalid(G), cmap=cm, vmin=0, vmax=1, aspect="equal", interpolation="nearest", rasterized=True, extent=(-.5, nf - .5, nf - .5, -.5))
    for k in range(nf + 1):
        ax.plot([-.5, nf - .5], [k - .5] * 2, color="#E3E8EC", lw=.2, zorder=1); ax.plot([k - .5] * 2, [-.5, nf - .5], color="#E3E8EC", lw=.2, zorder=1)
    for i, f in enumerate(fl):
        g_ = max(.30, .90 - .055 * (sizes[f] - 2)); ax.add_patch(Rectangle((i - .5, i - .5), 1, 1, fc=(g_, g_, g_), ec="#7a7a7a", lw=.2, zorder=3))
        if sizes[f] >= size_min_label and tick_fs >= 3:
            ax.text(i, i + .03, str(sizes[f]), fontsize=tick_fs, color="white" if g_ < .62 else S.NAVY, ha="center", va="center", zorder=4, fontweight="bold")
    SN = sub[sub.nominated]
    cells = []
    for r in SN.itertuples():
        i, j = pos[r.unit_a], pos[r.unit_b]; lo, hi = min(i, j), max(i, j)
        ax.add_patch(Rectangle((lo - .5, hi - .5), 1, 1, fill=False, ec="#1B1B2F", lw=.9, zorder=5))
        st, ev = link_key(r)
        if st in ("supported", "Gate-1 only"):
            ax.add_patch(Rectangle((hi - .5, lo - .5), 1, 1, fc=S.G1_FILL, ec=S.GREEN, lw=.7, zorder=5))
        col, mk, sz_, ec = symbol(st, ev, sym_scale)
        ax.scatter([hi], [lo], s=sz_, c=col, marker=mk, ec=ec if ec else "none", lw=.5, zorder=6)
        cells.append((lo, hi))
    if 0 < n_comp < nf:
        if full_dotted:
            ax.axvline(n_comp - .5, color="#6b6b6b", ls=":", lw=.6); ax.axhline(n_comp - .5, color="#6b6b6b", ls=":", lw=.6)
        else:
            ax.plot([-.5, n_comp - .5], [n_comp - .5] * 2, color="#6b6b6b", ls=":", lw=.6, zorder=2); ax.plot([n_comp - .5] * 2, [n_comp - .5, nf - .5], color="#6b6b6b", ls=":", lw=.6, zorder=2)
    if tick_fs >= 3:
        ax.set_xticks(range(nf)); ax.set_xticklabels(fl, rotation=90, fontsize=tick_fs); ax.set_yticks(range(nf)); ax.set_yticklabels(fl, fontsize=tick_fs)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(length=1.2, pad=.8); ax.set_xlim(-.5, nf - .5); ax.set_ylim(nf - .5, -.5)
    S.frame(ax, .4)
    return SN, cells, cm


def _counts(SN):
    keys = [link_key(r) for r in SN.itertuples()]
    out = {}
    for k in keys:
        out[k] = out.get(k, 0) + 1
    return out


def _failed_names(SN, limit=2):
    bad = [f"{r.unit_a} × {r.unit_b}" for r in SN.itertuples() if link_key(r)[0] == "rejected"]
    return (": " + ", ".join(bad)) if 0 < len(bad) <= limit else ""


def _key_A(kax, cm, KW, KH, kf, SN, n_nom, n_g1, aud, q_max, word, sym_scale):
    """Key of panel A, drawn in a 34.6 x 17.8 mm coordinate box (the v13 inset)."""
    from matplotlib.patches import Rectangle
    kax.set_xlim(0, KW); kax.set_ylim(KH, 0); kax.set_xticks([]); kax.set_yticks([]); kax.set_facecolor("white")
    for sp in kax.spines.values():
        sp.set_visible(True); sp.set_linewidth(.3); sp.set_color("#BBBBBB")
    kt = dict(fontsize=kf, va="center", ha="left"); KY = [v * KH / 20.7 for v in (1.9, 4.4, 6.9, 9.7, 12.2, 14.7, 17.0, 19.3)]
    kax.text(1.2, KY[0], "Below diagonal: nomination", fontweight="bold", **kt)
    kax.imshow(np.linspace(0, 1, 64)[None, :], cmap=cm, aspect="auto", extent=(2.8, 7.8, KY[1] + .8, KY[1] - .8), zorder=2)
    kax.add_patch(Rectangle((2.8, KY[1] - .8), 5.0, 1.6, fill=False, ec=S.NAVY, lw=.3, zorder=3)); kax.set_xlim(0, KW); kax.set_ylim(KH, 0)
    kax.text(1.9, KY[1], "0", fontsize=kf, va="center", ha="center"); kax.text(8.7, KY[1], "1", fontsize=kf, va="center", ha="center"); kax.text(10.2, KY[1], "pair fraction at 2nd-degree θ", **kt)
    kax.add_patch(Rectangle((1.4, KY[2] - 1.0), 2.0, 2.0, fc=cm(.85), ec="#1B1B2F", lw=.9)); kax.text(4.6, KY[2], f"nominated, BH $\\mathit{{q}}$ < {q_max:g} ({n_nom})", **kt)
    cnt = _counts(SN)
    if aud:
        kax.text(1.2, KY[3], "Above diagonal: audit", fontweight="bold", **kt)
        kax.add_patch(Rectangle((1.4, KY[4] - 1.0), 2.0, 2.0, fc=S.G1_FILL, ec=S.GREEN, lw=.7)); kax.text(4.6, KY[4], f"Gate 1 pass ({n_g1})", **kt)
        n_fail = sum(v for (st, _e), v in cnt.items() if st == "rejected")
        kax.scatter([21.6], [KY[4]], s=16 * sym_scale, c=S.RED, marker="X", lw=0); kax.text(23.4, KY[4], f"fail ({n_fail}){_failed_names(SN)}", **kt)
        kax.text(1.2, KY[5], "SV support:", **kt)
        sup = sorted({e for (st, e) in cnt if st == "supported"}, key=lambda e: (-e.count("+"), e))
        items = [(e, ("supported", e)) for e in sup] + [("none", ("Gate-1 only", ""))]
        slots = [(13.4, KY[5]), (2.4, KY[6]), (13.4, KY[6]), (24.4, KY[6])]
        for q, (lab, key) in enumerate(items[:4]):
            col, mk, sz_, ec = symbol(key[0], key[1], sym_scale); xx, yy = slots[q] if len(items) <= 4 else slots[q]
            kax.scatter([xx], [yy], s=sz_ * .85, c=col, marker=mk, ec=ec if ec else "none", lw=.5); kax.text(xx + 1.6, yy, f"{lab} ({cnt.get(key, 0)})", fontsize=kf, va="center", ha="left")
    else:
        kax.text(1.2, KY[3], "Above diagonal: no audit was run", fontweight="bold", **kt)
        col, mk, sz_, ec = symbol("nominated, not audited", "", sym_scale)
        kax.scatter([2.4], [KY[4]], s=sz_, c=col, marker=mk, ec=ec, lw=.5); kax.text(4.6, KY[4], f"nominated, not audited ({n_nom})", **kt)
    kax.add_patch(Rectangle((1.4, KY[7] - 1.0), 2.0, 2.0, fc=(.52, .52, .52), ec="#7a7a7a", lw=.2)); kax.text(2.4, KY[7] + .03, "n", fontsize=5.5, color="white", fontweight="bold", ha="center", va="center")
    kax.text(4.6, KY[7], f"diagonal: {word} per family", **kt)


# ------------------------------------------------------------------ grid panels B-G
def panel_score(ax, run, sz):
    E = edges(run); q_max = run.cfg["nomination"]["q_max"]
    sig = E.nominated.astype(bool); fr = E.fraction; null_ = float(fr[~sig].mean()) if (~sig).any() else float("nan")
    n_units = len(set(E.unit_a) | set(E.unit_b)); bb = np.linspace(0, 1, 26)
    ax.hist(fr, bins=bb, color="#C9D1D8", lw=0, label=f"all {len(E):,} pairs of {n_units} units")
    ax.hist(fr[sig], bins=bb, color="#C2601C", lw=0, label=f"nominated ({int(sig.sum())}, BH $\\mathit{{q}}$ < {q_max:g})")
    ax.set_yscale("log"); ax.set_ylim(.7, max(30000, len(E) * 10))
    if np.isfinite(null_):
        ax.axvline(null_, color=S.RED, ls="--", lw=.8)
        ax.text(null_ + .03, .78, f"null mean {null_:.3f}", transform=ax.get_xaxis_transform(), fontsize=sz["small"] - .4, color=S.RED, fontweight="bold", bbox=dict(fc="white", ec="none", alpha=.9, pad=.8), zorder=10)
    ax.set_xlabel("fraction of cross-unit pairs at 2nd-degree θ"); ax.set_ylabel("unit pairs")
    ax.legend(loc="upper right", fontsize=sz["small"] - .3, handlelength=1, labelspacing=.15, borderaxespad=.1)


def _aud(run):
    E = edges(run)
    return E[E.auditable.astype(bool) & E.LLR_one.notna()]


def panel_gate1(ax, run, sz):
    aud = _aud(run)
    bb = np.unique(np.sort(np.append(np.linspace(aud.LLR_one.min() - .02, aud.LLR_one.max() + .02, 36), 0)))
    ax.hist(aud.LLR_one[aud.LLR_one <= 0], bins=bb, color="#B9BFC9", lw=0, label=f"rejected ({int((aud.LLR_one <= 0).sum())})")
    ax.hist(aud.LLR_one[aud.LLR_one > 0], bins=bb, color=S.GREEN, lw=0, label=f"pass ({int((aud.LLR_one > 0).sum())})")
    ax.set_yscale("log"); ax.axvline(0, color=S.RED, ls="--", lw=.9); ax.set_xlabel("LLR per marker: one shared parent"); ax.set_ylabel("family pairs")
    ax.legend(loc="upper left", fontsize=sz["small"] - .3, handlelength=1, labelspacing=.15, borderaxespad=.1)


def gate2_classes(run):
    E = edges(run)
    return [c[len("share_"):] for c in E.columns if c.startswith("share_") and E[c].notna().any()] if audited(run) else []


def panel_gate2(ax, run, sz, k):
    E = edges(run); col = f"share_{k}"; base = k.split("_")[0]
    vc = S.class_colour(base); vdark = S.mixw(vc, -.3)
    ov, nullv = matched_null(E, col, run.cfg["audit"]["k_match"])
    name = {"DEL": "deletion", "DUP": "duplication"}.get(base, base)
    if len(nullv):
        ax.hist(nullv, bins=18, density=True, color=vc, alpha=.5, lw=0)
    for v in ov:
        ax.axvline(v, color=vdark, lw=.8, alpha=.75)
    txt = ""
    if len(ov):
        om = float(np.median(ov)); ax.axvline(om, color=vdark, lw=1.7); txt = f"observed median {om:.3f}"
        if len(nullv):
            nm = float(np.median(nullv)); ax.axvline(nm, color=S.RED, ls="--", lw=.9)
            P_ = (1 + (nullv >= om).sum()) / (1 + len(nullv))
            txt += f"\nmatched null    {nm:.3f}\n" + S.fmt_p(P_)
    ax.set_xlabel(f"{name} sharing between families"); ax.set_ylabel("density"); ax.set_ylim(0, ax.get_ylim()[1] * 1.55)
    if txt:
        ax.text(.03, .97, txt, transform=ax.transAxes, fontsize=sz["note"], color=S.NOTE_GREY, va="top", bbox=dict(fc="white", ec="none", alpha=1, pad=1.4), zorder=10)


def panel_compat(ax, run, sz):
    from matplotlib.patches import Rectangle
    aud = _aud(run); st = aud.apply(lambda r: status_of(r.link_class) if r.nominated else "", axis=1) if len(aud) else aud.LLR_one
    obs = aud[st.isin(["supported", "Gate-1 only"])]; bad = aud[st == "rejected"]; oth = aud[~aud.index.isin(obs.index) & ~aud.index.isin(bad.index)]
    ax.axvline(0, color="#555555", lw=.6); ax.axhline(0, color="#555555", lw=.6)
    ax.scatter(oth.LLR_one, oth.LLR_both, s=7, marker="s", c="#B9BFC9", ec="#6b6b6b", lw=.2, label=f"other family pairs ({len(oth) + len(bad)})")
    ax.scatter(obs.LLR_one, obs.LLR_both, s=22, c=S.GREEN, ec="#143d24", lw=.4, zorder=4, label=f"half-sib linked ({len(obs)})")
    if len(bad):
        ax.scatter(bad.LLR_one, bad.LLR_both, s=22, marker="X", c=S.RED, zorder=5, label="failed Gate 1")
    ax.set_xlabel("LLR: one shared parent"); ax.set_ylabel("LLR: both parents")
    lo = float(np.nanmin(aud.LLR_both)) if len(aud) else -1.0
    span = max(-lo, 1e-6); ax.set_ylim(-span * 1.5, span * .05)
    ax.legend(loc="lower right", fontsize=sz["small"] - .3, handletextpad=.1, labelspacing=.1, borderaxespad=.2, frameon=True, facecolor="white", edgecolor="none", framealpha=1).set_zorder(12)
    xl_ = ax.get_xlim(); ax.add_patch(Rectangle((0, -span), max(xl_[1], 0), span, fc="#EAF3EA", ec="none", lw=0, zorder=0)); ax.set_xlim(xl_)
    ax.text(.03, .95, "shaded = half-sib:\none parent yes, both no", transform=ax.transAxes, fontsize=sz["small"] - .4, color=S.GREEN, fontweight="bold", va="top", ha="left", zorder=6)


def latent_parents(run):
    LP = run.t.get("network", {}).get("latent_parents")
    if LP is None or not len(LP):
        return None
    t = int(run.cfg["network"]["figure_min_family_size"])
    LP = LP[LP.min_family_size.astype(int) == t]
    return LP if len(LP) else None


def panel_incidence(ax, run, sz, GH):
    from matplotlib.patches import Rectangle
    LP = latent_parents(run)
    fams = sorted({f for x in LP.families for f in str(x).split(",")}, key=S.fam_key); pars = list(LP.latent_parent)
    inc = {f: 0 for f in fams}
    for j, r in enumerate(LP.itertuples()):
        for f in str(r.families).split(","):
            inc[f] += 1; ax.scatter([fams.index(f)], [j], s=14, c=CCOL[(int(r.component) - 1) % len(CCOL)], ec="#1B1B2F", lw=.35, zorder=3)
    for f, n_ in inc.items():
        if n_ > 2:
            i_ = fams.index(f); ax.add_patch(Rectangle((i_ - .46, -.62), .92, len(pars) + .24, fill=False, ec=S.RED, lw=.9, zorder=4))
            ax.text(i_, -1.15, f"{n_} > 2", fontsize=sz["small"] - .3, color=S.RED, fontweight="bold", ha="center", va="bottom", zorder=5)
    fs_x = 5.5 if len(fams) <= 32 else max(3.0, 5.5 * 32 / len(fams)); fs_y = sz["small"] - .3 if len(pars) <= 16 else max(3.0, (sz["small"] - .3) * 16 / len(pars))
    ax.set_yticks(range(len(pars))); ax.set_yticklabels(pars, fontsize=fs_y); ax.set_xticks(range(len(fams))); ax.set_xticklabels(fams, rotation=90, fontsize=fs_x)
    ax.set_ylim(len(pars) - .5, -1.0); ax.set_xlim(-.7, len(fams) - .3); ax.tick_params(length=1, pad=.8); ax.grid(True, lw=.25, color="#E3E8EC"); ax.set_axisbelow(True)
    viol = run.s["network"].get("slot_violations", {})
    tiers = sorted({int(k) for k in viol}, reverse=True)
    get = lambda t: list(viol.get(t, viol.get(str(t), [])))
    _sl = [f"n ≥ {t}: {len(get(t))}" + (f" ({', '.join(get(t)[:6])}{', ...' if len(get(t)) > 6 else ''})" if get(t) else "") for t in tiers]
    lines = ["slot-limit violations (> 2 parents)"] + (["  ·  ".join(_sl[:2])] + _sl[2:3] if len(_sl) > 2 else ["  ·  ".join(_sl)])
    ax.text(.5, -8.1 / GH, "\n".join(lines), fontsize=5.5, va="center", ha="center", color="black", transform=ax.transAxes, linespacing=1.0)


# ------------------------------------------------------------------ H, I
ST_TAB = [("closed", "closed", "#71A29C"), ("not closed", "not closed", "#F1E492"), ("imported", "imported label", "#95BADD"), ("unassigned", "unassigned", "#A7AEB3")]
RL_TAB = [("core", "core member", "#71A29C"), ("attached", "attached", "#95BADD"), ("unresolved member", "unresolved member", "#D08C6E"), ("imported", "imported label", "#C9D1D8"),
          ("unassigned", "unassigned", "#A7AEB3")]


def panel_completeness(cv, run, x, yax, hax, sz):
    """Stacked bars: every individual by block status of its family (left) and by membership class (right)."""
    from matplotlib.patches import Rectangle
    M = run.t["families"]["membership"].copy(); F = run.t["families"]["families"]
    stat = dict(zip(F.family.astype(str), F.status.astype(str)))
    M["bs"] = [stat.get(f, "unassigned") if f != "unassigned" else "unassigned" for f in M.family.astype(str)]
    M["rl"] = M.membership.astype(str)
    ST = [t for t in ST_TAB if (M.bs == t[0]).any()]; RL = [t for t in RL_TAB if (M.rl == t[0]).any()]
    M["o1"] = M.bs.map({k: i for i, (k, _, _) in enumerate(ST)}); M["o2"] = M.rl.map({k: i for i, (k, _, _) in enumerate(RL)})
    M = M.sort_values(["o1", "o2"], kind="stable").reset_index(drop=True); n = len(M)
    hw = round(24 * hax / 42, 1)
    ax = cv.ax(x + 11, yax, hw, hax); ax.set_xlim(0, 28); ax.set_ylim(n, 0)
    for key, tab, x0_ in (("bs", ST, 0), ("rl", RL, 15)):
        cm_ = {k: c for k, _, c in tab}
        for v_, g_ in M.groupby(key, sort=False):
            ax.add_patch(Rectangle((x0_, g_.index.min()), 13, len(g_), fc=cm_[v_], ec="white", lw=.25))
    ax.text(6.5, -.02 * n, "block\nstatus", fontsize=sz["small"] - .3, ha="center", va="bottom"); ax.text(21.5, -.02 * n, f"{noun(run, 1)}\nrole", fontsize=sz["small"] - .3, ha="center", va="bottom")
    for sp_ in ("top", "right", "bottom"):
        ax.spines[sp_].set_visible(False)
    ax.spines["left"].set_position(("data", -1.5)); ax.spines["left"].set_linewidth(.4); ax.set_xticks([])
    step = max(1, int(round(n / 4 / 10 ** np.floor(np.log10(max(n / 4, 1)))) * 10 ** np.floor(np.log10(max(n / 4, 1)))))
    tk = [t for t in range(0, n, step) if n - t > step * .3] + [n]
    ax.set_yticks(tk); ax.set_yticklabels([str(t) for t in tk], fontsize=sz["small"] - .5); ax.tick_params(axis="y", length=1.5, pad=.8)
    ax.set_ylabel(f"number of {noun(run)}", fontsize=sz["small"] - .3, labelpad=1)
    x0k = x + 11 + hw + 2.5; DYK = 2.9; yy = yax + DYK / 2
    for head, tab, col_ in (("BLOCK STATUS", ST, "bs"), (f"{noun(run, 1).upper()} ROLE", RL, "rl")):
        cv.text(x0k, yy, head, fs=sz["small"] - .2, bold=True); yy += DYK
        for k, lab, c in tab:
            n_ = int((M[col_] == k).sum())
            cv.fig.add_artist(Rectangle((x0k / cv.W, 1 - (yy + .9) / cv.H), 1.8 / cv.W, 1.8 / cv.H, transform=cv.fig.transFigure, fc=c, ec="none"))
            cv.text(x0k + 2.8, yy, f"{lab}  {n_}", fs=sz["small"] - .2); yy += DYK
        yy += 1.0
    cv.text(x + 11, yax + hax + 1.5, f"Same {n} {noun(run)} classified two ways:\nblock status (left) and {noun(run, 1)} role (right).", fs=5.5, color=S.NOTE_GREY, va="top")


OH_SPEC = {"R0": (["within-family (full-sib)", "cross-family, same shared-parent group", "second degree or closer, no link", "third degree", "unrelated"],
                  ["full-sib", "shared parent", "2nd degree", "3rd degree", "unrelated"], [S.ROH_BASE["brick"], S.EXTRA["terra"], "#F1E492", S.ROH_BASE["sky"], "#C9D1D8"], "class_R0", "biSNP R$_0$"),
           "DEL_opposite_homozygote_rate": (["within-family (full-sib)", "directly linked families", "same group, not directly linked", "no link"],
                                            ["full-sib", "linked", "same group", "no link"], [S.ROH_BASE["brick"], S.EXTRA["terra"], "#F1E492", "#C9D1D8"], "class_DEL", "DEL OH rate")}


def oh_available(run):
    P = run.t.get("network", {}).get("opposite_pairs")
    if P is None or not len(P):
        return []
    return [c for c in OH_SPEC if c in P and P[c].notna().any() and P.loc[P[OH_SPEC[c][3]] == OH_SPEC[c][0][-1], c].notna().sum() >= 3]


def panel_oh(axc, run, col, sz):
    P = run.t["network"]["opposite_pairs"]; keys, labs, cols, ccol, name = OH_SPEC[col]
    base = P.loc[P[ccol] == keys[-1], col].median(); rg = np.random.default_rng(3)
    data, kk = [], []
    for k, lab, c in zip(keys, labs, cols):
        v = P.loc[P[ccol] == k, col].dropna().to_numpy() / base
        if len(v) >= 3:
            data.append(v); kk.append((lab, c))
    vp = axc.violinplot(data, positions=range(1, len(data) + 1), widths=.85, showextrema=False, showmedians=False)
    for body, (lab, c) in zip(vp["bodies"], kk):
        body.set_facecolor(c); body.set_edgecolor(S.NAVY); body.set_linewidth(.3); body.set_alpha(.75)
    for q, v in enumerate(data, start=1):
        sub = rg.choice(v, size=min(len(v), 160), replace=False); axc.scatter(q + rg.uniform(-.18, .18, len(sub)), sub, s=.9, c=S.NAVY, lw=0, alpha=.16, zorder=2, rasterized=True)
        axc.plot([q - .28, q + .28], [np.median(v)] * 2, color=S.NAVY, lw=.9, zorder=4)
        axc.text(q, np.percentile(v, 99.5) + .06, f"{np.median(v):.2f}", fontsize=sz["small"] - .4, ha="center", va="bottom", color=S.NOTE_GREY, bbox=dict(fc="white", ec="none", alpha=.9, pad=.5), zorder=8)
    axc.set_xticks(range(1, len(data) + 1)); axc.set_xticklabels([l for l, _ in kk], fontsize=sz["small"] - .3, rotation=40, ha="right")
    axc.set_ylabel(f"{name} (× {labs[-1]})", fontsize=sz["small"]); axc.axhline(1, color=S.NOTE_GREY, ls=":", lw=.4); axc.set_xlim(.4, len(data) + .6)
    return len(data)


# ------------------------------------------------------------------ J: network
def network_data(run, max_components=12):
    """Families and kept links of the network at the figure tier, per component (largest first)."""
    net = run.t.get("network", {})
    if "links" not in net or not len(net["links"]):
        return None
    t = int(run.cfg["network"]["figure_min_family_size"])
    L = net["links"]; L = L[(L.min_n.astype(int) >= t) & (L.status != "rejected")]
    if not len(L):
        return None
    sizes = family_sizes(run)
    import networkx as nx
    G = nx.Graph()
    for r in L.itertuples():
        st, ev = link_key(r); G.add_edge(str(r.unit_a), str(r.unit_b), status=st, evidence=ev)
    comps = sorted((sorted(c, key=S.fam_key) for c in nx.connected_components(G)), key=lambda c: (-len(c), -sum(sizes.get(f, 0) for f in c), S.fam_key(c[0])))
    return {"G": G, "comps": comps[:max_components], "n_hidden": max(0, len(comps) - max_components), "sizes": sizes, "tier": t}


def _template(G, comp, seed):
    """Abstract layout of one component: pairs on a diagonal, chains as zigzags in link order, a
    triangle flat, anything else by a force layout (Kamada-Kawai, deterministic)."""
    import networkx as nx
    sub = G.subgraph(comp); n = len(comp)
    if n == 2:
        a, b = comp; return {a: (0, 0), b: (1, 1)}, set()
    deg = dict(sub.degree())
    if sub.number_of_edges() == n - 1 and max(deg.values()) <= 2:
        end = sorted([f for f in comp if deg[f] == 1], key=S.fam_key)[0]; path = [end]
        while len(path) < n:
            path.append(next(x for x in sorted(sub.neighbors(path[-1]), key=S.fam_key) if x not in path))
        tpl = {f: (i, 1 - i % 2) for i, f in enumerate(path)}
        return tpl, {f for f, (x_, y_) in tpl.items() if y_ == 1}
    if n == 3:
        a, b, c = comp; tpl = {a: (0, 1), b: (1, 0), c: (2, 1)}
        return tpl, {a, c}
    p = nx.kamada_kawai_layout(sub)
    return {f: (float(p[f][0]), float(p[f][1])) for f in comp}, set()


def draw_network(cv, run, x0, y_top, aw, ah, sz, kj=None, key=True, nd=None):
    """Audited network in a rectangle of aw x ah mm whose top-left corner is (x0, y_top): a
    two-line key, then one rounded box per component (CC1 full height on the left, the others in
    two rows), nodes sized by family size, links coloured by audit outcome."""
    import matplotlib.patheffects as pe
    from matplotlib.patches import FancyBboxPatch
    nd = nd or network_data(run)
    G, comps, sizes = nd["G"], nd["comps"], nd["sizes"]
    KJ = 16.48 / 13.64 if kj is None else kj; PT = 72 / 25.4
    F_NODE, F_N, F_HEAD, F_KEY = max(5.5, sz["small"] * KJ), 5.5, 5.5, 5.5; KL = F_N / 4.0
    key_h = 6.9 if key else 0.0
    ax = cv.ax(x0, y_top + key_h, aw, ah - key_h); ah_ = ah - key_h; ax.set_xlim(0, aw); ax.set_ylim(0, ah_); ax.axis("off")
    NS = lambda n: (34 + 13 * n) * KJ ** 2
    RAD = {f: np.sqrt(NS(sizes.get(f, 1))) / 2 for c in comps for f in c}; NB = {f: list(G.neighbors(f)) for f in RAD}
    PADB, HDH, GAPB, LH = 1.3, 4.6, 1.6, 2.1
    lab_w = lambda f: (.56 * (len(f"n = {sizes.get(f, 1)}") - 2) + .556) * F_N / PT
    TPL, HIGH = {}, set()
    for k, c in enumerate(comps):
        TPL[k], hi = _template(G, c, run.seed); HIGH |= hi

    def nlab(f, pos_):
        x0_, y0_ = pos_[f]; rad = RAD[f] + 1.0 * KL; best = None
        order = ((90, "center", "bottom"), (270, "center", "top")) if f in HIGH else ((270, "center", "top"), (90, "center", "bottom"))
        for rank, (ang, ha, va) in enumerate(order + ((0, "left", "center"), (180, "right", "center"))):
            dx, dy = np.cos(np.radians(ang)), np.sin(np.radians(ang)); cx, cy = dx * (rad + (4.5 * KL if abs(dx) > .5 else 2.0 * KL)), dy * (rad + 2.0 * KL)
            sep = min([abs((np.degrees(np.arctan2((pos_[b][1] - y0_), (pos_[b][0] - x0_))) - ang + 180) % 360 - 180) for b in NB.get(f, []) if b in pos_] or [180])
            clear = min([np.hypot((pos_[g][0] - x0_) * PT - cx, (pos_[g][1] - y0_) * PT - cy) - RAD[g] for g in pos_ if g != f] or [99])
            ok = sep >= 40 and clear >= 4.0 * KL; score = (ok, min(sep / 40, 1.5) + min(clear / (4.0 * KL), 1.5) - .01 * rank)
            if best is None or score > best[0]:
                best = (score, (dx * rad if abs(dx) > .5 else 0.0, dy * rad if abs(dy) > .5 else 0.0), ha, va)
            if ok and rank == 0:
                break
        return best[1:]

    def extents(pos_):
        xs, ys, labs = [], [], {}
        for f, (x_, y_) in pos_.items():
            r_ = RAD[f] / PT; xs += [x_ - r_, x_ + r_]; ys += [y_ - r_, y_ + r_]; (ox, oy), ha, va = nlab(f, pos_); axx, ayy = x_ + ox / PT, y_ + oy / PT; w_ = lab_w(f); labs[f] = ((ox, oy), ha, va)
            lx0 = axx - w_ / 2 if ha == "center" else (axx if ha == "left" else axx - w_); ly0 = ayy - LH / 2 if va == "center" else (ayy if va == "bottom" else ayy - LH)
            xs += [lx0, lx0 + w_]; ys += [ly0, ly0 + LH]
        return min(xs), max(xs), min(ys), max(ys), labs

    def fit(tpl, bx, by, bw, bh, max_aniso=None):
        ux = [v[0] for v in tpl.values()]; uy = [v[1] for v in tpl.values()]; um = ((max(ux) + min(ux)) / 2, (max(uy) + min(uy)) / 2)
        tx0, tx1, ty0, ty1 = bx + PADB, bx + bw - PADB, by + PADB, by + bh - PADB
        sx = (tx1 - tx0) / max(max(ux) - min(ux), 1e-9) * .6; sy = (ty1 - ty0) / max(max(uy) - min(uy), 1e-9) * .5; cx, cy = bx + bw / 2, by + bh / 2
        best = None
        for it in range(60):
            pos_ = {f: (cx + (u[0] - um[0]) * sx, cy + (u[1] - um[1]) * sy) for f, u in tpl.items()}; x0_, x1_, y0_, y1_, labs = extents(pos_)
            err = max(abs(x0_ - tx0), abs(x1_ - tx1), abs(y0_ - ty0), abs(y1_ - ty1)) + 5 * max(0, tx0 - x0_, x1_ - tx1, ty0 - y0_, y1_ - ty1)
            if it >= 15 and (best is None or err < best[0]):
                best = (err, pos_, labs)
            sx *= ((tx1 - tx0) / max(x1_ - x0_, 1e-6)) ** .5; sy *= ((ty1 - ty0) / max(y1_ - y0_, 1e-6)) ** .5
            sx = max(sx, 1e-3); sy = max(sy, 1e-3)
            if max_aniso:
                sx, sy = min(sx, sy * max_aniso), min(sy, sx * max_aniso)
            cx += ((tx0 + tx1) - (x0_ + x1_)) / 2; cy += ((ty0 + ty1) - (y0_ + y1_)) / 2
        return best[1:]
    # ---- boxes: CC1 on the left over the full height, the others in two rows (one row when few)
    nC = len(comps); box = {}
    if nC == 1:
        box[0] = (0, 0, aw, ah_ - HDH)
    else:
        rest = list(range(1, nC)); ntot = sum(len(c) for c in comps)
        W1 = float(np.clip(aw * (len(comps[0]) / ntot + .04), aw * .22, aw * .60)); xr = W1 + GAPB; WR = aw - xr
        box[0] = (0, 0, W1, ah_ - HDH)
        two = len(rest) >= 3 and ah_ >= 2 * (HDH + 14)
        rows = [rest[:len(rest) // 2], rest[len(rest) // 2:]] if two else [rest]
        hb = (ah_ - (len(rows) - 1) * GAPB - len(rows) * HDH) / len(rows)
        for ri, row in enumerate(rows):
            wts = np.array([len(comps[k]) + 1.5 for k in row], float); ws = (WR - GAPB * (len(row) - 1)) * wts / wts.sum(); xx = xr
            by = (len(rows) - 1 - ri) * (hb + HDH + GAPB)
            for k, w_ in zip(row, ws):
                box[k] = (xx, by, w_, hb); xx += w_ + GAPB
    posn, LAB = {}, {}
    for k in range(nC):
        p_, l_ = fit(TPL[k], *box[k], max_aniso=1.25 if (k == 0 and len(comps[0]) > 3) else None); posn.update(p_); LAB.update(l_)
    STY = {"Gate-1 only": (S.NOTE_GREY, 1.1, "-")}
    used = {}
    for a, b, d in G.edges(data=True):
        st, ev = d["status"], d["evidence"]
        if st == "supported":
            col, lw, ls = evidence_style(ev), (2.0 if " + " in ev else 1.6), "-"; used[("sup", ev)] = (col, lw)
        elif st == "Gate-1 only":
            col, lw, ls = STY[st]; used[("g1", "")] = (col, lw)
        else:
            col, lw, ls = PROV; used[("prov", st)] = (col, lw)
        ax.plot(*zip(posn[a], posn[b]), color=col, lw=lw * KJ, ls=ls, zorder=1)
    fcol = S.family_colours(sorted(sizes, key=S.fam_key))
    for f in posn:
        ax.scatter([posn[f][0]], [posn[f][1]], s=NS(sizes.get(f, 1)), c=fcol.get(f, "#C9D1D8"), ec=S.NAVY, lw=.4 * KJ, zorder=3)
        ax.text(posn[f][0], posn[f][1], f, fontsize=F_NODE if len(f) <= 3 else max(4.0, F_NODE * 3.3 / len(f)), ha="center", va="center", zorder=4, color="white", fontweight="bold",
                path_effects=[pe.withStroke(linewidth=1.2 * KJ, foreground="#1B1B2F")])
        o_, ha_, va_ = LAB[f]; ax.annotate(f"n = {sizes.get(f, 1)}", posn[f], xytext=o_, textcoords="offset points", fontsize=F_N, ha=ha_, va=va_, color="black", zorder=6)
    for k, c in enumerate(comps):
        bx, by, bw_, bh_ = box[k]
        ax.add_patch(FancyBboxPatch((bx, by), bw_, bh_, boxstyle="round,pad=0,rounding_size=2.2", fc="none", ec="#8A949C", lw=.6 * KJ, zorder=0, clip_on=False))
        e_ = [d["status"] for a, b, d in G.subgraph(c).edges(data=True)]
        n_sup = sum(s in ("supported", "Gate-1 only") for s in e_); n_prov = len(e_) - n_sup
        head = f"CC{k + 1} · {len(c)} fam · {sum(sizes.get(f, 0) for f in c)} {noun(run)}"
        if lab_w("x") and cv.text_width(head, F_HEAD, True) > bw_:
            head = f"CC{k + 1} · {len(c)} fam"
        ax.text(bx + .5, by + bh_ + .45, f"{head}\n{n_sup} sup + {n_prov} prov", fontsize=F_HEAD, ha="left", va="bottom", color=S.NOTE_GREY, fontweight="bold", linespacing=1.05)
    if not key:
        return used
    kx = cv.ax(x0, y_top + .2, aw, 6.4); kx.set_xlim(0, aw); kx.set_ylim(6.4, 0); kx.axis("off"); jt = dict(fontsize=F_KEY, va="center", ha="left", color=S.NOTE_GREY)
    SAMP, GAP_S, GAP_I = 3.6, 0.9, 2.6

    def put(x, y, s, **kw):
        kx.text(x, y, s, **{**jt, **kw}); return x + cv.text_width(s, F_KEY, kw.get("fontweight") == "bold")

    def item(x, y, lab, col, lw, ls="-"):
        kx.plot([x, x + SAMP], [y, y], color=col, lw=lw, ls=ls, solid_capstyle="butt"); return put(x + SAMP + GAP_S, y, lab) + GAP_I
    sup_items = [(ev, v) for (kind, ev), v in sorted(used.items(), key=lambda kv: (-kv[0][1].count("+"), kv[0][1])) if kind == "sup"]
    x_ = put(0, 1.5, "Supported link:", fontweight="bold") + 1.6
    for ev, (col, lw) in sup_items:
        x_ = item(x_, 1.5, ev, col, lw * KJ)
    if ("g1", "") in used:
        x_ = item(x_, 1.5, "Gate 1 only", S.NOTE_GREY, 1.1 * KJ)
    if not sup_items and ("g1", "") not in used:
        x_ = put(x_, 1.5, "none") + GAP_I
    put(x_ + .6, 1.5, "node area ~ n")
    x_ = put(0, 4.7, "Provisional link:", fontweight="bold") + 1.6
    prov = sorted({st for (kind, st) in used if kind == "prov"})
    mn = run.cfg["audit"]["min_family_size"]
    txt = ("nominated, not audited" if prov == ["nominated, not audited"] else f"nominated, not validated (a family < {mn} {noun(run)})")
    CCDEF = "CC = component (sup / prov = supported / provisional links)"
    if x_ + SAMP + GAP_S + cv.text_width(txt, F_KEY) + GAP_I + .6 + cv.text_width(CCDEF, F_KEY) > aw:
        txt = "nominated, not audited" if prov == ["nominated, not audited"] else "nominated, not validated"
    xe = item(x_, 4.7, txt, PROV[0], (PROV[1] + .2) * KJ, PROV[2])
    if xe + .6 + cv.text_width(CCDEF, F_KEY) <= aw:
        put(xe + .6, 4.7, CCDEF)
    else:
        put(xe + .6, 4.7, "CC = component")
    return used


# ------------------------------------------------------------------ figures
def fig_shared_parent(run, plt):
    """Port of the published shared-parent figure (v13): A matrix of the auditable families,
    B nomination score, C Gate 1, D / E Gate 2 per SV class, F both Gate-1 scores, G incidence,
    H assignment completeness, I opposite homozygotes, J audited network."""
    if "pairs" not in run.t.get("nomination", {}):
        return None
    E = edges(run); sizes = family_sizes(run); mn = int(run.cfg["audit"]["min_family_size"]); aud = audited(run)
    big = [f for f in sorted(sizes, key=S.fam_key) if sizes[f] >= mn]
    sz = dict(S.SZ_MAIN); W = 240.0; notes = []
    Sb = E[E.unit_a.isin(big) & E.unit_b.isin(big)]
    SNb = Sb[Sb.nominated.astype(bool)]
    fl, n_comp = order_linked_first(big, sizes, list(zip(SNb.unit_a, SNb.unit_b))) if big else ([], 0)
    nf = len(fl)
    # ---- which panels exist
    grid = [("score", "Nomination score distribution")]
    g2 = gate2_classes(run)
    has_llr = aud and len(_aud(run)) > 0
    if has_llr:
        grid.append(("gate1", "Gate 1: shared parent"))
        for q, k in enumerate(g2[:2]):
            nm = {"DEL": "deletions", "DUP": "duplications"}.get(k.split("_")[0], k)
            grid.append((("gate2", k), f"Gate 2{'ab'[q] if len(g2) > 1 else ''}: {nm}"))
        if len(g2) > 2:
            notes.append(f"Gate 2 panels drawn for {g2[0]} and {g2[1]} only ({len(g2)} SV classes in the audit)")
        if not g2:
            notes.append("no Gate 2 panel: no SV class in the audit set")
        grid.append(("compat", "Latent-parent compatibility"))
    else:
        notes.append("panels of Gate 1 and Gate 2 omitted: no audit was run")
    if latent_parents(run) is not None:
        grid.append(("incidence", "Family × latent-parent incidence"))
    else:
        notes.append("incidence panel omitted: no latent parent")
    oh = oh_available(run)
    if "R0" not in oh:
        notes.append("biSNP R0 panel omitted: no R0 in the relatedness table and no SNP input")
    if "DEL_opposite_homozygote_rate" not in oh:
        notes.append("deletion opposite-homozygote panel omitted: no deletion input")
    nd = network_data(run)
    if nd is None:
        notes.append("network panel omitted: no link between families at the figure tier")
    elif nd["n_hidden"]:
        notes.append(f"network panel: the {len(nd['comps'])} largest components are drawn, {nd['n_hidden']} smaller ones are not")
    letters = iter("ABCDEFGHIJKLMN")
    # ---- geometry (mm): v13 positions; the key of A moves below the matrix when it would cover a cell
    AX0, AY0, GA = 6.5, 10.0, 68.0
    C = [85.5, 140.0, 194.5]; GW, GH = 43.0, 24.0; R1, R2 = 10.0, 54.0
    KX0, KY0, KCW, KCH = .3357 * max(nf, 1), -.30 * max(nf, 1) / 23, .6374 * max(nf, 1), .3283 * max(nf, 1)
    pos = {f: i for i, f in enumerate(fl)}
    inside = nf >= 12 and not any(min(pos[a], pos[b]) <= KY0 + KCH + .5 and max(pos[a], pos[b]) >= KX0 - .5 for a, b in zip(SNb.unit_a, SNb.unit_b)) \
        and not any(KX0 - .5 <= i <= KY0 + KCH + .5 for i in range(nf))
    key_below = (not inside) and nf > 0
    rows_grid = int(np.ceil(len(grid) / 3))
    top_block = max(AY0 + GA + 14.0 + (.3283 * 68 + 2.5 if key_below else 0.0), (R1 if rows_grid < 2 else R2) + GH + 14.0)
    RY = max(92.0, top_block) if (nf or rows_grid > 1) else R1 + GH + 16.0
    HAX = 30.0; BOT = 1.5; YAX = RY + 10.5; H = round(YAX + HAX + 10.7 + BOT, 2)
    cv = S.canvas(plt, W, H, sz); fig = cv.fig
    # ---- A
    if nf:
        cv.title(next(letters), "Family-pair nomination and audit", 4.0, 1.0)
        ax = cv.ax(AX0, AY0, GA, GA)
        tick_fs = 5.5 if nf <= 32 else 5.5 * 32 / nf
        FA = GA / 63
        SN, cells, cm = _matrix_core(ax, E, fl, sizes, n_comp, 0, FA ** 2 * float(np.clip(23 / nf, .25, 2.0)), tick_fs, False)
        n_pairs, n_nom = len(Sb), len(SN)
        n_g1 = sum(link_key(r)[0] in ("supported", "Gate-1 only") for r in SN.itertuples())
        FK = GA / 59; KW, KH = .6374 * 59, .3283 * 59
        if inside:
            kax = ax.inset_axes([(KX0 + .5) / nf, 1 - (KY0 + KCH + .5) / nf, KCW / nf, KCH / nf], zorder=8)
        else:
            kax = cv.ax(AX0, AY0 + GA + 10.5, KW * FK, KH * FK)
        _key_A(kax, cm, KW, KH, 5.5 * FK, SN, n_nom, n_g1, aud, run.cfg["nomination"]["q_max"], noun(run), FK ** 2)
        line = f"{n_pairs} family pairs  ·  {n_nom} nominated" + (f"  ·  {n_g1} pass Gate 1" if aud else "")
        fig.text(AX0 / W, 1 - (AY0 + GA + 7.2) / H, line, fontsize=6.5, fontweight="bold", va="center", ha="left")
    else:
        notes.append(f"matrix panel omitted: no family of at least {mn} individuals")
    # ---- B-G grid
    gaxes = []
    for q, (key, ttl) in enumerate(grid):
        x = C[q % 3]; y = R1 if q < 3 else R2
        cv.title(next(letters), ttl, y - 6, x - 9.5)
        ax = cv.ax(x, y, GW, GH); gaxes.append(ax)
        if key == "score":
            panel_score(ax, run, sz)
        elif key == "gate1":
            panel_gate1(ax, run, sz)
        elif key == "compat":
            panel_compat(ax, run, sz)
        elif key == "incidence":
            panel_incidence(ax, run, sz, GH)
        else:
            panel_gate2(ax, run, sz, key[1])
        if key != "incidence":
            ax.xaxis.set_label_coords(.5, -7.2 / GH); ax.yaxis.set_label_coords(-8.0 / GW, .5)
    # ---- bottom row: H, I, J
    cv.title(next(letters), "Family assignment completeness", RY, 1.0)
    panel_completeness(cv, run, 1.0, YAX, HAX, sz)
    xcur = 59.5
    if oh:
        cv.title(next(letters), "Opposite-homozygote gate", RY, xcur)
        cv.text(xcur, RY + 2.5, f"Ratio to unrelated pairs (dotted line = 1). Lower = fewer opposite\nhomozygotes than unrelated {noun(run)}, as expected when a parent is\nshared; full-sib is the reference.",
                fs=5.5, color=S.NOTE_GREY, va="top", linespacing=1.08)
        xa = xcur + 8.0
        for col in oh:
            w_ = round((35 if col == "R0" else 25) * HAX / 36, 1)
            a_ = cv.ax(xa, YAX, w_, HAX); panel_oh(a_, run, col, sz)
            tk_ = [float(v) for v in a_.get_yticks() if a_.get_ylim()[0] - 1e-9 <= v <= 2.5 + 1e-9]; yl_ = a_.get_ylim()
            if tk_ and abs(max(tk_) - 2.5) < 1e-9 and yl_[1] > 2.5:
                a_.set_yticks(tk_); a_.set_ylim(yl_); a_.spines["left"].set_bounds(min(tk_), 2.5)
            xa += w_ + 11.0
        xcur = max(xa - 11.0 + 2.0, xcur + 71.0 if len(oh) == 2 else xcur + 42.0)
    if nd is not None:
        JX0 = max(xcur, 130.5) if len(oh) == 2 else xcur + 2.0
        JX1 = 237.5
        cv.title(next(letters), "Audited shared-parent family network" if aud else "Nominated family network (not audited)", RY, JX0)
        draw_network(cv, run, JX0, RY + 3.1, JX1 - JX0, H - BOT - (RY + 3.1), sz, nd=nd)
    mnw = f"families of at least {mn} {noun(run)}"
    desc = (f"Port of the published shared-parent figure. Matrix: {mnw}; below the diagonal the fraction of cross-family pairs at second-degree kinship or closer with the nominated pairs boxed, "
            "above the diagonal the audit outcome. Then the nomination score of every unit pair, the Gate 1 score, Gate 2 sharing against its matched null per SV class, both Gate 1 scores, "
            "the family x latent-parent incidence with slot-limit violations, assignment completeness, opposite homozygotes by relationship class, and the network by component. "
            "Latent parents are bookkeeping labels; links are candidates.")
    return fig, "Shared-parent links: nomination, audit gates and network", desc, "; ".join(notes)


def fig_matrix_all_units(run, plt):
    """Port of the supplementary matrix of all units: every family in the matrix, links to
    unassigned individuals counted in the right margin."""
    from matplotlib.patches import Rectangle
    if "pairs" not in run.t.get("nomination", {}):
        return None
    E = edges(run); sizes = family_sizes(run); aud = audited(run); mn = int(run.cfg["audit"]["min_family_size"])
    fams = sorted(sizes, key=S.fam_key)
    if len(fams) < 2:
        return None
    M = run.t["families"]["membership"]; n_unass = int((M.family == "unassigned").sum())
    isF = E.unit_a.isin(sizes) & E.unit_b.isin(sizes); EF = E[isF]; Nm = EF[EF.nominated.astype(bool)]; single = E[E.nominated.astype(bool) & ~isF]
    fl, n_comp = order_linked_first(fams, sizes, list(zip(Nm.unit_a, Nm.unit_b)), optimise=False)
    nf = len(fl); pos = {f: i for i, f in enumerate(fl)}
    sz = dict(S.SZ_MAIN); W = 183.0; GA = 150.0; AX0, AY0 = 13.0, 10.0
    KW, KH = 67.0, 37.0
    # the key sits in the upper triangle when no audit cell is under it, otherwise below the matrix
    kx0, ky1 = .47 * nf - .5, (1 - .705) * nf - .5
    inside = not any(min(pos[a], pos[b]) <= ky1 + .5 and max(pos[a], pos[b]) >= kx0 - .5 for a, b in zip(Nm.unit_a, Nm.unit_b)) and not any(kx0 - .5 <= i <= ky1 + .5 for i in range(nf))
    H = AY0 + GA + 19 + (0 if inside else KH + 3)
    cv = S.canvas(plt, W, H, sz); fig = cv.fig
    cv.title("", f"Family-pair nomination and audit: all {nf} families and the unassigned {noun(run)}", 4.0, 1.0)
    ax = cv.ax(AX0, AY0, GA, GA)
    tick_fs = 5.5 if nf <= 60 else 5.5 * 60 / nf
    SN, cells, cm = _matrix_core(ax, E, fl, sizes, n_comp, 6, float(np.clip(51 / nf, .25, 2.5)), tick_fs, True)
    for t in ax.texts:
        t.set_color("white"); t.set_fontsize(min(5.0, tick_fs))
    sc = single.assign(f=np.where(single.unit_a.isin(sizes), single.unit_a, single.unit_b)).groupby("f").size() if len(single) else {}
    if len(single):
        ax.text(nf - .3, -1.3 * nf / 51, f"links to\nunassigned {noun(run)}", fontsize=5.5, ha="left", va="center", color="black", clip_on=False, zorder=20, linespacing=1.0)
        for f, c_ in sc.items():
            if f in pos:
                ax.text(nf + .3 * nf / 51, pos[f], str(c_), fontsize=min(5.5, tick_fs), ha="left", va="center", color="black", clip_on=False)
    lax = ax.inset_axes([.47, .705, KW / GA, KH / GA]) if inside else cv.ax(AX0, AY0 + GA + 17, KW, KH)
    lax.set_xlim(0, KW); lax.set_ylim(KH, 0); lax.set_xticks([]); lax.set_yticks([])
    for sp in lax.spines.values():
        sp.set_visible(True); sp.set_linewidth(.3); sp.set_color("#BBBBBB")
    kf = 5.8; kt = dict(fontsize=kf, va="center", ha="left"); n_ff = len(Nm); q_max = run.cfg["nomination"]["q_max"]
    lax.text(2, 2.6, "Lower triangle: nomination", fontweight="bold", **kt)
    lax.imshow(np.linspace(0, 1, 64)[None, :], cmap=cm, aspect="auto", extent=(2.2, 12.2, 6.6, 4.9), zorder=2); lax.add_patch(Rectangle((2.2, 4.9), 10, 1.7, fill=False, ec=S.NAVY, lw=.3, zorder=3)); lax.set_xlim(0, KW); lax.set_ylim(KH, 0)
    lax.text(2.2, 8.0, "0", fontsize=kf - .6, va="center", ha="center"); lax.text(12.2, 8.0, "1", fontsize=kf - .6, va="center", ha="center")
    lax.text(14, 5.75, "fraction of cross-family pairs at 2nd-degree θ", **kt)
    lax.add_patch(Rectangle((2.2, 9.6), 2.2, 2.2, fc=cm(.85), ec="#1B1B2F", lw=1.0)); lax.text(6.2, 10.7, f"nominated, BH $\\mathit{{q}}$ < {q_max:g} ({n_ff} family pairs)", **kt)
    cnt = _counts(SN); n_g1 = sum(v for (st, _e), v in cnt.items() if st in ("supported", "Gate-1 only"))
    lax.text(2, 14.6, "Upper triangle: audit outcome" if aud else "Upper triangle: no audit was run", fontweight="bold", **kt)
    items = []
    if aud:
        lax.add_patch(Rectangle((2.2, 16.6), 2.2, 2.2, fc=S.G1_FILL, ec=S.GREEN, lw=.7)); lax.text(6.2, 17.7, f"Gate 1 pass ({n_g1}); symbol = SV support", **kt)
        for ev in sorted({e for (st, e) in cnt if st == "supported"}, key=lambda e: (-e.count("+"), e)):
            items.append((f"{ev} support" + ("" if " + " in ev else " only"), ("supported", ev)))
        items += [("no SV support", ("Gate-1 only", "")), (f"Gate 1 fail{_failed_names(SN)}", ("rejected", "")), (f"not testable (a family < {mn} {noun(run)})", ("provisional", ""))]
    else:
        items.append(("nominated, not audited", ("nominated, not audited", "")))
    for q, (lab, key) in enumerate(items[:6]):
        col, mk, sz_, ec = symbol(*key); x0 = 3.3 if q % 2 == 0 else 29.5; yy = 21.2 + (q // 2) * 3.1
        lax.scatter([x0], [yy], s=sz_ * .9, c=col, marker=mk, ec=ec if ec else "none", lw=.5); lax.text(x0 + 2.4, yy, f"{lab} ({cnt.get(key, 0)})", **kt)
    lax.add_patch(Rectangle((2.2, 31.2), 2.2, 2.2, fc=(.52, .52, .52), ec="#7a7a7a", lw=.2)); lax.text(6.2, 32.3, f"diagonal: family, shaded by size (number = {noun(run)} if ≥ 6)", **kt)
    n_units = len(set(E.unit_a) | set(E.unit_b)); n_all = len(E); n_all_nom = int(E.nominated.sum()); n_fs = len(single); n_big = sum(v >= mn for v in sizes.values())
    fig.text(1.0 / W, 1 - (AY0 + GA + 9.5) / H,
             f"{nf} families and {n_unass} unassigned {noun(run)} = {n_units} units; {n_all:,} comparisons, {n_all_nom} nominated (BH $\\mathit{{q}}$ < {q_max:g}): {n_ff} between two families (boxed cells) and {n_fs} between a family\n"
             f"and an unassigned {noun(run, 1)} (counts in the right margin). The dotted lines separate the {n_comp} families with at least one nominated family link from the other {nf - n_comp}. The main-figure panel\n"
             f"shows the {n_big} families with ≥ {mn} {noun(run)} ({n_big * (n_big - 1) // 2} pairs), the only pairs that can be audited.", fontsize=5.8, color=S.NOTE_GREY, va="top", linespacing=1.15)
    desc = ("Port of the published supplementary matrix. Every family; lower triangle = fraction of cross-family pairs at second-degree kinship or closer, nominated pairs boxed; "
            "upper triangle = audit outcome; right margin = number of nominated links between the family and unassigned individuals.")
    return fig, "Family-pair nomination matrix, all units", desc, "" if inside else "key drawn below the matrix (it would have covered a cell inside)"
