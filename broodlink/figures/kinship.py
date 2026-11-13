"""Kinship per marker class: ports of the published Figure "SV genotypes follow the full-sib
families and kinship recovers them" (`build_Figure_4.py`, panels A-H) and of the Extended Data
figure on cross-marker concordance (`build_ED_Fig_H4.py`, panels A-D), and the overview figure
of the README (kinship matrices + network).

All numbers come from step 8 (family effect: `family_effect_per_marker.tsv`,
`family_effect_summary.tsv`) and step 9 (`cross_marker_concordance.tsv`, `kinship_grm.npz`,
`kinship_grm_individuals.tsv`).
"""
from __future__ import annotations
import numpy as np

from ..io import allele_freq
from . import style as S
from .shared_parent import draw_network, network_data, noun

SRC_T = "09_cross_marker_concordance/cross_marker_concordance.tsv"
SRC_G = "09_cross_marker_concordance/kinship_grm.npz + kinship_grm_individuals.tsv"
SRC_E = "08_family_effect/family_effect_summary.tsv"
SRC_V = "08_family_effect/family_effect_per_marker.tsv"
L = S.LED


def hap_label(run, cap=False, short=False):
    w = str(run.cfg["figures"].get("haplotype_label", "structural haplotypes"))
    if short and w == "structural haplotypes":
        w = "struct. haplotypes"
    return w[0].upper() + w[1:] if cap else w


def klabel(run, klass, cap=False):
    """Label of a marker class as printed: biSNPs, DELs, DUPs, the configured name of the structural haplotypes."""
    k = str(klass)
    if k == "HAPLOTYPE":
        return hap_label(run, cap)
    if k.startswith("SNP_matched_"):
        return f"biSNPs ({klabel_short(run, k[12:])}-matched)"
    return S.class_label(k)


def klabel_short(run, klass):
    k = str(klass)
    return {"HAPLOTYPE": "hap."}.get(k, S.class_label(k, plural=False)) if run.cfg["figures"].get("haplotype_label", "structural haplotypes") == "structural haplotypes" or k != "HAPLOTYPE" else "INV"


def kin_classes(run, masked=True):
    """Marker classes of step 9 that have a GRM, in table order: dicts with key, klass, label, colour, K, row."""
    c = run.t.get("concordance", {})
    if "summary" not in c or "grm_order" not in c:
        return [], None
    T = c["summary"]; grm = run.get_grm(); O = c["grm_order"]
    out = []
    for r in T.itertuples():
        if str(r.grm_key) in grm and (masked or r.marker_class != "SNP_MASKED"):
            out.append({"key": str(r.grm_key), "klass": str(r.marker_class), "label": klabel(run, r.marker_class), "cap": klabel(run, r.marker_class, True),
                        "colour": S.class_colour(r.marker_class), "K": grm[str(r.grm_key)], "row": r})
    return out, O


def _heat(cv, ax, K, n_asg, cmap):
    A = np.array(K, float); np.fill_diagonal(A, np.nan)
    ax.imshow(np.clip(A, 0, .6), cmap=cmap, vmin=0, vmax=.6, interpolation="nearest", aspect="auto", rasterized=True)
    if 0 < n_asg < len(A):
        ax.axvline(n_asg - .5, color=S.NAVY, lw=.3, ls=":"); ax.axhline(n_asg - .5, color=S.NAVY, lw=.3, ls=":")
    ax.set_xticks([]); ax.set_yticks([]); S.frame(ax, .35)


def _colourbar(cv, x, y, h, sz):
    cax = cv.ax(x, y, 1.4, h); cax.imshow(np.linspace(1, 0, 64)[:, None], cmap=S.ros_seq(), aspect="auto"); cax.set_xticks([])
    cax.set_yticks([0, 63]); cax.set_yticklabels(["≥ 0.6", "0"], fontsize=sz["small"]); cax.yaxis.tick_right(); cax.tick_params(length=1.2, pad=1)
    for s in cax.spines.values():
        s.set_visible(False)


def kinship_row(cv, run, KC, O, top, letters, sz, x0=1.0, pitch=35.3, side=30.0, fig_name=""):
    """One row of kinship heat maps (panels D-H of the published figure): title, matrix, AUC and marker count."""
    asg = int((O.family != "unassigned").sum()); cm = S.ros_seq()
    for q, kc in enumerate(KC):
        x = x0 + q * pitch; let = next(letters) if letters is not None else ""
        cv.title(let, kc["cap"], top, x, dx=4)
        ax = cv.ax(x + 1, top + 3.5, side, side); _heat(cv, ax, kc["K"], asg, cm)
        r = kc["row"]
        if np.isfinite(r.auc_within_vs_between):
            au = f"{L(let or kc['key'], 'AUC ' + kc['label'], float(r.auc_within_vs_between), SRC_T + ' auc_within_vs_between'):.4f}"
            lab = f"AUC {au}\n{L(let or kc['key'], 'markers ' + kc['label'], int(r.n_markers), SRC_T + ' n_markers'):,} variants"
        else:
            lab = f"{int(r.n_markers):,} variants"
        ax.set_xlabel(lab, fontsize=sz["small"], labelpad=1.5, linespacing=1.15)
    _colourbar(cv, x0 + (len(KC) - 1) * pitch + 1 + side + 1.2, top + 3.5, side, sz)
    return asg


def _eta_tables(run):
    fe = run.t.get("family_effect", {})
    if "per_marker" not in fe or "summary" not in fe or not len(fe["per_marker"]):
        return None, None
    return fe["per_marker"], fe["summary"].drop_duplicates("marker_class").set_index("marker_class")


def _genotype_blocks(run, order_asg):
    """Genotype matrices for panel A (markers x assigned individuals in family order): all structural
    haplotypes; 40 random deletions and 30 random duplications with allele frequency 0.2-0.8
    called in at least 90% of the individuals (170 of 188 in the original), seed 3."""
    rng = np.random.default_rng(3); blocks = []
    hap, _i = run.hap_set()

    def reorder(ms):
        ix = {s: i for i, s in enumerate(ms.samples)}
        cols = [ix.get(s, -1) for s in order_asg]
        G = np.full((ms.n, len(cols)), -1, np.int8); ok = np.array(cols) >= 0
        G[:, ok] = ms.G[:, np.array(cols)[ok]]
        return G
    if hap is not None and hap.n:
        G = reorder(hap)
        cap = 60
        blocks.append((f"{hap_label(run, True)} (all {hap.n})" if hap.n <= cap else f"{hap_label(run, True)} (first {cap} of {hap.n})", G[:cap], S.HAP_C))
    for klass, k in (("DEL", 40), ("DUP", 30)):
        cand = [m for m in run.sets if m.klass == klass]
        if not cand:
            continue
        G = reorder(max(cand, key=lambda m: m.n))
        af = allele_freq(G); ok = (af > .2) & (af < .8) & ((G >= 0).sum(1) >= .9 * G.shape[1])
        G = G[ok]
        if not len(G):
            continue
        sel = rng.choice(len(G), min(k, len(G)), replace=False)
        blocks.append((f"{S.class_label(klass)} ({len(sel)} random, allele frequency 0.2-0.8)", G[sel], S.class_colour(klass)))
    return blocks


def fig_sv_kinship(run, plt):
    """Port of the published figure: A genotypes in family order, B eta squared against allele
    frequency, C eta squared above the null, D-H kinship matrix per marker class."""
    from matplotlib.colors import ListedColormap
    KC, O = kin_classes(run)
    if not KC:
        return None
    notes = []
    if len(KC) > 5:
        notes.append(f"kinship matrices drawn for the first 5 of {len(KC)} marker classes"); KC = KC[:5]
    if not any(k["klass"] == "SNP" for k in KC):
        notes.append("no biSNP kinship matrix: no SNP input")
    if not any(k["klass"] == "SNP_MASKED" for k in KC):
        notes.append("no SV-masked biSNP matrix: it needs SNPs together with SV calls or intervals")
    V2, S2 = _eta_tables(run)
    sz = dict(S.SZ_ED); W = 183.0
    fam = O.family.to_numpy(); order = O.individual.astype(str).tolist(); asg = fam != "unassigned"
    order_asg = [s for s, a in zip(order, asg) if a]; fam_asg = fam[asg]
    fl = list(dict.fromkeys(fam_asg)); size = {f: int((fam_asg == f).sum()) for f in fl}
    top = V2 is not None and len(order_asg) > 0
    blocks = _genotype_blocks(run, order_asg) if top else []
    yA = 17.5
    for _lab, G, _c in blocks:
        yA += .64 * len(G) + 1 + 4.6
    T2 = max(117.0, yA + 6.0) if top else 0.0
    H = T2 + 55.0
    cv = S.canvas(plt, W, H, sz); fig = cv.fig; letters = iter("ABCDEFGHIJKL")
    if top:
        sv = [k for k in ("HAPLOTYPE", "DEL", "DUP") if k in S2.index] + [k for k in S2.index if k not in ("HAPLOTYPE", "DEL", "DUP", "SNP", "SNP_INDEL", "INDEL") and not str(k).startswith("SNP_matched_")]
        f2 = lambda c: f"{L('banner 1 / B key', 'median eta2 ' + c, float(S2.loc[c, 'eta2_obs_median']), SRC_E + ' eta2_obs_median'):.2f}"
        nf_ = lambda c: f"{L('banner 1 / B key', 'null median ' + c, float(S2.loc[c, 'null_median_of_medians']), SRC_E + ' null_median_of_medians'):.2f}"
        if sv:
            nulls = sorted(float(S2.loc[c, "null_median_of_medians"]) for c in sv)
            cv.banner(0, "1", "SVs co-segregate with families: median η² " + ", ".join(f"{f2(c)} {klabel(run, c)}" for c in sv[:3]) + " between families; null "
                      + (f"{nulls[0]:.2f}-{nulls[-1]:.2f}" if nulls[-1] - nulls[0] >= .005 else f"{nulls[0]:.2f}"))
        else:
            cv.banner(0, "1", "Genotypes follow the families")
        # ---- A
        cv.title(next(letters), f"Genotypes of {L('A', 'individuals in families', len(order_asg), '02_families/family_membership.tsv')} {noun(run)}, in family order", 10, 0.5, dx=4)
        y0 = 17.5; bnd = np.cumsum([size[f] for f in fl]); ax = None; nA = len(order_asg)
        for lab, G, cc in blocks:
            h = .64 * len(G) + 1; ax = cv.ax(1, y0, 118, h); pal = ["#F1F3F5", S.mixw(cc, .55), S.mixw(cc, -.2)]
            ax.imshow(np.ma.masked_less(G, 0), cmap=ListedColormap(pal), vmin=0, vmax=2, aspect="auto", interpolation="nearest", rasterized=True)
            for b_ in bnd[:-1]:
                ax.axvline(b_ - .5, color="white", lw=.3)
            ax.set_xticks([]); ax.set_yticks([]); S.frame(ax, .35)
            cv.text(1, y0 - 1.6, lab, fs=sz["tick"], bold=True)
            for q, (cp, lb) in enumerate(zip(pal, ("hom-ref", "het", "hom-alt"))):
                cv.swatch(78 + q * 14, y0 - 1.6, cp, lb, ec=S.BACKBONE)
            y0 += h + 4.6
        if ax is not None:
            ax.set_xticks(bnd - .5); ax.set_xticklabels([]); ax.tick_params(length=1.2, width=.3)
            st_, last = 0, -9
            for f_ in fl:
                k_ = size[f_]; xm = (st_ + k_ / 2 - .5) * 118 / nA
                if xm - last >= 2.9:
                    d = "".join(ch for ch in str(f_) if ch.isdigit())
                    ax.text(st_ + k_ / 2 - .5, -.05, d if d and len(str(f_)) - len(d) <= 1 else str(f_), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=sz["small"]); last = xm
                st_ += k_
            cv.text(1, y0 + .3, f"full-sib family, largest first (ticks: all {len(fl)} families; number printed where it fits); missing genotype = white", fs=sz["small"], color=S.INK)
        else:
            notes.append("genotype heat maps omitted: no deletion, duplication or structural-haplotype genotypes")
        # ---- B
        cv.title(next(letters), "η² against allele frequency", 10, 124, dx=4)
        edges = np.arange(.05, .96, .1); rb = np.random.default_rng(5)
        ref = next((c for c in ("DEL", "DUP") if c in sv), sv[0] if sv else None)
        lines = [(c, S.class_colour(c), "-") for c in ("DEL", "DUP", "HAPLOTYPE") if c in sv]
        if ref is not None and f"SNP_matched_{ref}" in S2.index:
            lines.append((f"SNP_matched_{ref}", S.BACKBONE, "--"))
        elif not any(str(k).startswith("SNP_matched_") for k in S2.index):
            notes.append("matched biSNPs not drawn: no SNP input")
        rl = klabel(run, ref) if ref else ""
        items = [(col, f"{klabel(run, key) if not key.startswith('SNP_matched') else 'matched biSNPs'} {f2(key)}", ls) for key, col, ls in lines]
        if ref is not None:
            items.append((S.INK, f"permutation null, {rl} {nf_(ref)}", ":"))
        xk, yk = 126.0, 15.0
        for c, l, ls in items:
            w_ = 4.1 + cv.text_width(l, sz["small"]) + 2.5
            if xk + w_ > 183 and xk > 126.0:
                xk, yk = 126.0, yk + 4.0
            cv.line(xk, xk + 3.4, yk, c, .9, ls); cv.text(xk + 4.1, yk, l, fs=sz["small"]); xk += w_
        cv.text(126, yk + 4.0, "values = overall medians; bands = 95% bootstrap per bin", fs=sz["small"], color=S.INK)
        btop = yk + 7.0
        axb = cv.ax(132, btop, 50, 54.0 - btop)
        nbin = {}
        for key, col, ls in lines:
            d = V2[V2.marker_class == key]; b_ = np.digitize(d.af, edges); xm, ym, lo_, hi_ = [], [], [], []
            for j in range(len(edges) + 1):
                v_ = d.eta2_obs[b_ == j].to_numpy()
                if len(v_) >= (5 if key != "HAPLOTYPE" else 3):
                    bs = [np.median(rb.choice(v_, len(v_))) for _ in range(300)]
                    xm.append(float(d.af[b_ == j].mean())); ym.append(float(np.median(v_))); lo_.append(float(np.percentile(bs, 2.5))); hi_.append(float(np.percentile(bs, 97.5)))
                    if key == ref:
                        nbin[j] = (xm[-1], len(v_))
            axb.fill_between(xm, lo_, hi_, color=col, alpha=.22, lw=0); axb.plot(xm, ym, color=col, lw=.8, ls=ls, marker="o", ms=1.5)
        if ref is not None:
            nul = V2[V2.marker_class == ref]; b_ = np.digitize(nul.af, edges); js = [j for j in range(len(edges) + 1) if (b_ == j).sum() >= 5]
            axb.plot([nul.af[b_ == j].mean() for j in js], [nul.null_median[b_ == j].median() for j in js], color=S.INK, lw=.7, ls=":")
            for j, (xx, nn) in nbin.items():
                t = f"{nn / 1000:.1f}k" if nn >= 1000 else str(nn)
                axb.text(xx, .02, t, fontsize=sz["small"], ha="center", color=S.INK); L("B", f"{klabel(run, ref)} in allele-frequency bin {j}", nn, SRC_V + " af")
        axb.set_xlabel(f"allele frequency (numbers = {rl} per bin)", fontsize=sz["tick"], labelpad=1.5); axb.set_ylabel("median η²", fontsize=sz["tick"]); axb.set_ylim(0, .8); axb.set_xlim(.05, .95); axb.tick_params(labelsize=sz["tick"])
        # ---- C
        cv.title(next(letters), "η² above the null of each variant", 66, 124, dx=4)
        KEYS = []
        for c in ("DEL", "DUP", "HAPLOTYPE"):
            if c in sv:
                KEYS.append((c, klabel(run, c), S.class_colour(c)))
                if f"SNP_matched_{c}" in S2.index:
                    KEYS.append((f"SNP_matched_{c}", klabel(run, f"SNP_matched_{c}"), S.BACKBONE))
        if not KEYS:
            KEYS = [(c, klabel(run, c), S.class_colour(c)) for c in list(S2.index)[:6]]
        axc = cv.ax(148, 70, 34, 36)
        for q_, (key, nm, col) in enumerate(KEYS):
            d = V2[V2.marker_class == key]; ex = (d.eta2_obs - d.null_median).to_numpy(); y = len(KEYS) - 1 - q_
            if len(ex) > 2 and np.ptp(ex) > 0:
                vp = axc.violinplot([ex], positions=[y], vert=False, widths=.82, showextrema=False); vp["bodies"][0].set_facecolor(col); vp["bodies"][0].set_alpha(.4); vp["bodies"][0].set_edgecolor("none")
            else:
                axc.scatter(ex, np.full(len(ex), y), s=4, color=col, lw=0)
            if len(ex):
                axc.plot([np.median(ex)] * 2, [y - .34, y + .34], color=S.NAVY, lw=.7)
                axc.text(.84, y, f"{L('C', 'median excess ' + key, float(np.median(ex)), SRC_V + ' eta2_obs - null_median'):.2f}", fontsize=sz["small"], va="center", ha="right", color=S.NAVY)
        axc.axvline(0, color=S.RED, lw=.6, ls="--"); axc.set_yticks(range(len(KEYS))); axc.set_yticklabels([k[1] for k in KEYS][::-1], fontsize=sz["small"]); axc.set_xlim(-.15, .86)
        axc.set_xlabel("η² − null median (dashed: none)", fontsize=sz["tick"], labelpad=1.5); axc.tick_params(axis="x", labelsize=sz["tick"]); axc.set_ylim(-.6, len(KEYS) - .4)
    else:
        notes.append("genotype and eta-squared panels omitted: the family effect did not run")
    # ---- kinship row
    a4 = lambda kc: f"{float(kc['row'].auc_within_vs_between):.4f}"
    svk = [kc for kc in KC if kc["klass"] not in ("SNP", "SNP_MASKED", "SNP_INDEL", "INDEL") and np.isfinite(kc["row"].auc_within_vs_between)] or [kc for kc in KC if np.isfinite(kc["row"].auc_within_vs_between)]
    src_fam = "biSNP families" if str(getattr(run, "rel_mode", "")) == "given" or any(k["klass"] == "SNP" and bool(k["row"].is_source_of_families) for k in KC) else "families"
    cv.banner(T2, "2" if top else "1", f"Kinship from each {'SV' if svk and svk[0]['klass'] not in ('SNP', 'SNP_MASKED') else 'marker'} class agrees with the {src_fam}: AUC " + ", ".join(f"{a4(kc)} {kc['label']}" for kc in svk[:4]))
    n_asg = kinship_row(cv, run, KC, O, T2 + 10, letters, sz)
    n = len(O); nfam = len(fl)
    cv.text(1, T2 + 51.5, f"VanRaden kinship, {L('note', 'individuals', n, SRC_G)} {noun(run)} in family order; dotted lines separate the {L('note', 'assigned', n_asg, SRC_G)} {noun(run)} in {nfam} families from the {L('note', 'unassigned', n - n_asg, SRC_G)} unassigned",
            fs=sz["small"], color=S.INK)
    desc = ("Port of the published figure. A genotypes of the individuals in families, ordered by family (largest first): every structural haplotype, and random deletions and duplications at allele frequency 0.2-0.8 (seed 3). "
            "B median eta squared (variance explained by family) per allele-frequency bin with a bootstrap band, against frequency-matched biSNPs and the permutation null. C eta squared minus the null median of each variant. "
            "Kinship row: VanRaden GRM per marker class, clipped at 0 and 0.6, with the AUC of within- against between-family pairs. Families come from the relatedness source: concordance, not pedigree truth.")
    return fig, "SV genotypes follow the families and kinship recovers them", desc, "; ".join(notes)


def fig_concordance(run, plt):
    """Port of the Extended Data figure: A AUC with its family-level interval, B kinship of within-
    and between-family pairs per class, C SV kinship against masked-biSNP kinship, D pairs recovered."""
    KC, O = kin_classes(run)
    KC = [k for k in KC if np.isfinite(k["row"].auc_within_vs_between)]
    if not KC:
        return None
    notes = []
    if len(KC) > 5:
        notes.append(f"the first 5 of {len(KC)} marker classes are drawn"); KC = KC[:5]
    sz = dict(S.SZ_ED); W, H = 183.0, 124.0
    fam = O.family.to_numpy(); asg = fam != "unassigned"; n = len(O)
    iu = np.triu_indices(n, 1); both = asg[iu[0]] & asg[iu[1]]; within = (fam[iu[0]] == fam[iu[1]]) & both; between = both & ~within
    V = {k["key"]: k["K"][iu] for k in KC}
    nw, nb = int(within.sum()), int(between.sum())
    masked = next((k for k in KC if k["klass"] == "SNP_MASKED"), None)
    svc = [k for k in KC if k["klass"] not in ("SNP", "SNP_MASKED", "SNP_INDEL", "INDEL")][:3]
    has_c = masked is not None and bool(svc)
    cv = S.canvas(plt, W, H if True else H, sz, thin=False)
    COL = {k["key"]: k["colour"] for k in KC}; SHORT = {k["key"]: k["cap"] if k["klass"] == "HAPLOTYPE" else k["label"] for k in KC}
    a4 = lambda k: f"{float(k['row'].auc_within_vs_between):.4f}"
    hi_k = max(KC, key=lambda k: k["row"].auc_within_vs_between); lo_k = min(KC, key=lambda k: k["row"].auc_within_vs_between)
    L("banner 1", "highest AUC " + hi_k["label"], float(hi_k["row"].auc_within_vs_between), SRC_T + " auc_within_vs_between"); L("banner 1", "lowest AUC " + lo_k["label"], float(lo_k["row"].auc_within_vs_between), SRC_T + " auc_within_vs_between")
    cv.banner(0, "1", f"Every marker class separates sibling pairs from other pairs, from AUC {a4(hi_k)} ({hi_k['label']}) to {a4(lo_k)} ({lo_k['label']})" if len(KC) > 1
              else f"{hi_k['cap']} separate sibling pairs from other pairs: AUC {a4(hi_k)}")
    # ---- A
    cv.title("A", "AUC, family-level 95% CI", 10, .5, dx=4)
    ax = cv.ax(25, 16, 34, 37); nk = len(KC)
    lo_all = min(float(k["row"].auc_ci95_low) if np.isfinite(k["row"].auc_ci95_low) else float(k["row"].auc_within_vs_between) for k in KC)
    x_lo = min(.97, np.floor((lo_all - .002) * 100) / 100)
    for q, k in enumerate(KC):
        y = nk - 1 - q; r = k["row"]; a, lo, hi = float(r.auc_within_vs_between), float(r.auc_ci95_low), float(r.auc_ci95_high)
        if np.isfinite(lo):
            ax.plot([lo, hi], [y, y], color=S.NAVY, lw=.7, zorder=2)
        ax.scatter([a], [y], s=9, color=COL[k["key"]], ec=S.NAVY, lw=.3, zorder=3)
        t = f"{a:.4f} [{lo:.4f}, {hi:.4f}]" if np.isfinite(lo) else f"{a:.4f}"
        ax.text(1.0, y + .36, t, fontsize=sz["small"], ha="right", va="center"); L("A", "AUC [CI] " + k["label"], t, SRC_T + " auc_within_vs_between, auc_ci95_low, auc_ci95_high")
    tk = np.linspace(x_lo, 1.0, 4)
    ax.set_yticks(range(nk)); ax.set_yticklabels([SHORT[k["key"]] for k in KC[::-1]], fontsize=sz["small"]); ax.set_ylim(-.5, nk - .25); ax.set_xlim(x_lo, 1.0 + (1 - x_lo) * .0167)
    ax.set_xticks(tk); ax.set_xticklabels([f"{t:.2f}" for t in tk], fontsize=sz["tick"]); ax.set_xlabel("AUC, within vs between families", fontsize=sz["tick"], labelpad=1.5)
    # ---- B
    cv.title("B", f"Kinship of {L('B', 'within-family pairs', nw, SRC_G):,} within-family and {L('B', 'between-family pairs', nb, SRC_G):,} between-family pairs", 10, 63, dx=4)
    x0b, x1b = -.35, 1.0                      # range of the original; values outside are drawn in the edge bins
    bw = np.linspace(x0b, x1b, 80)
    for q, k in enumerate(KC):
        x0 = 66 + (q % 3) * 39.5; y0 = 19.5 + (q // 3) * 19.5; v = V[k["key"]]
        ax = cv.ax(x0, y0, 35, 12.5); vb, vw = v[between & np.isfinite(v)], v[within & np.isfinite(v)]
        ax.hist(np.clip(vb, x0b, x1b), bins=bw, density=True, color="#C9D1D8", lw=0); ax.hist(np.clip(vw, x0b, x1b), bins=bw, density=True, color=S.ROH_BASE["brick"], alpha=.85, lw=0)
        ax.set_xlim(x0b, x1b); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.tick_params(labelsize=sz["tick"], length=1.5); ax.set_xticks([0, .5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
        mw, mb = float(np.median(vw)), float(np.median(vb))
        ax.axvline(mb, color=S.NAVY, lw=.4, ls=":"); ax.axvline(mw, color=S.NAVY, lw=.4)
        t = f"{mb:.2f} | {mw:.2f}"; ax.text(.99, 1.0, t, transform=ax.transAxes, fontsize=sz["small"], ha="right", va="bottom"); L("B", "median between | within " + k["label"], t, SRC_G)
        ax.text(0, 1.0, SHORT[k["key"]], transform=ax.transAxes, fontsize=sz["small"], fontweight="bold", ha="left", va="bottom")
        if q >= len(KC) - 3 or q >= 2:
            ax.set_xlabel("kinship (GRM)", fontsize=sz["tick"], labelpad=1)
    qk = len(KC); xk = 66 + (qk % 3) * 39.5 + 1.0 if qk < 6 else 146.0; yk = 19.5 + (qk // 3) * 19.5 + 2.0
    if qk % 3 == 0 and qk < 6:
        xk, yk = 146.0, 41.0
    if qk in (4, 5):
        xk, yk = 146.0, 41.0
    cv.swatch(xk, yk, "#C9D1D8", "between-family"); cv.swatch(xk, yk + 3.3, S.ROH_BASE["brick"], "within-family")
    cv.text(xk, yk + 7.5, "top right: median\nbetween | within\n(dotted | solid line)", fs=sz["small"], color=S.INK, linespacing=1.1)
    # ---- banner 2, C, D
    T2 = 62.0
    rec = lambda k: f"{float(k['row'].recovery_pct):.0f}"
    lo_r = min(KC, key=lambda k: k["row"].recovery_pct); hi_r = max(KC, key=lambda k: k["row"].recovery_pct)
    if has_c:
        XS = V[masked["key"]]; rr = {}
        for k in svc:
            ok = np.isfinite(XS) & np.isfinite(V[k["key"]]); rr[k["key"]] = float(np.corrcoef(XS[ok], V[k["key"]][ok])[0, 1])
        best = max(svc, key=lambda k: rr[k["key"]])
        L("banner 2", "Pearson r " + best["label"], rr[best["key"]], SRC_G); L("banner 2", "recovery " + lo_r["label"], float(lo_r["row"].recovery_pct), SRC_T + " recovery_pct")
        b2 = f"{klabel_short(run, best['klass'])} kinship tracks masked-biSNP kinship (r = {rr[best['key']]:.2f})"
        b2 += f"; {lo_r['label']} discriminate but recover only {rec(lo_r)}% of pairs at a strict cut" if float(lo_r["row"].recovery_pct) < 50 else f"; at a strict cut {lo_r['label']} recover {rec(lo_r)}% of pairs"
        cv.banner(T2, "2", b2)
        npair = int(np.isfinite(XS).sum())
        cv.title("C", f"SV kinship against masked-biSNP kinship, all {L('C', 'pairs', npair, SRC_G):,} pairs", T2 + 10, .5, dx=4)
        for q, k in enumerate(svc):
            ax = cv.ax(9 + q * 40.5, T2 + 16, 35, 35); b = ~within; y = V[k["key"]]
            ax.scatter(XS[b], y[b], s=.5, c=COL[k["key"]], lw=0, alpha=.35, rasterized=True); ax.scatter(XS[within], y[within], s=1.2, c=S.ROH_BASE["brick"], lw=0, alpha=.7, rasterized=True)
            ax.set_xlabel("masked-biSNP kinship", fontsize=sz["tick"], labelpad=1.5); ax.tick_params(labelsize=sz["tick"], length=1.5)
            if q == 0:
                ax.set_ylabel("SV kinship", fontsize=sz["tick"], labelpad=1.5)
            t = f"r = {rr[k['key']]:.2f}"; ax.text(.03, .97, f"{SHORT[k['key']]}\n{t}", transform=ax.transAxes, fontsize=sz["small"], va="top", fontweight="bold", linespacing=1.1); L("C", "Pearson r " + k["label"], t, SRC_G)
        cv.swatch(86, T2 + 13.0, S.ROH_BASE["brick"], "within-family pairs; other colour = all other pairs", fs=sz["small"])
        dx0, dlet = 131.0, "D"
    else:
        notes.append("panel of SV kinship against masked-biSNP kinship omitted: it needs SNPs and SV calls together")
        L("banner 2", "recovery range", f"{rec(lo_r)}-{rec(hi_r)}", SRC_T + " recovery_pct")
        cv.banner(T2, "2", f"At a strict cut, recovery ranges from {rec(lo_r)}% ({lo_r['label']}) to {rec(hi_r)}% ({hi_r['label']})" if len(KC) > 1 else f"At a strict cut {hi_r['label']} recover {rec(hi_r)}% of within-family pairs")
        dx0, dlet = .5, "C"
    cv.title(dlet, "Pairs recovered", T2 + 10, dx0, dx=4)
    ax = cv.ax(dx0 + 19 if has_c else 25, T2 + 16, 31 if has_c else 60, 37)
    for q, k in enumerate(KC):
        y = nk - 1 - q; r = k["row"]; v, lo, hi = float(r.recovery_pct), float(r.recovery_ci95_low), float(r.recovery_ci95_high)
        ax.barh(y, v, color=COL[k["key"]], height=.62)
        if np.isfinite(lo):
            ax.errorbar(v, y, xerr=[[max(v - lo, 0)], [max(hi - v, 0)]], color=S.NAVY, lw=.6, capsize=1.2)
        t = f"{v:.0f}% [{lo:.0f}, {hi:.0f}]" if np.isfinite(lo) else f"{v:.0f}%"
        ax.text(1.5, y + .44, t, fontsize=sz["small"], va="center"); L(dlet, "recovery [CI] " + k["label"], t, SRC_T + " recovery_pct, recovery_ci95_low, recovery_ci95_high")
    ax.set_yticks(range(nk)); ax.set_yticklabels([SHORT[k["key"]] for k in KC[::-1]], fontsize=sz["small"]); ax.set_xlim(0, 105); ax.set_ylim(-.5, nk - .15); ax.set_xticks([0, 50, 100]); ax.tick_params(axis="x", labelsize=sz["tick"], length=1.5)
    ax.set_xlabel(f"% of {nw:,} within-family pairs above\nthe 99.9th pct of between-family pairs", fontsize=sz["small"], labelpad=1.5, linespacing=1.1)
    desc = ("Port of the published Extended Data figure. A rank AUC of within- against between-family pairs per marker class with the 95% interval of a bootstrap over families. B GRM values of within- and between-family pairs "
            "(medians printed). C GRM of each SV class against the GRM of the biSNPs outside SV regions, all pairs, with the Pearson correlation. Last panel: share of within-family pairs above the 99.9th percentile of the "
            "between-family pairs, with its bootstrap interval. The families come from the relatedness source: concordance between marker classes, not pedigree truth.")
    return cv.fig, "Cross-marker concordance of kinship", desc, "; ".join(notes)


def fig_overview(run, plt):
    """Overview for the front page: one kinship matrix per marker class (the kinship row of the
    published figure, without the masked SNP class) and the network of families below it."""
    KC, O = kin_classes(run, masked=False)
    KC = KC[:4]; nd = network_data(run)
    if not KC and nd is None:
        return None
    sz = dict(S.SZ_ED); sz.update(title=8.5, small=6.5, letter=10.0); W = 183.0; notes = []
    side = 40.0; pitch = 44.2
    hk = 62.0 if KC else 0.0
    nh = 62.0 if nd is not None else 0.0
    cv = S.canvas(plt, W, hk + nh + 1.0, sz); letters = None        # front-page figure: titles without panel letters
    if KC:
        x0 = 1.0 + (4 - len(KC)) * pitch / 2
        n_asg = kinship_row(cv, run, KC, O, 5.0, letters, sz, x0=x0, pitch=pitch, side=side)
        n = len(O); nfam = O.family[O.family != "unassigned"].nunique()
        cv.text(1, 58.5, f"VanRaden kinship per marker class, {n} {noun(run)} in family order (largest family first); dotted lines separate the {n_asg} {noun(run)} in {nfam} families from the {n - n_asg} unassigned",
                fs=6.0, color=S.INK)
    else:
        notes.append("kinship matrices omitted: no genotype input")
    for k in ("SNP", "DEL", "DUP", "HAPLOTYPE"):
        if KC and not any(c["klass"] == k for c in KC):
            notes.append(f"no {S.class_label(k) if k != 'HAPLOTYPE' else 'structural-haplotype'} matrix: that input is absent")
    if nd is not None:
        aud = bool(run.s.get("audit", {}).get("ran"))
        cv.title("", "Audited shared-parent family network" if aud else "Nominated family network (not audited)", hk + 5.0, 1.0)
        draw_network(cv, run, 1.0, hk + 8.5, W - 2.0, nh - 9.5, sz, nd=nd)
    else:
        notes.append("network omitted: no link between families")
    desc = ("Top: VanRaden kinship matrix of every marker class, individuals in family order, with the AUC of within- against between-family pairs and the number of markers. "
            "Bottom: network of families (nodes sized by family size) and links coloured by audit outcome, one box per component.")
    return cv.fig, "Overview: kinship per marker class and the network of families", desc, "; ".join(notes)
