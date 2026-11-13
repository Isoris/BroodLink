"""Parameter sweeps. Every sweep reports a DISCRIMINATOR (a measure that tells good settings from
bad: agreement with the reference families, recovery of the reference links or of the planted
intervals, stability of the nominated set across seeds), never just "the output changed", and a
PROVENANCE note saying where the default and each tested value come from.

`broodlink sweep` runs on the synthetic cohort, where the reference is the simulated truth. The
same functions run on any data set given as a dictionary (see `load_synthetic`); on real data the
reference is a reference run, so the discriminator measures agreement, not correctness.
"""
from __future__ import annotations
import copy
import itertools
import os
import time

import numpy as np
import pandas as pd

from . import audit, barcode, families, inversions, network, nomination, relatedness
from .evaluate import family_map, family_scores, link_scores, link_set
from .io import log, read_haplotypes, read_intervals_bed, read_relatedness, read_variants, read_fai
from .mask import apply_mask, build_mask

PROVENANCE = {
    "families_theta": {"parameter": "families.theta_min", "default": "0.177: original toolkit common.py line 12 (THETA_FD) and get_FS.py line 50; Methods",
                       "values": "0.10 to 0.25 in steps: chosen here to bracket the default; 0.25 is the expectation for full sibs"},
    "families_k2": {"parameter": "families.k2_min", "default": "0.15: common.py line 15 (K2_FS); Methods",
                    "values": "0.12: older value in sibling_build.py line 7; 0 (gate off), 0.05, 0.10, 0.20, 0.25: chosen here"},
    "nomination_theta": {"parameter": "nomination.theta_min / theta_max", "default": "theta_min 0.0884, no upper bound: sp_lib.py line 40 (THETA_HS, common.py line 13); Methods ('second degree or closer')",
                         "values": "0.0442: third-degree bound, common.py line 14; upper bound 0.177: the BroodLink task description (not the original code); 0.0625, 0.10, 0.125: chosen here"},
    "audit_min_family_size": {"parameter": "audit.min_family_size", "default": "4: run_audit.py lines 49-50; Methods, 'Eligibility'",
                              "values": "3 and 2: the extension and provisional tiers of the original; 5 and 6: chosen here"},
    "permutations": {"parameter": "nomination.n_permutations (and seed)", "default": "20,000, seed 1: the primary run of the original analysis (Methods)",
                     "values": "1,500: run_nomination.py and the original seed sweep (nomination_seed_sweep.tsv); 200 to 5,000 and seeds 1-5: chosen here"},
    "barcode_posterior": {"parameter": "barcode.posterior_call", "default": "0.90: get_INV_barcode_SP_get_LP.py line 27 (--post) and lp01_profile.py line 31",
                          "values": "0.95: the stricter count reported in the Methods; 0.50 to 0.99 otherwise: chosen here"},
    "grm_kinship": {"parameter": "families.theta_min applied to GRM kinship", "default": "0.177 on GRM / 2: grm_verify.py line 23 and relatedness_atlas.py cluster_sibships (fs_lo); "
                    "the GRM benchmark itself fixed no threshold", "values": "0.10 to 0.225: chosen here"},
    "scan_snp": {"parameter": "inversions.snp_scan.window_snps / z_outlier", "default": "100 SNPs, |z| >= 3: Methods (local PCA in the style of lostruct); the scan script was not archived, the scan is reimplemented from Methods",
                 "values": "15, 25, 50 SNPs and |z| 2.0, 2.5: chosen here (the synthetic genome has 750 SNPs per chromosome)"},
    "scan_del": {"parameter": "inversions.del_scan.window_bp / hit_quantile", "default": "1 Mb windows, 90th percentile: step2_inversions.py",
                 "values": "0.5 and 2 Mb, 85th and 95th percentile: chosen here"},
}
DISCRIMINATOR = {
    "families_theta": "families identical to the reference; pair precision and recall",
    "families_k2": "families identical to the reference; pair precision and recall",
    "nomination_theta": "reference links nominated (of the reference links) and family pairs nominated that are not reference links",
    "audit_min_family_size": "reference links that can be audited; supported links that are reference links / that are not",
    "permutations": "stability of the nominated set across seeds (mean pairwise Jaccard index) and reference links nominated",
    "barcode_posterior": "share of called states equal to the reference parent genotype, against the share of states called",
    "grm_kinship": "families identical to the reference; pair precision and recall, per marker class",
    "scan_snp": "planted intervals overlapped by a scan candidate and candidates that overlap no planted interval",
    "scan_del": "planted intervals overlapped by a scan candidate and candidates that overlap no planted interval",
}


def load_synthetic(data, cfg):
    j = lambda f: os.path.join(data, f)
    genome = read_fai(j("genome.fai")); mask = build_mask([j("callable.bed")], None, genome)
    R, _ = read_relatedness(j("relatedness.tsv"))
    T = pd.read_csv(j("truth/truth_families.tsv"), sep="\t", dtype=str).fillna("")
    L = pd.read_csv(j("truth/truth_links.tsv"), sep="\t", dtype=str)
    snp = apply_mask(read_variants(j("snps.beagle.gz"), "bisnp", cfg)[0], mask, cfg)[0]
    svs = [apply_mask(m, mask, cfg)[0] for m in read_variants(j("sv.vcf.gz"), "gt", cfg)]
    ids = sorted(T.individual)
    return {"name": "synthetic cohort (reference = simulated truth)", "ids": ids, "R": R, "ref_fam": dict(zip(T.individual, T.family)),
            "ref_links": {frozenset((a, b)) for a, b in zip(L.family_a, L.family_b)}, "snp": snp, "svs": svs,
            "hap": read_haplotypes(j("haplotypes.tsv"), set(ids)), "planted": read_intervals_bed(j("inversions.bed")),
            "parent_truth": pd.read_csv(j("truth/truth_parent_haplotypes.tsv"), sep="\t", dtype={"parent": str}).set_index("parent"),
            "link_parent": {frozenset((a, b)): p for a, b, p in zip(L.family_a, L.family_b, L.shared_parent)},
            "grm_sets": {"SNP (hard calls)": snp, "DEL": [m for m in svs if m.klass == "DEL"][0]}}


def _cfg(cfg, **kv):
    c = copy.deepcopy(cfg)
    for k, v in kv.items():
        a, b = k.split("__"); c[a][b] = v
    return c


def _fam(D, cfg, R=None):
    return families.reconstruct(D["R"] if R is None else R, D["ids"], cfg)[0]


def _nom(D, M, cfg, seed=1, n_perm=2000):
    return nomination.nominate(D["R"], dict(zip(M.individual, M.unit)), cfg, seed, n_perm=n_perm)


def sweep_families(D, cfg, quick):
    rows = []
    for name, key, vals, dflt in (("families_theta", "theta_min", [0.10, 0.125, 0.15, 0.177, 0.20, 0.225, 0.25], 0.177),
                                  ("families_k2", "k2_min", [0.0, 0.05, 0.10, 0.12, 0.15, 0.20, 0.25], 0.15)):
        for v in vals:
            M = _fam(D, _cfg(cfg, **{"families__" + key: v}))
            rows.append({"sweep": name, "setting": f"{key} = {v}", "value": v, "is_default": v == dflt, **family_scores(M, D["ref_fam"])})
    return rows


def sweep_nomination_theta(D, cfg, quick):
    rows = []; M = _fam(D, cfg); fm = family_map(M, D["ref_fam"]); npm = 1000 if quick else 5000
    for tmax in (None, 0.177):
        for v in [0.0442, 0.0625, 0.0884, 0.10, 0.125]:
            nom = _nom(D, M, _cfg(cfg, nomination__theta_min=v, nomination__theta_max=tmax), n_perm=npm)
            sc = link_scores(link_set(nom, fm, nom.nominated & (nom.unassigned_sides == 0)), D["ref_links"])
            rows.append({"sweep": "nomination_theta", "setting": f"theta_min = {v}, theta_max = {'none' if tmax is None else tmax}", "value": v,
                         "is_default": v == 0.0884 and tmax is None, "n_nominated_all_units": int(nom.nominated.sum()), **sc, "n_permutations": npm})
    return rows


def sweep_audit_size(D, cfg, quick):
    rows = []; M = _fam(D, cfg); fm = family_map(M, D["ref_fam"]); nom = _nom(D, M, cfg, n_perm=1000 if quick else 5000)
    unit_of = dict(zip(M.individual, M.unit))
    for v in [2, 3, 4, 5, 6]:
        E, s = audit.run_audit(nom, unit_of, D["svs"], _cfg(cfg, audit__min_family_size=v), 1)
        aud = E.auditable
        sup = link_set(E, fm, aud & E.link_class.str.startswith("supported"))
        g1 = link_set(E, fm, aud & E.nominated & (E.gate1 == 1))
        rej = link_set(E, fm, aud & E.link_class.str.startswith("rejected"))
        rows.append({"sweep": "audit_min_family_size", "setting": f"min_family_size = {v}", "value": v, "is_default": v == 4, "auditable_pairs": int(aud.sum()),
                     "reference_links_auditable": len(link_set(E, fm, aud) & D["ref_links"]), "links_reference": len(D["ref_links"]),
                     "gate1_pass_true": len(g1 & D["ref_links"]), "gate1_pass_false": len(g1 - D["ref_links"]), "rejected_true": len(rej & D["ref_links"]),
                     "supported_true": len(sup & D["ref_links"]), "supported_false": len(sup - D["ref_links"])})
    return rows


def sweep_permutations(D, cfg, quick, counts=None, seeds=None):
    rows = []; M = _fam(D, cfg); fm = family_map(M, D["ref_fam"])
    counts = counts or ([200, 1000, 5000] if quick else [200, 500, 1000, 1500, 2000, 5000, 20000])
    seeds = seeds or ([1, 2, 3] if quick else [1, 2, 3, 4, 5])
    for n in counts:
        sets, fsets, t0 = [], [], time.time()
        for sd in seeds:
            nom = _nom(D, M, cfg, seed=sd, n_perm=n)
            sets.append({(a, b) for a, b in zip(nom.unit_a[nom.nominated], nom.unit_b[nom.nominated])})
            fsets.append(link_set(nom, fm, nom.nominated & (nom.unassigned_sides == 0)))
        jac = [len(a & b) / len(a | b) if (a | b) else 1.0 for a, b in itertools.combinations(sets, 2)]
        rows.append({"sweep": "permutations", "setting": f"n_permutations = {n}", "value": n, "is_default": n == 20000, "n_seeds": len(seeds),
                     "nominated_min": min(map(len, sets)), "nominated_max": max(map(len, sets)), "jaccard_across_seeds_mean": round(float(np.mean(jac)), 4),
                     "jaccard_across_seeds_min": round(float(np.min(jac)), 4), "links_true_min": min(len(f & D["ref_links"]) for f in fsets),
                     "links_false_max": max(len(f - D["ref_links"]) for f in fsets), "links_reference": len(D["ref_links"]),
                     "seconds_per_run": round((time.time() - t0) / len(seeds), 2)})
    return rows


def sweep_barcode(D, cfg, quick):
    if D.get("hap") is None or D.get("parent_truth") is None:
        return []
    M = _fam(D, cfg); fm = family_map(M, D["ref_fam"]); members = {}
    for s, u in zip(M.individual, M.unit):
        members.setdefault(u, []).append(s)
    inv = {}
    for f, t in fm.items():
        inv.setdefault(t, f)
    links = [(inv[a], inv[b], "reference link") for a, b in map(sorted, D["ref_links"]) if a in inv and b in inv]
    B = barcode.barcode(D["hap"].subset_samples(D["ids"]), members, links, {}, cfg)
    code = {"0/0": 0, "0/1": 1, "1/1": 2}
    truth = [int(D["parent_truth"].loc[D["link_parent"][frozenset((fm[r.families.split(",")[0]], fm[r.families.split(",")[1]]))], r.interval]) for r in B.itertuples()]
    ok = np.array([code[g] for g in B.best_genotype]) == np.array(truth)
    rows = []
    for v in [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]:
        c = (B.posterior >= v).to_numpy()
        rows.append({"sweep": "barcode_posterior", "setting": f"posterior_call = {v}", "value": v, "is_default": v == 0.9, "n_states": len(B), "n_called": int(c.sum()),
                     "share_called": round(float(c.mean()), 4), "accuracy_of_called": round(float(ok[c].mean()), 4) if c.any() else np.nan, "n_called_wrong": int((~ok[c]).sum())})
    return rows


def sweep_grm(D, cfg, quick):
    rows = []
    for label, ms in (D.get("grm_sets") or {}).items():
        if isinstance(ms, pd.DataFrame):
            R, used = ms, np.nan                                 # a precomputed GRM pair table
        else:
            R, used = relatedness.relatedness_from_markers(ms.subset_samples(D["ids"]), cfg)
        for v in [0.10, 0.125, 0.15, 0.177, 0.20, 0.225]:
            M = _fam(D, _cfg(cfg, families__theta_min=v), R=R)
            rows.append({"sweep": "grm_kinship", "setting": f"{label}: kinship >= {v}", "value": v, "is_default": v == 0.177, "marker_class": label, "n_markers": used,
                         **family_scores(M, D["ref_fam"])})
    return rows


def _scan_scores(C, planted):
    ov = inversions.overlap_table(planted, C) if len(C) else None
    rec = int((ov.category == "supplied interval recovered by the scan").sum()) if ov is not None else 0
    extra = int((ov.category == "extra candidate interval found by the scan").sum()) if ov is not None else 0
    return {"planted_recovered": rec, "planted_total": len(planted), "extra_candidates": extra, "n_candidates": int(len(C))}


def sweep_scan(D, cfg, quick):
    rows = []
    if D.get("planted") is None:
        return rows
    if D.get("snp") is not None:
        for w in ([25, 100] if quick else [15, 25, 50, 100]):
            for z in ([2.5, 3.0] if quick else [2.0, 2.5, 3.0]):
                c = copy.deepcopy(cfg); c["inversions"]["snp_scan"].update({"window_snps": w, "z_outlier": z})
                C, W, notes, sec, _st = inversions.scan_snps(D["snp"], c)
                rows.append({"sweep": "scan_snp", "setting": f"window_snps = {w}, z_outlier = {z}", "value": w, "is_default": w == 100 and z == 3.0,
                             **_scan_scores(C, D["planted"]), "chromosomes_skipped": len(notes), "seconds": round(sec, 1)})
    dels = [m for m in D.get("svs", []) if m.klass == "DEL"]
    if dels:
        for w in [500000, 1000000, 2000000]:
            for q in [0.85, 0.90, 0.95]:
                c = copy.deepcopy(cfg); c["inversions"]["del_scan"].update({"window_bp": w, "hit_quantile": q, "core_quantile": max(q, 0.95)})
                C, W, notes, sec = inversions.scan_dels(dels[0], c)
                rows.append({"sweep": "scan_del", "setting": f"window_bp = {w}, hit_quantile = {q}", "value": w, "is_default": w == 1000000 and q == 0.90,
                             **_scan_scores(C, D["planted"]), "windows_tested": int(len(W)), "seconds": round(sec, 1)})
    return rows


SWEEPS = {"families": sweep_families, "nomination_theta": sweep_nomination_theta, "audit_min_family_size": sweep_audit_size, "permutations": sweep_permutations,
          "barcode_posterior": sweep_barcode, "grm_kinship": sweep_grm, "scan": sweep_scan}


def run_sweeps(data, out, cfg, quick=False, only=None, D=None):
    os.makedirs(out, exist_ok=True)
    D = D or load_synthetic(data, cfg)
    allrows = []
    for name, fn in SWEEPS.items():
        if only and name not in only:
            continue
        t0 = time.time(); rows = fn(D, cfg, quick)
        log(f"sweep {name}: {len(rows)} settings in {time.time() - t0:.0f} s")
        allrows += rows
        for sw in sorted({r["sweep"] for r in rows}):                # partial results are saved after every sweep
            T = pd.DataFrame([r for r in rows if r["sweep"] == sw])
            T["discriminator"] = DISCRIMINATOR[sw]
            T.to_csv(os.path.join(out, f"sweep_{sw}.tsv"), sep="\t", index=False, lineterminator="\n")
    P = pd.DataFrame([{"sweep": k, **v, "discriminator": DISCRIMINATOR[k]} for k, v in PROVENANCE.items()])
    P.to_csv(os.path.join(out, "sweep_provenance.tsv"), sep="\t", index=False, lineterminator="\n")
    with open(os.path.join(out, "sweep_report.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(report_md(allrows, D["name"]))
    try:
        figure(allrows, os.path.join(out, "sweep_summary"), D["name"])
    except Exception as exc:                                         # the tables are the result; the figure is a view of them
        log(f"WARNING: sweep figure failed ({type(exc).__name__}: {exc})")
    log(f"sweeps written to {out}")
    return allrows


def report_md(rows, name):
    L = [f"Data: {name}.", ""]
    for sw in dict.fromkeys(r["sweep"] for r in rows):
        T = pd.DataFrame([r for r in rows if r["sweep"] == sw]).drop(columns=["sweep", "value"])
        T["is_default"] = np.where(T.is_default, "default", "")
        T = T.rename(columns={"is_default": ""})
        p = PROVENANCE[sw]
        L += [f"### {sw} ({p['parameter']})", "", f"Discriminator: {DISCRIMINATOR[sw]}.", "", f"Provenance. Default: {p['default']}. Tested values: {p['values']}.", "",
              "| " + " | ".join(T.columns) + " |", "|" + "---|" * T.shape[1]]
        L += ["| " + " | ".join("" if (isinstance(x, float) and np.isnan(x)) else str(x) for x in r) + " |" for r in T.itertuples(index=False)]
        L.append("")
    return "\n".join(L) + "\n"


def figure(rows, stem, name):
    from .figures import BLUE, GREY, INK, ORANGE, RED, _plt
    plt = _plt()
    A = pd.DataFrame(rows)
    panels = [("families_theta", [("families_identical", BLUE, "families identical to the reference")], "families_reference"),
              ("families_k2", [("families_identical", BLUE, "families identical to the reference")], "families_reference"),
              ("nomination_theta", [("links_true", BLUE, "reference links nominated"), ("links_false", RED, "other family pairs nominated")], "links_reference"),
              ("audit_min_family_size", [("reference_links_auditable", GREY, "reference links auditable"), ("supported_true", BLUE, "supported, reference link"), ("supported_false", RED, "supported, not a reference link")], "links_reference"),
              ("permutations", [("jaccard_across_seeds_mean", BLUE, "nominated set: mean Jaccard across seeds"), ("jaccard_across_seeds_min", ORANGE, "minimum Jaccard")], None),
              ("barcode_posterior", [("accuracy_of_called", BLUE, "called states equal to the reference"), ("share_called", GREY, "share of states called")], None),
              ("grm_kinship", [("pair_recall", BLUE, "pair recall"), ("pair_precision", ORANGE, "pair precision")], None),
              ("scan_snp", [("planted_recovered", BLUE, "planted intervals recovered"), ("extra_candidates", RED, "extra candidates")], "planted_total"),
              ("scan_del", [("planted_recovered", BLUE, "planted intervals recovered"), ("extra_candidates", RED, "extra candidates")], "planted_total")]
    panels = [p for p in panels if len(A) and (A.sweep == p[0]).any()]
    nc = 3; nr = int(np.ceil(len(panels) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.6 * nc, 3.6 * nr), constrained_layout=True, squeeze=False)
    for ax in axs.ravel()[len(panels):]:
        ax.set_axis_off()
    for ax, (sw, series, total) in zip(axs.ravel(), panels):
        T = A[A.sweep == sw].reset_index(drop=True); x = np.arange(len(T))
        mk = ["o", "s", "^"]
        for k, (col, c, lab) in enumerate(series):
            ax.plot(x, T[col], marker=mk[k], ms=4, lw=1.2 if sw not in ("grm_kinship", "nomination_theta", "scan_snp", "scan_del") else 0, color=c, label=lab)
        if total and total in T:
            ax.axhline(T[total].iloc[0], color=INK, lw=0.6, ls=":", label=f"reference total ({int(T[total].iloc[0])})")
        for i in np.flatnonzero(T.is_default.to_numpy()):
            ax.axvspan(i - 0.4, i + 0.4, color="#EEEEEE", zorder=0, label="default" if i == np.flatnonzero(T.is_default.to_numpy())[0] else None)
        ax.set_xticks(x); ax.set_xticklabels([s.replace(", ", "\n").replace(": ", "\n") for s in T.setting], rotation=90, fontsize=5.5)
        ax.set_title(f"{sw}\n{PROVENANCE[sw]['parameter']}", fontsize=8); ax.legend(fontsize=6)
        ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
    fig.suptitle(f"Parameter sweeps with their discriminators - {name}", fontsize=10)
    fig.canvas.draw(); fig.set_layout_engine("none")
    fig.savefig(stem + ".png", dpi=200); fig.savefig(stem + ".pdf"); plt.close(fig)
