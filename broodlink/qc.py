"""Step 0: quality control of every input. QC flags; it removes nothing unless asked.

REIMPLEMENTED FROM METHODS (no QC script of the original analysis was archived). Rules that are
new in BroodLink are marked as such in the rule text printed with every flag.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import kruskal

from .io import allele_freq

POINT = ("SNP", "INDEL", "SNP_INDEL")


def per_sample(sets, mask, depth):
    rows = []
    span_mb = (mask.span() / 1e6) if mask.active else ((mask.assembly() or 0) / 1e6 or np.nan)
    for ms in sets:
        G = ms.G; called = G >= 0
        n = called.sum(0); het = (G == 1).sum(0); hom = (G == 2).sum(0)
        carriers = (G > 0).sum(1)
        single = ((G > 0) & (carriers == 1)[:, None]).sum(0)
        length = np.maximum(ms.end - ms.pos, 0)
        for i, s in enumerate(ms.samples):
            r = {"sample": s, "input": ms.name, "marker_class": ms.klass, "n_markers": ms.n, "n_genotyped": int(n[i]),
                 "call_rate": n[i] / max(ms.n, 1), "heterozygosity": het[i] / max(n[i], 1),
                 "het_hom_alt_ratio": het[i] / hom[i] if hom[i] else np.nan, "n_singletons": int(single[i]),
                 "het_sites_per_callable_Mb": het[i] / span_mb if span_mb and ms.klass in POINT else np.nan}
            if ms.klass not in POINT and ms.klass != "HAPLOTYPE":
                r["n_sv_calls_carried"] = int((G[:, i] > 0).sum()); r["sv_span_carried_bp"] = int(length[G[:, i] > 0].sum())
            rows.append(r)
    T = pd.DataFrame(rows)
    if depth is not None and len(T):
        d = depth.rename(columns={"id": "sample"})
        T = T.merge(d, on="sample", how="left")
    return T


def per_input(sets, genome):
    rows, chrom_rows = [], []
    for ms in sets:
        af = allele_freq(ms.G); maf = np.minimum(af, 1 - af)
        r = {"input": ms.name, "marker_class": ms.klass, "flag": "--" + ms.flag, "file": ms.source, "n_samples": len(ms.samples), "n_markers": ms.n,
             "filters": "; ".join(f"{s}: {b} -> {a}" for s, b, a in ms.filters) or "none",
             "af_median": float(np.nanmedian(af)) if ms.n else np.nan, "share_maf_below_0.05": float(np.nanmean(maf < 0.05)) if ms.n else np.nan,
             "share_monomorphic": float(np.nanmean(maf == 0)) if ms.n else np.nan, "mean_call_rate": float((ms.G >= 0).mean()) if ms.n else np.nan,
             "notes": "; ".join(ms.notes)}
        rows.append(r)
        for c, k in pd.Series(ms.chrom).value_counts(sort=False).items():
            L = genome.get(c) if genome else None
            chrom_rows.append({"input": ms.name, "chrom": c, "n_markers": int(k), "chrom_length_bp": L if L else np.nan,
                               "markers_per_Mb": k / (L / 1e6) if L else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(chrom_rows)


def _concordance(ms, max_markers=20000, seed=1):
    G = ms.G
    if G.shape[0] > max_markers:
        G = G[np.sort(np.random.default_rng(seed).choice(G.shape[0], max_markers, replace=False))]
    C = (G >= 0).astype(np.float32); both = C.T @ C
    match = sum(((G == g).astype(np.float32)).T @ (G == g).astype(np.float32) for g in (0, 1, 2))
    return match / np.maximum(both, 1), both


def sample_checks(sets, R, PS, samples_meta, all_samples, cfg):
    """Returns (flag table, identical-pair table, group-effect table)."""
    q = cfg["qc"]; flags, eff = [], []
    for ms in sets:                                              # samples missing from an input
        for s in sorted(set(all_samples) - set(ms.samples)):
            flags.append({"sample": s, "reason": "missing from an input", "rule": "sample absent from the input's sample list", "value": ms.name, "auto_exclude": False})
    for name, g in (PS.groupby("input", sort=False) if len(PS) else []):     # low coverage
        med = g.call_rate.median(); thr = q["low_coverage_factor"] * med
        for r in g[g.call_rate < thr].itertuples():
            flags.append({"sample": r.sample, "reason": "low coverage", "rule": f"call rate below {q['low_coverage_factor']} x cohort median ({thr:.3f}) in {name}",
                          "value": f"{r.call_rate:.3f}", "auto_exclude": True})
        h = g.heterozygosity; mad = 1.4826 * (h - h.median()).abs().median()
        if mad > 0 and g.marker_class.iloc[0] in POINT:
            for r in g[h > h.median() + q["contamination_mad"] * mad].itertuples():
                flags.append({"sample": r.sample, "reason": "possible contamination", "value": f"{r.heterozygosity:.4f}", "auto_exclude": False,
                              "rule": f"heterozygosity above cohort median + {q['contamination_mad']} MAD in {name} (new in BroodLink; the original fixed no number)"})
    ident = []
    if R is not None and len(R):
        long = pd.concat([R[["id1", "theta"]].rename(columns={"id1": "s"}), R[["id2", "theta"]].rename(columns={"id2": "s"})])
        medrel = long.groupby("s").theta.median(); mad = 1.4826 * (medrel - medrel.median()).abs().median()
        if mad > 0:
            for s, v in medrel[medrel > medrel.median() + q["contamination_mad"] * mad].items():
                flags.append({"sample": s, "reason": "possible contamination", "value": f"{v:.4f}", "auto_exclude": False,
                              "rule": f"median relatedness to the cohort above cohort median + {q['contamination_mad']} MAD (new in BroodLink)"})
        for r in R[R.theta >= q["identical_theta"]].itertuples():
            ident.append({"sample_1": r.id1, "sample_2": r.id2, "theta": r.theta, "genotype_concordance": np.nan, "n_markers_compared": np.nan,
                          "rule": f"theta >= {q['identical_theta']}"})
    dense = [m for m in sets if m.klass != "HAPLOTYPE" and m.D is None]
    if dense:                                                    # hard genotypes only: likelihood calls are not comparable at low coverage
        ms = max(dense, key=lambda m: m.n)
        conc, both = _concordance(ms)
        iu, ju = np.triu_indices(len(ms.samples), 1)
        hit = (conc[iu, ju] >= q["identical_concordance"]) & (both[iu, ju] >= 100)
        have = {(d["sample_1"], d["sample_2"]) for d in ident}
        for i, j in zip(iu[hit], ju[hit]):
            a, b = sorted((ms.samples[i], ms.samples[j]))
            if (a, b) in have:
                for d in ident:
                    if (d["sample_1"], d["sample_2"]) == (a, b):
                        d["genotype_concordance"] = float(conc[i, j]); d["n_markers_compared"] = int(both[i, j])
                continue
            ident.append({"sample_1": a, "sample_2": b, "theta": np.nan, "genotype_concordance": float(conc[i, j]), "n_markers_compared": int(both[i, j]),
                          "rule": f"genotype concordance >= {q['identical_concordance']} in {ms.name} (new in BroodLink)"})
    ID = pd.DataFrame(ident, columns=["sample_1", "sample_2", "theta", "genotype_concordance", "n_markers_compared", "rule"])
    cr = PS.groupby("sample").call_rate.mean() if len(PS) else pd.Series(dtype=float)
    for r in ID.itertuples():                                    # pairs with identical genotypes: flag both, propose the lower call rate
        lo = min((r.sample_1, r.sample_2), key=lambda s: (cr.get(s, 0), s))
        for s in (r.sample_1, r.sample_2):
            flags.append({"sample": s, "reason": "pair of samples with identical genotypes", "rule": r.rule + "; with --auto-exclude the sample of the pair with the lower call rate is removed",
                          "value": f"pair with {r.sample_2 if s == r.sample_1 else r.sample_1}", "auto_exclude": s == lo})
    if samples_meta is not None and len(PS):                     # sex / batch effects (new in BroodLink)
        for col in [c for c in samples_meta.columns if c in ("sex", "batch", "group")]:
            lab = samples_meta.set_index("id")[col]
            for name, g in PS.groupby("input", sort=False):
                g = g.assign(lab=g["sample"].map(lab)).dropna(subset=["lab"])
                groups = [x for _, x in g.groupby("lab") if len(x) >= 3]
                for stat in ("call_rate", "heterozygosity"):
                    row = {"variable": col, "input": name, "statistic": stat, "n_groups": len(groups), "test": "Kruskal-Wallis (new in BroodLink)", "p_value": np.nan, "flag": ""}
                    if len(groups) >= 2 and g[stat].nunique() > 1:
                        row["p_value"] = float(kruskal(*[x[stat].to_numpy() for x in groups]).pvalue)
                        row["flag"] = "possible effect (P < 0.01)" if row["p_value"] < 0.01 else ""
                    eff.append(row)
    FL = pd.DataFrame(flags, columns=["sample", "reason", "rule", "value", "auto_exclude"])
    EF = pd.DataFrame(eff, columns=["variable", "input", "statistic", "n_groups", "test", "p_value", "flag"])
    return FL, ID, EF


def exclusions(FL, cfg, auto):
    """Exclusion table: flagged samples with the action taken. Default: flag, do not remove."""
    conf = set(map(str, cfg["samples"].get("exclude_confirmed") or []))
    rows = []
    for r in FL.itertuples():
        if r.sample in conf:
            act = "removed (confirmed in the config)"
        elif auto and r.auto_exclude:
            act = "removed (--auto-exclude)"
        else:
            act = "flagged, kept"
        rows.append({"sample": r.sample, "reason": r.reason, "rule": r.rule, "value": r.value, "action": act})
    for s in sorted(conf - set(FL["sample"])):
        rows.append({"sample": s, "reason": "confirmed by the user", "rule": "samples.exclude_confirmed in the config", "value": "", "action": "removed (confirmed in the config)"})
    X = pd.DataFrame(rows, columns=["sample", "reason", "rule", "value", "action"])
    return X, sorted(set(X[X.action.str.startswith("removed")]["sample"]))


def relatedness_qc(R, classes, cfg):
    q = cfg["qc"]
    counts = pd.Series(classes).value_counts().rename_axis("class").reset_index(name="n_pairs")
    has_k = R.k1.notna().any()
    po = int(((R.k1 >= q["po_k1_min"]) & (R.k2 <= q["po_k2_max"])).sum()) if has_k else None
    desc = R.theta.describe(percentiles=[0.5, 0.9, 0.99])
    summ = {"n_pairs": int(len(R)), "theta_median": float(desc["50%"]), "theta_p99": float(desc["99%"]), "theta_max": float(desc["max"]),
            "has_k2": bool(R.k2.notna().any()), "k2_median_first_degree": float(R[R.theta >= cfg["relatedness"]["theta_first"]].k2.median()) if R.k2.notna().any() else None,
            "parent_offspring_pairs": po,
            "parent_offspring_rule": f"k1 >= {q['po_k1_min']} and k2 <= {q['po_k2_max']}" if has_k else "not evaluated: no k1 / k2 (GRM-based relatedness)"}
    return counts, summ
