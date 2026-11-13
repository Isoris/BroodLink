"""Tables (TSV + XLSX with a definitions sheet), statistics.tsv, qc_report.md, run_report.md with the
capability table, and methods_text.md (a Methods paragraph generated from the config and the steps
that actually ran)."""
from __future__ import annotations
import datetime
import os

import pandas as pd

from . import __version__, families
from .io import log, write_tsv
from .pipeline import FOLDER, STEPS

STEP_TITLE = {"qc": "0. Quality control", "relatedness": "1. Pairwise relatedness", "families": "2. Full-sib families",
              "nomination": "3. Half-sib nomination", "audit": "4. Audit", "network": "5. Network, latent parents, slot limit",
              "haplotypes": "6a. Structural haplotypes (genotyping and scan)", "barcode": "6b. Barcode of the shared parent",
              "family_effect": "7. Family effect (eta squared)", "concordance": "8. Cross-marker concordance"}

# (step, key, table name, sheet name)
TABLES = [("qc", "flags", "T01_qc_exclusion_table", "T01 QC flags"), ("families", "families", "T02_families", "T02 families"),
          ("families", "membership", "T03_family_membership", "T03 membership"), ("network", "links", "T04_links", "T04 links"),
          ("audit", "gate2", "T05_audit_gate2", "T05 Gate 2"), ("network", "latent_parents", "T06_latent_parents", "T06 latent parents"),
          ("network", "slot", "T07_slot_limit", "T07 slot limit"), ("haplotypes", "summary", "T08_structural_haplotypes", "T08 haplotypes"),
          ("haplotypes", "tags", "T09_arrangement_tag_variants", "T09 tag variants"), ("barcode", "barcode", "T10_barcode", "T10 barcode"),
          ("family_effect", "summary", "T11_family_effect", "T11 family effect"), ("concordance", "summary", "T12_cross_marker_concordance", "T12 concordance")]

COLUMN_DEFS = {
    "theta": "kinship coefficient of the pair (0.25 full sibs, 0.125 half sibs in expectation)",
    "k0": "probability that the pair shares no allele identical by descent", "k1": "probability of sharing one allele identical by descent",
    "k2": "probability of sharing two alleles identical by descent (separates full sibs from parent-offspring)",
    "unit": "family label, or S_<id> for an unassigned individual", "family": "family label (relatedness-based group, not a pedigree record)",
    "membership": "core = in the largest complete set of full-sib pairs of the family; attached = at least two full-sib pairs into the core; unresolved member = other member of the block",
    "closure": "share of within-family pairs that pass the full-sib thresholds", "n": "number of individuals",
    "fraction": "share of cross-unit pairs of individuals at second-degree kinship or closer", "p_perm": "permutation P of the fraction (unit labels permuted among individuals)",
    "q_value": "Benjamini-Hochberg q over all unit pairs", "nominated": "q below the configured threshold: a candidate link (screening only)",
    "min_n": "size of the smaller unit of the pair", "auditable": "both units are families of at least the minimum size",
    "LLR_one": "mean per-marker log-likelihood ratio, one shared parent against none", "LLR_both": "mean per-marker log-likelihood ratio, both parents shared against none",
    "gate1": "1 when LLR_one > 0 and LLR_both < 0 (one shared parent is not excluded)", "link_class": "supported / Gate-1 only / provisional / rejected / nominated, not audited",
    "evidence": "SV classes whose sharing exceeds the Gate-2 bar", "status": "short form of link_class (or closed / not closed for families)",
    "component": "connected component of the network of kept links", "latent_parent": "bookkeeping label of a maximal set of mutually linked families; not a reconstructed parent",
    "min_family_size": "tier: only families of at least this size are considered", "interval": "structural-haplotype interval",
    "tier": "high_confidence when all five criteria pass, otherwise candidate", "genotyped_from": "SNPs, or 'from deletions only'",
    "posterior": "posterior probability of the best genotype", "call": "genotype when the posterior reaches the call threshold, else unresolved",
    "best_genotype": "most probable genotype of the shared parent (0/0, 0/1, 1/1 = arrangement classes; labels carry no ancestral meaning)",
    "eta2_obs_median": "median over markers of the variance explained by family", "null_median_of_medians": "median over markers of the median eta squared after permuting family labels",
    "auc_within_vs_between": "rank AUC of within-family against between-family kinship", "recovery_pct": "percentage of within-family pairs above the 99.9th percentile of between-family pairs",
    "sample": "sample identifier", "reason": "why the sample was flagged", "rule": "rule and threshold that fired", "action": "flagged, kept / removed",
    "difference": "absolute difference of carrier frequency between the two homokaryotype classes", "fully_private": "carried by every called individual of one class and by none of the other",
    "marker_class": "marker class of the row (never mixed)", "n_markers": "number of markers behind the row"}


def statistics(run):
    rows = []

    def add(step, name, value, definition):
        rows.append({"step": FOLDER.get(step, step), "statistic": name, "value": value, "definition": definition})
    s = run.s
    q = s.get("qc", {})
    add("qc", "individuals_input", q.get("n_samples_input"), "individuals found in the inputs (or in --samples)")
    add("qc", "individuals_analysis_set", q.get("n_analysis_set"), "individuals every later step works on")
    add("qc", "individuals_removed", q.get("n_removed"), "removed by QC (confirmed in the config or --auto-exclude)")
    add("qc", "individuals_flagged_kept", q.get("n_flagged_kept"), "flagged by QC and kept")
    add("qc", "pairs_identical_genotypes", q.get("n_pairs_identical_genotypes"), "pairs of samples with identical genotypes")
    add("qc", "callable_span_bp", q.get("callable_span_bp") if q.get("callable_span_bp") is not None else "unmasked", "span of the callable mask, the denominator of per-base quantities")
    add("qc", "assembly_bp", q.get("assembly_bp"), "assembly length from --genome")
    if "median_depth_over_mask" in q:
        add("qc", "median_depth_over_mask", q["median_depth_over_mask"], "cohort median of the per-sample depth over the callable mask")
        add("qc", "median_depth_over_assembly", q["median_depth_over_assembly"], "cohort median of the per-sample depth over the assembly")
    r = s.get("relatedness", {})
    add("relatedness", "relatedness_source", r.get("source"), "where theta comes from")
    add("relatedness", "relatedness_markers", r.get("n_markers"), "markers in the GRM (empty when --relatedness is used as is)")
    add("relatedness", "pairs", r.get("n_pairs"), "pairs of individuals of the analysis set")
    add("relatedness", "parent_offspring_pairs", r.get("parent_offspring_pairs"), r.get("parent_offspring_rule"))
    add("relatedness", "unrelated_subset", r.get("n_unrelated_subset"), "greedy set without first-degree pairs")
    f = s.get("families", {})
    add("families", "families", f.get("n_families"), "relatedness-based groups: " + str(f.get("group_rule")) + "; edge: " + str(f.get("edge_rule")))
    add("families", "individuals_in_families", f.get("n_individuals_in_families"), "individuals assigned to a family")
    add("families", "individuals_unassigned", f.get("n_unassigned"), "individuals without a full-sib pair")
    add("families", "units", f.get("n_units"), "families plus unassigned individuals")
    add("families", "families_closed", f.get("n_closed"), "families in which every pair passes the full-sib thresholds")
    add("families", "families_not_closed", f.get("n_not_closed"), "families with at least one pair failing a threshold (connected blocks, not complete cliques)")
    add("families", "within_family_pairs_passing", f"{f.get('n_within_family_pairs_full_sib')} of {f.get('n_within_family_pairs')}", "closure over all families")
    n = s.get("nomination", {})
    add("nomination", "unit_pairs", n.get("n_unit_pairs"), "pairs of units tested")
    add("nomination", "nominated", n.get("n_nominated"), f"unit pairs with q < {n.get('q_max')}; {n.get('statistic')}; {n.get('n_permutations')} permutations")
    add("nomination", "nominated_family_pairs", n.get("n_nominated_family_pairs"), "nominated pairs with a family on both sides")
    a = s.get("audit", {})
    if a.get("ran"):
        add("audit", "auditable_pairs", a.get("n_auditable_pairs"), f"pairs of families with at least {a.get('min_family_size')} individuals each")
        add("audit", "nominated_auditable", a.get("n_nominated_auditable"), "nominated pairs among them")
        add("audit", "gate1_pass", a.get("n_gate1_pass"), "nominated auditable pairs with LLR_one > 0 and LLR_both < 0")
        add("audit", "rejected", a.get("n_rejected"), "nominated auditable pairs failing Gate 1")
        add("audit", "supported_links", a.get("n_supported"), "Gate-1-compatible pairs whose SV sharing exceeds the Gate-2 bar in at least one class")
        for k, v in (a.get("gate2") or {}).items():
            add("audit", f"gate2_{k}_median_sharing", v.get("threshold_median_of_gate1_pass"), f"median sharing of the Gate-1-compatible nominated pairs, {v.get('marker_set')}")
            add("audit", f"gate2_{k}_null_median", v.get("null_median"), "median sharing of the matched non-nominated pairs")
            add("audit", f"gate2_{k}_p", v.get("p_one_sided"), "one-sided P of the observed median against the matched null")
    w = s.get("network", {})
    if w.get("n_links_kept") is not None and "n_components" in w:
        add("network", "links_kept", w.get("n_links_kept"), "nominated links between families that were not rejected")
        add("network", "families_on_supported_links", w.get("n_families_on_supported_links"), "families joined by supported links")
        add("network", "components_supported", w.get("n_components_supported"), "connected components of the supported links")
        add("network", "latent_parents", w.get("n_latent_parents"), f"maximal sets of mutually linked families of at least {w.get('latent_parent_tier')} individuals")
        for t, v in (w.get("slot_violations") or {}).items():
            add("network", f"slot_limit_violations_min{t}", len(v), f"families of at least {t} individuals whose partners cannot be split into two groups of mutually linked families")
    h = s.get("haplotypes", {})
    if h:
        add("haplotypes", "intervals_supplied", h.get("n_intervals_supplied"), "intervals given with --inversions")
        add("haplotypes", "intervals_genotyped", h.get("n_intervals_genotyped"), "intervals genotyped by local PCA (supplied and scanned)")
        add("haplotypes", "intervals_high_confidence", h.get("n_high_confidence"), "intervals passing all five criteria")
        if h.get("scan_ran"):
            add("haplotypes", "scan_candidates", h.get("n_scan_candidates"), "candidate intervals found by the scan (detected automatically, not validated)")
            if "n_supplied_recovered" in h:
                add("haplotypes", "supplied_recovered_by_scan", h.get("n_supplied_recovered"), "supplied intervals overlapped by a scan candidate")
        if "n_tag_variants" in h:
            add("haplotypes", "arrangement_tag_variants", h.get("n_tag_variants"), "DEL / DUP with a carrier-frequency difference between homokaryotype classes at or above the threshold")
    b = s.get("barcode", {})
    if b.get("n_states"):
        add("barcode", "barcode_states", b.get("n_states"), "shared parent x interval states evaluated (links and junctions)")
        add("barcode", "barcode_called", b.get("n_called"), f"states at posterior >= {b.get('posterior_call')}")
        add("barcode", "barcode_called_strict", b.get("n_called_strict"), f"states at posterior >= {b.get('posterior_strict')}")
    for k, v in (s.get("family_effect", {}).get("classes") or {}).items():
        add("family_effect", f"eta2_median [{k}]", v["eta2_median"], f"median variance explained by family over {v['n_tested']} markers; null median {v['null_median']:.4g}")
    for k, v in (s.get("concordance", {}).get("classes") or {}).items():
        add("concordance", f"auc_within_vs_between [{k}]", v["auc"], f"{v['n_markers']} markers; recovery {v['recovery_pct']:.1f}% of within-family pairs above the 99.9th percentile of between-family pairs")
    T = pd.DataFrame(rows)
    T["value"] = ["" if v is None else v for v in T["value"]]
    return T


def analysis_line(run):
    return (f"**Analysis set: {len(run.analysis)} individuals** ({len(run.removed)} removed by QC, {run.n_flagged} flagged and kept). "
            f"{families.HEADER} {run.mask.label()}.")


def qc_report(run):
    s = run.s.get("qc", {}); t = run.t.get("qc", {})
    L = ["# BroodLink QC report", "", analysis_line(run), "", "QC flags samples and continues; it removes a sample only when the sample is confirmed in the config "
         "(`samples.exclude_confirmed`) or the run uses `--auto-exclude`.", "", "## Rules", ""]
    for k, v in (s.get("rules") or {}).items():
        L.append(f"- **{k}**: {v}")
    L += ["", "## Inputs", "", _md(t.get("per_input"), ["input", "marker_class", "n_samples", "n_markers", "filters", "mean_call_rate", "af_median"]), "",
          "## Mask", "", run.mask.label() + ".", "", _md(t.get("mask")), "", "## Flagged samples", ""]
    X = t.get("flags")
    L.append(_md(X) if X is not None and len(X) else "No sample was flagged.")
    L += ["", "## Pairs of samples with identical genotypes", ""]
    ID = t.get("identical")
    L.append(_md(ID) if ID is not None and len(ID) else "None found.")
    if "median_depth_over_mask" in s:
        L += ["", "## Depth", "", f"Cohort median depth over the mask: {s['median_depth_over_mask']:.2f}x; over the assembly: {s['median_depth_over_assembly']:.2f}x "
              "(two different denominators, reported in separate columns of `qc_per_sample.tsv`)."]
    E = t.get("effects")
    if E is not None and len(E):
        L += ["", "## Sex / batch / group effects", "", _md(E)]
    r = run.s.get("relatedness")
    if r:
        L += ["", "## Relatedness", "", f"Source: {r.get('source')}.", "", _md(run.t["relatedness"].get("classes")), "",
              f"Parent-offspring pairs ({r.get('parent_offspring_rule')}): {r.get('parent_offspring_pairs')}."]
    return "\n".join(L) + "\n"


def _md(df, cols=None, maxrows=60):
    if df is None or not len(df):
        return "_(empty)_"
    d = df[[c for c in (cols or df.columns) if c in df.columns]]
    out = ["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * d.shape[1]]
    for r in d.head(maxrows).itertuples(index=False):
        out.append("| " + " | ".join((f"{x:.4g}" if isinstance(x, float) else str(x)).replace("|", "/") for x in r) + " |")
    if len(d) > maxrows:
        out.append(f"\n_{len(d) - maxrows} more rows in the table file._")
    return "\n".join(out)


def capability(run):
    rows = []
    for st in STEPS:
        r = run.ran.get(st)
        rows.append({"step": STEP_TITLE[st], "ran": r[0] if r else "no", "inputs used": (r[1] if r else "") or "-",
                     "confidence": r[2] if r else "skipped", "why": r[3] if r else "step not run"})
    return pd.DataFrame(rows)


def run_report(run, fig_manifest):
    s = run.s
    L = [f"# BroodLink run report", "", f"BroodLink {__version__}, {datetime.date.today().isoformat()}, config hash {run.cfg_hash}, seed {run.seed}.", "", analysis_line(run), "",
         "Links between families are candidates. Nomination is screening; a link is 'supported' only after both audit gates, and the support bar of Gate 2 is relative.", "",
         "## Inputs detected and what each is used for", ""]
    for m in run.sets:
        L.append(f"- **{m.name}** (`--{m.flag}`, {m.source}): class {m.klass}, {m.n:,} markers, {len(m.samples)} individuals. Used for: {run.used.get(m.name, '')}.")
    L.append(f"- **Relatedness**: {run.rel_source}.")
    if run.hap_given is not None:
        L.append(f"- **--haplotypes**: {run.hap_given.n} intervals.")
    if run.intervals is not None:
        L.append(f"- **--inversions**: {len(run.intervals)} intervals.")
    L.append(f"- **Mask**: {run.mask.description}; {run.mask.label()}.")
    if run.notes:
        L += ["", "## Notes and warnings", ""] + [f"- {n}" for n in dict.fromkeys(run.notes)]
    L += ["", "## Headline numbers", "", _md(statistics(run)[["step", "statistic", "value"]], maxrows=200), "",
          "Definitions of every number: `statistics.tsv`."]
    L += ["", "## Run time", ""]
    for k, v in run.timing.items():
        L.append(f"- {FOLDER.get(k, k)}: {v:.1f} s")
    if fig_manifest is not None and len(fig_manifest):
        L += ["", "## Figures", "", _md(fig_manifest[["file", "title", "status"]])]
    L += ["", "## Capability table", "", "Confidence: *benchmarked* = the step ran as in the analysis the code was extracted from, on inputs of the kind it was benchmarked on; "
          "*provisional* = it ran in a mode that was not benchmarked (single input, too few markers, intervals from the scan, links not audited); "
          "*reimplemented from Methods* = no original code existed; *skipped* = not run.", "", _md(capability(run)), ""]
    return "\n".join(L) + "\n"


def methods_text(run):
    cfg, s = run.cfg, run.s
    P = []
    q = s.get("qc", {})
    P.append(f"Analyses were run with BroodLink {__version__} (seed {run.seed}). Input quality control flagged samples with a call rate below "
             f"{cfg['qc']['low_coverage_factor']} times the cohort median, pairs of samples with identical genotypes (theta >= {cfg['qc']['identical_theta']}) and outliers in "
             f"heterozygosity or median relatedness; {q.get('n_removed', 0)} samples were removed and {q.get('n_analysis_set')} individuals were analysed. "
             + (f"Markers were restricted to the callable mask ({run.mask.label().replace('callable span: ', '')}). " if run.mask.active else "No callable mask was applied. "))
    r = s.get("relatedness", {})
    if run.rel_mode == "given":
        P.append(f"Pairwise relatedness (theta, k0, k1, k2) was read from a precomputed table ({run.rel_info['format']}).")
    else:
        g = cfg["relatedness"]["grm"]
        P.append(f"Pairwise kinship was estimated as half of a VanRaden genomic relationship matrix with per-pair denominators, from {r.get('n_markers'):,} markers of {run.primary.name} "
                 f"(allele frequency {g['af_min']}-{g['af_max']}); no k2 can be estimated from this matrix.")
    f = s.get("families", {})
    if f:
        P.append(f"Full-sib families were defined as {f.get('group_rule')} ({f.get('edge_rule')}); they are relatedness-based groups, not pedigree records. "
                 f"This gave {f.get('n_families')} families holding {f.get('n_individuals_in_families')} individuals, {f.get('n_unassigned')} individuals unassigned; "
                 f"{f.get('n_closed')} families were closed (every pair passing) and {f.get('n_not_closed')} were not.")
    n = s.get("nomination", {})
    if n.get("n_unit_pairs"):
        P.append(f"Candidate half-sib links were nominated between units (families and unassigned individuals) from the {n.get('statistic')}, compared with a size-matched null of "
                 f"{n.get('n_permutations'):,} permutations of unit labels (P = (exceedances + 1) / (permutations + 1); Benjamini-Hochberg q < {n.get('q_max')}): "
                 f"{n.get('n_nominated')} of {n.get('n_unit_pairs'):,} unit pairs were nominated. Nomination is screening only."
                 + (f" Because a single input was available, its markers were split by {r.get('single_input_split')} into a nomination set and an audit set." if run.single_input else ""))
    a = s.get("audit", {})
    if a.get("ran"):
        au = cfg["audit"]
        P.append(f"Nominated pairs of families with at least {au['min_family_size']} individuals each were audited on markers not used for nomination ({'; '.join(a.get('audit_classes', []))}; "
                 f"allele frequency {au['af_min']}-{au['af_max']}). Gate 1 compared the per-marker likelihood of one shared parent and of two shared parents with that of unrelated parents "
                 f"(genotype error {au['genotype_error']}); a pair was compatible when the first log ratio was positive and the second negative. "
                 + ("Gate 2 compared the sharing of inherited structural variants with the " + str(au["k_match"]) + " non-nominated pairs nearest in the Gate-1 score, once per variant class; "
                    "a link was supported when its sharing exceeded the median of the Gate-1-compatible pairs in at least one class. " if a.get("gate2") else "No structural-variant class was available for Gate 2. ")
                 + f"Of {a.get('n_nominated_auditable')} nominated auditable pairs, {a.get('n_gate1_pass')} passed Gate 1, {a.get('n_rejected')} were rejected and {a.get('n_supported')} were supported."
                 + (" Links are labelled single-input, provisional." if run.single_input else ""))
    elif n.get("n_unit_pairs"):
        P.append("No audit was run; links are nominated, not audited.")
    w = s.get("network", {})
    if "n_components" in w:
        P.append(f"Families and kept links formed a network with {w.get('n_components')} components; supported links connected {w.get('n_families_on_supported_links')} families in "
                 f"{w.get('n_components_supported')} components. Latent parents (bookkeeping labels for maximal sets of mutually linked families) and violations of the two-parent slot limit "
                 f"were reported at minimum family sizes {', '.join(map(str, cfg['network']['slot_tiers']))}.")
    h = s.get("haplotypes", {})
    if h.get("n_intervals_genotyped"):
        t = cfg["inversions"]["tier"]
        P.append(f"Large structural haplotypes were genotyped by local PCA (standardised dosage, three clusters on PC1 by one-dimensional k-means; individuals with a class posterior below "
                 f"{cfg['inversions']['ambiguous_posterior']} set to missing) and tiered on five criteria (separation >= {t['separation_sd']} SD, heterozygote excess in the middle cluster, "
                 f"concordance of subset and full-cohort calls >= {t['concordance']}, at least {t['min_minor_copies']} minor-arrangement copies, family eta squared above its permutation null): "
                 f"{h.get('n_intervals_genotyped')} intervals genotyped, {h.get('n_high_confidence')} high confidence."
                 + (f" An automatic scan (local PCA in windows of {cfg['inversions']['snp_scan']['window_snps']} SNPs, reimplemented from the published description; and linkage disequilibrium among deletions in "
                    f"{cfg['inversions']['del_scan']['window_bp']:,}-bp windows) proposed {h.get('n_scan_candidates')} candidate intervals, which are not validated." if h.get("scan_ran") else ""))
    b = s.get("barcode", {})
    if b.get("n_states"):
        P.append(f"The genotype of each shared parent at each structural haplotype was inferred from the descendants' genotypes (Hardy-Weinberg prior, Mendelian transmission with genotype error "
                 f"{cfg['barcode']['genotype_error']}); {b.get('n_called')} of {b.get('n_states')} states reached a posterior of {b.get('posterior_call')}.")
    e = s.get("family_effect", {}).get("classes")
    if e:
        P.append("The variance explained by family (eta squared) was computed per marker class among individuals assigned to a family and compared with "
                 f"{cfg['family_effect']['n_permutations']:,} permutations of family labels: " + "; ".join(f"{k}: median {v['eta2_median']:.3f} (null {v['null_median']:.3f})" for k, v in e.items()) + ".")
    return "# Methods text generated by BroodLink\n\nGenerated from the configuration and the steps that ran. Edit before use.\n\n" + "\n\n".join(P) + "\n"


def write_tables(run):
    d = os.path.join(run.out, "tables"); os.makedirs(d, exist_ok=True)
    made, defs = [], []
    for step, key, name, sheet in TABLES:
        df = run.t.get(step, {}).get(key)
        if df is None:
            continue
        desc = run.files[step][key][1]
        write_tsv(df, os.path.join(d, name + ".tsv"), run.header(desc))
        made.append((name, sheet, df, desc))
        for c in df.columns:
            defs.append({"table": name, "column": c, "definition": COLUMN_DEFS.get(c, "see the description of the table")})
    try:
        import openpyxl  # noqa: F401
        with pd.ExcelWriter(os.path.join(d, "broodlink_tables.xlsx"), engine="openpyxl") as xw:
            pd.DataFrame({"table": [m[0] for m in made], "sheet": [m[1] for m in made], "rows": [len(m[2]) for m in made], "description": [m[3] for m in made]}
                         ).to_excel(xw, sheet_name="tables", index=False)
            pd.DataFrame(defs).to_excel(xw, sheet_name="definitions", index=False)
            pd.DataFrame({"header": run.header()}).to_excel(xw, sheet_name="run", index=False)
            for name, sheet, df, _d in made:
                df.head(1000000).to_excel(xw, sheet_name=sheet[:31], index=False)
    except ImportError:
        run.note("openpyxl is not installed: tables were written as TSV only (pip install openpyxl for the XLSX workbook)")
    with open(os.path.join(d, "README.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"BroodLink {__version__} - tables\n{families.HEADER}\nThe same tables as TSV and as sheets of broodlink_tables.xlsx (sheet 'definitions' explains the columns).\n\n")
        for name, _s, df, desc in made:
            fh.write(f"{name}.tsv\n  {desc}\n  rows: {len(df)}\n\n")
    return made


def write_all(run, figures=True):
    out = run.out
    with open(os.path.join(out, FOLDER["qc"], "qc_report.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(qc_report(run))
    write_tsv(statistics(run), os.path.join(out, "statistics.tsv"), run.header("Every headline number of the run with its definition."))
    write_tables(run)
    from .tracks import write_tracks
    write_tracks(run)
    man = None
    if figures:
        from .figures import make_figures
        man = make_figures(run)
    with open(os.path.join(out, "run_report.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(run_report(run, man))
    with open(os.path.join(out, "methods_text.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(methods_text(run))
    from .config import dump_config
    dump_config(run.cfg, os.path.join(out, "config_used.yaml"))
    log(f"done: {out} (run_report.md, statistics.tsv, methods_text.md, tables/, figures/, tracks/)")
