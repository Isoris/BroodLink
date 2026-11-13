"""Tables (TSV + XLSX with a definitions sheet), statistics.tsv, qc_report.md, run_report.md with the
capability table, and methods_text.md (a Methods paragraph generated from the config and the steps
that actually ran)."""
from __future__ import annotations
import datetime
import os

import pandas as pd

from . import __version__, families
from .io import detail, log, write_tsv
from .pipeline import FOLDER, STEPS

STEP_TITLE = {"qc": "0. Quality control", "relatedness": "1. Pairwise relatedness", "families": "2. Full-sib families",
              "nomination": "3. Half-sib nomination", "audit": "4. Audit", "network": "5. Network, latent parents, slot limit",
              "haplotypes": "6a. Structural haplotypes (genotyping and scan)", "barcode": "6b. Barcode of the shared parent",
              "family_effect": "7. Family effect (eta squared)", "concordance": "8. Cross-marker concordance",
              "pcr": "9. PCR markers of the structural haplotypes (in-silico design)"}

# (step, key, table name, sheet name)
TABLES = [("qc", "flags", "T01_qc_exclusion_table", "T01 QC flags"), ("families", "families", "T02_families", "T02 families"),
          ("families", "membership", "T03_family_membership", "T03 membership"), ("network", "links", "T04_links", "T04 links"),
          ("audit", "gate2", "T05_audit_gate2", "T05 Gate 2"), ("network", "latent_parents", "T06_latent_parents", "T06 latent parents"),
          ("network", "slot", "T07_slot_limit", "T07 slot limit"), ("haplotypes", "summary", "T08_structural_haplotypes", "T08 haplotypes"),
          ("haplotypes", "tags", "T09_arrangement_tag_variants", "T09 tag variants"), ("barcode", "barcode", "T10_barcode", "T10 barcode"),
          ("family_effect", "summary", "T11_family_effect", "T11 family effect"), ("concordance", "summary", "T12_cross_marker_concordance", "T12 concordance"),
          ("pcr", "candidates", "T13_pcr_candidate_deletions", "T13 PCR candidates"), ("pcr", "primers", "T14_pcr_primers", "T14 PCR primers"),
          ("pcr", "per_interval", "T15_pcr_markers_per_interval", "T15 PCR per interval"), ("pcr", "order", "T16_pcr_primer_order_list", "T16 primer order list")]

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
    "marker_class": "marker class of the row (never mixed)", "n_markers": "number of markers behind the row",
    # PCR markers (step 10)
    "chrom": "sequence (chromosome) name", "start": "start position (for a deletion: POS of the VCF, the last base before it)", "end": "end position (for a deletion: END, its last base)",
    "interval_start": "start of the interval", "interval_end": "end of the interval",
    "marker": "label of the primer pair (M01, M02, ... in interval order)", "variant_id": "identifier of the deletion in the input",
    "length_bp": "length of the deletion (END - POS)", "arrangement_with_deletion": "homokaryotype class in which the deletion is the frequent allele (class 0 or class 2; labels carry no ancestral meaning)",
    "deletion_frequency_class_0": "allele frequency of the deletion among the called class-0 homokaryotypes", "deletion_frequency_class_2": "allele frequency of the deletion among the called class-2 homokaryotypes",
    "deletion_frequency_heterokaryotypes": "allele frequency of the deletion among the called heterokaryotypes (about 0.5 for a marker of the arrangement)",
    "frequency_difference": "absolute difference of the deletion allele frequency between the two homokaryotype classes",
    "carrier_frequency_class_0": "share of called class-0 homokaryotypes carrying the deletion", "carrier_frequency_class_2": "share of called class-2 homokaryotypes carrying the deletion",
    "n_called_class_0": "class-0 homokaryotypes with a genotype at the deletion", "n_called_class_2": "class-2 homokaryotypes with a genotype at the deletion", "n_called_heterokaryotypes": "heterokaryotypes with a genotype at the deletion",
    "screen_band_deleted_bp": "short band of the generic screening pair (two flanks)", "screen_band_undeleted_bp": "long band of the generic screening pair (two flanks + deletion length)",
    "diagnostic": "frequency difference at or above the threshold and heterokaryotype frequency inside its window", "gel_resolvable": "deletion length inside the gel window of the screen",
    "step06_tag": "the deletion is also in the tag-variant table of step 06 (carrier-frequency rule)", "sequence_source": "sequence the design was tried on", "rank_in_interval": "rank among the screened deletions of the interval, largest frequency difference first",
    "design_status": "result of the primer design where it was tried (ok or the reason)", "selected": "a primer pair was designed on this deletion",
    "forward_primer": "forward primer, 5'-3', + strand of the design sequence", "reverse_primer": "reverse primer, 5'-3', - strand of the design sequence",
    "length_forward": "length of the forward primer (nt)", "length_reverse": "length of the reverse primer (nt)", "tm_forward": "melting temperature of the forward primer (model and conditions in the table header)",
    "tm_reverse": "melting temperature of the reverse primer", "gc_forward_percent": "GC content of the forward primer", "gc_reverse_percent": "GC content of the reverse primer",
    "product_with_deletion_bp": "designed product on the allele with the deletion (short band)", "product_without_deletion_bp": "designed product on the allele without the deletion (short band + deletion length)",
    "product_difference_bp": "difference of the two products (the deletion length)", "design_source": "sequence the pair was designed on: split-read consensus of the SV VCF, or flanks cut from --reference",
    "engine": "design engine that produced the pair", "design_sequence": "the sequence across the deletion junction the pair was designed on", "breakpoint_index": "0-based index of the first base after the junction in the design sequence",
    "forward_start_index": "0-based index of the first base of the forward primer in the design sequence", "reverse_end_index": "0-based index of the last base of the reverse primer site in the design sequence",
    "check_primers_exact_in_design_sequence": "both primers are exact substrings of the design sequence", "check_primers_flank_breakpoint": "forward primer ends before the breakpoint and the reverse site starts after it",
    "check_products_recomputed": "product sizes recomputed from the primer positions and the deletion length equal the designed ones", "check_frequency_difference": "frequency difference at or above the threshold",
    "check_gel_resolvable_design": "designed pair: products differ by >= 60 bp, short product 60-1,200 bp, long product <= 6,000 bp (rule of the deletion toolkit)",
    "tm_difference_pair": "absolute Tm difference between the two primers", "min_gc_percent": "lower GC content of the two primers", "forward_3prime_gc": "the forward primer ends with G or C",
    "reverse_3prime_gc": "the reverse primer ends with G or C", "longest_homopolymer": "longest run of one base in either primer",
    "in_silico_products_at_locus": "products predicted in --reference that overlap the deletion locus (empty without --reference)", "in_silico_product_at_locus_bp": "size of the product at the locus (the reference carries the undeleted allele unless it has the deletion)",
    "in_silico_product_minus_design_bp": "product at the locus minus the designed product without the deletion", "in_silico_products_off_target": "products predicted elsewhere in --reference",
    "in_silico_off_target_products": "the first off-target products (sequence:start-end(size, mismatches))", "check_in_silico_product_within_5bp": "a product at the locus within 5 bp of the designed product without the deletion",
    "hairpin_forward_tm": "new in BroodLink: Tm of the most stable hairpin of the forward primer (primer3-py; empty when none or without primer3-py)", "hairpin_forward_dg_kcal": "new in BroodLink: its free energy (kcal/mol)",
    "hairpin_reverse_tm": "new in BroodLink: hairpin Tm of the reverse primer", "hairpin_reverse_dg_kcal": "new in BroodLink: its free energy (kcal/mol)",
    "self_dimer_forward_tm": "new in BroodLink: Tm of the forward primer self-dimer", "self_dimer_forward_dg_kcal": "new in BroodLink: its free energy (kcal/mol)",
    "self_dimer_reverse_tm": "new in BroodLink: Tm of the reverse primer self-dimer", "self_dimer_reverse_dg_kcal": "new in BroodLink: its free energy (kcal/mol)",
    "cross_dimer_tm": "new in BroodLink: Tm of the forward-reverse dimer", "cross_dimer_dg_kcal": "new in BroodLink: its free energy (kcal/mol)",
    "n_class_0": "individuals of homokaryotype class 0 at the interval", "n_class_2": "individuals of homokaryotype class 2", "n_heterokaryotypes": "heterokaryotypes at the interval",
    "n_deletions_screened": "deletions inside the interval genotyped in enough individuals of both classes", "n_diagnostic": "of these, diagnostic", "n_diagnostic_gel_resolvable": "diagnostic and gel-resolvable",
    "n_deletion_on_class_0": "diagnostic gel-resolvable deletions frequent on class 0", "n_deletion_on_class_2": "diagnostic gel-resolvable deletions frequent on class 2", "n_primer_pairs": "primer pairs designed for the interval",
    "reason_no_marker": "why the interval has no primer pair", "deletion_start": "POS of the selected deletion", "deletion_end": "END of the selected deletion",
    "name": "primer name for the order: marker_F or marker_R", "sequence": "primer sequence 5'-3'", "length": "primer length (nt)", "tm": "melting temperature (model and conditions in the table header)"}


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
    p = s.get("pcr", {})
    if p.get("ran"):
        add("pcr", "pcr_primer_engine", p.get("engine"), "engine that designed the pairs (primer3 = the original procedure; builtin = fallback routine); " + str(p.get("tm_conditions")))
        add("pcr", "pcr_deletions_screened", p.get("n_deletions_screened"), "deletions inside the intervals genotyped in enough individuals of both homokaryotype classes")
        add("pcr", "pcr_diagnostic_gel_resolvable", p.get("n_diagnostic_gel_resolvable"), "screened deletions that are diagnostic and gel-resolvable")
        add("pcr", "pcr_primer_pairs", p.get("n_primer_pairs"), "primer pairs designed (in-silico design, not tested in a laboratory); sources: " + ("; ".join(f"{k}: {v}" for k, v in (p.get("design_sources") or {}).items()) or "none"))
        add("pcr", "pcr_intervals_with_marker", f"{p.get('n_intervals_with_marker')} of {p.get('n_intervals')}", "intervals with at least one primer pair")
        if p.get("in_silico_ran"):
            add("pcr", "pcr_in_silico_product_at_locus", p.get("n_in_silico_product_at_locus"), "pairs with a product predicted at the deletion locus in --reference")
            add("pcr", "pcr_in_silico_off_target", p.get("n_in_silico_off_target"), "pairs with products predicted elsewhere in --reference")
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
    L += pcr_section(run)
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


def pcr_section(run):
    """Run-report section of step 10: engine, Tm model, the sequence each pair was designed on, intervals without a marker."""
    p = run.s.get("pcr", {})
    if not p.get("ran"):
        return []
    P = run.t.get("pcr", {}).get("primers")
    L = ["", "## PCR markers (in-silico design; no laboratory test)", "",
         f"- Engine: **{'Primer3 (' + 'primer3-py ' + str(p.get('primer3_version')) + '), the procedure of the original study' if p.get('engine') == 'primer3' else 'built-in routine (fallback; it does not return the primers Primer3 returns)'}**"
         + ("" if p.get("engine") == p.get("engine_requested") else f"; requested: {p.get('engine_requested')}, not available because primer3-py is not installed") + ".",
         f"- {p.get('tm_conditions')}.",
         f"- {p.get('n_deletions_screened')} deletions screened in {p.get('n_intervals')} intervals; {p.get('n_diagnostic')} diagnostic; {p.get('n_diagnostic_gel_resolvable')} diagnostic and gel-resolvable; "
         f"{p.get('n_primer_pairs')} primer pairs for {p.get('n_intervals_with_marker')} intervals" + (" (deletions taken from --pcr-picks)" if p.get("picks_supplied") else "") + ".",
         "- Design sequence per pair: " + ("; ".join(f"{k}: {v} pairs" for k, v in (p.get("design_sources") or {}).items()) or "none") + " (column `design_source` of `pcr_primers.tsv`)."]
    if p.get("in_silico_ran"):
        L.append(f"- In-silico PCR against --reference: {p.get('n_in_silico_product_at_locus')} pairs with a product at the locus ({p.get('n_in_silico_within_5bp')} within 5 bp of the designed long product), {p.get('n_in_silico_off_target')} with off-target products.")
    else:
        L.append("- In-silico PCR: not run (no --reference).")
    if P is not None and len(P):
        L.append(f"- Checks: primers exact substrings of the design sequence {p.get('n_checks_exact')} of {len(P)}; flanking the breakpoint {p.get('n_checks_flank')}; product sizes recomputed {p.get('n_checks_products')}.")
        L += ["", _md(P, ["marker", "interval", "variant_id", "arrangement_with_deletion", "forward_primer", "reverse_primer", "tm_forward", "tm_reverse", "product_with_deletion_bp", "product_without_deletion_bp", "design_source"], maxrows=100)]
    wo = p.get("intervals_without_marker") or {}
    if wo:
        L += ["", "Intervals without a marker:", ""] + [f"- {k}: {v}" for k, v in wo.items()]
    L += ["", "Simulated gels in `figures/` are drawn from the designed product sizes; they are not photographs."]
    return L


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
    c = s.get("concordance", {}).get("classes")
    if c:
        P.append("Kinship was computed from every marker class as a VanRaden genomic relationship matrix. For each class the rank AUC of within-family against between-family pairs of assigned individuals "
                 "and the share of within-family pairs above the 99.9th percentile of the between-family pairs (recovery) were computed, with 95% intervals from a bootstrap over families "
                 f"({cfg['concordance']['n_bootstrap']:,} resamples): " + "; ".join(f"{k}: AUC {v['auc']:.4f} [{v['auc_ci95'][0]:.4f}, {v['auc_ci95'][1]:.4f}], recovery {v['recovery_pct']:.1f}%" for k, v in c.items())
                 + ". The families were built from one of the inputs, so this is concordance between marker classes, not pedigree truth.")
    p = s.get("pcr", {})
    if p.get("ran"):
        c = cfg["pcr"]; sel = c["selection"]; p3 = c["primer3"]; bi = c["builtin"]
        eng = (f"with Primer3 (primer3-py {p.get('primer3_version')}; optimum Tm {p3['opt_tm']} C, range {p3['min_tm']}-{p3['max_tm']} C, primer length {p3['min_size']}-{p3['max_size']} nt with optimum {p3['opt_size']}, "
               f"GC at least {p3['min_gc']}%, product {p3['product_min']}-{p3['product_max']} bp, the breakpoint as target)" if p.get("engine") == "primer3" else
               f"with the built-in routine of BroodLink (the window of {bi['len_min']}-{bi['len_max']} nt nearest to a Tm of {bi['tm_target']} C within {bi['tm_lo']}-{bi['tm_hi']} C on each side, 3' end {bi['gap_min']}-{bi['gap_max']} nt from the breakpoint, "
               f"GC {bi['gc_lo']}-{bi['gc_hi']}%, no run above {bi['max_homopolymer']})")
        P.append(f"PCR markers of the structural haplotypes were designed in silico. Deletions inside each interval, genotyped in at least {sel['min_called_per_class']} individuals of each homokaryotype class, were screened: a deletion was "
                 f"diagnostic when its allele frequency differed by at least {sel['min_frequency_difference']} between the two homokaryotype classes with the heterokaryotypes at {sel['het_frequency_min']}-{sel['het_frequency_max']}, "
                 f"and gel-resolvable when {sel['min_deletion_bp']}-{sel['max_deletion_bp']} bp long. For each interval the diagnostic gel-resolvable deletions were tried in order of decreasing frequency difference until "
                 f"{c['markers_per_interval']} flanking primer pair(s) could be designed {eng} on the sequence across the deletion junction ({'; '.join(f'{k}: {v} pairs' for k, v in (p.get('design_sources') or {}).items()) or 'none'}). "
                 f"{p.get('tm_conditions')}. The allele with the deletion gives the short product and the allele without it the short product plus the deletion length. "
                 f"{p.get('n_deletions_screened')} deletions were screened, {p.get('n_diagnostic_gel_resolvable')} were diagnostic and gel-resolvable, and {p.get('n_primer_pairs')} primer pairs were designed for "
                 f"{p.get('n_intervals_with_marker')} of {p.get('n_intervals')} intervals."
                 + (f" Each pair was searched in the reference assembly (at most {c['in_silico']['max_mismatches']} mismatches per primer, matching 3' base, products of {c['in_silico']['product_min']}-{c['in_silico']['product_max']} bp): "
                    f"{p.get('n_in_silico_product_at_locus')} pairs gave a product at the deletion locus and {p.get('n_in_silico_off_target')} gave off-target products." if p.get("in_silico_ran") else "")
                 + " No primer was tested in the laboratory; the gels shown are simulated from the designed product sizes.")
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
    log("outputs: writing reports, tables, tracks" + (" and figures" if figures else ""))
    with open(os.path.join(out, FOLDER["qc"], "qc_report.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(qc_report(run))
    detail("00_qc/qc_report.md")
    write_tsv(statistics(run), os.path.join(out, "statistics.tsv"), run.header("Every headline number of the run with its definition."))
    detail("statistics.tsv")
    write_tables(run)
    detail("tables/ (TSV + XLSX with a definitions sheet)")
    from .tracks import write_tracks
    write_tracks(run)
    detail("tracks/ (GFF3 + BED)")
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
    log(f"start with {os.path.join(out, 'run_report.md')}")
