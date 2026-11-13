"""Configuration: every threshold and permutation count lives here (and in config.yaml).

Each default is traced to its source in docs/EXTRACTED_METHOD.md. Keys marked
"new in BroodLink" are not rules of the original analysis.
"""
from __future__ import annotations
import copy
import os

DEFAULTS = {
    "seed": 1,
    "samples": {"exclude_confirmed": []},          # sample ids confirmed for removal by the user
    "genotypes": {
        "gl_hard_call": 0.90,                      # hard call when the best genotype likelihood exceeds this
        "bisnp_step": 1,                           # keep every n-th SNP of --bisnp (original analyses used 3)
        "max_markers_snp": 2000000,                # safety cap after thinning
    },
    "mask": {"sv_min_callable_fraction": 0.5, "sv_class_loss_warning": 0.5},
    "qc": {
        "low_coverage_factor": 0.5,                # call rate below factor x cohort median
        "identical_theta": 0.354,                  # duplicate / MZ class
        "identical_concordance": 0.99,             # new in BroodLink
        "contamination_mad": 5.0,                  # new in BroodLink (the original fixed no number)
        "po_k1_min": 0.8, "po_k2_max": 0.05,       # parent-offspring region
    },
    "relatedness": {
        "theta_duplicate": 0.354, "theta_first": 0.177, "theta_second": 0.0884, "theta_third": 0.0442,
        "ngsrelate_binary": "ngsRelate",
        "method": "auto",                          # auto | ml | grm (auto: ml, the likelihood estimator; new in BroodLink)
        "ml": {
            "max_sites": 50000,                    # evenly spaced sites kept for the estimator (run time grows with it)
            "maf_min": 0.05,
            "min_informative_fraction": 0.5,       # share of individuals with an informative likelihood at the site
            "genotype_error": 0.01,                # for inputs with hard genotypes only
            "em_iterations": 200, "em_tolerance": 1e-5,
        },
        "grm": {
            "af_min": 0.05, "af_max": 0.95,
            "min_called_fraction_snp": 0.6944,     # 150 of 216 in the original
            "min_called_fraction_sv": 0.4629,      # 100 of 216 in the original
            "markers_benchmarked_high": 17276,     # benchmark: recovery 95.9% (deletions)
            "markers_benchmarked_low": 1178,       # benchmark: recovery 74.1% (duplications)
        },
    },
    "families": {"theta_min": 0.177, "k2_min": 0.15, "rule": "components"},   # or "clique"
    "nomination": {
        "theta_min": 0.0884, "theta_max": None,    # original: second degree or closer (no upper bound)
        "n_permutations": 20000, "q_max": 0.05,
        "single_input_split": "chromosome_parity",  # or "blocks"
        "block_size_bp": 5000000,
    },
    "audit": {
        "min_family_size": 4, "af_min": 0.02, "af_max": 0.98, "genotype_error": 0.03,
        "k_match": 40, "support_rule": "median_of_gate1_pass",
        "gate1_max_markers": 50000,                # new in BroodLink: seeded subsample above this
    },
    "network": {"slot_tiers": [4, 3, 2], "figure_min_family_size": 3},   # original post-audit network: families >= 3
    "inversions": {
        "scan": True,                              # default on; --no-inversion-scan turns it off
        "min_unrelated": 30,                       # new in BroodLink (original subset: 78)
        "min_snps": 50,                            # per interval / window after filters
        "min_dels": 8,                             # arrangement_classes(min_block=8)
        "del_af_min": 0.1, "del_af_max": 0.9, "del_link_cut": 0.3,
        "ambiguous_posterior": 0.95,
        "tier": {"separation_sd": 2.0, "het_p": 0.05, "concordance": 0.90, "min_minor_copies": 5,
                 "eta2_null_quantile": 0.95, "eta2_permutations": 200},
        "barcode_min_tier": "all",                 # supplied intervals: all | high_confidence
        "barcode_scan_criteria": ["separation", "heterozygote_excess"],
        "snp_scan": {"window_snps": 100, "n_pcs": 2, "mds_coords": 20, "z_outlier": 3.0,
                     "test_coords": 2, "merge_gap_windows": 3, "min_windows": 3,
                     "thin_to": 200000, "min_windows_per_chrom": 10, "min_size_bp": 0},
        "del_scan": {"window_bp": 1000000, "min_dels": 5, "hit_quantile": 0.90,
                     "core_quantile": 0.95, "merge_gap_bp": 2000000},
        "subclusters": {"min_class": 12, "min_component": 5, "k_max": 4, "restarts": 12,
                        "merge_sd": 2.0, "del_af_min": 0.05, "del_af_max": 0.95,
                        "diagnostic_delta": 1.0, "n_permutations": 200},
        "tags": {"min_called_per_class": 5, "min_delta": 0.6},
    },
    "barcode": {"posterior_call": 0.90, "posterior_strict": 0.95, "genotype_error": 0.03,
                # tag variants the shared parent must carry (compute_augmentation_relabelled_v1.py: called in >= 2 per family, k = 2, delta >= 0.8)
                # (allele frequency strictly inside af_min-af_max among the individuals in families: function prep of that script)
                "tag_evidence": {"min_called_per_family": 2, "min_tags": 2, "min_delta": 0.8, "af_min": 0.02, "af_max": 0.98},
                # sweep of the evidence rule and its permutation check (same script: k = 1, 2, 3; difference >= 0.6, 0.8, 0.9; 500 permutations)
                "augmentation": {"k_sweep": [1, 2, 3], "delta_sweep": [0.6, 0.8, 0.9], "n_permutations": 500},
                # recurrence of a junction parent across two links (refine_LP.py: rare deletions, pseudo-chain null)
                "recurrence": {"af_min": 0.01, "af_max": 0.10, "min_call_rate": 0.80, "n_null_max": 200, "p_max": 0.05}},
    # deletion genotypes along linked families (build_supp_chromosome_DEL_heatmaps_v2.py: fit_family, block2)
    "deletion_paths": {
        "low_call": 0.30,                          # individuals with more missing genotypes on the chromosome are not used for the fit
        "genotype_error": 0.02, "fit_p_min": 0.05, "fit_draws": 300, "fit_seed": 1, "fit_min_called": 3,
        "min_call": 0.90,                          # call rate among the well-called individuals of every family
        "thin_kb": 100, "thin_concordance": 0.90, "thin_min_called": 3,
        "max_columns": 160,                        # deletions drawn (even subsample of the eligible ones)
        "max_chain_families": 6,                   # new in BroodLink
        "chain_links": "supported",                # new in BroodLink: supported = chains around junction families (two or more supported links); kept = also families on two or more kept links
    },
    "family_effect": {"af_min": 0.05, "af_max": 0.95, "min_called": 20, "n_permutations": 1000,
                      "max_markers_per_class": 20000},
    "concordance": {
        "n_bootstrap": 1000,                       # resamples of the families (build_extdata_SV_validation_v2.py, boot_ci)
        "bootstrap_seed": 7,
        "sv_mask_flank_bp": 100000,                # SNPs "outside SV regions": structural-haplotype intervals + this flank, and every DEL / DUP span
    },
    # step 10, PCR markers of the structural haplotypes (docs/EXTRACTED_METHOD.md, section "PCR markers")
    "pcr": {
        "engine": "primer3",                       # primer3 = the original procedure (optional package primer3-py; falls back to builtin with a note) | builtin
        "markers_per_interval": 1,                 # the original designed one pair per interval
        "reference_flank_bp": 150,                 # new in BroodLink: flank cut on each side of a deletion when the design sequence comes from --reference
        "selection": {                             # the rules of the published screening table
            "min_called_per_class": 5,             # called individuals in each homokaryotype class (smallest value in the published table)
            "min_frequency_difference": 0.90,      # deletion allele frequency, class 2 against class 0
            "het_frequency_min": 0.30, "het_frequency_max": 0.70,   # heterokaryotypes near one half
            "screen_flank_bp": 100,                # generic pair of the screen: 2 x 100 bp of flank (bands 200 and 200 + length)
            "min_deletion_bp": 60, "max_deletion_bp": 1800,         # gel-resolvable: long band of the screen <= 2 kb; band difference >= 60 bp
        },
        "tm": {"mv_mM": 50.0, "dv_mM": 1.5, "dntp_mM": 0.6, "dna_nM": 50.0},   # conditions of primer3-py calc_tm (its defaults)
        "primer3": {"opt_tm": 56.0, "min_tm": 54.0, "max_tm": 57.0, "opt_size": 20, "min_size": 18, "max_size": 27,
                    "min_gc": 30.0, "max_gc": 80.0, "max_poly_x": 5, "product_min": 60, "product_max": 300,
                    "target_left": 1, "target_length": 2,   # target: the base before and the base after the breakpoint
                    "extra": {}},                  # any other Primer3 global argument, e.g. {PRIMER_GC_CLAMP: 1}
        # fallback routine (deltools.pick_flanking_primers of the deletion toolkit, its defaults)
        "builtin": {"tm_lo": 54.5, "tm_hi": 56.6, "tm_target": 55.5, "len_min": 18, "len_max": 26, "gc_lo": 30.0, "gc_hi": 60.0,
                    "gap_min": 40, "gap_max": 160, "max_homopolymer": 4},
        "in_silico": {"run": True, "max_mismatches": 2, "product_min": 60, "product_max": 6000, "locus_window_bp": 2500},
        "gels_per_page": 9,                        # simulated gels per page (the published sheets held 18 and 17; half of that)
        "ladder_bp": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000],
    },
    "figures": {
        "dpi": 450,                                # raster resolution of the PNG files
        "individual_noun": "individuals",          # word printed in figures for the individuals of the cohort (e.g. fish)
        "haplotype_label": "structural haplotypes",  # name of the large structural haplotypes in figure text (e.g. candidate polyINVs)
        "max_chain_figures": 4,                    # deletion-genotype figures: at most this many chains / link panels
        "gene_level_max_intervals": 4,             # gene-level barcode figure (needs --genes): number of intervals opened
        "gene_level_intervals": [],                # interval ids to open instead of the automatic choice (most called parent states)
    },
}


def _merge(a, b):
    for k, v in (b or {}).items():
        if isinstance(v, dict) and isinstance(a.get(k), dict):
            _merge(a[k], v)
        else:
            a[k] = v
    return a


def load_config(path=None):
    cfg = copy.deepcopy(DEFAULTS)
    if path:
        import yaml
        if not os.path.exists(path):
            raise SystemExit(f"ERROR: config file not found: {path}")
        with open(path, encoding="utf-8") as fh:
            _merge(cfg, yaml.safe_load(fh) or {})
    return cfg


def dump_config(cfg, path):
    import yaml
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
