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
    "barcode": {"posterior_call": 0.90, "posterior_strict": 0.95, "genotype_error": 0.03},
    "family_effect": {"af_min": 0.05, "af_max": 0.95, "min_called": 20, "n_permutations": 1000,
                      "max_markers_per_class": 20000},
    "figures": {"dpi": 450},                    # raster resolution of the PNG files (author rule: at least 450)
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
