"""The run: load inputs, decide what each marker class is used for, run the steps in order.

Every step writes its tables into a numbered folder and a summary.json next to them, so a
single step can be re-run later (`broodlink <step>`): the steps before it are then read back
from their tables instead of being recomputed.
"""
from __future__ import annotations
import hashlib
import json
import os
import re
import time

import numpy as np
import pandas as pd

from . import __version__, audit, barcode, concordance, families, family_effect, inversions, network, nomination, opposite, qc, relatedness
from .io import (InputError, MarkerSet, detail, expand, log, read_beagle, read_depth, read_fai, read_haplotypes,
                 read_intervals_bed, read_relatedness, read_samples, read_variants, sniff, write_tsv)
from .mask import apply_mask, build_mask

STEPS = ["qc", "relatedness", "families", "nomination", "audit", "network", "haplotypes", "barcode",
         "family_effect", "concordance"]
FOLDER = {"qc": "00_qc", "relatedness": "01_relatedness", "families": "02_families", "nomination": "03_nomination",
          "audit": "04_audit", "network": "05_network", "haplotypes": "06_structural_haplotypes", "barcode": "07_barcode",
          "family_effect": "08_family_effect", "concordance": "09_cross_marker_concordance"}
ID_COLS = ("individual", "sample", "id1", "id2", "sample_1", "sample_2", "unit_a", "unit_b", "family", "unit", "interval",
           "parent", "families", "supplied_interval", "scan_interval", "variant_id", "id", "chrom", "latent_parent")
POINT = ("SNP", "INDEL", "SNP_INDEL")
WHAT = {"qc": "quality control of every input (flag, do not remove)",
        "relatedness": "pairwise relatedness (theta, k0, k1, k2) and relationship classes",
        "families": "full-sib families from theta and k2 (relatedness-based groups)",
        "nomination": "candidate half-sib links between units (size-matched permutation test, BH q)",
        "audit": "audit of nominated family pairs on markers not used for nomination (Gate 1, Gate 2)",
        "network": "network of families, components, latent parents, slot limit",
        "haplotypes": "large structural haplotypes: local PCA genotyping and scan",
        "barcode": "posterior genotype of each shared parent at the structural haplotypes",
        "family_effect": "variance explained by family (eta squared) per marker class",
        "concordance": "kinship of every marker class within and between families"}


def _result_line(step, s):
    """One line of headline numbers printed when a step finishes."""
    g = s.get
    try:
        if step == "qc":
            return f"{g('n_analysis_set')} of {g('n_samples_input')} individuals in the analysis set; {g('n_flagged_kept')} flagged and kept, {g('n_removed')} removed; {g('n_pairs_identical_genotypes')} pairs with identical genotypes"
        if step == "relatedness":
            m = f"; {g('n_markers'):,} markers" if g("n_markers") else ""
            return f"source: {g('mode')}{m}; {g('n_unrelated_subset')} individuals in the unrelated subset"
        if step == "families":
            return f"{g('n_families')} families holding {g('n_individuals_in_families')} individuals; {g('n_unassigned')} unassigned; {g('n_units')} units; {g('n_not_closed')} families not fully closed"
        if step == "nomination":
            return f"{g('n_unit_pairs', 0):,} unit pairs tested, {g('n_nominated', 0)} nominated ({g('n_nominated_family_pairs', 0)} between families)"
        if step == "audit":
            if not g("ran"):
                return "not run: links stay 'nominated, not audited'"
            return f"{g('n_auditable_pairs')} auditable pairs, {g('n_nominated_auditable')} nominated, {g('n_gate1_pass')} pass Gate 1, {g('n_rejected')} rejected; {g('n_supported')} supported, {g('n_gate1_only')} Gate-1 only"
        if step == "network":
            ks = [k for k in s if "component" in k or "latent" in k or "kept" in k]
            return "; ".join(f"{k.replace('_', ' ')}: {s[k]}" for k in ks if not isinstance(s[k], (dict, list)))[:300] or "built"
        if step == "haplotypes":
            r = f"{g('n_intervals_genotyped', 0)} intervals genotyped ({g('n_high_confidence', 0)} high confidence, {g('n_candidate', 0)} candidate); {g('n_barcode_intervals', 0)} used for the barcode"
            if g("scan_ran"):
                r += f"; scan: {g('n_scan_candidates', 0)} candidate intervals"
            return r
        if step == "barcode":
            return f"{g('n_states', 0)} parent x interval states; {g('n_called', 0)} called at posterior >= {g('posterior_call')}, {g('n_called_strict', 0)} at >= {g('posterior_strict')}" if g("n_states") else "not run"
        if step == "family_effect":
            return "; ".join(f"{k}: median eta2 {v['eta2_median']:.3f} (null {v['null_median']:.3f}, {v['n_tested']:,} markers)" for k, v in (g("classes") or {}).items()) or "not run"
        if step == "concordance":
            return "; ".join(f"{k}: AUC {v['auc']:.3f}" for k, v in (g("classes") or {}).items()) or "not run"
    except Exception:                                           # a progress line never stops a run
        return ""
    return ""


def _json(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


class Run:
    def __init__(self, args, cfg):
        self.a, self.cfg = args, cfg
        self.out = args.out
        self.seed = int(cfg["seed"])
        self.t = {}            # step -> {key: DataFrame}
        self.s = {}            # step -> summary dict
        self.files = {}        # step -> {key: (file name, description)}
        self.ran = {}          # step -> (status, inputs used, confidence, why)
        self.notes = []        # messages repeated in the run report
        self.timing = {}
        self.cfg_hash = hashlib.sha256(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()[:10]
        self.analysis = []
        self.removed = []
        self.n_flagged = 0

    # ------------------------------------------------------------------ helpers
    def note(self, msg):
        log(msg); self.notes.append(msg)

    def header(self, description=""):
        classes = "; ".join(f"{m.name}: {m.n:,} markers" for m in self.sets) or "none"
        h = [f"BroodLink {__version__} | config {self.cfg_hash} | seed {self.seed}",
             f"analysis set: {len(self.analysis)} individuals ({len(self.removed)} removed by QC; {self.n_flagged} flagged and kept)",
             f"marker classes: {classes}", f"relatedness source: {getattr(self, 'rel_source', 'not yet determined')}",
             families.HEADER, self.mask.label()]
        return h + ([description] if description else [])

    def put(self, step, key, df, fname, description):
        self.t.setdefault(step, {})[key] = df
        self.files.setdefault(step, {})[key] = (fname, description)
        write_tsv(df, os.path.join(self.out, FOLDER[step], fname), self.header(description))

    def finish(self, step, summary, t0):
        self.s[step] = summary
        self.timing[step] = time.time() - t0
        st = self.ran.get(step) or ("", "", "", "")
        res = _result_line(step, summary)
        if st[0] == "no":
            log(f"  -> skipped ({self.timing[step]:.1f} s): {st[3]}")
        else:
            log(f"  -> done in {self.timing[step]:.1f} s" + (f" [{st[2]}]" if st[2] else "") + (f": {res}" if res else ""))
        d = os.path.join(self.out, FOLDER[step]); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "summary.json"), "w", encoding="utf-8", newline="\n") as fh:
            json.dump({"summary": summary, "files": self.files.get(step, {}), "status": self.ran.get(step)}, fh, indent=1, default=_json, sort_keys=True)
        with open(os.path.join(d, "README.txt"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(f"BroodLink {__version__} - {FOLDER[step]}\n{families.HEADER}\nEvery table starts with '#' header lines (version, config hash, analysis set, marker classes, callable span).\n\n")
            for key, (fname, desc) in self.files.get(step, {}).items():
                df = self.t[step][key]
                fh.write(f"{fname}\n  {desc}\n  rows: {len(df)}; columns: {', '.join(map(str, df.columns))}\n\n")
            fh.write("summary.json\n  headline numbers of this step, the list of files, and the status of the step.\n")

    def load(self, step):
        """Read a step back from its tables (used when a later step is re-run alone)."""
        d = os.path.join(self.out, FOLDER[step]); p = os.path.join(d, "summary.json")
        if not os.path.exists(p):
            raise InputError(f"step '{step}' has not been run in {self.out}: run `broodlink run` first (missing {FOLDER[step]}/summary.json)")
        with open(p, encoding="utf-8") as fh:
            j = json.load(fh)
        self.s[step] = j["summary"]; self.files[step] = {k: tuple(v) for k, v in j["files"].items()}
        self.ran[step] = tuple(j["status"]) if j.get("status") else None
        self.t[step] = {}
        for key, (fname, _d) in self.files[step].items():
            df = pd.read_csv(os.path.join(d, fname), sep="\t", comment="#", dtype={c: str for c in ID_COLS}, keep_default_na=True)
            for c in df.columns:
                if c in ID_COLS or not (pd.api.types.is_numeric_dtype(df[c]) or pd.api.types.is_bool_dtype(df[c])):
                    df[c] = df[c].fillna("").astype(str)
            self.t[step][key] = df
        self.timing[step] = 0.0

    def status(self, step, ran, inputs, confidence, why):
        self.ran[step] = (ran, inputs, confidence, why)

    # ------------------------------------------------------------------ inputs
    def load_inputs(self):
        a, cfg = self.a, self.cfg
        t0 = time.time()
        if not any([a.bisnp, a.gt, a.snp_indel, a.vcf, a.relatedness]):
            raise InputError("at least one of --bisnp, --gt, --snp_indel, --vcf, --relatedness is required")
        log("reading inputs")
        self.genome = read_fai(a.genome) if a.genome else {}
        if a.genome:
            detail(f"--genome: {len(self.genome)} sequences, {sum(self.genome.values()):,} bp")
        self.samples_meta = read_samples(a.samples) if a.samples else None
        if self.samples_meta is not None:
            detail(f"--samples: {len(self.samples_meta)} individuals")
        self.depth = read_depth(a.depth) if a.depth else None
        if self.depth is not None:
            detail(f"--depth: {len(self.depth)} samples")
        sets, bg = [], []
        if a.bisnp:
            files = expand(a.bisnp)
            bg = [f for f in files if sniff(f) == "beagle"]
            if bg:
                detail(f"--bisnp: reading {len(bg)} Beagle file(s) of genotype likelihoods (every {cfg['genotypes']['bisnp_step']} SNP kept)")
                ms = read_beagle(bg, None, cfg["genotypes"]["bisnp_step"], cfg["genotypes"]["gl_hard_call"], "bisnp", cfg["genotypes"]["max_markers_snp"])[0]
                if all(re.fullmatch(r"Ind\d+", s) for s in ms.samples):
                    if self.samples_meta is not None and len(self.samples_meta) == len(ms.samples):
                        ms.samples = self.samples_meta.id.tolist()
                        self.note("--bisnp: Beagle columns are named Ind0..IndN; they were mapped to the rows of --samples in order")
                    else:
                        self.note("--bisnp: Beagle columns are named Ind0..IndN and no --samples file of the same length was given; the names are kept")
                sets.append(ms)
            for f in files:
                if f not in bg:
                    sets += [m for m in read_variants(f, "bisnp", cfg)]
        for flag, paths in (("gt", a.gt), ("snp_indel", [a.snp_indel] if a.snp_indel else []), ("vcf", a.vcf)):
            for p in expand(paths or []):
                detail(f"--{flag}: reading {os.path.basename(p)}")
                new = read_variants(p, flag, cfg)
                for m in new:
                    detail(f"  {m.name}: {m.n:,} markers x {len(m.samples)} samples")
                sets += new
        self.R_given, self.rel_info = (None, None)
        if a.relatedness:
            names = self.samples_meta.id.tolist() if self.samples_meta is not None else None
            self.R_given, self.rel_info = read_relatedness(a.relatedness, names)
            detail(f"--relatedness: {self.rel_info['format']}, {self.rel_info['n_samples']} samples, {self.rel_info['n_pairs']:,} pairs")
        elif a.bisnp and bg:
            # the original analysis used the external tool NGSrelate; this branch needs the binary
            # and was never exercised on the development machine (tool not installed)
            Rn, msg = relatedness.try_ngsrelate(bg, len(sets[0].samples), sets[0].samples, cfg)
            self.note("--bisnp without --relatedness: " + msg)
            if Rn is not None:
                ids_ = sorted(set(Rn.id1) | set(Rn.id2))
                self.R_given = Rn
                self.rel_info = {"format": "NGSrelate run by BroodLink", "n_files": 1, "n_samples": len(ids_), "n_pairs": len(Rn),
                                 "n_pairs_expected": len(ids_) * (len(ids_) - 1) // 2, "complete": True}
        self.hap_given = None
        known = set().union(*[set(m.samples) for m in sets]) if sets else set()
        if self.R_given is not None:
            known |= set(self.R_given.id1) | set(self.R_given.id2)
        if a.haplotypes:
            self.hap_given = read_haplotypes(a.haplotypes, known)
        self.intervals = read_intervals_bed(a.inversions) if a.inversions else None
        # ---- sample matching: the only place where QC stops a run
        if self.samples_meta is not None:
            universe = self.samples_meta.id.tolist()
        else:
            universe = sorted(known | (set(self.hap_given.samples) if self.hap_given is not None and not known else set()))
        uni = set(universe)
        named = [(m.name, m.samples) for m in sets]
        if self.R_given is not None:
            named.append(("--relatedness", sorted(set(self.R_given.id1) | set(self.R_given.id2))))
        if self.hap_given is not None:
            named.append(("--haplotypes", self.hap_given.samples))
        ref_name, ref = named[0]
        for name, smp in named:
            ov = set(smp) & (uni if self.samples_meta is not None else set(ref))
            if not ov:
                other = "--samples" if self.samples_meta is not None else ref_name
                raise InputError(f"sample names of {name} cannot be matched to {other}: none of its {len(smp)} samples is found there "
                                 f"(unmatched, first 10: {', '.join(list(smp)[:10])})")
            extra = sorted(set(smp) - uni)
            if extra and self.samples_meta is not None:
                self.note(f"{name}: {len(extra)} samples are not in --samples and are ignored")
        self.all_samples = universe
        # ---- mask
        self.mask = build_mask(a.mask, a.mask_exclude, self.genome)
        rows, masked = [], []
        for m in sets:
            m2, row = apply_mask(m.subset_samples([s for s in universe if s in set(m.samples)]), self.mask, cfg)
            rows.append(row); masked.append(m2)
        self.sets = masked
        self.mask_report = pd.DataFrame(rows, columns=["input", "marker_class", "markers_total", "markers_inside_mask", "markers_outside_mask", "rule", "warning"])
        for w in self.mask_report.warning:
            if w:
                self.notes.append("WARNING: " + w)
        if self.hap_given is not None:
            self.hap_given = self.hap_given.subset_samples([s for s in universe if s in set(self.hap_given.samples)])
        self.analysis = list(universe)
        self.plan()
        self.timing["read inputs"] = time.time() - t0

    def plan(self):
        """Decide the relatedness source and what each marker class is used for; print it."""
        a, cfg = self.a, self.cfg
        geno = self.sets
        self.primary, self.single_input, self.audit_note = None, False, ""
        bis = [m for m in geno if m.flag == "bisnp"]
        if self.R_given is not None:
            self.rel_source = f"--relatedness, used as is ({self.rel_info['format']}, {self.rel_info['n_samples']} samples, {self.rel_info['n_pairs']} pairs)"
            self.rel_mode = "given"
            self.audit_plan = [m for m in geno if m.flag != "bisnp"]
            if bis and not self.audit_plan:
                self.audit_note = ("the only genotype input is --bisnp, which is assumed to be the source of the supplied relatedness: "
                                   "auditing on it would be circular, so no audit is run")
            elif bis:
                self.audit_note = "--bisnp is assumed to be the source of the supplied relatedness and is kept out of the audit"
            elif not geno:
                self.audit_note = "no genotype input: links are nominated, not audited"
        else:
            self.rel_mode = "grm"
            self.primary = bis[0] if bis else max(geno, key=lambda m: m.n)
            self.rel_source = f"genomic relationship matrix (VanRaden, kinship = GRM / 2) from {self.primary.name}; no k2 can be estimated"
            if bis:
                self.note("--bisnp given without --relatedness: BroodLink has no likelihood-based estimator of its own; "
                          "relatedness comes from a GRM of hard-called genotypes (likelihood > %s)" % cfg["genotypes"]["gl_hard_call"])
            if len(geno) == 1:
                self.single_input = True
                self.audit_plan = []
                self.audit_note = "single input: its markers are split into a nomination set and an audit set; links are labelled 'single-input, provisional'"
            else:
                self.audit_plan = [m for m in geno if m is not self.primary]
        used = {}
        for m in geno:
            u = []
            if m is self.primary:
                u.append("relatedness, families and nomination" + (" (nomination half) and audit (audit half)" if self.single_input else ""))
            if m in self.audit_plan:
                u.append("audit Gate 1" + (" and Gate 2" if m.klass in audit.SV_LIKE else ""))
            if m.klass in POINT or m.klass == "DEL":
                u.append("structural-haplotype genotyping / scan")
            u += ["family effect", "cross-marker concordance", "QC"]
            used[m.name] = ", ".join(u)
        self.used = used
        log(f"BroodLink {__version__}: inputs detected")
        for m in geno:
            log(f"  {m.name}: class {m.klass}, {m.n:,} markers, {len(m.samples)} samples -> {used[m.name]}")
        log(f"  relatedness: {self.rel_source}")
        if self.hap_given is not None:
            log(f"  --haplotypes: {self.hap_given.n} intervals, {len(self.hap_given.samples)} samples -> barcode, family effect")
        if self.intervals is not None:
            log(f"  --inversions: {len(self.intervals)} intervals -> genotyped by local PCA")
        log(f"  mask: {self.mask.description}; {self.mask.label()}")
        if self.audit_note:
            self.note(self.audit_note)

    # ------------------------------------------------------------------ relatedness helpers
    def _grm_pairs(self, ms):
        R, used = relatedness.relatedness_from_markers(ms, self.cfg, f"GRM:{ms.klass}")
        return R, used

    def _marker_confidence(self, used):
        g = self.cfg["relatedness"]["grm"]
        if used >= g["markers_benchmarked_high"]:
            return f"{used:,} markers: at or above the larger benchmarked marker set ({g['markers_benchmarked_high']:,} markers; 95.9% of within-family pairs recovered in the benchmark)"
        if used >= g["markers_benchmarked_low"]:
            return (f"{used:,} markers: between the two benchmarked marker sets ({g['markers_benchmarked_low']:,} markers: 74.1% of within-family pairs recovered; "
                    f"{g['markers_benchmarked_high']:,}: 95.9%). The benchmark fixed no minimum marker count; these are observations")
        return (f"{used:,} markers: BELOW the smallest benchmarked marker set ({g['markers_benchmarked_low']:,} markers, 74.1% recovery); families are provisional")

    # ------------------------------------------------------------------ step 0
    def step_qc(self):
        t0 = time.time(); cfg = self.cfg
        if self.rel_mode == "given":
            R0 = self.R_given[self.R_given.id1.isin(self.all_samples) & self.R_given.id2.isin(self.all_samples)]
        else:
            R0, _ = self._grm_pairs(self.primary)
        PS = qc.per_sample(self.sets, self.mask, self.depth)
        PI, PC = qc.per_input(self.sets, self.genome)
        FL, ID, EF = qc.sample_checks(self.sets, R0, PS, self.samples_meta, self.all_samples, cfg)
        X, removed = qc.exclusions(FL, cfg, self.a.auto_exclude)
        self.removed = removed
        self.n_flagged = int(X[~X.action.str.startswith("removed")]["sample"].nunique())
        self.analysis = [s for s in self.all_samples if s not in set(removed)]
        self.sets = [m.subset_samples([s for s in self.analysis if s in set(m.samples)]) for m in self.sets]
        self._rebind_plan()
        self.put("qc", "per_sample", PS, "qc_per_sample.tsv", "Per sample and input: markers genotyped, call rate, heterozygosity (share of called genotypes that are heterozygous), het / hom-alt ratio, singletons, SV calls carried and their span (bp); depth over the mask and over the assembly in separate columns when --depth is given.")
        self.put("qc", "per_input", PI, "qc_per_input.tsv", "Per input: markers, filters applied (before -> after), allele-frequency summary, mean call rate.")
        self.put("qc", "per_chrom", PC, "qc_markers_per_chromosome.tsv", "Markers per chromosome against chromosome length (needs --genome).")
        self.put("qc", "mask", self.mask_report, "qc_mask_report.tsv", "Markers inside and outside the callable mask per input, with the rule applied.")
        self.put("qc", "flags", X, "qc_exclusion_table.tsv", "Every flagged sample: reason, rule, value and the action taken. Default: flag, do not remove.")
        self.put("qc", "identical", ID, "qc_pairs_identical_genotypes.tsv", "Pairs of samples with identical genotypes (not dropped automatically).")
        self.put("qc", "effects", EF, "qc_group_effects.tsv", "Sex / batch / group effects on call rate and heterozygosity (Kruskal-Wallis).")
        self.put("qc", "analysis_set", pd.DataFrame({"individual": self.analysis}), "analysis_set.tsv", "The individuals every later step works on.")
        dm = {}
        if self.depth is not None and len(self.depth):
            d = self.depth[self.depth.id.isin(self.analysis)]
            dm = {"median_depth_over_mask": float(d.depth_mask.median()), "median_depth_over_assembly": float(d.depth_assembly.median())}
        summ = {"n_samples_input": len(self.all_samples), "n_analysis_set": len(self.analysis), "n_removed": len(removed), "removed": removed,
                "n_flagged_kept": self.n_flagged, "flags_by_reason": X.groupby("reason")["sample"].nunique().to_dict(),
                "n_pairs_identical_genotypes": int(len(ID)), "callable_span_bp": self.mask.span(), "assembly_bp": self.mask.assembly(), **dm,
                "rules": {"low coverage": f"call rate below {cfg['qc']['low_coverage_factor']} x cohort median, per input",
                          "pairs of samples with identical genotypes": f"theta >= {cfg['qc']['identical_theta']} or genotype concordance >= {cfg['qc']['identical_concordance']} (hard genotypes)",
                          "possible contamination": f"median relatedness or heterozygosity above cohort median + {cfg['qc']['contamination_mad']} MAD (new in BroodLink)",
                          "removal": "only samples confirmed in the config (samples.exclude_confirmed) or, with --auto-exclude, low-coverage samples and one sample of each identical pair"}}
        self.status("qc", "yes", "all inputs", "reimplemented from Methods", "always runs")
        self.finish("qc", summ, t0)

    def _rebind_plan(self):
        """After samples were subset the MarkerSet objects are new: re-point primary / audit plan by name."""
        by = {m.name: m for m in self.sets}
        if self.primary is not None:
            self.primary = by[self.primary.name]
        self.audit_plan = [by[m.name] for m in self.audit_plan]

    def restore_analysis(self):
        A = self.t["qc"]["analysis_set"].individual.astype(str).tolist()
        self.removed = [s for s in self.all_samples if s not in set(A)]
        self.n_flagged = int(self.s["qc"].get("n_flagged_kept", 0))
        self.analysis = A
        self.sets = [m.subset_samples([s for s in A if s in set(m.samples)]) for m in self.sets]
        self._rebind_plan()

    # ------------------------------------------------------------------ step 1
    def step_relatedness(self):
        t0 = time.time(); cfg = self.cfg; A = set(self.analysis)
        summ = {"source": self.rel_source, "mode": self.rel_mode}
        Rn = None
        if self.rel_mode == "given":
            R = self.R_given[self.R_given.id1.isin(A) & self.R_given.id2.isin(A)].reset_index(drop=True)
            miss = sorted(A - (set(R.id1) | set(R.id2)))
            if miss:
                self.note(f"--relatedness: {len(miss)} individuals of the analysis set have no pair in the relatedness table; they stay unassigned")
            conf, cls = "benchmarked", "relatedness table"
            summ["n_markers"] = None
        else:
            R, used = self._grm_pairs(self.primary)
            summ["n_markers"] = used; summ["marker_class"] = self.primary.name
            summ["marker_count_check"] = self._marker_confidence(used)
            self.note("relatedness: " + summ["marker_count_check"])
            conf = "benchmarked" if used >= cfg["relatedness"]["grm"]["markers_benchmarked_low"] else "provisional"
            cls = self.primary.name
            if self.single_input:
                m0, m1, how = nomination.split_markers(self.primary, cfg)
                self.split = (m0, m1, how)
                Rn, usedn = self._grm_pairs(self.primary.subset_markers(m0))
                summ["single_input_split"] = how; summ["n_markers_nomination_half"] = usedn
                self.note(f"single input: nomination on {usedn:,} markers, audit on the other half ({how})")
        R = R.copy(); R["class"] = relatedness.classify_pairs(R, cfg)
        counts, rq = qc.relatedness_qc(R, R["class"], cfg)
        if rq["parent_offspring_pairs"]:
            self.note(f"relatedness QC: {rq['parent_offspring_pairs']} pairs fall in the parent-offspring region ({rq['parent_offspring_rule']}); families assume one generation")
        unrel = relatedness.unrelated_subset(R, sorted(A), cfg)
        self.put("relatedness", "pairs", R, "relatedness_pairs.tsv", "One row per pair: theta (kinship), k0, k1, k2 (empty when relatedness comes from a GRM), source and relationship class.")
        if Rn is not None:
            self.put("relatedness", "pairs_nomination", Rn, "relatedness_pairs_nomination_half.tsv", "Single-input mode: kinship from the nomination half of the markers only.")
        self.put("relatedness", "classes", counts, "relatedness_class_counts.tsv", "Number of pairs per relationship class (thresholds in the config).")
        self.put("relatedness", "unrelated", pd.DataFrame({"individual": unrel}), "unrelated_subset.tsv", "Greedy set of individuals without a first-degree pair among them (the original analysis used an exact maximum set).")
        summ.update(rq); summ["n_unrelated_subset"] = len(unrel)
        summ["thresholds"] = {k: cfg["relatedness"][k] for k in ("theta_duplicate", "theta_first", "theta_second", "theta_third")}
        self.status("relatedness", "yes", cls, conf, self.rel_source)
        self.finish("relatedness", summ, t0)

    # ------------------------------------------------------------------ step 2
    def step_families(self):
        t0 = time.time(); cfg = self.cfg
        R = self.t["relatedness"]["pairs"]
        if getattr(self.a, "families", None):
            M = families.read_families(self.a.families, sorted(self.analysis))
            g = M[M.family != "unassigned"].groupby("family").size()
            if len(R) and len(g):
                # imported labels are kept; closure, status and membership class are described with the rules of the reconstruction
                given = {f: M.individual[M.family == f].astype(str).tolist() for f in sorted(g.index, key=lambda x: (len(x), x))}
                M, F, info = families.reconstruct(R, sorted(self.analysis), cfg, given=given)
                info = {**info, "edge_rule": "imported with --families; closure described with " + info["edge_rule"], "group_rule": "imported"}
            else:
                F = pd.DataFrame({"family": g.index, "n": g.values}); F["status"] = "imported"; F["closure"] = np.nan
                info = {"edge_rule": "imported with --families", "group_rule": "imported", "has_k2": False}
        else:
            M, F, info = families.reconstruct(R, sorted(self.analysis), cfg)
        self.put("families", "membership", M, "family_membership.tsv", "One row per individual: unit (family, or S_<id> for an unassigned individual), family, membership class (core / attached / unresolved member / unassigned), family size.")
        self.put("families", "families", F, "families.tsv", "One row per family: size, number of within-family pairs, pairs passing the full-sib thresholds, closure (their ratio), status, core size.")
        nf = len(F)
        summ = {"n_families": nf, "n_individuals_in_families": int(F.n.sum()) if nf else 0, "n_unassigned": int((M.family == "unassigned").sum()),
                "n_units": int(M.unit.nunique()), "n_closed": int((F.status == "closed").sum()) if nf else 0,
                "n_not_closed": int((F.status == "not closed").sum()) if nf else 0,
                "n_within_family_pairs": int(F.n_pairs.sum()) if nf and "n_pairs" in F else None,
                "n_within_family_pairs_full_sib": int(F.n_pairs_full_sib.sum()) if nf and "n_pairs_full_sib" in F else None,
                "n_families_at_least_min_size": int((F.n >= cfg["audit"]["min_family_size"]).sum()) if nf else 0, **info}
        conf = self.ran["relatedness"][2] if self.ran.get("relatedness") else "benchmarked"
        why = info["edge_rule"] + "; " + info["group_rule"]
        if self.rel_mode == "grm":
            why += "; families rest on the GRM kinship threshold alone"
        self.status("families", "yes", self.ran["relatedness"][1] if self.ran.get("relatedness") else "", conf, why)
        self.finish("families", summ, t0)

    def units(self):
        M = self.t["families"]["membership"]
        unit_of = dict(zip(M.individual.astype(str), M.unit.astype(str)))
        fam_of = dict(zip(M.individual.astype(str), M.family.astype(str)))
        members = {}
        for s, u in unit_of.items():
            members.setdefault(u, []).append(s)
        return unit_of, fam_of, members

    # ------------------------------------------------------------------ step 3
    def step_nomination(self):
        t0 = time.time(); cfg = self.cfg
        unit_of, _f, _m = self.units()
        R = self.t["relatedness"].get("pairs_nomination", self.t["relatedness"]["pairs"])
        if len(set(unit_of.values())) < 2:
            self.status("nomination", "no", "", "skipped", "fewer than two units")
            self.finish("nomination", {"n_units": len(set(unit_of.values())), "n_nominated": 0}, t0); return
        cls = "relatedness table" if self.rel_mode == "given" else self.primary.name + (" (nomination half)" if self.single_input else "")
        nm = self.s["relatedness"].get("n_markers_nomination_half") or self.s["relatedness"].get("n_markers") or np.nan
        nom = nomination.nominate(R, unit_of, cfg, self.seed, marker_class=cls, n_markers=nm)
        c = cfg["nomination"]
        self.put("nomination", "pairs", nom, "unit_pair_nomination.tsv", "Every pair of units: sizes, cross-unit pairs, pairs at second-degree kinship or closer, their fraction, permutation P, Benjamini-Hochberg q, nominated (screening only).")
        summ = {"n_units": int(len(set(unit_of.values()))), "n_unit_pairs": int(len(nom)), "n_nominated": int(nom.nominated.sum()),
                "n_nominated_family_pairs": int((nom.nominated & (nom.unassigned_sides == 0)).sum()),
                "n_nominated_by_smaller_unit": {k: int(((nom.nominated) & (nom.min_n >= lo) & (nom.min_n <= hi)).sum()) for k, lo, hi in (("at least 4", 4, 10 ** 9), ("3", 3, 3), ("2", 2, 2), ("1 (unassigned individual)", 1, 1))},
                "statistic": f"fraction of cross-unit pairs with theta >= {c['theta_min']}" + (f" and < {c['theta_max']}" if c["theta_max"] is not None else " (no upper bound, as in the original code)"),
                "n_permutations": c["n_permutations"], "q_max": c["q_max"]}
        self.status("nomination", "yes", cls, "provisional" if self.single_input or self.ran["relatedness"][2] != "benchmarked" else "benchmarked",
                    "screening only; " + ("single-input, provisional" if self.single_input else "nomination is not confirmation"))
        self.finish("nomination", summ, t0)

    # ------------------------------------------------------------------ step 4
    def step_audit(self):
        t0 = time.time(); cfg = self.cfg
        if "pairs" not in self.t.get("nomination", {}):
            self.status("audit", "no", "", "skipped", "no nomination table")
            self.finish("audit", {"ran": False}, t0); return
        nom = self.t["nomination"]["pairs"]
        unit_of, _f, _m = self.units()
        sets = list(self.audit_plan)
        if self.single_input:
            if not hasattr(self, "split"):
                self.split = nomination.split_markers(self.primary, cfg)
            sets = [self.primary.subset_markers(self.split[1], "audit half")]
        E, summ = audit.run_audit(nom, unit_of, sets, cfg, self.seed, single_input=self.single_input)
        g2 = pd.DataFrame([{"class": k, **v} for k, v in summ["gate2"].items()])
        keep = [c for c in E.columns]
        self.put("audit", "edges", E[keep], "audit_edges.tsv", "Every pair of units with the audit columns: auditable (both families at least the minimum size), Gate-1 scores (LLR_one, LLR_both), sharing per SV class, evidence and link class.")
        self.put("audit", "gate2", g2, "audit_gate2_summary.tsv", "Gate 2 per SV class: threshold (median of the Gate-1-compatible nominated pairs), matched-null median, one-sided P.")
        aud = E[E.auditable]
        lc = E[E.nominated].link_class.value_counts().to_dict()
        summ.update({"n_auditable_pairs": int(len(aud)), "n_nominated_auditable": int(aud.nominated.sum()),
                     "n_gate1_pass": int((aud.nominated & (aud.gate1 == 1)).sum()) if summ["ran"] else None,
                     "n_rejected": int(aud.link_class.str.startswith("rejected").sum()), "n_supported": int(aud.link_class.str.startswith("supported").sum()),
                     "n_gate1_only": int(aud.link_class.str.startswith("Gate-1 only").sum()), "link_classes": lc})
        if not summ["ran"]:
            self.status("audit", "no", "", "skipped", self.audit_note or "no marker class left for the audit: links are 'nominated, not audited'")
        else:
            self.status("audit", "yes", "; ".join(summ["audit_classes"]), "provisional" if self.single_input else "benchmarked",
                        ("single input split by " + self.split[2] + "; " if self.single_input else "") + summ.get("gate2_note", "Gate 1 and Gate 2; the support bar of Gate 2 is relative (median of the Gate-1-compatible pairs)"))
        self.finish("audit", summ, t0)

    # ------------------------------------------------------------------ step 5
    def step_network(self):
        t0 = time.time()
        if "edges" not in self.t.get("audit", {}):
            self.status("network", "no", "", "skipped", "no nominated links")
            self.finish("network", {"n_links_kept": 0}, t0); return
        E = self.t["audit"]["edges"]
        F = self.t["families"]["families"]
        sizes = dict(zip(F.family.astype(str), F.n.astype(int)))
        net = network.build(E, self.cfg, sizes)
        self.net = net
        s = net["summary"]
        self.put("network", "links", net["edges"], "network_links.tsv", "Nominated links between families with their status (supported / Gate-1 only / provisional / rejected / nominated, not audited) and component.")
        self.put("network", "nodes", net["nodes"], "network_families.tsv", "Families in the network: size, component, links kept, supported links, number of latent parents.")
        self.put("network", "latent_parents", net["latent_parents"], "latent_parents.tsv", "Latent parents (bookkeeping labels, not reconstructed parents): maximal sets of mutually linked families, once per minimum family size.")
        v = pd.DataFrame([{"min_family_size": t, "n_families_violating": len(x), "families": ",".join(x)} for t, x in s["slot_violations"].items()])
        self.put("network", "slot", v, "slot_limit_violations.tsv", "Families whose link partners cannot be split into two groups of mutually linked families (a family has two parents).")
        # opposite homozygotes by relationship class (biSNP R0 and deletions)
        _u, fam_of, _m = self.units()
        kept = net["edges"][net["edges"].status != "rejected"]
        comp_of = dict(zip(net["nodes"].family.astype(str), net["nodes"].component.astype(int)))
        snps = [m for m in self.sets if m.klass == "SNP"]; dls = [m for m in self.sets if m.klass == "DEL"]
        OP, OS, oinfo = opposite.opposite_homozygotes(self.t["relatedness"]["pairs"], fam_of, kept, comp_of, max(snps, key=lambda m: (m.flag == "bisnp", m.n)) if snps else None,
                                                      max(dls, key=lambda m: m.n) if dls else None, self.cfg)
        if len(OS):
            self.put("network", "opposite_pairs", OP, "opposite_homozygotes_pairs.tsv", "Every pair of individuals (no identifiers): kinship, relationship class for the biSNP R0 comparison and for the deletion comparison, R0, opposite-homozygote rate of the deletions and the number of deletions genotyped in both.")
            self.put("network", "opposite_summary", OS, "opposite_homozygotes_by_class.tsv", "Median R0 of the biSNPs and median opposite-homozygote rate of the deletions per relationship class, and the ratio to the reference class (unrelated pairs / pairs without a link).")
            s["opposite_homozygotes"] = {**oinfo, "classes": {f"{r.quantity}: {r['class']}": {"n_pairs": int(r.n_pairs), "ratio": float(r.ratio_to_reference_class)} for _i, r in OS.iterrows()}}
        aud = bool(self.s.get("audit", {}).get("ran"))
        self.status("network", "yes", self.ran["nomination"][1] if self.ran.get("nomination") else "", "benchmarked" if aud and not self.single_input else "provisional",
                    "links audited" if aud else "links are nominated, not audited")
        self.finish("network", {k: v for k, v in s.items()}, t0)

    # ------------------------------------------------------------------ step 6
    def step_haplotypes(self):
        t0 = time.time(); cfg = self.cfg; c = cfg["inversions"]
        snps = [m for m in self.sets if m.klass == "SNP"] or [m for m in self.sets if m.klass == "SNP_INDEL"]
        snp = max(snps, key=lambda m: (m.flag == "bisnp", m.n)) if snps else None
        dls = [m for m in self.sets if m.klass == "DEL"]
        dels = max(dls, key=lambda m: m.n) if dls else None
        _u, fam_of, _m = self.units()
        unrel = self.t["relatedness"]["unrelated"].individual.astype(str).tolist()
        A = self.analysis
        summ = {"scan_ran": False, "notes": []}
        S_parts, I_parts, geno = [], [], {}
        if self.intervals is not None and (snp is not None or dels is not None):
            iv = self.intervals.copy(); iv["source"] = "supplied (--inversions)"
            S, I, g = inversions.genotype_intervals(iv, snp, dels, A, unrel, fam_of, cfg, self.seed)
            S_parts.append(S); I_parts.append(I); geno.update(g)
        scan_on = bool(c["scan"]) and not getattr(self.a, "no_inversion_scan", False) and (snp is not None or dels is not None)
        cand = pd.DataFrame()
        if scan_on:
            parts = []
            if snp is not None:
                C1, W1, n1, sec, step = inversions.scan_snps(snp, cfg)
                summ["snp_scan_seconds"] = round(sec, 1); summ["notes"] += n1
                self.timing["haplotypes: SNP scan"] = sec
                if len(W1):
                    self.put("haplotypes", "scan_windows_snp", W1, "scan_windows_snp.tsv", "SNP scan (local PCA in windows of SNPs, reimplemented from Methods): one row per window with the scaling coordinates, z scores, PC1 share of variance and the heterozygosity ratio.")
                parts.append(C1)
            if dels is not None:
                C2, W2, n2, sec = inversions.scan_dels(dels, cfg)
                summ["del_scan_seconds"] = round(sec, 1); summ["notes"] += n2
                self.timing["haplotypes: deletion scan"] = sec
                if len(W2):
                    self.put("haplotypes", "scan_windows_del", W2, "scan_windows_del.tsv", "Deletion scan: fixed windows, linkage disequilibrium (mean |r|) among deletion dosages.")
                parts.append(C2)
            parts = [p for p in parts if len(p)]
            cand = inversions.combine_scans(parts)
            summ["scan_ran"] = True; summ["n_scan_candidates"] = int(len(cand))
            if len(cand):
                cand["source"] = "scan"
                S, I, g = inversions.genotype_intervals(cand, snp, dels, A, unrel, fam_of, cfg, self.seed + 1000)
                if len(S):
                    S = S.merge(cand[["interval", "scan"]], on="interval", how="left")
                    S["label"] = inversions.SCAN_LABEL
                S_parts.append(S); I_parts.append(I)
                geno_scan = g
            else:
                geno_scan = {}
            if self.intervals is not None:
                ov = inversions.overlap_table(self.intervals, cand)
                self.put("haplotypes", "overlap", ov, "scan_overlap_with_supplied.tsv", "Supplied intervals recovered by the scan, not recovered, and extra candidate intervals found by the scan.")
                summ["n_supplied_recovered"] = int((ov.category == "supplied interval recovered by the scan").sum())
                summ["n_supplied_not_recovered"] = int((ov.category == "supplied interval not recovered").sum())
                summ["n_extra_candidates"] = int((ov.category == "extra candidate interval found by the scan").sum())
        else:
            geno_scan = {}
            why = "turned off (--no-inversion-scan or inversions.scan: false)" if (snp is not None or dels is not None) else "no SNP or deletion input"
            summ["notes"].append("scan not run: " + why)
        S = pd.concat([p for p in S_parts if len(p)], ignore_index=True) if any(len(p) for p in S_parts) else pd.DataFrame(columns=["interval", "chrom", "start", "end", "source", "tier"])
        I = pd.concat([p for p in I_parts if len(p)], ignore_index=True) if any(len(p) for p in I_parts) else pd.DataFrame(columns=["interval", "markers", "individual", "genotype", "pc1", "pc2"])
        # ---- the genotype matrix used by the barcode
        final, info = {}, {}
        if self.intervals is not None:
            for r in S[S.source.str.startswith("supplied")].itertuples():
                if r.interval in geno and (c["barcode_min_tier"] == "all" or r.tier == "high_confidence"):
                    final[r.interval] = geno[r.interval]
                    info[r.interval] = {"source": "supplied interval, genotyped by local PCA", "tier": r.tier, "label": "", "markers": r.genotyped_from}
        elif len(S):
            crit = c["barcode_scan_criteria"]
            for r in S.itertuples():
                if r.interval in geno_scan and all(getattr(r, "pass_" + k, "fail") == "pass" for k in crit):
                    final[r.interval] = geno_scan[r.interval]
                    info[r.interval] = {"source": "scan", "tier": r.tier, "label": inversions.SCAN_LABEL, "markers": r.genotyped_from}
        if self.hap_given is not None:
            six = {s: i for i, s in enumerate(self.hap_given.samples)}
            conc = {}
            for j, name in enumerate(self.hap_given.ids):
                v = np.array([self.hap_given.G[j, six[s]] if s in six else -1 for s in A], np.int8)
                for src in (geno, geno_scan):
                    if name in src:
                        both = (v >= 0) & (src[name] >= 0)
                        if both.sum():
                            conc[name] = float(max((v[both] == src[name][both]).mean(), (v[both] == 2 - src[name][both]).mean()))
                final[name] = v
                info[name] = {"source": "supplied genotypes (--haplotypes)", "tier": "", "label": "", "markers": "supplied"}
            if conc and len(S):
                S["concordance_with_supplied_genotypes"] = S.interval.map(conc)
                summ["median_concordance_with_supplied_genotypes"] = float(np.median(list(conc.values())))
        names = list(final)
        W = pd.DataFrame({"individual": A, **{n: final[n] for n in names}})
        self.put("haplotypes", "summary", S, "structural_haplotypes.tsv", "One row per interval (supplied or found by the scan): markers used, cluster sizes, the five criteria, tier. Genotype labels 0 / 2 are the two homokaryotype clusters (0 = the more frequent); they carry no ancestral meaning. Intervals found by the scan are candidates, detected automatically, not validated.")
        self.put("haplotypes", "per_individual", I, "structural_haplotype_calls_long.tsv", "One row per interval, marker type and individual: genotype (-1 = ambiguous, set to missing), forced genotype, posterior, PC1, PC2, heterozygosity inside the interval.")
        self.put("haplotypes", "genotypes", W, "structural_haplotype_genotypes.tsv", "Individuals x intervals genotype matrix used by the barcode and the family effect (0 / 1 / 2, -1 missing). Can be passed back as --haplotypes.")
        self.put("haplotypes", "interval_info", pd.DataFrame([{"interval": k, **v} for k, v in info.items()], columns=["interval", "source", "tier", "label", "markers"]),
                 "barcode_interval_sources.tsv", "Where the genotypes of each barcode interval come from.")
        gl = {**geno_scan, **geno, **{k: final[k] for k in final}}
        sv_sets = [m for m in self.sets if m.klass in ("DEL", "DUP")]
        allint = pd.concat([S[["interval", "chrom", "start", "end"]]] + ([pd.DataFrame({"interval": self.hap_given.ids, "chrom": self.hap_given.chrom, "start": self.hap_given.pos, "end": self.hap_given.end})] if self.hap_given is not None else []),
                           ignore_index=True).drop_duplicates("interval") if len(S) or self.hap_given is not None else pd.DataFrame(columns=["interval", "chrom", "start", "end"])
        allint = allint[allint.end.astype(int) > 0] if len(allint) else allint
        if len(allint) and sv_sets:
            TG = inversions.tag_variants(allint, gl, A, sv_sets, cfg)
            self.put("haplotypes", "tags", TG, "arrangement_tag_variants.tsv", "Deletions / duplications whose carrier frequency differs between the two homokaryotype classes by at least the configured difference; fully_private = carried by every called individual of one class and none of the other.")
            summ["n_tag_variants"] = int(len(TG)); summ["n_tag_variants_fully_private"] = int(TG.fully_private.sum()) if len(TG) else 0
        if len(I) and len(S):
            SC, SA = inversions.subclusters(S[S.tier != "not genotyped"], I, dels, cfg, self.seed)
            self.put("haplotypes", "subclusters", SC, "haplotype_subclusters.tsv", "Sub-clusters inside each genotype class (Gaussian mixture on PC1 and PC2) and whether deletion dosage supports each split (label-permutation null).")
            self.put("haplotypes", "subcluster_assignments", SA, "haplotype_subcluster_assignments.tsv", "Sub-cluster of every individual, per interval and genotype class.")
            if len(SC) and "split_supported_by_deletions" in SC:
                sp = SC[SC.n_subclusters > 1]
                summ["n_classes_split"] = int(len(sp)); summ["n_splits_supported_by_deletions"] = int((sp.split_supported_by_deletions == True).sum())  # noqa: E712
        g_ok = S[S.tier != "not genotyped"] if len(S) else S
        summ.update({"n_intervals_supplied": 0 if self.intervals is None else int(len(self.intervals)), "n_intervals_genotyped": int(len(g_ok)),
                     "n_high_confidence": int((g_ok.tier == "high_confidence").sum()) if len(g_ok) else 0,
                     "n_candidate": int((g_ok.tier == "candidate").sum()) if len(g_ok) else 0,
                     "n_barcode_intervals": len(names), "genotyped_from": sorted(set(g_ok.genotyped_from)) if len(g_ok) else [],
                     "haplotypes_supplied": self.hap_given is not None})
        for n in summ["notes"]:
            self.note("structural haplotypes: " + n)
        if not len(names) and not len(S):
            self.status("haplotypes", "no", "", "skipped", "; ".join(summ["notes"]) or "no --inversions, no --haplotypes and nothing found by the scan")
        else:
            used = ", ".join(x for x in [snp.name if snp is not None else "", dels.name if dels is not None else "", "--haplotypes" if self.hap_given is not None else ""] if x)
            conf = "benchmarked" if (self.intervals is not None or self.hap_given is not None) else "provisional"
            self.status("haplotypes", "yes", used, conf, ("genotyping extracted from the original code; " if self.intervals is not None else "") +
                        ("scan: SNP part reimplemented from Methods, candidates not validated" if scan_on else "scan not run"))
        self.finish("haplotypes", summ, t0)

    def hap_set(self):
        W = self.t.get("haplotypes", {}).get("genotypes")
        if W is None or W.shape[1] < 2:
            return None, {}
        names = [c for c in W.columns if c != "individual"]
        G = W[names].to_numpy(np.int8).T
        chrom, pos, end = [], [], []
        for x in names:
            m = re.match(r"(?:scan_)?(.+):(\d+)-(\d+)$", x)
            chrom.append(m.group(1) if m else "NA"); pos.append(int(m.group(2)) if m else 0); end.append(int(m.group(3)) if m else 0)
        ms = MarkerSet("HAPLOTYPE (structural haplotypes)", "HAPLOTYPE", "haplotypes", "structural_haplotype_genotypes.tsv", W.individual.astype(str).tolist(),
                       np.array(chrom), np.array(pos, np.int64), np.array(end, np.int64), np.array(names), G)
        II = self.t["haplotypes"]["interval_info"]
        info = {r.interval: {"source": r.source, "tier": r.tier, "label": r.label, "markers": r.markers} for r in II.itertuples()}
        return ms, info

    # ------------------------------------------------------------------ step 7
    def step_barcode(self):
        t0 = time.time(); cfg = self.cfg
        hap, info = self.hap_set()
        links_t = self.t.get("network", {}).get("links")
        if hap is None:
            self.status("barcode", "no", "", "skipped", "no structural-haplotype genotypes (give --haplotypes or --inversions, or let the scan find intervals)")
            self.finish("barcode", {"n_states": 0}, t0); return
        if links_t is None or not len(links_t):
            self.status("barcode", "no", "", "skipped", "no link between families")
            self.finish("barcode", {"n_states": 0}, t0); return
        mn = cfg["audit"]["min_family_size"]
        L = links_t[(links_t.status != "rejected") & (links_t.min_n >= mn)]
        links = [(r.unit_a, r.unit_b, r.link_class) for r in L.itertuples()]
        _u, _f, members = self.units()
        junc = self.s["network"].get("junction_families", {})
        B = barcode.barcode(hap, members, links, junc, cfg, info)
        self.put("barcode", "barcode", B, "shared_parent_barcode.tsv", "Posterior genotype of the shared parent at every structural haplotype, per link and per junction family (same parent assumed on every supported link: candidate, not tested). call = genotype at posterior >= the configured threshold, else unresolved.")
        pc, ps = cfg["barcode"]["posterior_call"], cfg["barcode"]["posterior_strict"]
        per = {p: {"n_states": int(len(g)), "n_called": int((g.posterior >= pc).sum()), "n_called_strict": int((g.posterior >= ps).sum())} for p, g in B.groupby("parent")} if len(B) else {}
        summ = {"n_intervals": int(hap.n), "n_links": len(links), "n_junctions": len(junc), "n_states": int(len(B)),
                "n_called": int((B.posterior >= pc).sum()) if len(B) else 0, "n_called_strict": int((B.posterior >= ps).sum()) if len(B) else 0,
                "posterior_call": pc, "posterior_strict": ps, "per_parent": per}
        src = sorted({v["source"] for v in info.values()})
        scan_only = src == ["scan"]
        self.status("barcode", "yes", "; ".join(src), "provisional" if scan_only or self.single_input or not self.s.get("audit", {}).get("ran") else "benchmarked",
                    "intervals found by the scan: candidate, detected automatically, not validated" if scan_only else "posterior of the shared parent from the descendants' genotypes")
        self.finish("barcode", summ, t0)

    # ------------------------------------------------------------------ step 8
    def step_family_effect(self):
        t0 = time.time()
        _u, fam_of, _m = self.units()
        hap, _info = self.hap_set()
        rows, parts = [], []
        for ms in self.sets + ([hap] if hap is not None else []):
            d, row = family_effect.family_effect(ms, fam_of, self.cfg, self.seed)
            rows.append(row)
            if d is not None:
                parts.append(d)
        # frequency-matched SNPs as the reference level of every non-SNP class (when SNPs are given)
        snps = [m for m in self.sets if m.klass == "SNP"]
        if snps and parts:
            snp = max(snps, key=lambda m: (m.flag == "bisnp", m.n))
            targets = {}
            for d in parts:
                k = str(d.marker_class.iloc[0])
                if k not in POINT and k not in targets:
                    targets[k] = d.af.to_numpy()
            for k, m in family_effect.matched_snps(snp, fam_of, self.cfg, self.seed, targets).items():
                d, row = family_effect.family_effect(m, fam_of, self.cfg, self.seed + 1)
                row["note"] = (row.get("note", "") + "; " if row.get("note") else "") + f"one SNP per tested {k} marker, drawn from the same allele-frequency bin (reference level, not an independent null)"
                rows.append(row)
                if d is not None:
                    parts.append(d)
        T = pd.DataFrame(rows)
        self.put("family_effect", "summary", T, "family_effect_summary.tsv", "Per marker class: number of markers tested, median eta squared (variance explained by family), median of the per-marker null medians (family labels permuted), share above the null median / 95th percentile / every permutation, Benjamini-Hochberg count.")
        if parts:
            self.put("family_effect", "per_marker", pd.concat(parts, ignore_index=True), "family_effect_per_marker.tsv", "Per marker: observed eta squared, null median, null 95th percentile, null maximum, empirical P.")
        ok = T[T.n_tested > 0] if len(T) else T
        summ = {"classes": {r.marker_set: {"marker_class": r.marker_class, "n_tested": int(r.n_tested), "eta2_median": float(r.eta2_obs_median), "null_median": float(r.null_median_of_medians)} for r in ok.itertuples()} if len(ok) else {},
                "n_permutations": self.cfg["family_effect"]["n_permutations"]}
        if len(ok):
            self.status("family_effect", "yes", "; ".join(ok.marker_set), "benchmarked", "once per marker class; individuals assigned to a family only")
        else:
            self.status("family_effect", "no", "", "skipped", "; ".join(str(x) for x in T.get("note", [])) or "no family")
        self.finish("family_effect", summ, t0)

    # ------------------------------------------------------------------ step 9
    def step_concordance(self):
        """Kinship of every marker class (VanRaden GRM) within and between families: AUC, recovery
        and their family-level bootstrap intervals (cross-marker concordance, not pedigree truth:
        the families were built from one of the inputs). The GRM of every class is written out."""
        t0 = time.time()
        _u, fam_of, _m = self.units()
        hap, _i = self.hap_set()
        A = list(self.analysis)
        sets, masked_note = concordance.marker_sets(self.sets, hap, self.intervals, self.cfg)
        T, pairs, grm, order = concordance.run(sets, fam_of, A, self.cfg, self.primary.name if self.primary is not None else None)
        if masked_note:
            self.note("cross-marker concordance: " + masked_note)
        self.grm = grm
        self.put("concordance", "summary", T, "cross_marker_concordance.tsv", "Per marker class: markers used by the GRM, rank AUC of within-family against between-family pairs of assigned individuals, recovery_pct (share of within-family pairs above the 99.9th percentile of the between-family pairs), their 95% intervals from a bootstrap over families, median GRM value within and between, and the Pearson correlation of the GRM with the GRM of the SNPs outside SV regions (all pairs). Concordance between marker classes, not pedigree truth.")
        if len(pairs):
            self.put("concordance", "pairs", pairs, "cross_marker_kinship_pairs.tsv", "Kinship (GRM / 2) of every pair of assigned individuals per marker class, with the same-family flag (no identifiers).")
        if grm:
            d = os.path.join(self.out, FOLDER["concordance"]); os.makedirs(d, exist_ok=True)
            np.savez_compressed(os.path.join(d, "kinship_grm.npz"), **{k: v.astype(np.float32) for k, v in grm.items()})
            self.put("concordance", "grm_order", order, "kinship_grm_individuals.tsv", "Row / column order of the matrices in kinship_grm.npz (compressed NumPy archive, one array per marker class, named as in the column grm_key of cross_marker_concordance.tsv; VanRaden GRM, individuals x individuals, NaN for an individual absent from that input), with the family of each individual. Individuals are in family order (largest family first), unassigned individuals last.")
        ok = T[T.n_markers > 0] if len(T) else T
        summ = {"classes": {r.marker_set: {"n_markers": int(r.n_markers), "auc": float(r.auc_within_vs_between), "recovery_pct": float(r.recovery_pct),
                                           "auc_ci95": [float(r.auc_ci95_low), float(r.auc_ci95_high)], "recovery_ci95": [float(r.recovery_ci95_low), float(r.recovery_ci95_high)]}
                            for r in ok.itertuples() if "auc_within_vs_between" in T and pd.notna(r.auc_within_vs_between)},
                "n_bootstrap": self.cfg["concordance"]["n_bootstrap"], "grm_file": "kinship_grm.npz" if grm else None}
        if summ["classes"]:
            self.status("concordance", "yes", "; ".join(summ["classes"]), "benchmarked", "cross-marker concordance against the reconstructed families, not pedigree truth")
        else:
            self.status("concordance", "no", "", "skipped", "no genotype input or no family")
        self.finish("concordance", summ, t0)

    def get_grm(self):
        """GRM per marker class of step 9 (from memory, or read back from kinship_grm.npz)."""
        if getattr(self, "grm", None):
            return self.grm
        p = os.path.join(self.out, FOLDER["concordance"], "kinship_grm.npz")
        if not os.path.exists(p):
            return {}
        with np.load(p) as z:
            self.grm = {k: z[k].astype(float) for k in z.files}
        return self.grm

    # ------------------------------------------------------------------ driver
    def execute(self, compute=None):
        """compute: set of steps to compute; the steps before the first of them are loaded from tables."""
        compute = set(STEPS) if compute is None else set(compute)
        first = min(STEPS.index(s) for s in compute) if compute else len(STEPS)
        t_all = time.time()
        self.load_inputs()
        for i, step in enumerate(STEPS):
            if step in compute:
                log(f"[{i + 1}/{len(STEPS)}] {FOLDER[step]}: {WHAT[step]}")
                getattr(self, "step_" + step)()
            else:
                try:
                    self.load(step)
                except InputError:
                    if i < first:
                        raise
                    continue
                if step == "qc":
                    self.restore_analysis()
        self.timing["total"] = time.time() - t_all
        return self
