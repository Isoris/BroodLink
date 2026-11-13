"""A quick benchmark on the synthetic cohort: each input mode against the simulated truth, with
run time. Small and clearly scoped: one simulated cohort of 80 individuals, one seed. It shows
that every mode runs and what each mode recovers on easy data; it is not a performance study.
"""
from __future__ import annotations
import gzip
import json
import os
import time

import numpy as np
import pandas as pd

from .evaluate import family_map, family_scores, link_scores, link_set
from .io import log, read_beagle


def _snp_vcf(data, out):
    """Hard-call SNP VCF made from the Beagle likelihoods (stands in for a small-variant cohort VCF)."""
    ms = read_beagle([os.path.join(data, "snps.beagle.gz")])[0]
    p = os.path.join(out, "small_variants_from_beagle.vcf.gz")
    code = np.array(["0/0", "0/1", "1/1", "./."])
    with gzip.open(p, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ms.samples) + "\n")
        for j in range(ms.n):
            fh.write(f"{ms.chrom[j]}\t{ms.pos[j]}\t.\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(code[ms.G[j]]) + "\n")
    return p


def _read(out, rel, **kw):
    return pd.read_csv(os.path.join(out, rel), sep="\t", comment="#", dtype={"individual": str, "parent": str, "families": str}, **kw)


def score_run(out, data):
    T = pd.read_csv(os.path.join(data, "truth/truth_families.tsv"), sep="\t", dtype=str).fillna("")
    L = pd.read_csv(os.path.join(data, "truth/truth_links.tsv"), sep="\t", dtype=str)
    ref_fam = dict(zip(T.individual, T.family)); ref_links = {frozenset((a, b)) for a, b in zip(L.family_a, L.family_b)}
    M = _read(out, "02_families/family_membership.tsv")
    row = {**family_scores(M, ref_fam)}
    fm = family_map(M, ref_fam)
    p = os.path.join(out, "05_network", "network_links.tsv")
    if os.path.exists(p):
        N = _read(out, "05_network/network_links.tsv")
        kept = link_set(N, fm, N.status != "rejected"); sup = link_set(N, fm, N.status == "supported")
        row.update(link_scores(kept, ref_links))
        row.update({"supported_true": len(sup & ref_links), "supported_false": len(sup - ref_links),
                    "link_labels": "; ".join(sorted(set(N.link_class)))})
    p = os.path.join(out, "07_barcode", "shared_parent_barcode.tsv")
    if os.path.exists(p):
        B = _read(out, "07_barcode/shared_parent_barcode.tsv"); B = B[B.kind == "link"]
        P = pd.read_csv(os.path.join(data, "truth/truth_parent_haplotypes.tsv"), sep="\t", dtype={"parent": str}).set_index("parent")
        K = pd.read_csv(os.path.join(data, "truth/truth_karyotypes.tsv"), sep="\t", dtype={"individual": str}).set_index("individual")
        W = _read(out, "06_structural_haplotypes/structural_haplotype_genotypes.tsv").set_index("individual")
        parent = {frozenset((a, b)): q for a, b, q in zip(L.family_a, L.family_b, L.shared_parent)}
        flip = {}
        for c in W.columns:                                       # labels 0 / 2 are arbitrary: orient each interval on the truth
            if c in K:
                g = W[c].reindex(K.index); ok = g >= 0
                flip[c] = bool((g[ok] == 2 - K[c][ok]).mean() > (g[ok] == K[c][ok]).mean())
        code = {"0/0": 0, "0/1": 1, "1/1": 2}; hit = tot = n = 0
        for r in B.itertuples():
            a, b = r.families.split(",")
            key = frozenset((fm.get(a), fm.get(b)))
            if key not in parent or r.interval not in flip:
                continue
            n += 1
            if r.call != "unresolved":
                g = code[r.call]; g = 2 - g if flip[r.interval] else g
                tot += 1; hit += int(g == int(P.loc[parent[key], r.interval]))
        row.update({"barcode_states": n, "barcode_called": tot, "barcode_called_correct": hit})
        if len(flip):
            acc = []
            for c in flip:
                g = W[c].reindex(K.index); ok = g >= 0
                acc.append(max((g[ok] == K[c][ok]).mean(), (g[ok] == 2 - K[c][ok]).mean()))
            row["karyotype_accuracy_min"] = round(float(min(acc)), 4)
    return row


def run_benchmark(data, out, cfg):
    from .cli import main
    os.makedirs(out, exist_ok=True)
    j = lambda f: os.path.join(data, f)
    cfgp = os.path.join(out, "benchmark_config.yaml")
    c = {"seed": cfg["seed"], "nomination": {"n_permutations": cfg["nomination"]["n_permutations"]}, "family_effect": {"n_permutations": 100},
         "inversions": {"min_unrelated": 20, "snp_scan": {"window_snps": 25, "z_outlier": 2.5}}}
    import yaml
    with open(cfgp, "w", encoding="utf-8") as fh:
        yaml.safe_dump(c, fh)
    vcf = _snp_vcf(data, out)
    common = ["--mask", j("callable.bed"), "--genome", j("genome.fai")]
    modes = [("full run (relatedness + bisnp + gt + inversions BED)", ["--relatedness", j("relatedness.tsv"), "--bisnp", j("snps.beagle.gz"), "--gt", j("sv.vcf.gz"), "--inversions", j("inversions.bed")]),
             ("relatedness + gt + haplotypes (supplied genotypes)", ["--relatedness", j("relatedness.tsv"), "--gt", j("sv.vcf.gz"), "--haplotypes", j("haplotypes.tsv"), "--no-inversion-scan"]),
             ("bisnp + gt, no relatedness table (GRM from SNPs)", ["--bisnp", j("snps.beagle.gz"), "--gt", j("sv.vcf.gz"), "--inversions", j("inversions.bed")]),
             ("--relatedness only", ["--relatedness", j("relatedness.tsv")]),
             ("--bisnp only", ["--bisnp", j("snps.beagle.gz"), "--no-inversion-scan"]),
             ("--gt only (DEL + DUP in one VCF)", ["--gt", j("sv.vcf.gz"), "--no-inversion-scan"]),
             ("--snp_indel only (hard-call VCF)", ["--snp_indel", vcf, "--no-inversion-scan"])]
    rows = []
    for k, (name, args) in enumerate(modes):
        o = os.path.join(out, f"mode_{k + 1}")
        t0 = time.time()
        main(["run", *args, *common, "--config", cfgp, "--out", o, "--no-figures"])
        sec = time.time() - t0
        with open(os.path.join(o, "01_relatedness", "summary.json"), encoding="utf-8") as fh:
            rel = json.load(fh)
        row = {"mode": name, "relatedness": "table" if rel["summary"]["mode"] == "given" else f"GRM, {rel['summary'].get('n_markers')} markers", **score_run(o, data), "seconds": round(sec, 1)}
        rows.append(row)
        pd.DataFrame(rows).to_csv(os.path.join(out, "benchmark_synthetic.tsv"), sep="\t", index=False, lineterminator="\n")       # partial results are kept
    T = pd.DataFrame(rows)
    with open(os.path.join(out, "benchmark_report.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(report_md(T, cfg))
    log(f"benchmark written to {out}")
    return T


def report_md(T, cfg):
    def cell(x):
        if isinstance(x, float) and np.isnan(x):
            return "-"
        if isinstance(x, float) and float(x).is_integer():
            return str(int(x))
        return str(x)
    cols = ["mode", "relatedness", "families_identical", "families_reference", "pair_precision", "pair_recall", "links_true", "links_false", "links_reference",
            "supported_true", "supported_false", "barcode_called", "barcode_called_correct", "barcode_states", "seconds"]
    cols = [c for c in cols if c in T.columns]
    L = [f"Synthetic cohort, seed {cfg['seed']}, {cfg['nomination']['n_permutations']} permutations, figures off. Seconds are wall-clock time of the whole run on the machine that produced this table.", "",
         "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    L += ["| " + " | ".join(cell(x) for x in r) + " |" for r in T[cols].itertuples(index=False)]
    return "\n".join(L) + "\n"
