"""Command line: `broodlink run ...`, one subcommand per step, `report`, `sweep`, `benchmark`."""
from __future__ import annotations
import argparse
import sys

from . import __version__
from .config import load_config
from .pipeline import FOLDER, STEPS, Run

STEP_CMD = {"qc": "qc", "relatedness": "relatedness", "families": "families", "nominate": "nomination", "audit": "audit",
            "network": "network", "haplotypes": "haplotypes", "barcode": "barcode", "family-effect": "family_effect",
            "concordance": "concordance"}


def add_inputs(p):
    g = p.add_argument_group("inputs (at least one of --bisnp, --gt, --snp_indel, --vcf, --relatedness)")
    g.add_argument("--bisnp", nargs="+", metavar="FILE", help="biallelic SNPs from ANGSD: Beagle genotype likelihoods (file, files or folder), or BCF/VCF with GL/PL; low coverage accepted")
    g.add_argument("--gt", action="append", metavar="FILE", help="genotype VCF with GT (Manta, DELLY, ...; repeatable; class DEL / DUP / INV / INS read from SVTYPE), or a genotype matrix TSV")
    g.add_argument("--snp_indel", metavar="FILE", help="small-variant VCF or gVCF with GT (e.g. a Clair3 cohort)")
    g.add_argument("--vcf", action="append", metavar="FILE", help="any other variant VCF with GT (repeatable)")
    g.add_argument("--relatedness", metavar="FILE_OR_DIR", help="precomputed relatedness: NGSrelate .res (or a folder of them), KING .kin/.kin0, or a TSV (id1, id2, theta or kinship, k0, k1, k2)")
    g.add_argument("--haplotypes", metavar="FILE", help="structural-haplotype genotypes per individual (interval id, 0/1/2), for the barcode; overrides pipeline calls for the intervals it covers")
    g.add_argument("--inversions", metavar="BED", help="intervals of large structural haplotypes to genotype by local PCA (chrom, start, end, name)")
    g.add_argument("--no-inversion-scan", action="store_true", help="do not scan for large structural haplotypes (the scan is on by default when variants are given)")
    g.add_argument("--mask", action="append", metavar="FILE", help="callable regions to KEEP (BED) or a soft-masked FASTA (uppercase = callable); several masks are intersected")
    g.add_argument("--mask-exclude", action="append", metavar="FILE", help="regions to drop (BED)")
    g.add_argument("--samples", metavar="FILE", help="sample list / metadata TSV (id, optional sex, batch, group)")
    g.add_argument("--depth", metavar="FILE_OR_DIR", help="per-sample depth table or a folder of mosdepth / samtools-depth summaries (QC only)")
    g.add_argument("--genome", metavar="FAI", help="FASTA index of the assembly (chromosome lengths)")
    g.add_argument("--families", metavar="FILE", help="import family labels (individual, family) instead of reconstructing them")
    g.add_argument("--genes", metavar="FILE", help="gene models of the assembly, GFF3 (gene features, CDS when present) or BED (chrom, start, end, name, score, strand); optional: genes of each barcode interval and the gene-level barcode figure")
    o = p.add_argument_group("run")
    o.add_argument("--auto-exclude", action="store_true", help="remove low-coverage samples and one sample of each pair with identical genotypes (default: flag only)")
    o.add_argument("--out", default="broodlink_out", metavar="DIR")
    o.add_argument("--config", metavar="config.yaml")
    o.add_argument("--seed", type=int, help="overrides the seed of the config")
    o.add_argument("--permutations", type=int, metavar="N", help="overrides nomination.n_permutations (for quick tests)")
    o.add_argument("--no-figures", action="store_true", help="skip the figures")
    o.add_argument("-q", "--quiet", action="store_true", help="print only the step headers and their results (no details, no progress of long loops)")


def build_parser():
    ap = argparse.ArgumentParser(prog="broodlink", description="BroodLink: full-sib families, candidate shared-parent links and structural-haplotype barcodes "
                                 "from one cohort of related individuals without a pedigree. Families are relatedness-based groups, not pedigree records.")
    ap.add_argument("--version", action="version", version=f"broodlink {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    add_inputs(sub.add_parser("run", help="run every step, then the tables, statistics, figures, tracks and reports"))
    for name, step in STEP_CMD.items():
        add_inputs(sub.add_parser(name, help=f"re-run step {FOLDER[step]} alone (earlier steps are read from --out; later steps are not re-run)"))
    add_inputs(sub.add_parser("report", help="rebuild figures, tables, tracks and reports from the tables in --out"))
    sp = sub.add_parser("sweep", help="parameter sweeps with a discriminator for every setting (synthetic data by default)")
    sp.add_argument("--data", default="tests/data", help="folder made by tests/data/simulate_pedigree.py")
    sp.add_argument("--out", default="sweep_out")
    sp.add_argument("--config")
    sp.add_argument("--quick", action="store_true", help="fewer values and permutations")
    sp.add_argument("--only", nargs="+", help="names of the sweeps to run")
    bp = sub.add_parser("benchmark", help="a quick benchmark on the synthetic data: each input mode against the truth, with run time")
    bp.add_argument("--data", default="tests/data")
    bp.add_argument("--out", default="benchmark_out")
    bp.add_argument("--config")
    bp.add_argument("--permutations", type=int, default=2000)
    return ap


def make_cfg(a):
    cfg = load_config(getattr(a, "config", None))
    if getattr(a, "seed", None) is not None:
        cfg["seed"] = a.seed
    if getattr(a, "permutations", None):
        cfg["nomination"]["n_permutations"] = a.permutations
    return cfg


def main(argv=None):
    a = build_parser().parse_args(argv)
    from .io import set_verbosity
    # sweeps and the benchmark run the pipeline many times: keep their output to one line per run
    set_verbosity(1 if getattr(a, "quiet", False) or a.cmd in ("sweep", "benchmark") else 2)
    cfg = make_cfg(a)
    if a.cmd == "sweep":
        from .sweep import run_sweeps
        run_sweeps(a.data, a.out, cfg, quick=a.quick, only=a.only); return 0
    if a.cmd == "benchmark":
        from .benchmark import run_benchmark
        run_benchmark(a.data, a.out, cfg); return 0
    from .report import write_all
    run = Run(a, cfg)
    if a.cmd == "run":
        run.execute()
    elif a.cmd == "report":
        run.execute(compute=[])
    else:
        run.execute(compute=[STEP_CMD[a.cmd]])
        run.note(f"only step {FOLDER[STEP_CMD[a.cmd]]} was re-run; later steps were read from their existing tables and may be out of date")
    write_all(run, figures=not a.no_figures)
    return 0


if __name__ == "__main__":
    sys.exit(main())
