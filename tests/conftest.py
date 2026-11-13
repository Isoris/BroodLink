"""Shared fixtures: the synthetic cohort in tests/data and a few pipeline runs (one per session)."""
import gzip
import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, "tests", "data")
EXAMPLE_CFG = os.path.join(ROOT, "example", "config.example.yaml")

from broodlink.cli import main  # noqa: E402

FULL = ["--relatedness", os.path.join(DATA, "relatedness.tsv"), "--bisnp", os.path.join(DATA, "snps.beagle.gz"), "--gt", os.path.join(DATA, "sv.vcf.gz"),
        "--inversions", os.path.join(DATA, "inversions.bed"), "--mask", os.path.join(DATA, "callable.bed"), "--genome", os.path.join(DATA, "genome.fai"),
        "--samples", os.path.join(DATA, "samples.tsv"), "--depth", os.path.join(DATA, "depth.tsv"), "--config", EXAMPLE_CFG]

QUICK_CFG = """seed: 1
nomination: {n_permutations: 1000}
family_effect: {n_permutations: 50}
inversions: {min_unrelated: 20, snp_scan: {window_snps: 25, z_outlier: 2.5}, tier: {eta2_permutations: 50}, subclusters: {restarts: 2, n_permutations: 20}}
"""


def read(out, rel):
    return pd.read_csv(os.path.join(out, rel), sep="\t", comment="#", dtype={"individual": str, "sample": str})


def summary(out, folder):
    import json
    with open(os.path.join(out, folder, "summary.json"), encoding="utf-8") as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def truth():
    T = pd.read_csv(os.path.join(DATA, "truth", "truth_families.tsv"), sep="\t", dtype=str).fillna("")
    L = pd.read_csv(os.path.join(DATA, "truth", "truth_links.tsv"), sep="\t", dtype=str)
    P = pd.read_csv(os.path.join(DATA, "truth", "truth_parent_haplotypes.tsv"), sep="\t", dtype={"parent": str}).set_index("parent")
    K = pd.read_csv(os.path.join(DATA, "truth", "truth_karyotypes.tsv"), sep="\t", dtype={"individual": str}).set_index("individual")
    Q = pd.read_csv(os.path.join(DATA, "truth", "truth_qc.tsv"), sep="\t", dtype=str)
    return {"families": T, "links": L, "parents": P, "karyotypes": K, "qc": Q}


@pytest.fixture(scope="session")
def quick_cfg(tmp_path_factory):
    p = tmp_path_factory.mktemp("cfg") / "quick.yaml"
    p.write_text(QUICK_CFG, encoding="utf-8")
    return str(p)


@pytest.fixture(scope="session")
def full_run(tmp_path_factory):
    """The README quick-start command (example config, all figures)."""
    out = str(tmp_path_factory.mktemp("full") / "out")
    main(["run", *FULL, "--out", out])
    return out


@pytest.fixture(scope="session")
def snp_vcf(tmp_path_factory):
    """A hard-call SNP VCF made from the Beagle likelihoods, standing in for a small-variant cohort VCF."""
    from broodlink.config import load_config
    from broodlink.io import read_beagle
    ms = read_beagle([os.path.join(DATA, "snps.beagle.gz")])[0]
    p = str(tmp_path_factory.mktemp("vcf") / "small_variants.vcf.gz")
    code = np.array(["0/0", "0/1", "1/1", "./."])
    with gzip.open(p, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ms.samples) + "\n")
        for j in range(ms.n):
            fh.write(f"{ms.chrom[j]}\t{ms.pos[j]}\t.\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(code[ms.G[j]]) + "\n")
    return p


def run_mode(tmp_path_factory, name, args, quick_cfg):
    out = str(tmp_path_factory.mktemp(name) / "out")
    main(["run", *args, "--config", quick_cfg, "--out", out, "--no-figures", "--no-inversion-scan"])
    return out


@pytest.fixture(scope="session")
def run_relatedness_only(tmp_path_factory, quick_cfg):
    return run_mode(tmp_path_factory, "rel", ["--relatedness", os.path.join(DATA, "relatedness.tsv")], quick_cfg)


@pytest.fixture(scope="session")
def run_bisnp_only(tmp_path_factory, quick_cfg):
    return run_mode(tmp_path_factory, "bisnp", ["--bisnp", os.path.join(DATA, "snps.beagle.gz")], quick_cfg)


@pytest.fixture(scope="session")
def run_gt_only(tmp_path_factory, quick_cfg):
    return run_mode(tmp_path_factory, "gt", ["--gt", os.path.join(DATA, "sv.vcf.gz")], quick_cfg)


@pytest.fixture(scope="session")
def run_snp_indel_only(tmp_path_factory, quick_cfg, snp_vcf):
    return run_mode(tmp_path_factory, "snpindel", ["--snp_indel", snp_vcf], quick_cfg)


@pytest.fixture(scope="session")
def run_with_haplotypes(tmp_path_factory, quick_cfg):
    return run_mode(tmp_path_factory, "hap", ["--relatedness", os.path.join(DATA, "relatedness.tsv"), "--gt", os.path.join(DATA, "sv.vcf.gz"),
                                               "--haplotypes", os.path.join(DATA, "haplotypes.tsv")], quick_cfg)
