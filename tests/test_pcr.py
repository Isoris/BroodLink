"""Step 10, PCR markers: Tm model, design engines, FASTA reader, the step on the example cohort,
skipping without sequence, the built-in engine, --reference with in-silico PCR, figure paging."""
import gzip
import importlib.util
import os
import random
import shutil
import sys

import numpy as np
import pandas as pd
import pytest

from conftest import DATA, FULL, main, read, summary
from broodlink import pcr
from broodlink.config import load_config

HAVE_P3 = importlib.util.find_spec("primer3") is not None
needs_primer3 = pytest.mark.skipif(not HAVE_P3, reason="optional package primer3-py is not installed (pip install broodlink[primers])")


def _random_seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


# ------------------------------------------------------------------ units
def test_tm_formula_known_value_and_conditions():
    # value of primer3-py calc_tm (2.4.0) for this oligo: 59.758
    assert abs(pcr.tm_formula("ACGTTGCATGCAAGTCCATG") - 59.75825) < 1e-4
    assert pcr.tm_formula("ACGTNACGTACGT") is None
    assert pcr.tm_toolkit_fallback("ACGTTGCATGCAAGTCCATG") < pcr.tm_formula("ACGTTGCATGCAAGTCCATG") - 8      # the toolkit fallback runs about 10 C cooler
    assert "50 mM monovalent" in pcr.tm_conditions() and "50 nM primer" in pcr.tm_conditions()


@needs_primer3
def test_tm_formula_equals_primer3_calc_tm():
    import primer3
    rnd = random.Random(3)
    for _ in range(500):
        s = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(16, 30)))
        assert abs(pcr.tm_formula(s) - primer3.calc_tm(s)) < 1e-6, s
    assert abs(pcr.tm("ACGTTGCATGCAAGTCCATG") - primer3.calc_tm("ACGTTGCATGCAAGTCCATG")) < 1e-9


def _check_pair(res, template, bp, c_len, c_tm):
    f, r = res["fwd"], res["rev"]
    fs = template.find(f); rs = template.find(pcr.revcomp(r))
    assert fs >= 0 and rs >= 0 and fs + len(f) <= bp <= rs                       # both primers in the template, on either side of the breakpoint
    assert res["band_deleted"] == rs + len(r) - fs and res["fpos"] == fs and res["rpos"] == rs + len(r) - 1
    for p, t in ((f, res["Tm_F"]), (r, res["Tm_R"])):
        assert c_len[0] <= len(p) <= c_len[1] and c_tm[0] - 0.051 <= t <= c_tm[1] + 0.051, (p, t)


def test_builtin_engine_constraints_and_failures():
    cfg = load_config(); b = cfg["pcr"]["builtin"]; n_ok = 0
    for seed in range(12):
        t = _random_seq(420, seed); res = pcr.design_builtin(t, 210, 500, b)
        if res["status"] != "ok":
            continue
        n_ok += 1
        _check_pair(res, t, 210, (b["len_min"], b["len_max"]), (b["tm_lo"], b["tm_hi"]))
        for p in (res["fwd"], res["rev"]):
            assert b["gc_lo"] <= pcr.gc(p) <= b["gc_hi"] and pcr.max_run(p) <= b["max_homopolymer"]
        # 3' ends lie gap_min..gap_max from the breakpoint
        assert b["gap_min"] <= 210 - (res["fpos"] + res["len_F"]) <= b["gap_max"] and b["gap_min"] <= (res["rpos"] - res["len_R"] + 1) - 210 <= b["gap_max"]
        assert res["band_undeleted"] == res["band_deleted"] + 500
    assert n_ok >= 6
    assert pcr.design_builtin("ACGT" * 5, 10, 100)["status"] == "no_consensus"
    assert pcr.design_builtin(_random_seq(200, 1), 20, 100)["status"] == "breakpoint_too_close_to_end"
    assert pcr.design_builtin("A" * 300, 150, 100)["status"] == "no_primer_in_window"


@needs_primer3
def test_primer3_engine_constraints():
    cfg = load_config(); c = cfg["pcr"]["primer3"]; n_ok = 0
    for seed in range(8):
        t = _random_seq(200, 100 + seed); res = pcr.design_primer3(t, 95, 700, c, cfg["pcr"]["tm"])
        if res["status"] != "ok":
            continue
        n_ok += 1
        _check_pair(res, t, 95, (c["min_size"], c["max_size"]), (c["min_tm"], c["max_tm"]))
        assert c["product_min"] <= res["band_deleted"] <= c["product_max"] and res["band_undeleted"] == res["band_deleted"] + 700
        assert res["engine"].startswith("Primer3") and min(res["GC_F"], res["GC_R"]) >= c["min_gc"] - 0.05
        th = pcr.thermo(res["fwd"], res["rev"])
        assert set(th) == set(pcr.THERMO_COLS)
    assert n_ok >= 5
    assert pcr.design_primer3("A" * 200, 100, 500, c)["status"] == "no_primer_pair"


def test_engine_falls_back_without_primer3(monkeypatch):
    cfg = load_config()
    assert cfg["pcr"]["engine"] == "primer3"                                       # the original procedure is the default
    monkeypatch.setattr(pcr, "HAVE_PRIMER3", False)
    eng, msg = pcr.resolve_engine(cfg)
    assert eng == "builtin" and "primer3-py is not installed" in msg and "pip install broodlink[primers]" in msg
    assert pcr.resolve_engine(cfg, "builtin") == ("builtin", "")
    assert pcr.design_primer3(_random_seq(200, 1), 100, 500, cfg["pcr"]["primer3"])["status"] == "primer3_not_installed"
    assert pcr.thermo("ACGTTGCATGCAAGTCCATG", "TTGCATGCAAGTCCATGACC") == {}


def test_fasta_reader_with_and_without_index(tmp_path):
    seqs = {"a": _random_seq(205, 1), "b": _random_seq(61, 2), "c": _random_seq(60, 3)}
    p = tmp_path / "x.fa"
    with open(p, "w", newline="\n") as fh:
        for k, s in seqs.items():
            fh.write(f">{k} description\n" + "\n".join(s[i:i + 50] for i in range(0, len(s), 50)) + "\n")
    fa = pcr.Fasta(str(p))
    assert fa.indexed_from == "scan of the file" and not os.path.exists(str(p) + ".fai")      # nothing written next to the input
    off = 0; lines = []
    for k, s in seqs.items():
        off += len(f">{k} description\n"); lines.append(f"{k}\t{len(s)}\t{off}\t50\t51"); off += len(s) + (len(s) + 49) // 50
    with open(str(p) + ".fai", "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    fb = pcr.Fasta(str(p))
    assert fb.indexed_from == ".fai"
    for f in (fa, fb):
        for k, s in seqs.items():
            assert f.length(k) == len(s) and f.chromosome(k) == s
            for a, b in ((0, 10), (45, 55), (49, 51), (100, 205), (190, 400), (-5, 3)):
                assert f.fetch(k, a, b) == s[max(a, 0):b], (k, a, b)
        assert "a" in f and "zz" not in f
        f.close()
    t, bp = pcr.junction_from_reference(pcr.Fasta(str(p)), "a", 100, 140, 40)
    assert t == seqs["a"][60:100] + seqs["a"][140:180] and bp == 40
    with pytest.raises(BaseException):
        pcr.Fasta(str(tmp_path / "x.fa.gz"))


def test_figure_paging_rules():
    from broodlink.figures import pcr as fp
    per, ncol, nrow = fp.gel_layout(9)[:3]
    assert (per, ncol, nrow) == (9, 3, 3)
    assert int(np.ceil(35 / per)) == 4                       # the 35 pairs of the original cohort: four consecutive pages instead of two
    assert fp.gel_layout(18)[:3] == (18, 6, 3)
    assert fp.gel_layout(500)[0] <= 36                      # never more rows than fit an A4 page
    assert fp.assay_pages(4) == [4] and fp.assay_pages(35) == [35]
    pg = fp.assay_pages(300)
    assert sum(pg) == 300 and len(pg) > 1 and pg[0] < pg[1]


# ------------------------------------------------------------------ the step on the example cohort
def test_step_runs_on_the_example(full_run):
    s = summary(full_run, "10_pcr_markers")
    assert s["status"][0] == "yes" and s["summary"]["ran"]
    S = read(full_run, "10_pcr_markers/pcr_candidate_deletions.tsv"); P = read(full_run, "10_pcr_markers/pcr_primers.tsv")
    I = read(full_run, "10_pcr_markers/pcr_markers_per_interval.tsv"); O = read(full_run, "10_pcr_markers/pcr_primer_order_list.tsv")
    cfg = load_config(); sel = cfg["pcr"]["selection"]
    assert len(I) == 4 and len(P) >= 3 and len(O) == 2 * len(P) and int(I.n_primer_pairs.sum()) == len(P)
    assert (I[I.n_primer_pairs == 0].reason_no_marker.fillna("") != "").all()
    assert list(O.name[:2]) == ["M01_F", "M01_R"] and set(O.sequence) == set(P.forward_primer) | set(P.reverse_primer)
    # screening rules
    d = S[S.diagnostic]
    assert (d.frequency_difference >= sel["min_frequency_difference"] - 5e-4).all() and (S.n_called_class_0 >= sel["min_called_per_class"]).all()
    assert (S[S.gel_resolvable].length_bp <= sel["max_deletion_bp"]).all() and (S.screen_band_undeleted_bp - S.screen_band_deleted_bp == S.length_bp).all()
    assert S.selected.sum() == len(P) and (S[S.selected].diagnostic & S[S.selected].gel_resolvable).all()
    for _iv, g in S.groupby("interval", sort=False):
        assert (np.diff(g.frequency_difference.to_numpy()) <= 1e-9).all()        # ranked by decreasing frequency difference
    # every pair: constraints of the engine that ran, primers flank the deletion, product sizes
    eng = s["summary"]["engine"]
    lim = ((cfg["pcr"]["primer3"]["min_size"], cfg["pcr"]["primer3"]["max_size"], cfg["pcr"]["primer3"]["min_tm"], cfg["pcr"]["primer3"]["max_tm"]) if eng == "primer3" else
           (cfg["pcr"]["builtin"]["len_min"], cfg["pcr"]["builtin"]["len_max"], cfg["pcr"]["builtin"]["tm_lo"], cfg["pcr"]["builtin"]["tm_hi"]))
    for r in P.itertuples():
        t = r.design_sequence; fs = t.find(r.forward_primer); rs = t.find(pcr.revcomp(r.reverse_primer))
        assert fs >= 0 and rs >= 0 and fs + len(r.forward_primer) <= r.breakpoint_index <= rs, r.marker
        assert r.product_with_deletion_bp == rs + len(r.reverse_primer) - fs and r.product_without_deletion_bp == r.product_with_deletion_bp + r.length_bp
        assert r.product_without_deletion_bp - r.product_with_deletion_bp == r.end - r.start == r.product_difference_bp
        for p_, tm_, gc_ in ((r.forward_primer, r.tm_forward, r.gc_forward_percent), (r.reverse_primer, r.tm_reverse, r.gc_reverse_percent)):
            assert lim[0] <= len(p_) <= lim[1] and lim[2] - 0.051 <= tm_ <= lim[3] + 0.051 and abs(pcr.gc(p_) - gc_) < 0.06
            assert abs(pcr.tm_formula(p_) - tm_) < 0.051                          # the Tm of the table is the stated model
        assert r.check_primers_exact_in_design_sequence and r.check_primers_flank_breakpoint and r.check_products_recomputed and r.check_frequency_difference
        assert r.design_source == "split-read consensus (SV VCF)"
        assert r.arrangement_with_deletion == ("class 2" if r.deletion_frequency_class_2 > r.deletion_frequency_class_0 else "class 0")
    # the combined table carries the primers of each interval
    m = I.merge(P, on="interval", suffixes=("", "_p"))
    assert (m.forward_primer == m.forward_primer_p).all() and (m.tm_forward == m.tm_forward_p).all() and (m.product_without_deletion_bp == m.product_without_deletion_bp_p).all()
    # reports and tables
    rep = open(os.path.join(full_run, "run_report.md"), encoding="utf-8").read()
    assert "## PCR markers (in-silico design; no laboratory test)" in rep and "split-read consensus (SV VCF)" in rep and "9. PCR markers" in rep
    assert "No primer was tested in the laboratory" in open(os.path.join(full_run, "methods_text.md"), encoding="utf-8").read()
    head = open(os.path.join(full_run, "10_pcr_markers", "pcr_primers.tsv"), encoding="utf-8").read(3000)
    assert "Tm: nearest neighbour, SantaLucia (1998)" in head and "50 mM monovalent" in head
    for f in ("T13_pcr_candidate_deletions", "T14_pcr_primers", "T15_pcr_markers_per_interval", "T16_pcr_primer_order_list"):
        assert os.path.getsize(os.path.join(full_run, "tables", f + ".tsv")) > 0
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(full_run, "tables", "broodlink_tables.xlsx"), read_only=True)
    assert {"T14 PCR primers", "T15 PCR per interval", "T16 primer order list"} <= set(wb.sheetnames)
    defs = pd.read_excel(os.path.join(full_run, "tables", "broodlink_tables.xlsx"), sheet_name="definitions")
    d14 = defs[defs.table == "T14_pcr_primers"]
    assert len(d14) == P.shape[1] and not (d14.definition == "see the description of the table").any()


def test_pcr_figures_drawn_and_in_the_ledger(full_run):
    M = read(full_run, "figures/MANIFEST.tsv").fillna("")
    row = lambda name: M[M.file.str.contains(name)]
    assert row("fig16_pcr_assay").status.iloc[0].startswith("drawn") and row("fig18_pcr_primer_schematic").status.iloc[0].startswith("drawn")
    g = row("fig17_pcr_simulated_gels")
    assert len(g) >= 1 and g.status.str.startswith("drawn").all() and list(g.file) == [f"fig17_pcr_simulated_gels_page{i + 1}of{len(g)}" for i in range(len(g))]
    assert all(f"page {i + 1} of {len(g)}" in n for i, n in enumerate(g.note))
    assert "no laboratory test" in g.description.iloc[0] and "not gel images" in row("fig16_pcr_assay").description.iloc[0]
    LG = read(full_run, "figures/LEDGER.tsv"); P = read(full_run, "10_pcr_markers/pcr_primers.tsv")
    a = LG[LG.figure == "fig16_pcr_assay"]
    assert len(a) >= 3 * len(P) and (LG.figure == "fig17_pcr_simulated_gels").sum() >= 2 * len(P)
    assert a.value.astype(str).str.contains(P.forward_primer.iloc[0]).any() or a.value.astype(str).str.contains(P.forward_primer[P.arrangement_with_deletion == "class 2"].iloc[0]).any()


@pytest.fixture(scope="module")
def vcf_without_consensus(tmp_path_factory):
    p = str(tmp_path_factory.mktemp("nocons") / "sv_no_consensus.vcf.gz")
    import re
    with gzip.open(os.path.join(DATA, "sv.vcf.gz"), "rt") as fi, gzip.open(p, "wt", newline="\n") as fo:
        for line in fi:
            fo.write(re.sub(r";CONSENSUS=[ACGT]+;CONSBP=\d+", "", line))
    return p


def test_skipped_without_sequence(tmp_path, quick_cfg, vcf_without_consensus):
    out = str(tmp_path / "out")
    main(["run", "--relatedness", os.path.join(DATA, "relatedness.tsv"), "--gt", vcf_without_consensus, "--haplotypes", os.path.join(DATA, "haplotypes.tsv"),
          "--config", quick_cfg, "--out", out, "--no-figures", "--no-inversion-scan", "--quiet"])
    s = summary(out, "10_pcr_markers")
    assert s["status"][0] == "no" and "no sequence to design primers on" in s["status"][3] and "--reference" in s["status"][3]
    assert not os.path.exists(os.path.join(out, "10_pcr_markers", "pcr_primers.tsv"))
    rep = open(os.path.join(out, "run_report.md"), encoding="utf-8").read()
    assert "9. PCR markers" in rep and "no sequence to design primers on" in rep


def test_skipped_without_structural_haplotypes(run_relatedness_only):
    s = summary(run_relatedness_only, "10_pcr_markers")
    assert s["status"][0] == "no" and "no structural-haplotype genotypes" in s["status"][3]


def test_builtin_engine_rerun_alone(full_run, tmp_path):
    """`broodlink pcr-markers` re-runs the step alone; --primer-engine builtin uses the fallback routine."""
    out = str(tmp_path / "copy"); shutil.copytree(full_run, out, ignore=shutil.ignore_patterns("figures"))
    main(["pcr-markers", *FULL, "--out", out, "--no-figures", "--primer-engine", "builtin", "--quiet"])
    s = summary(out, "10_pcr_markers"); b = load_config()["pcr"]["builtin"]
    assert s["summary"]["engine"] == "builtin" and s["status"][2] == "provisional"
    P = read(out, "10_pcr_markers/pcr_primers.tsv")
    assert len(P) >= 1 and (P.engine == "built-in routine").all()
    for r in P.itertuples():
        for p_, tm_ in ((r.forward_primer, r.tm_forward), (r.reverse_primer, r.tm_reverse)):
            assert b["len_min"] <= len(p_) <= b["len_max"] and b["tm_lo"] - 0.051 <= tm_ <= b["tm_hi"] + 0.051 and b["gc_lo"] <= pcr.gc(p_) <= b["gc_hi"]
        assert r.check_primers_exact_in_design_sequence and r.check_primers_flank_breakpoint and r.check_products_recomputed
    assert "only step 10_pcr_markers was re-run" in open(os.path.join(out, "run_report.md"), encoding="utf-8").read()


@needs_primer3
def test_primer3_engine_used_by_default(full_run):
    s = summary(full_run, "10_pcr_markers")
    assert s["summary"]["engine"] == "primer3" and s["summary"]["primer3_installed"] and s["status"][2] == "benchmarked"
    P = read(full_run, "10_pcr_markers/pcr_primers.tsv")
    assert P.engine.str.startswith("Primer3 (primer3-py").all() and P.hairpin_forward_dg_kcal.notna().any() or P.cross_dimer_tm.notna().any()


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    sys.path.insert(0, DATA)
    import simulate_pedigree as sim
    return sim.write_reference(str(tmp_path_factory.mktemp("ref")))


def test_reference_flanks_and_in_silico_pcr(full_run, reference, tmp_path, vcf_without_consensus):
    """--reference: the consensus of the simulated VCF is the reference junction; without a consensus the pairs are
    designed on reference flanks; every pair gives its long product at the locus in the reference."""
    from broodlink.io import read_consensus
    fa = pcr.Fasta(reference)
    assert fa.indexed_from == ".fai" and len(fa.index) == 8
    cons = read_consensus(os.path.join(DATA, "sv.vcf.gz"))
    assert len(cons) > 500
    V = pd.read_csv(os.path.join(DATA, "sv.vcf.gz"), sep="\t", comment="#", header=None, usecols=[0, 1, 2, 7])
    end = {(r[0], r[1], r[2]): int(r[3].split("END=")[1].split(";")[0]) for r in V.itertuples(index=False)}
    for k in list(cons)[:40]:
        seq, bp = cons[k]
        t, b2 = pcr.junction_from_reference(fa, k[0], k[1], end[k], 150)
        assert t[b2 - bp:b2] == seq[:bp] and t[b2:b2 + len(seq) - bp] == seq[bp:], k
    fa.close()
    out = str(tmp_path / "copy"); shutil.copytree(full_run, out, ignore=shutil.ignore_patterns("figures"))
    args = [a if a != os.path.join(DATA, "sv.vcf.gz") else vcf_without_consensus for a in FULL]
    main(["pcr-markers", *args, "--reference", reference, "--out", out, "--no-figures", "--quiet"])
    s = summary(out, "10_pcr_markers")["summary"]
    P = read(out, "10_pcr_markers/pcr_primers.tsv")
    assert s["reference_given"] and s["in_silico_ran"] and len(P) >= 3
    assert P.design_source.str.startswith("reference flanks (150 bp").all()
    assert (P.in_silico_products_at_locus >= 1).all() and (P.in_silico_product_minus_design_bp == 0).all() and P.check_in_silico_product_within_5bp.all()
    assert (P.in_silico_products_off_target == 0).all()                           # random sequence: no second site
    assert s["n_in_silico_product_at_locus"] == len(P)
    rep = open(os.path.join(out, "run_report.md"), encoding="utf-8").read()
    assert "reference flanks (150 bp on each side)" in rep and "In-silico PCR against --reference" in rep
