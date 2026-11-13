"""Fast unit tests of the building blocks (no pipeline run)."""
import gzip
import os

import numpy as np
import pandas as pd
import pytest
import yaml

from conftest import DATA, ROOT
from broodlink import audit, config, families, io, mask, network, nomination, relatedness
from broodlink.cli import build_parser


def test_config_yaml_equals_defaults():
    with open(os.path.join(ROOT, "config.yaml"), encoding="utf-8") as fh:
        assert yaml.safe_load(fh) == config.DEFAULTS
    cfg = config.load_config(os.path.join(ROOT, "example", "config.example.yaml"))
    assert cfg["nomination"]["n_permutations"] == 5000 and cfg["families"]["k2_min"] == 0.15


def test_defaults_keep_the_original_rules():
    c = config.DEFAULTS
    assert c["nomination"]["theta_max"] is None            # original: second degree or closer, no upper bound
    assert c["families"]["rule"] == "components"           # connected blocks, not complete cliques
    assert c["families"]["theta_min"] == 0.177 and c["nomination"]["theta_min"] == 0.0884 and c["nomination"]["n_permutations"] == 20000
    assert c["barcode"]["posterior_call"] == 0.90 and c["audit"]["min_family_size"] == 4 and c["inversions"]["scan"] is True


def test_cli_has_run_and_one_subcommand_per_step():
    ap = build_parser()
    sub = [a for a in ap._actions if a.dest == "cmd"][0].choices
    assert {"run", "qc", "relatedness", "families", "nominate", "audit", "network", "haplotypes", "barcode", "family-effect", "concordance", "report", "sweep", "benchmark"} <= set(sub)
    a = ap.parse_args(["run", "--bisnp", "x", "--gt", "a", "--gt", "b", "--mask", "m", "--mask-exclude", "e", "--auto-exclude"])
    assert a.gt == ["a", "b"] and a.auto_exclude and a.out == "broodlink_out"


def test_bh():
    q = nomination.bh([0.01, 0.02, 0.03, 0.5])
    assert np.allclose(q, [0.04, 0.04, 0.04, 0.5])


def _pairs(rows):
    return pd.DataFrame(rows, columns=["id1", "id2", "theta", "k0", "k1", "k2"])


def test_families_are_blocks_with_closure():
    cfg = config.load_config()
    R = _pairs([("a", "b", .25, .25, .5, .25), ("b", "c", .25, .25, .5, .25), ("a", "c", .20, .4, .5, .10),      # a-c fails k2
                ("d", "e", .26, .2, .5, .3), ("a", "d", .0, 1, 0, 0), ("a", "e", .0, 1, 0, 0), ("b", "d", .0, 1, 0, 0),
                ("b", "e", .0, 1, 0, 0), ("c", "d", .0, 1, 0, 0), ("c", "e", .0, 1, 0, 0), ("a", "f", .0, 1, 0, 0)])
    M, F, info = families.reconstruct(R, list("abcdef"), cfg)
    f1 = F.iloc[0]
    assert f1.n == 3 and f1.status == "not closed" and np.isclose(f1.closure, 2 / 3)
    assert dict(zip(M.individual, M.family))["f"] == "unassigned" and dict(zip(M.individual, M.unit))["f"] == "S_f"
    cfg["families"]["rule"] = "clique"
    M2, F2, _ = families.reconstruct(R, list("abcdef"), cfg)
    assert (F2.status == "closed").all() and F2.n.max() == 2
    R["k2"] = np.nan                                           # GRM mode: theta alone
    M3, F3, info3 = families.reconstruct(R, list("abcdef"), config.load_config())
    assert not info3["has_k2"] and "alone" in info3["edge_rule"] and F3.iloc[0].status == "closed"


def test_nomination_theta_max_is_optional():
    cfg = config.load_config(); cfg["nomination"]["n_permutations"] = 200
    ids = [f"i{k}" for k in range(12)]
    unit = {s: ("F1" if k < 4 else "F2" if k < 8 else f"S_{s}") for k, s in enumerate(ids)}
    rows = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            th = 0.25 if unit[a] == unit[b] else (0.20 if {unit[a], unit[b]} == {"F1", "F2"} else 0.0)
            rows.append((a, b, th, np.nan, np.nan, np.nan))
    R = _pairs(rows)
    g = nomination.nominate(R, unit, cfg, seed=1)
    r = g[(g.unit_a == "F1") & (g.unit_b == "F2")].iloc[0]
    assert r.fraction == 1.0 and r.p_perm == g.p_perm.min() and r.n_cross_pairs == 16
    cfg["nomination"]["theta_max"] = 0.177                    # the upper bound removes these first-degree pairs
    g2 = nomination.nominate(R, unit, cfg, seed=1)
    assert g2[(g2.unit_a == "F1") & (g2.unit_b == "F2")].iloc[0].fraction == 0.0
    assert g.equals(nomination.nominate(R, unit, config.load_config() | {"nomination": {**cfg["nomination"], "theta_max": None}}, seed=1))


def test_split_markers_never_random():
    ms = io.MarkerSet("x", "SNP", "bisnp", "x", ["a"], np.array(["c1"] * 4 + ["c2"] * 4 + ["c3"] * 2), np.arange(10) * 3_000_000, np.arange(10), np.arange(10).astype(str), np.zeros((10, 1), np.int8))
    a, b, how = nomination.split_markers(ms, config.load_config())
    assert "parity" in how and list(a) == [True] * 4 + [False] * 4 + [True] * 2 and not (a & b).any()
    one = ms.subset_markers(np.arange(4))
    a, b, how = nomination.split_markers(one, config.load_config())
    assert "blocks" in how and a.sum() == 2 and b.sum() == 2


def test_slot_limit():
    E = pd.DataFrame({"unit_a": ["F1", "F1", "F1"], "unit_b": ["F2", "F3", "F4"], "min_n": [4, 4, 4]})
    assert network.slot_violations(E, 4) == ["F1"]            # three partners, none linked to another: needs three parents
    E2 = pd.concat([E, pd.DataFrame({"unit_a": ["F2"], "unit_b": ["F3"], "min_n": [4]})])
    assert network.slot_violations(E2, 4) == []               # F2 and F3 can share one parent of F1


def test_transmission_tensor_and_gate1_direction():
    T = audit.transmission_tensor(0.0)
    assert np.allclose(T.sum(2), 1) and np.isclose(T[1, 1, 1], 0.5) and np.isclose(T[0, 2, 1], 1.0)
    rng = np.random.default_rng(0); n = 600
    af = rng.uniform(0.1, 0.9, n)
    par = lambda: (rng.random(n) < af).astype(int) + (rng.random(n) < af).astype(int)
    kid = lambda p, q: rng.binomial(1, p / 2) + rng.binomial(1, q / 2)
    P, Q1, Q2, Q3, Q4 = par(), par(), par(), par(), par()
    G = np.array([kid(P, Q1) for _ in range(5)] + [kid(P, Q2) for _ in range(5)] + [kid(Q3, Q4) for _ in range(5)] + [par() for _ in range(15)]).T.astype(np.int8)
    s = [f"s{i}" for i in range(30)]
    ms = io.MarkerSet("DEL", "DEL", "gt", "x", s, np.array(["c"] * n), np.arange(n), np.arange(n), np.arange(n).astype(str), G)
    fam = {"A": s[:5], "B": s[5:10], "C": s[10:15]}
    st, _n = audit.gate1_scores(ms, fam, [("A", "B"), ("A", "C")], config.load_config(), 1)
    ab, ac = st.iloc[0], st.iloc[1]
    assert ab.LLR_one > 0 and ab.LLR_both < 0 and ac.LLR_one < 0


def test_grm_separates_sibs():
    cfg = config.load_config()
    ms = io.read_variants(os.path.join(DATA, "sv.vcf.gz"), "gt", cfg)
    assert [m.klass for m in ms] == ["DEL", "DUP"] or {m.klass for m in ms} == {"DEL", "DUP"}
    d = [m for m in ms if m.klass == "DEL"][0]
    R, used = relatedness.relatedness_from_markers(d, cfg)
    T = pd.read_csv(os.path.join(DATA, "truth", "truth_families.tsv"), sep="\t").set_index("individual").family
    same = (R.id1.map(T) == R.id2.map(T)) & (R.id1.map(T) != "unassigned")
    assert used > 500 and R[same].theta.median() > 0.18 and abs(R[~same].theta.median()) < 0.03 and R.k2.isna().all()


def test_mask_rules(tmp_path):
    bed = tmp_path / "keep.bed"; bed.write_text("c1\t0\t100\nc1\t200\t300\nc2\t0\t50\n")
    ex = tmp_path / "drop.bed"; ex.write_text("c1\t50\t60\n")
    m = mask.build_mask([str(bed)], [str(ex)], {"c1": 1000, "c2": 100})
    assert m.span() == 240 and "24.00% of the 1,000-bp" not in m.label() and "21.82%" in m.label()
    assert list(m.contains(np.array(["c1"] * 4), np.array([1, 55, 100, 101]))) == [True, False, True, False]
    assert list(m.covered(np.array(["c1", "c1"]), np.array([1, 91]), np.array([100, 210]))) == [90, 20]
    cfg = config.load_config()
    sv = io.MarkerSet("DEL", "DEL", "gt", "x", ["a"], np.array(["c1", "c1", "c1", "c3"]), np.array([10, 90, 210, 5]), np.array([40, 210, 290, 9]),
                      np.array(list("wxyz")), np.zeros((4, 1), np.int8))
    out, row = mask.apply_mask(sv, m, cfg)
    assert list(out.ids) == ["w", "y"] and row["markers_outside_mask"] == 2 and row["warning"] == ""
    sv2 = sv.subset_markers(np.array([False, True, False, True]))
    _o, row2 = mask.apply_mask(sv2, m, cfg)
    assert "removes 100% of the DEL calls" in row2["warning"]                 # never mask SVs away silently
    fa = tmp_path / "g.fa"; fa.write_text(">c1 x\nACGTacgtNN\nAC\n>c2\nnnAA\n")
    k = mask.read_softmasked_fasta(str(fa))
    assert list(k["c1"][0]) == [0, 10] and list(k["c1"][1]) == [4, 12] and list(k["c2"][0]) == [2]
    assert mask.Mask(None, {"c1": 10}).label().startswith("callable span: unmasked")


def test_vcf_reader_classes_and_gvcf(tmp_path):
    p = tmp_path / "x.vcf"
    p.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
                 "c1\t10\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:5\t1|1:3\n"
                 "c1\t20\t.\tAT\tA\t.\t.\t.\tGT\t0/0\t./.\n"
                 "c1\t25\t.\tA\t<NON_REF>\t.\t.\tEND=30\tGT\t0/0\t0/0\n"
                 "c1\t40\t.\tA\tC,G\t.\t.\t.\tGT\t1/2\t0/1\n"
                 "c1\t50\tdel1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=150\tGT\t0/1\t0/0\n"
                 "c1\t60\tdup1\tN\t<DUP:TANDEM>\t.\t.\tSVTYPE=DUP:TANDEM;END=90\tGT\t1/1\t0/1\n")
    cfg = config.load_config()
    a = {m.klass: m for m in io.read_variants(str(p), "vcf", cfg)}
    assert set(a) == {"SNP", "INDEL", "DEL", "DUP"} and a["DEL"].end[0] == 150 and list(a["SNP"].G[0]) == [1, 2] and a["INDEL"].G[0, 1] == -1
    assert any("multi-allelic" in n for n in a["SNP"].notes) and any("reference-only" in n for n in a["SNP"].notes)
    b = io.read_variants(str(p), "snp_indel", cfg)
    assert {m.klass for m in b} == {"SNP_INDEL", "DEL", "DUP"}
    q = tmp_path / "gl.vcf"
    q.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
                 "c1\t10\t.\tA\tC\t.\t.\t.\tGT:PL\t0/1:40,0,40\t0/0:0,3,30\t./.:0,0,0\n")
    g = io.read_variants(str(q), "bisnp", cfg)[0]
    assert list(g.G[0]) == [1, -1, -1] and g.D[0, 0] == 100 and g.D[0, 2] == 255


def test_beagle_reader_and_sniff():
    assert io.sniff(os.path.join(DATA, "snps.beagle.gz")) == "beagle" and io.sniff(os.path.join(DATA, "sv.vcf.gz")) == "vcf"
    ms = io.read_beagle([os.path.join(DATA, "snps.beagle.gz")], step=3)[0]
    assert ms.n == 2000 and len(ms.samples) == 80 and ms.D is not None and ms.chrom[0] == "chr1" and (ms.G >= -1).all()
    d = ms.dosage()
    assert np.nanmax(d) <= 2.0 and np.isnan(d).any()


def test_relatedness_formats_and_folder_checks(tmp_path):
    res = tmp_path / "a.res"
    res.write_text("a\tb\tida\tidb\tnSites\tJ9\tJ8\tJ7\ttheta\n0\t1\ts1\ts2\t10\t0.25\t0.5\t0.25\t0.25\n0\t2\ts1\ts3\t10\t1\t0\t0\t0\n1\t2\ts2\ts3\t10\t1\t0\t0\t0\n")
    R, info = io.read_relatedness(str(res))
    assert info["format"] == "NGSrelate .res" and info["complete"] and R.k2.iloc[0] == 0.25 and R.k0.iloc[0] == 0.25
    kin = tmp_path / "x.kin0"
    kin.write_text("FID1 ID1 FID2 ID2 NSNP HETHET IBS0 KINSHIP\nf s1 f s2 10 0.1 0.0 0.24\n")
    K, ki = io.read_relatedness(str(kin))
    assert ki["format"] == "KING .kin0" and K.theta.iloc[0] == 0.24 and K.k2.isna().all()
    d = tmp_path / "folder"; d.mkdir()
    (d / "p1.res").write_text("a\tb\tida\tidb\tJ9\tJ8\tJ7\ttheta\n0\t1\ts1\ts2\t0.25\t0.5\t0.25\t0.25\n")
    (d / "p2.res").write_text("a\tb\tida\tidb\tJ9\tJ8\tJ7\ttheta\n0\t1\ts1\ts3\t1\t0\t0\t0\n")
    with pytest.raises(SystemExit) as e:
        io.read_relatedness(str(d))
    assert "complete pair set" in str(e.value)
    (d / "p3.res").write_text("a\tb\tida\tidb\tJ9\tJ8\tJ7\ttheta\n0\t1\ts2\ts3\t1\t0\t0\t0\n")
    R3, i3 = io.read_relatedness(str(d))
    assert i3["n_files"] == 3 and i3["complete"]
    (d / "p4.res").write_text("a\tb\tida\tidb\tJ9\tJ8\tJ7\ttheta\n0\t1\ts2\ts3\t1\t0\t0\t0.2\n")
    with pytest.raises(SystemExit) as e:
        io.read_relatedness(str(d))
    assert "inconsistent" in str(e.value)


def test_haplotype_table_orientations(tmp_path):
    w = pd.read_csv(os.path.join(DATA, "haplotypes.tsv"), sep="\t", dtype=str)
    a = io.read_haplotypes(os.path.join(DATA, "haplotypes.tsv"), set(w.individual))
    assert a.n == 4 and len(a.samples) == 80 and a.chrom[0] == "chr1" and a.end[0] == 2800000
    p = tmp_path / "t.tsv"; w.set_index("individual").T.to_csv(p, sep="\t", index_label="interval")
    b = io.read_haplotypes(str(p), set(w.individual))
    assert (a.G == b.G).all()
    long = w.melt(id_vars="individual", var_name="interval", value_name="genotype")[["interval", "individual", "genotype"]]
    q = tmp_path / "l.tsv"; long.to_csv(q, sep="\t", index=False)
    c = io.read_haplotypes(str(q), set(w.individual))
    assert c.n == 4 and (c.subset_samples(a.samples).G == a.G).all()


def test_simulator_is_deterministic(tmp_path):
    import hashlib
    import subprocess
    import sys
    subprocess.run([sys.executable, os.path.join(DATA, "simulate_pedigree.py"), "--out", str(tmp_path)], check=True, capture_output=True)
    for f in ["relatedness.tsv", "sv.vcf.gz", "snps.beagle.gz", "haplotypes.tsv", "callable.bed", "truth/truth_families.tsv"]:
        h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
        assert h(os.path.join(tmp_path, f)) == h(os.path.join(DATA, f)), f
