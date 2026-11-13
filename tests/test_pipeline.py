"""End-to-end tests on the synthetic cohort: recovery of the truth (with stated tolerances),
single-input modes, absent-input paths, QC flags, outputs, determinism, expected outputs."""
import itertools
import os

import numpy as np
import pandas as pd
import pytest

from conftest import DATA, FULL, ROOT, read, summary
from broodlink.cli import main

pytestmark = pytest.mark.slow


def partition(M):
    out = {}
    for s, f in zip(M.individual, M.family):
        if f != "unassigned":
            out.setdefault(f, set()).add(s)
    return {frozenset(v) for v in out.values()}


def truth_partition(truth):
    T = truth["families"]
    return {frozenset(g.individual) for f, g in T.groupby("family") if f != "unassigned"}


def pair_scores(A, T):
    pa = {frozenset(p) for v in A for p in itertools.combinations(sorted(v), 2)}
    pt = {frozenset(p) for v in T for p in itertools.combinations(sorted(v), 2)}
    return len(pa & pt) / max(len(pa), 1), len(pa & pt) / max(len(pt), 1)


def family_map(out, truth):
    """called family label -> true family label (majority)."""
    M = read(out, "02_families/family_membership.tsv")
    t = dict(zip(truth["families"].individual, truth["families"].family))
    M = M[M.family != "unassigned"].assign(true=lambda d: d.individual.map(t))
    return M.groupby("family").true.agg(lambda x: x.value_counts().index[0]).to_dict()


def called_links(out, truth, statuses=None):
    L = read(out, "05_network/network_links.tsv")
    if statuses is not None:
        L = L[L.status.isin(statuses)]
    fm = family_map(out, truth)
    return {frozenset((fm[a], fm[b])) for a, b in zip(L.unit_a, L.unit_b)}


def true_links(truth):
    return {frozenset((a, b)) for a, b in zip(truth["links"].family_a, truth["links"].family_b)}


# ------------------------------------------------------------------ recovery of the truth
def test_families_recovered_exactly_from_relatedness(full_run, truth):
    """Tolerance: every true family recovered as the same set of individuals."""
    A = partition(read(full_run, "02_families/family_membership.tsv")); T = truth_partition(truth)
    assert A == T
    s = summary(full_run, "02_families")["summary"]
    assert s["n_families"] == 15 and s["n_unassigned"] == 8
    assert s["n_not_closed"] >= 1, "the synthetic families include blocks that are not complete cliques"


def test_links_nominated_and_audited(full_run, truth):
    """Tolerance: all 7 true links nominated and none rejected; no false link kept;
    at least 3 'supported' (the Gate-2 bar is the median of the compatible pairs, so about half pass)."""
    T = true_links(truth)
    kept = called_links(full_run, truth, ["supported", "Gate-1 only"])
    assert kept == T
    sup = called_links(full_run, truth, ["supported"])
    assert sup <= T and len(sup) >= 3
    assert not called_links(full_run, truth, ["rejected"])


def test_network_and_slot_limit(full_run):
    s = summary(full_run, "05_network")["summary"]
    assert s["n_components"] == 4                      # TF01-02, TF03-04, TF05-06-07, TF08-09-10
    assert all(len(v) == 0 for v in s["slot_violations"].values())
    LP = read(full_run, "05_network/latent_parents.tsv")
    assert (LP[LP.min_family_size == 4].n_families == 3).sum() == 1      # the parent with three families


def test_karyotypes_recovered(full_run, truth):
    """Tolerance: at least 95% of called genotypes equal the truth, per supplied interval, up to the
    arbitrary swap of labels 0 and 2."""
    W = read(full_run, "06_structural_haplotypes/structural_haplotype_genotypes.tsv").set_index("individual")
    K = truth["karyotypes"]
    assert set(K.columns) <= set(W.columns)
    for c in K.columns:
        g = W[c].reindex(K.index); ok = g >= 0
        assert ok.mean() >= 0.9
        acc = max((g[ok] == K[c][ok]).mean(), (g[ok] == 2 - K[c][ok]).mean())
        assert acc >= 0.95, (c, acc)
    S = read(full_run, "06_structural_haplotypes/structural_haplotypes.tsv")
    sup = S[S.source.str.startswith("supplied")]
    assert (sup.genotyped_from == "SNPs").all() and (sup.deletion_call_concordance >= 0.9).all()


def test_scan_recovers_some_intervals_and_labels_candidates(full_run):
    s = summary(full_run, "06_structural_haplotypes")["summary"]
    assert s["scan_ran"] and s["n_supplied_recovered"] >= 2          # tolerance: at least 2 of 4 planted intervals
    S = read(full_run, "06_structural_haplotypes/structural_haplotypes.tsv")
    assert (S[S.source == "scan"].label == "candidate, detected automatically, not validated").all()


def test_barcode_matches_true_parent(run_with_haplotypes, truth):
    """Tolerance: at least 90% of the states called at posterior >= 0.90 equal the true genotype of
    the shared parent, and at least half of the link states are called."""
    out = run_with_haplotypes
    B = read(out, "07_barcode/shared_parent_barcode.tsv"); B = B[B.kind == "link"]
    fm = family_map(out, truth)
    parent = {frozenset((a, b)): p for a, b, p in zip(truth["links"].family_a, truth["links"].family_b, truth["links"].shared_parent)}
    code = {"0/0": 0, "0/1": 1, "1/1": 2}
    hit = tot = 0
    for r in B.itertuples():
        a, b = r.families.split(",")
        p = parent[frozenset((fm[a], fm[b]))]
        if r.call != "unresolved":
            tot += 1; hit += int(code[r.call] == int(truth["parents"].loc[p, r.interval]))
    assert tot >= 0.5 * len(B) and hit / tot >= 0.9, (hit, tot, len(B))


def test_family_effect_above_null(full_run):
    T = read(full_run, "08_family_effect/family_effect_summary.tsv")
    # SNP, DEL, DUP, structural haplotypes, and frequency-matched SNPs for the three non-SNP classes
    assert len(T) == 7 and (T.eta2_obs_median > T.null_median_of_medians + 0.1).all()
    assert set(T.marker_class) == {"SNP", "DEL", "DUP", "HAPLOTYPE", "SNP_matched_DEL", "SNP_matched_DUP", "SNP_matched_HAPLOTYPE"}
    n = T.set_index("marker_class").n_tested
    assert n["SNP_matched_DEL"] == n["DEL"] and n["SNP_matched_DUP"] == n["DUP"]


# ------------------------------------------------------------------ single-input modes and absent inputs
def test_relatedness_alone(run_relatedness_only, truth):
    out = run_relatedness_only
    assert partition(read(out, "02_families/family_membership.tsv")) == truth_partition(truth)
    L = read(out, "05_network/network_links.tsv")
    assert (L.link_class == "nominated, not audited").all()
    assert called_links(out, truth) == true_links(truth)
    st = {k: summary(out, f)["status"] for k, f in (("audit", "04_audit"), ("barcode", "07_barcode"), ("hap", "06_structural_haplotypes"), ("eta", "08_family_effect"))}
    assert all(v[0] == "no" and v[2] == "skipped" for v in st.values())
    rep = open(os.path.join(out, "run_report.md"), encoding="utf-8").read()
    assert "nominated, not audited" in rep and "Capability table" in rep


def test_bisnp_alone_runs_to_network(run_bisnp_only, truth):
    """Hard-called low-coverage SNPs through the GRM: families exact; links provisional. The GRM of
    hard calls underestimates half-sib kinship, so the tolerance on links is: at least 4 of 7
    true links nominated and no false link kept."""
    out = run_bisnp_only
    A, T = partition(read(out, "02_families/family_membership.tsv")), truth_partition(truth)
    prec, rec = pair_scores(A, T)
    assert prec >= 0.98 and rec >= 0.95
    kept = called_links(out, truth, ["supported", "Gate-1 only"])
    assert kept <= true_links(truth) and len(kept) >= 4
    L = read(out, "05_network/network_links.tsv")
    assert L.link_class.str.contains("single-input, provisional").any()
    s = summary(out, "01_relatedness")["summary"]
    assert "single_input_split" in s and not s["has_k2"]
    assert summary(out, "04_audit")["status"][2] == "provisional"
    assert os.path.exists(os.path.join(out, "05_network", "summary.json"))


def test_gt_alone_runs_to_network(run_gt_only, truth):
    """One VCF with DEL and DUP: families and nomination on the denser class (DEL), audit on DUP.
    768 deletions are below the smallest benchmarked marker set, and the run says so: kinship is
    noisy enough that families sharing a parent can merge. Tolerance: pair precision >= 0.8,
    recall >= 0.8 (observed 0.84 / 1.0); the step is labelled provisional."""
    out = run_gt_only
    A, T = partition(read(out, "02_families/family_membership.tsv")), truth_partition(truth)
    prec, rec = pair_scores(A, T)
    assert prec >= 0.8 and rec >= 0.8, (prec, rec)
    s = summary(out, "01_relatedness")["summary"]
    assert "BELOW the smallest benchmarked marker set" in s["marker_count_check"]
    assert summary(out, "01_relatedness")["status"][2] == "provisional"
    assert os.path.exists(os.path.join(out, "05_network", "summary.json"))
    assert summary(out, "04_audit")["summary"]["audit_classes"] == ["DUP (gt:sv.vcf.gz)"]


def test_snp_indel_alone_runs_to_network(run_snp_indel_only, truth):
    out = run_snp_indel_only
    A, T = partition(read(out, "02_families/family_membership.tsv")), truth_partition(truth)
    prec, rec = pair_scores(A, T)
    assert prec >= 0.98 and rec >= 0.95
    a = summary(out, "04_audit")["summary"]
    assert a["ran"] and "gate2_note" in a                  # SNPs only: Gate 2 cannot run
    assert os.path.exists(os.path.join(out, "05_network", "summary.json"))


def test_at_least_one_input_required(tmp_path):
    with pytest.raises(SystemExit) as e:
        main(["run", "--out", str(tmp_path / "o")])
    assert "at least one of" in str(e.value)


# ------------------------------------------------------------------ QC
def test_qc_flags_fire(full_run, truth):
    X = read(full_run, "00_qc/qc_exclusion_table.tsv")
    q = truth["qc"]
    low = q[q.what == "low-coverage sample"]["sample"].iloc[0]
    pair = set(q[q.what.str.startswith("identical")]["sample"])
    assert low in set(X[X.reason == "low coverage"]["sample"])
    assert pair <= set(X[X.reason == "pair of samples with identical genotypes"]["sample"])
    assert (X.action == "flagged, kept").all()            # default: flag, do not remove
    ID = read(full_run, "00_qc/qc_pairs_identical_genotypes.tsv")
    assert {ID.sample_1.iloc[0], ID.sample_2.iloc[0]} == pair
    rep = open(os.path.join(full_run, "00_qc", "qc_report.md"), encoding="utf-8").read()
    assert "pairs of samples with identical genotypes" in rep.lower() and "duplicate" not in rep.lower()


def test_auto_exclude_removes_and_states_the_analysis_set(tmp_path, quick_cfg, truth):
    out = str(tmp_path / "o")
    main(["run", "--relatedness", os.path.join(DATA, "relatedness.tsv"), "--bisnp", os.path.join(DATA, "snps.beagle.gz"), "--gt", os.path.join(DATA, "sv.vcf.gz"), "--config", quick_cfg,
          "--out", out, "--no-figures", "--no-inversion-scan", "--auto-exclude"])
    s = summary(out, "00_qc")["summary"]
    assert s["n_removed"] == 2 and s["n_analysis_set"] == 78
    low = truth["qc"][truth["qc"].what == "low-coverage sample"]["sample"].iloc[0]
    assert low in s["removed"]
    head = open(os.path.join(out, "02_families", "family_membership.tsv"), encoding="utf-8").read().split("\n")[:6]
    assert any("analysis set: 78 individuals (2 removed" in h for h in head)
    assert any("not pedigree records" in h for h in head)


def test_unmatched_sample_names_stop_with_names(tmp_path, quick_cfg):
    R = pd.read_csv(os.path.join(DATA, "relatedness.tsv"), sep="\t")
    R["id1"] = "X" + R.id1; R["id2"] = "X" + R.id2
    p = str(tmp_path / "rel.tsv"); R.to_csv(p, sep="\t", index=False)
    with pytest.raises(SystemExit) as e:
        main(["run", "--relatedness", p, "--gt", os.path.join(DATA, "sv.vcf.gz"), "--config", quick_cfg, "--out", str(tmp_path / "o"), "--no-figures"])
    assert "cannot be matched" in str(e.value) and "XIND" in str(e.value)


def test_mask_reported_and_depth_in_two_columns(full_run):
    Mk = read(full_run, "00_qc/qc_mask_report.tsv")
    assert (Mk.markers_inside_mask + Mk.markers_outside_mask == Mk.markers_total).all() and (Mk.markers_outside_mask > 0).all()
    PS = read(full_run, "00_qc/qc_per_sample.tsv")
    assert {"depth_mask", "depth_assembly"} <= set(PS.columns)
    head = open(os.path.join(full_run, "statistics.tsv"), encoding="utf-8").read()
    assert "callable span: 47,620,925 bp (79.37% of the 60,000,000-bp assembly)" in head


# ------------------------------------------------------------------ outputs
def test_outputs_and_figures_exist(full_run):
    for d in ["00_qc", "01_relatedness", "02_families", "03_nomination", "04_audit", "05_network", "06_structural_haplotypes", "07_barcode",
              "08_family_effect", "09_cross_marker_concordance", "10_pcr_markers"]:
        assert os.path.exists(os.path.join(full_run, d, "README.txt")) and os.path.exists(os.path.join(full_run, d, "summary.json")), d
    for f in ["run_report.md", "statistics.tsv", "methods_text.md", "config_used.yaml", "00_qc/qc_report.md", "tables/broodlink_tables.xlsx", "tables/README.txt",
              "tracks/README.txt", "figures/MANIFEST.tsv"]:
        assert os.path.getsize(os.path.join(full_run, f)) > 0, f
    M = read(full_run, "figures/MANIFEST.tsv")
    assert len(M) >= 19 and M.status.str.startswith("drawn").all(), M[~M.status.str.startswith("drawn")]
    for stem in M.file:
        assert os.path.getsize(os.path.join(full_run, "figures", stem + ".png")) > 5000
        assert os.path.getsize(os.path.join(full_run, "figures", stem + ".pdf")) > 2000
    # the overview and every port of a published figure is drawn from the synthetic data
    for name in PORTS:
        assert M.file.str.contains(name).any(), name
    assert M.file.str.startswith("fig00_overview").any()
    # with SNPs every figure that has a biSNP panel shows it: no note about a missing SNP input
    assert not M.note.fillna("").str.contains("no SNP input|no biSNP").any(), M[M.note.fillna("").str.contains("no SNP input|no biSNP")]
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    for stem in M.file:
        with Image.open(os.path.join(full_run, "figures", stem + ".png")) as im:
            dpi = im.info.get("dpi", (0, 0))[0]
            assert dpi >= 449, (stem, dpi)                               # raster output at 450 dpi or more
            assert im.height / dpi * 25.4 <= 297.5, (stem, im.height / dpi * 25.4)   # nothing taller than an A4 page
    LG = read(full_run, "figures/LEDGER.tsv")
    assert len(LG) > 50 and set(LG.columns) == {"figure", "panel", "printed_item", "value", "source"} and LG.source.notna().all()
    assert set(LG.figure) <= set(stems_of(M))
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(full_run, "tables", "broodlink_tables.xlsx"), read_only=True)
    assert "definitions" in wb.sheetnames and len(wb.sheetnames) >= 10
    rep = open(os.path.join(full_run, "run_report.md"), encoding="utf-8").read()
    assert "## Capability table" in rep and "Analysis set: 80 individuals" in rep and "SNP scan" in rep
    meth = open(os.path.join(full_run, "methods_text.md"), encoding="utf-8").read()
    assert "relatedness-based groups, not pedigree records" in meth and "5,000 permutations" in meth


PORTS = ["pcr_assay", "pcr_simulated_gels", "pcr_primer_schematic", "shared_parent_nomination_audit_network", "nomination_matrix_all_units", "sv_genotypes_and_kinship", "cross_marker_concordance", "relatedness_and_families",
         "barcode_deletion_augmentation", "barcode_gene_level", "funnel_to_unsampled_parent", "deletion_genotypes_linked_families", "deletion_genotypes_other_links"]


def stems_of(M):
    from broodlink.figures import stems
    return stems()


def test_grm_outputs_and_bootstrap(full_run):
    """Step 9 writes the GRM of every marker class and the table with bootstrap intervals."""
    T = read(full_run, "09_cross_marker_concordance/cross_marker_concordance.tsv")
    O = read(full_run, "09_cross_marker_concordance/kinship_grm_individuals.tsv")
    assert list(T.marker_class) == ["SNP", "SNP_MASKED", "DEL", "DUP", "HAPLOTYPE"]
    with np.load(os.path.join(full_run, "09_cross_marker_concordance", "kinship_grm.npz")) as z:
        assert set(z.files) == set(T.grm_key)
        for k in z.files:
            K = z[k]
            assert K.shape == (len(O), len(O)) and np.allclose(K, K.T, equal_nan=True)
    assert (T.auc_ci95_low <= T.auc_within_vs_between + 1e-9).all() and (T.auc_within_vs_between <= T.auc_ci95_high + 1e-9).all()
    assert (T.recovery_ci95_low <= T.recovery_ci95_high).all()
    assert T.n_markers[T.marker_class == "SNP_MASKED"].iloc[0] < T.n_markers[T.marker_class == "SNP"].iloc[0]
    assert T.pearson_r_vs_snp_outside_sv[T.marker_class == "DEL"].iloc[0] > 0.5
    # family order: assigned individuals first, grouped by family
    f = O.family.tolist(); first_un = f.index("unassigned")
    assert all(x == "unassigned" for x in f[first_un:]) and len(set(f[:first_un])) == 15


def test_figures_without_snps_and_without_genes(run_no_snps_with_figures):
    """Optional inputs absent: the panels and figures that need them are left out, nothing fails."""
    out = run_no_snps_with_figures
    M = read(out, "figures/MANIFEST.tsv").fillna("")
    assert not M.status.str.startswith("failed").any(), M[M.status.str.startswith("failed")]
    row = lambda name: M[M.file.str.contains(name)].iloc[0]
    assert row("barcode_gene_level").status == "not drawn: no --genes input"
    for name in ("shared_parent_nomination_audit_network", "sv_genotypes_and_kinship", "cross_marker_concordance", "relatedness_and_families", "overview"):
        assert row(name).status.startswith("drawn"), (name, row(name).status)
    assert "biSNP R0 panel omitted" in row("shared_parent_nomination_audit_network").note
    assert "no biSNP kinship matrix" in row("sv_genotypes_and_kinship").note and "matched biSNPs not drawn" in row("sv_genotypes_and_kinship").note
    assert "needs SNPs and SV calls together" in row("cross_marker_concordance").note
    assert "no biSNPs matrix" in row("overview").note
    T = read(out, "09_cross_marker_concordance/cross_marker_concordance.tsv")
    assert "SNP" not in set(T.marker_class) and "SNP_MASKED" not in set(T.marker_class)
    assert not os.path.exists(os.path.join(out, "07_barcode", "barcode_interval_genes.tsv"))


def test_readme_images_exist():
    """The front page cannot ship with a broken image: every image the README shows is in the repository."""
    import re
    txt = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    imgs = re.findall(r"!\[[^\]]*\]\(([^)\s]+)\)", txt) + re.findall(r"<img[^>]+src=\"([^\"]+)\"", txt)
    assert imgs, "the README shows no image"
    for rel in imgs:
        if rel.startswith("http"):
            continue
        p_ = os.path.join(ROOT, rel.replace("/", os.sep))
        assert os.path.exists(p_) and os.path.getsize(p_) > 2000, rel
    assert "docs/figures/front_page.png" in imgs


def test_every_table_has_the_header(full_run):
    n = 0
    for root, _d, files in os.walk(full_run):
        for f in files:
            if f.endswith(".tsv"):
                head = open(os.path.join(root, f), encoding="utf-8").read(1500)
                assert head.startswith("# BroodLink 0.1.0 | config "), f
                assert "Families are relatedness-based groups, not pedigree records." in head, f
                assert "analysis set: 80 individuals" in head and "marker classes:" in head, f
                n += 1
    assert n > 40


def test_tracks_agree_with_tables(full_run):
    S = read(full_run, "06_structural_haplotypes/structural_haplotypes.tsv")
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(full_run, "tracks", "structural_haplotypes.gff3"), encoding="utf-8") if not l.startswith("#")]
    assert len(rows) == len(S)
    bed = [l.rstrip("\n").split("\t") for l in open(os.path.join(full_run, "tracks", "structural_haplotypes.bed"), encoding="utf-8") if not l.startswith("track")]
    for g, b, r in zip(rows, bed, S.itertuples()):
        assert len(g) == 9 and g[0] == r.chrom and g[2] == "structural_haplotype"
        assert int(g[3]) == r.start + 1 and int(g[4]) == r.end          # GFF3: 1-based, closed
        assert int(b[1]) == r.start and int(b[2]) == r.end              # BED: 0-based, half-open
        attr = dict(x.split("=", 1) for x in g[8].split(";"))
        assert attr["ID"] == r.interval and attr["tier"] == r.tier
    T = read(full_run, "06_structural_haplotypes/arrangement_tag_variants.tsv")
    tg = [l.split("\t") for l in open(os.path.join(full_run, "tracks", "arrangement_tag_variants.gff3"), encoding="utf-8") if not l.startswith("#")]
    assert len(tg) == len(T) and all(int(g[3]) == s and int(g[4]) == e for g, s, e in zip(tg, T.start, T.end))
    assert T.fully_private.any()


def _tables(out):
    res = {}
    for root, _d, files in os.walk(out):
        for f in files:
            if f.endswith(".tsv"):
                res[os.path.relpath(os.path.join(root, f), out).replace("\\", "/")] = open(os.path.join(root, f), encoding="utf-8").read()
    return res


def test_identical_results_for_a_fixed_seed(full_run, tmp_path):
    out = str(tmp_path / "again")
    main(["run", *FULL, "--out", out, "--no-figures"])
    a, b = _tables(full_run), _tables(out)
    a.pop("figures/MANIFEST.tsv", None); a.pop("figures/LEDGER.tsv", None)
    assert set(a) == set(b)
    diff = [k for k in a if a[k] != b[k]]
    assert not diff, diff


def test_step_can_be_rerun_alone(full_run, tmp_path):
    import shutil
    out = str(tmp_path / "copy"); shutil.copytree(full_run, out)
    before = open(os.path.join(out, "03_nomination", "unit_pair_nomination.tsv"), encoding="utf-8").read()
    os.remove(os.path.join(out, "03_nomination", "unit_pair_nomination.tsv"))
    main(["nominate", *FULL, "--out", out, "--no-figures"])
    assert open(os.path.join(out, "03_nomination", "unit_pair_nomination.tsv"), encoding="utf-8").read() == before
    assert "only step 03_nomination was re-run" in open(os.path.join(out, "run_report.md"), encoding="utf-8").read()


def test_matches_expected_example_outputs(full_run):
    """example/expected holds the outputs of the README command. Integers must be equal; floats agree to
    1e-6 relative (linear-algebra libraries differ in the last digits)."""
    exp = os.path.join(ROOT, "example", "expected")
    E = read(exp, "statistics.tsv"); G = read(full_run, "statistics.tsv")
    assert list(E.statistic) == list(G.statistic)
    import importlib.util
    have_p3 = importlib.util.find_spec("primer3") is not None   # the expected PCR-marker numbers were made with Primer3 (optional package)
    for s, e, g in zip(E.statistic, E.value, G.value):
        if str(s).startswith("pcr_") and not have_p3:
            continue
        try:
            assert np.isclose(float(e), float(g), rtol=1e-6, atol=1e-9, equal_nan=True), (s, e, g)
        except (TypeError, ValueError):
            assert str(e) == str(g), (s, e, g)
    for f in ["T02_families.tsv", "T03_family_membership.tsv", "T04_links.tsv", "T06_latent_parents.tsv", "T10_barcode.tsv"]:
        a = read(exp, "tables/" + f); b = read(full_run, "tables/" + f)
        assert list(a.columns) == list(b.columns) and len(a) == len(b), f
        for c in a.columns:
            if pd.api.types.is_float_dtype(a[c]):
                assert np.allclose(a[c].astype(float), b[c].astype(float), rtol=1e-5, atol=1e-8, equal_nan=True), (f, c)
            else:
                assert (a[c].fillna("").astype(str) == b[c].fillna("").astype(str)).all(), (f, c)
