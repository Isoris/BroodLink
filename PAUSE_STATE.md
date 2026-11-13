# BroodLink build: pause state (13 November 2026)

Paused on request so the external drive can be unplugged. On resume: read this file, check
that the drive with the sources is present, then continue with "Next three actions".
Remove this file before tagging v0.1.0.

Git on this drive (exFAT) needs `git -c safe.directory=* ...` on every call (do not edit
the global git config). Python: the user's Python 3.12, run with `-P`. `pytest` was
installed with `pip --user` (the only package installed). No PDF converter is available
locally (no pandoc / LibreOffice / reportlab): use python-docx for the Zenodo documents.
Heredocs with several files in one Bash call failed to parse: write files with the Write tool.

## State at the last milestone (13 November 2026, second session)

Run things with `PYTHONPATH=<repo>` or through pytest (tests/conftest.py adds the repo to the
path); the package is NOT pip-installed on this machine.

### Done and run
- 20,000-permutation validation on the original tables: `local/validate_original.py`, results in
  `local/reproduction/validation_20000.json` (git-ignored). All published numbers reproduced:
  51 families / 188; 79 units, 3,081 pairs, 90 nominated (with the original labels AND with
  BroodLink's own labels); 23 families >= 4, 253 pairs, 12 nominated, 11 Gate 1, 1 rejected;
  6 supported, 11 families, 5 components; tier 3: 21 kept links + 1 rejected = the 22 edges of the
  original table, 14 latent parents; barcode of the junction 25 (>= 0.90) / 24 (>= 0.95) of 41;
  eta squared 0.585 / 0.528 / 0.602 at 1,000 permutations.
- Package: config, io, mask, relatedness, families, nomination, audit, network (tiers), barcode,
  family_effect, inversions (debugged, runs), qc, pipeline (numbered folders, headers, README.txt
  and summary.json per step, step re-run from tables), cli (run, one subcommand per step, report,
  sweep, benchmark), report (statistics.tsv, run_report.md with capability table, methods_text.md,
  qc_report.md, tables TSV + XLSX), tracks (GFF3 + BED), figures (14 figures, PNG + PDF, MANIFEST).
- Synthetic cohort: `tests/data/simulate_pedigree.py` + data + truth (about 1.1 MB).
- `example/config.example.yaml`, `example/expected/` (outputs of the README command).
- Tests: `tests/test_units.py` (16 passed), `tests/test_pipeline.py` (end to end; see git log for
  whether it was green at the last commit).

### Not started / to do next, in order
1. `broodlink/sweep.py` and `broodlink/benchmark.py` (the CLI already imports them), docs
   `SWEEP.md`, `BENCHMARK.md`; run both on the synthetic data; original-cohort sweeps and
   single-input benchmark from existing intermediates into `local/` (GRM file: see Sources).
2. `docs/REPRODUCTION.md` (tables a-e), `docs/DISPLAY_ITEMS.md`; add display-item and track
   sections to `docs/EXTRACTED_METHOD.md`.
3. README (author: git config name), CITATION.cff, CHANGELOG, gallery PNGs in `docs/figures/`.
4. Zenodo folder `F:\dataroodlink_zenodo_v0.1.0\` (nothing there yet) and THESIS_INSERT.md.
5. Privacy grep, delete this file, tag v0.1.0, `git archive` tarball + SHA-256.
Not implemented so far (say so in the capability table and the final report unless done):
panel-by-panel ports of the published nomination figure (A-J), families figure and Extended
Data H3-H5 (fig14 is a compact composite only); reproduction table (e) from raw per-interval
dosage; reproduction (d) mask/depth from raw data.

## Pending scope (all coordinator messages, in the agreed priority order)
1. Core path with tests passing on the synthetic data (single-input modes, absent-input
   paths, QC flags, determinism, expected example outputs).
2. Per-step numbered output folders: `00_qc`, `01_relatedness`, `02_families`,
   `03_nomination`, `04_audit`, `05_network`, `06_structural_haplotypes`, `07_barcode`,
   `08_family_effect`, `09_cross_marker_concordance`, `figures`, `tables`, `tracks`;
   `run_report.md`, `statistics.tsv`, `methods_text.md` at the top; a README.txt per folder
   (files, columns, units); stable names; every table with a header comment (version,
   config hash, analysis set, marker classes, "families are relatedness-based groups").
3. Inversions: `--inversions BED`; scan ON by default whenever variants are given
   (`--no-inversion-scan` and config key to turn off; SNPs when available, deletions too,
   deletions alone with SV-only input); with a BED do both and write the overlap table;
   `--haplotypes` overrides calls for the intervals it covers; without a BED the barcode
   uses scanned intervals passing separation + heterozygote excess, labelled "candidate,
   detected automatically, not validated"; scan time printed in the run report; chromosomes
   with too few markers skipped with a note. `--detect-inversions` was replaced by this.
4. Local PCA on deletions: (a) karyotypes from deletions and concordance with SNP calls,
   "from deletions only" when no SNPs; (b) deletion window scan; (c) sub-clusters and
   deletion support of splits; (d) arrangement-tag deletions with "fully private".
5. Figures: QC panels, theta-k2 scatter, family sizes, nomination matrix with audit
   triangle, gate plots, latent-parent incidence, network, PC1 plot per interval; ports of
   the nomination figure A-J (v13 design), the families figure panels, Extended Data H3
   (SV calls co-segregate), H4 (cross-marker kinship concordance), H5 (arrangement-private
   DEL/DUP and the barcode); `figures/MANIFEST.tsv`; PNG 300 dpi + PDF; every colour keyed.
6. `docs/DISPLAY_ITEMS.md`: inventory of every figure panel and table of the implemented
   subsections with status implemented / planned / out of scope (diversity statistics, ROH
   and heterozygosity are out of scope).
7. Tracks: `structural_haplotypes.gff3` + BED, `arrangement_tag_variants.gff3` + BED,
   optional `family_informative_variants.bed`, `tracks/README.txt`; GFF3 1-based closed,
   BED 0-based half-open; a test that parses the GFF3 against the tables.
8. Tables (TSV + XLSX with a definitions sheet), `statistics.tsv`, `run_report.md` with
   the capability table, `methods_text.md`.
9. `docs/REPRODUCTION.md` tables (a) full run, (b) one input at a time, (c) relatedness
   alone, (d) mask and depth, (e) 41 intervals: karyotypes (SNP and deletions alone)
   against the original calls, barcode, sub-cluster support 32 / 21, scan recovery; plus
   the visual check of the regenerated figures against the original renders.
10. Zenodo folder: tarball + SHA-256, method description (thesis excerpt), GRM benchmark
    document, two figures with captions, example, REPRODUCTION.md, `.zenodo.json`, README,
    CHECKLIST.md (privacy grep for the sample-name prefix and drive paths), optional
    41-interval BED under `example/` (coordinates and names only).
11. `THESIS_INSERT.md` (Methods sentence, availability line with its location in the
    thesis, paragraph location, manuscript variant, placeholders).
Anything not finished goes into the capability table and the final report as
"not implemented", with what is missing.

## Sources located (do not search again)
Roots: `FAM` = `F:\data\Gariepinus_families_v1.0`, `HAT` = `F:\data\manuscript_ready\_hatchery`,
`POP` = `F:\data\fClaHyb_population_v1.0`.
- Toolkit: `FAM\_project_docs\code\` (`common.py`, `get_FS.py`, `get_HS_audit_get_SP.py`,
  `get_INV_barcode_SP_get_LP.py`, `refine_LP.py`).
- Canonical rerun: `FAM\_shared_parent_links\rerun_canonical188_2026-10-29\` (`code\sp_lib.py`,
  `run_audit.py`, `make_tables.py`, `lp01_profile.py`; `sp_edge_table_canonical188_rerun.tsv`,
  `postaudit_edges_canonical188_v02.tsv`, `rerun_summary.json`, `nomination_seed_sweep.tsv`,
  `README_RERUN.md`); older tables in `FAM\_shared_parent_links\tables\` (`latent_parents_v01.tsv`,
  `postaudit_network_layout_v01.tsv`); `docs\METHODS_shared_parent.md`.
- Relatedness: `FAM\_relatedness_ngsrelate\from_bundle_relatedness_firstdegree_v01\02_outputs\ngsrelate\catfish_226_relatedness.res`
  (226 individuals, has `ida`/`idb`); unrelated list `FAM\_relatedness_ngsrelate\tables\unrelated_set_v01.tsv`.
- Families: `FAM\_fs_families\tables\Supplementary_Table_FS_membership_216_v01.tsv` (analysis
  set of 216 = its `sample_id` column; `Fnum` gives the original family numbers),
  `family_membership_core_v02.tsv`, `Supplementary_Table_FS_closure_v01.tsv`.
- SV genotypes: `POP\_sv_catalogue_delly\catalogues_226\delly_DEL\catalog_226.DEL.GT_matrix.tsv`
  and `...\delly_DUP\catalog_226.DUP.GT_matrix.tsv` (+ `.vcf.gz`, `.bed`).
- Genotype likelihoods: `POP\_genotype_likelihoods_angsd_biSNP\` (28 beagle files, sample
  order `samples_order_226.tsv`, `*.pos.tsv`); reading every third SNP took about 11 min
  in the original run.
- GRM and benchmark: `FAM\_manuscript_figures_2026-10-28\code\compute_sv_kinship_v1.py`,
  `...\tables\sv_kinship_grm_216.npz` (K_snp_all, K_snp_masked, K_del, K_dup, K_inv: allows
  reproduction table (b) without recomputing), `sv_kinship_validation_summary_v1.tsv`,
  `...\upgrade_batch_v2_2026-10-29\tables\sv_kinship_validation_summary_v2.tsv`,
  `...\code\build_extdata_SV_validation_v2.py`, figure `ExtData_SV_kinship_recovers_families_v2.*`
  and caption; `FAM\_shared_parent_links\code\grm_verify.py`, `relatedness_atlas.py`
  (`cluster_sibships`); `FAM\_fs_families\from_bundle_MS_Inversions_DEL_vs_biSNP_fullsib_v01\README.md`.
- Eta squared: `FAM\_manuscript_figures_2026-10-28\code\compute_association_eta2_v1.py`,
  tables `association_eta2_*_canonical188.tsv`; ledger `HAT\LEDGER_SV_inheritance_2026-10-29.tsv` (V13).
- Inversions: `FAM\_inversions_localPCA_sv_discovery\unrelated78_rerun_2026-10-16\`
  (`inv41_unrelated_rerun.py`, `INV41_unrelated78_rerun_tiers.tsv`, README);
  `...\tables\candidate_position_class.tsv` (41 intervals: source of the BED),
  `...\tables\<LG>_<start>_<end>_karyotype_calls.tsv` (sample, geno_code, het_inside, pc1, pc2),
  `...\tables\INV41_eta2_v01.tsv`, `arrangement_private_DEL_catalogue.tsv`;
  per-interval dosage `*_data.npz` under `...\from_release_MS_inversions_release_v02_2026-07-28\`
  (dosage, pos, fish, het, geno, pc1); `FAM\_inv_vectors_haplotypes\tables\KARYOTYPE_VECTORS.tsv`
  (216 x 41) and `LP01_genotype_profile.tsv`; scans `...\code\snp_genome_scan.py`,
  `scan_large_inv.py`, `step2_inversions.py`, `arrangement_classes.py`;
  `FAM\_relatedness_ngsrelate\code\relatedness_inversion_pipeline.py` (`arrangement_classes`, line 297).
- Sub-clusters: `FAM\_manuscript_figures_2026-10-28\code\compute_subclusters_v1.py`,
  `build_extdata_haplotype_subclusters_v1.py`, tables `subclusters_41_assignments_sep2.0.tsv`.
- Audit note on intervals and labels: `HAT\per_inversion_first_pass\AUDIT_intervals_labels_barcode_2026-10-30.md`.
- Figure scripts: `HAT\Fig_sharedparent_nomination_upgrade_2026-11-12\code\build_fig3_sharedparent_v13.py`
  and `build_suppfig_matrix79_v1.py`, outputs in `...\out\` (v13 PNG/PDF, 79-unit matrix);
  `HAT\figures\main\code_sharedparent\build_fig2_halfsib_v8.py` (slot-limit function);
  `HAT\figures\main\code_relatedness\build_hatchery_fig_A_v6.py` (families figure);
  `HAT\figures\upgrades_2026-11-03\round_ed\ED_Fig_H3`, `ED_Fig_H4`, `ED_Fig_H5` (not opened yet).
- Text: hatchery part `HAT\Hatchery_MS_fClaHyb_v40g_10_NOV_2026_DEPTH_THETAPI_CLEAN.docx`
  (Methods headings identified; Data availability paragraph found); thesis
  `F:\data\thesis_ready\v15_2026-11-10\THESIS_2026-11-10_v15p.docx` (text dumped once to the
  session scratchpad as `thesis.txt`, which may be gone on resume: copy the docx and dump
  again; Chapter 5 sections not located yet). Registries:
  `F:\data\manuscript_ready\_REGISTRY_2026-11-12\REGISTRY_manuscript_items.tsv`, `lanes\LANE_3.tsv`,
  `F:\data\thesis_ready\v15_2026-11-10\REGISTRY_2026-11-12\REGISTRY_thesis_items.tsv`.
- Mask and depth: `POP\_callable_mask_denominator_2026-10\` (ledger, `callable_mask_per_LG.tsv`);
  per-sample depth table = Supplementary Data H57 (cohort coverage-QC table; path not located yet).

