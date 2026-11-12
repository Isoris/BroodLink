# BroodLink build: pause state (13 November 2026)

Paused on request so the external drive can be unplugged. On resume: read this file, check
that the drive with the sources is present, then continue with "Next three actions".
Remove this file before tagging v0.1.0.

Git on this drive (exFAT) needs `git -c safe.directory=* ...` on every call (do not edit
the global git config). Python: the user's Python 3.12, run with `-P`. `pytest` was
installed with `pip --user` (the only package installed). No PDF converter is available
locally (no pandoc / LibreOffice / reportlab): use python-docx for the Zenodo documents.
Heredocs with several files in one Bash call failed to parse: write files with the Write tool.

## Done
- `docs/EXTRACTED_METHOD.md` (committed): every step with source script, function, line,
  thresholds, nulls, cohort-specific hard-coding, and the "not found" list. Still to add:
  sections for the display items / figure scripts and for the tracks (nothing extracted yet).
- `broodlink/config.py`: all defaults (includes the inversion scan, sub-clusters, tags keys).
- `broodlink/io.py`: MarkerSet; readers for beagle, VCF/gVCF (pandas-free parser), genotype
  matrix TSV, haplotype tables, interval BED, relatedness (.res, folder, KING, TSV), samples,
  FAI, depth; `write_tsv`. Not yet run against real VCF/beagle files (only `read_relatedness`,
  `read_gt_matrix`, `read_haplotypes` were exercised).
- `broodlink/mask.py`: BED / soft-masked FASTA masks, intersect, subtract, SV rule. Untested.
- `broodlink/relatedness.py` (VanRaden GRM, NGSrelate call, pair matrices, unrelated subset),
  `families.py`, `nomination.py`, `audit.py`, `network.py`, `barcode.py`, `family_effect.py`.
- **Core validated on the original tables** (scratch script, 2,000 permutations): 51 families
  holding 188 (identical sets); 3,081 unit pairs; 253 pairs between families >= 4; 12
  nominated, 11 pass Gate 1, 1 rejected; 6 supported links, 11 families, 5 components;
  Gate-2 medians 0.779 (DEL) / 0.663 (DUP) vs null 0.701 / 0.532, P = 0.0023; barcode of the
  junction family (three families) 25 states at >= 0.90 and 24 at >= 0.95; eta squared
  medians 0.585 / 0.528 / 0.602 (200 permutations in that check). Total nominated was 82 at
  2,000 permutations; the 20,000-permutation run (expected 90 with the canonical labels) was
  stopped at the pause and must be rerun.
- `broodlink/inversions.py`: written in full, **never imported or run**: interval genotyping
  (SNP and deletion local PCA, five criteria, ambiguous -> missing), SNP scan (reimplemented
  from Methods), deletion scan, overlap table, tag variants, sub-clusters + deletion support.

## Half-done
- `broodlink/inversions.py`: untested; `tag_variants` contains a leftover no-op loop
  (`for cls, name in ((0, "0"), (2, "2")): pass`) to delete.
- `network.build`: latent parents are computed on all nominated family pairs. The original
  post-audit network (22 edges, 14 cliques) is restricted to families of at least 3
  individuals; add a tier (min family size 4 / 3 / 2) to the latent-parent table and use
  min size 3 for the network figure. The check of "22 edges / 14 cliques at tier 3" was
  running when paused (appended at the end of the scratch script; not confirmed).
- Family labels: BroodLink numbers families by size then first member, which differs from
  the original F numbers, so the Monte Carlo stream of the nomination differs (same n >= 4
  set expected). `families.read_families` exists to import labels (`--families`) for the
  exact reproduction of the 90 nominated pairs.

## Not started
- `broodlink/__init__.py` (empty file exists; needs `__version__ = "0.1.0"`), `cli.py`,
  `pipeline.py`, `qc.py`, `figures.py`, `report.py`, `tracks.py`, `pyproject.toml`,
  `config.yaml`, README, LICENSE, CITATION.cff, CHANGELOG, `.gitattributes` (`* text=auto eol=lf`).
- `tests/data/simulate_pedigree.py`, the synthetic data, all tests, `example/expected/`.
- docs: `REPRODUCTION.md`, `DISPLAY_ITEMS.md`.
- Zenodo folder (nothing written there yet), `THESIS_INSERT.md`.

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

## Next three actions on resume
1. Rerun the scratch validation at 20,000 permutations (families, nomination with the
   original labels through `read_families`, audit, network tiers, barcode, eta squared) and
   record the numbers; fix `network.build` to carry tiers.
2. Write `qc.py`, `pipeline.py` (numbered step folders, headers), `cli.py`, `report.py`,
   `pyproject.toml`, `config.yaml`; write `tests/data/simulate_pedigree.py` (families with
   shared parents, unassigned individuals, low-coverage likelihoods, DEL/DUP, two or three
   inversions with tag deletions, mask, identical pair, low-coverage sample) and get a
   first `broodlink run` on the synthetic data to finish; import and debug `inversions.py`.
3. Write the tests (core path first), then figures, tracks, reproduction tables, the Zenodo
   folder and the thesis insert, in the priority order above.
