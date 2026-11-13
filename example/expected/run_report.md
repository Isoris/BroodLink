# BroodLink run report

BroodLink 0.1.0, 2026-11-13, config hash afaa86ec32, seed 1.

**Analysis set: 80 individuals** (0 removed by QC, 3 flagged and kept). Families are relatedness-based groups, not pedigree records. callable span: 47,620,925 bp (79.37% of the 60,000,000-bp assembly).

Links between families are candidates. Nomination is screening; a link is 'supported' only after both audit gates, and the support bar of Gate 2 is relative.

## Inputs detected and what each is used for

- **SNP (bisnp)** (`--bisnp`, snps.beagle.gz): class SNP, 5,863 markers, 80 individuals. Used for: structural-haplotype genotyping / scan, family effect, cross-marker concordance, QC.
- **DEL (gt:sv.vcf.gz)** (`--gt`, sv.vcf.gz): class DEL, 768 markers, 80 individuals. Used for: audit Gate 1 and Gate 2, structural-haplotype genotyping / scan, family effect, cross-marker concordance, QC.
- **DUP (gt:sv.vcf.gz)** (`--gt`, sv.vcf.gz): class DUP, 203 markers, 80 individuals. Used for: audit Gate 1 and Gate 2, family effect, cross-marker concordance, QC.
- **Relatedness**: --relatedness, used as is (plain TSV, 80 samples, 3160 pairs).
- **--inversions**: 4 intervals.
- **Mask**: BED callable.bed; callable span: 47,620,925 bp (79.37% of the 60,000,000-bp assembly).

## Notes and warnings

- --bisnp is assumed to be the source of the supplied relatedness and is kept out of the audit

## Headline numbers

| step | statistic | value |
|---|---|---|
| 00_qc | individuals_input | 80 |
| 00_qc | individuals_analysis_set | 80 |
| 00_qc | individuals_removed | 0 |
| 00_qc | individuals_flagged_kept | 3 |
| 00_qc | pairs_identical_genotypes | 1 |
| 00_qc | callable_span_bp | 47620925 |
| 00_qc | assembly_bp | 60000000 |
| 00_qc | median_depth_over_mask | 5.985 |
| 00_qc | median_depth_over_assembly | 10.09 |
| 01_relatedness | relatedness_source | --relatedness, used as is (plain TSV, 80 samples, 3160 pairs) |
| 01_relatedness | relatedness_markers |  |
| 01_relatedness | pairs | 3160 |
| 01_relatedness | parent_offspring_pairs | 0 |
| 01_relatedness | unrelated_subset | 27 |
| 02_families | families | 15 |
| 02_families | individuals_in_families | 72 |
| 02_families | individuals_unassigned | 8 |
| 02_families | units | 23 |
| 02_families | families_closed | 7 |
| 02_families | families_not_closed | 8 |
| 02_families | within_family_pairs_passing | 144 of 159 |
| 03_nomination | unit_pairs | 253 |
| 03_nomination | nominated | 7 |
| 03_nomination | nominated_family_pairs | 7 |
| 04_audit | auditable_pairs | 66 |
| 04_audit | nominated_auditable | 7 |
| 04_audit | gate1_pass | 7 |
| 04_audit | rejected | 0 |
| 04_audit | supported_links | 5 |
| 04_audit | gate2_DEL_median_sharing | 0.7269 |
| 04_audit | gate2_DEL_null_median | 0.6186 |
| 04_audit | gate2_DEL_p | 0.003559 |
| 04_audit | gate2_DUP_median_sharing | 0.6885 |
| 04_audit | gate2_DUP_null_median | 0.5707 |
| 04_audit | gate2_DUP_p | 0.003559 |
| 05_network | links_kept | 7 |
| 05_network | families_on_supported_links | 7 |
| 05_network | components_supported | 3 |
| 05_network | latent_parents | 5 |
| 05_network | slot_limit_violations_min4 | 0 |
| 05_network | slot_limit_violations_min3 | 0 |
| 05_network | slot_limit_violations_min2 | 0 |
| 06_structural_haplotypes | intervals_supplied | 4 |
| 06_structural_haplotypes | intervals_genotyped | 7 |
| 06_structural_haplotypes | intervals_high_confidence | 5 |
| 06_structural_haplotypes | scan_candidates | 3 |
| 06_structural_haplotypes | supplied_recovered_by_scan | 3 |
| 06_structural_haplotypes | arrangement_tag_variants | 138 |
| 07_barcode | barcode_states | 40 |
| 07_barcode | barcode_called | 23 |
| 07_barcode | barcode_called_strict | 20 |
| 08_family_effect | eta2_median [SNP (bisnp)] | 0.5903 |
| 08_family_effect | eta2_median [DEL (gt:sv.vcf.gz)] | 0.561 |
| 08_family_effect | eta2_median [DUP (gt:sv.vcf.gz)] | 0.5821 |
| 08_family_effect | eta2_median [HAPLOTYPE (structural haplotypes)] | 0.6067 |
| 09_cross_marker_concordance | auc_within_vs_between [SNP (bisnp)] | 0.9999 |
| 09_cross_marker_concordance | auc_within_vs_between [DEL (gt:sv.vcf.gz)] | 0.9988 |
| 09_cross_marker_concordance | auc_within_vs_between [DUP (gt:sv.vcf.gz)] | 0.9904 |
| 09_cross_marker_concordance | auc_within_vs_between [HAPLOTYPE (structural haplotypes)] | 0.8021 |

Definitions of every number: `statistics.tsv`.

## Run time

- read inputs: 1.3 s
- 00_qc: 0.8 s
- 01_relatedness: 0.3 s
- 02_families: 0.2 s
- 03_nomination: 1.6 s
- 04_audit: 0.4 s
- 05_network: 0.2 s
- haplotypes: SNP scan: 3.2 s
- haplotypes: deletion scan: 0.2 s
- 06_structural_haplotypes: 20.6 s
- 07_barcode: 0.2 s
- 08_family_effect: 24.0 s
- 09_cross_marker_concordance: 0.8 s
- total: 50.6 s

## Figures

| file | title | status |
|---|---|---|
| fig01_qc_panels | QC panels | drawn (PNG + PDF) |
| fig02_relatedness_theta_k2 | Relatedness: theta and k2 with thresholds | drawn (PNG + PDF) |
| fig03_family_sizes_closure | Family sizes and closure | drawn (PNG + PDF) |
| fig04_nomination_matrix_audit | Nomination matrix with audit triangle | drawn (PNG + PDF) |
| fig05_audit_gates | Audit gates | drawn (PNG + PDF) |
| fig06_latent_parent_incidence | Latent-parent incidence | drawn (PNG + PDF) |
| fig07_network | Network of families and links | drawn (PNG + PDF) |
| fig08_structural_haplotypes_pc1 | Structural haplotypes: PC1 against heterozygosity per interval | drawn (PNG + PDF) |
| fig09_structural_haplotype_scan | Scan for large structural haplotypes | drawn (PNG + PDF) |
| fig10_arrangement_tag_variants | Arrangement-tag deletions and duplications | drawn (PNG + PDF) |
| fig11_barcode | Barcode of the shared parents | drawn (PNG + PDF) |
| fig12_family_effect | Family effect (eta squared) | drawn (PNG + PDF) |
| fig13_cross_marker_concordance | Cross-marker kinship concordance | drawn (PNG + PDF) |
| fig14_shared_parent_overview | Shared-parent overview | drawn (PNG + PDF) |

## Capability table

Confidence: *benchmarked* = the step ran as in the analysis the code was extracted from, on inputs of the kind it was benchmarked on; *provisional* = it ran in a mode that was not benchmarked (single input, too few markers, intervals from the scan, links not audited); *reimplemented from Methods* = no original code existed; *skipped* = not run.

| step | ran | inputs used | confidence | why |
|---|---|---|---|---|
| 0. Quality control | yes | all inputs | reimplemented from Methods | always runs |
| 1. Pairwise relatedness | yes | relatedness table | benchmarked | --relatedness, used as is (plain TSV, 80 samples, 3160 pairs) |
| 2. Full-sib families | yes | relatedness table | benchmarked | theta >= 0.177 and k2 >= 0.15; connected blocks of full-sib edges |
| 3. Half-sib nomination | yes | relatedness table | benchmarked | screening only; nomination is not confirmation |
| 4. Audit | yes | DEL (gt:sv.vcf.gz); DUP (gt:sv.vcf.gz) | benchmarked | Gate 1 and Gate 2; the support bar of Gate 2 is relative (median of the Gate-1-compatible pairs) |
| 5. Network, latent parents, slot limit | yes | relatedness table | benchmarked | links audited |
| 6a. Structural haplotypes (genotyping and scan) | yes | SNP (bisnp), DEL (gt:sv.vcf.gz) | benchmarked | genotyping extracted from the original code; scan: SNP part reimplemented from Methods, candidates not validated |
| 6b. Barcode of the shared parent | yes | supplied interval, genotyped by local PCA | benchmarked | posterior of the shared parent from the descendants' genotypes |
| 7. Family effect (eta squared) | yes | SNP (bisnp); DEL (gt:sv.vcf.gz); DUP (gt:sv.vcf.gz); HAPLOTYPE (structural haplotypes) | benchmarked | once per marker class; individuals assigned to a family only |
| 8. Cross-marker concordance | yes | SNP (bisnp); DEL (gt:sv.vcf.gz); DUP (gt:sv.vcf.gz); HAPLOTYPE (structural haplotypes) | benchmarked | cross-marker concordance against the reconstructed families, not pedigree truth |

