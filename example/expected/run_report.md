# BroodLink run report

BroodLink 0.1.0, 2026-11-14, config hash 64177d3d75, seed 1.

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
- cross-marker concordance: 712 of 5,863 SNPs (12.1%) lie in SV regions (structural-haplotype intervals with a flank of 100,000 bp, and the spans of 971 deletion / duplication calls) and are left out of the masked SNP class
- PCR markers: 4 primer pairs for 4 of 4 intervals, designed with Primer3 on: split-read consensus (SV VCF): 4 pairs. In-silico design only; no laboratory test.

## PCR markers (in-silico design; no laboratory test)

- Engine: **Primer3 (primer3-py 2.4.0), the procedure of the original study**.
- Tm: nearest neighbour, SantaLucia (1998) parameters and salt correction, 50 mM monovalent, 1.5 mM divalent, 0.6 mM dNTP, 50 nM primer (primer3-py 2.4.0 calc_tm).
- 105 deletions screened in 4 intervals; 63 diagnostic; 33 diagnostic and gel-resolvable; 4 primer pairs for 4 intervals.
- Design sequence per pair: split-read consensus (SV VCF): 4 pairs (column `design_source` of `pcr_primers.tsv`).
- In-silico PCR: not run (no --reference).
- Checks: primers exact substrings of the design sequence 4 of 4; flanking the breakpoint 4; product sizes recomputed 4.

| marker | interval | variant_id | arrangement_with_deletion | forward_primer | reverse_primer | tm_forward | tm_reverse | product_with_deletion_bp | product_without_deletion_bp | design_source |
|---|---|---|---|---|---|---|---|---|---|---|
| M01 | chr1:2000000-2800000 | DEL_chr1_2386926 | class 2 | ACCTAGCAGTAGAGTGTGAT | AACCATGATCGTTTGCTGTT | 55 | 56.2 | 92 | 1859 | split-read consensus (SV VCF) |
| M02 | chr2:4500000-5400000 | DEL_chr2_4542004 | class 0 | ACAGGATAGCGTCTCTGAAA | GTGGATCGACTCTCTGATCT | 56 | 55.9 | 88 | 515 | split-read consensus (SV VCF) |
| M03 | chr3:1000000-1900000 | DEL_chr3_1413242 | class 2 | GTCCGAAACCAAGCATTACT | TTTTGCGATTCTCCAGTCAC | 56.1 | 56.4 | 144 | 731 | split-read consensus (SV VCF) |
| M04 | chr5:5000000-5800000 | DEL_chr5_5485597 | class 0 | ACGCGTGTGTTTCATCAATA | GGACGTTAGCTAACTTCTCG | 56.1 | 55.7 | 106 | 1795 | split-read consensus (SV VCF) |

Simulated gels in `figures/` are drawn from the designed product sizes; they are not photographs.

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
| 08_family_effect | eta2_median [SNP matched to DEL] | 0.5839 |
| 08_family_effect | eta2_median [SNP matched to DUP] | 0.5837 |
| 08_family_effect | eta2_median [SNP matched to HAPLOTYPE] | 0.6579 |
| 09_cross_marker_concordance | auc_within_vs_between [SNP (bisnp)] | 0.9999 |
| 09_cross_marker_concordance | auc_within_vs_between [SNP, SV regions masked (bisnp)] | 0.9999 |
| 09_cross_marker_concordance | auc_within_vs_between [DEL (gt:sv.vcf.gz)] | 0.9988 |
| 09_cross_marker_concordance | auc_within_vs_between [DUP (gt:sv.vcf.gz)] | 0.9904 |
| 09_cross_marker_concordance | auc_within_vs_between [HAPLOTYPE (structural haplotypes)] | 0.8021 |
| 10_pcr_markers | pcr_primer_engine | primer3 |
| 10_pcr_markers | pcr_deletions_screened | 105 |
| 10_pcr_markers | pcr_diagnostic_gel_resolvable | 33 |
| 10_pcr_markers | pcr_primer_pairs | 4 |
| 10_pcr_markers | pcr_intervals_with_marker | 4 of 4 |

Definitions of every number: `statistics.tsv`.

## Run time

- read inputs: 1.4 s
- 00_qc: 0.6 s
- 01_relatedness: 0.2 s
- 02_families: 0.1 s
- 03_nomination: 1.4 s
- 04_audit: 0.4 s
- 05_network: 1.5 s
- haplotypes: SNP scan: 3.0 s
- haplotypes: deletion scan: 0.1 s
- 06_structural_haplotypes: 17.8 s
- 07_barcode: 2.8 s
- 08_family_effect: 12.2 s
- 09_cross_marker_concordance: 3.0 s
- 10_pcr_markers: 0.3 s
- total: 42.4 s

## Figures

| file | title | status |
|---|---|---|
| fig00_overview | Overview: kinship per marker class and the network of families | drawn (PNG + PDF) |
| fig01_qc_panels | QC panels | drawn (PNG + PDF) |
| fig02_structural_haplotypes_pc1 | Structural haplotypes: PC1 against heterozygosity per interval | drawn (PNG + PDF) |
| fig03_structural_haplotype_scan | Scan for large structural haplotypes | drawn (PNG + PDF) |
| fig04_arrangement_tag_variants | Arrangement-tag deletions and duplications | drawn (PNG + PDF) |
| fig05_barcode_all_parents | Barcode of the shared parents | drawn (PNG + PDF) |
| fig06_shared_parent_nomination_audit_network | Shared-parent links: nomination, audit gates and network | drawn (PNG + PDF) |
| fig07_nomination_matrix_all_units | Family-pair nomination matrix, all units | drawn (PNG + PDF) |
| fig08_sv_genotypes_and_kinship | SV genotypes follow the families and kinship recovers them | drawn (PNG + PDF) |
| fig09_cross_marker_concordance | Cross-marker concordance of kinship | drawn (PNG + PDF) |
| fig10_relatedness_and_families | Relatedness and the full-sib families | drawn (PNG + PDF) |
| fig11_barcode_deletion_augmentation | Barcode of a shared parent with DEL augmentation | drawn (PNG + PDF) |
| fig12_barcode_gene_level | Barcode at gene level | drawn (PNG + PDF) |
| fig13_funnel_to_unsampled_parent | Funnel from offspring relationships to an unsampled parent | drawn (PNG + PDF) |
| fig14_deletion_genotypes_linked_families_chain1_F7 | Deletion genotypes along linked families: F10 / F7 / F12 | drawn (PNG + PDF) |
| fig15_deletion_genotypes_other_links | Deletion genotypes, the other links | drawn (PNG + PDF) |
| fig16_pcr_assay | PCR markers of the structural haplotypes: assay, expected bands, all primer pairs | drawn (PNG + PDF) |
| fig17_pcr_simulated_gels_page1of1 | Simulated gels of the PCR markers, page 1 of 1 | drawn (PNG + PDF) |
| fig18_pcr_primer_schematic | Schematic of the flanking-pair deletion assay | drawn (PNG + PDF) |

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
| 7. Family effect (eta squared) | yes | SNP (bisnp); DEL (gt:sv.vcf.gz); DUP (gt:sv.vcf.gz); HAPLOTYPE (structural haplotypes); SNP matched to DEL; SNP matched to DUP; SNP matched to HAPLOTYPE | benchmarked | once per marker class; individuals assigned to a family only |
| 8. Cross-marker concordance | yes | SNP (bisnp); SNP, SV regions masked (bisnp); DEL (gt:sv.vcf.gz); DUP (gt:sv.vcf.gz); HAPLOTYPE (structural haplotypes) | benchmarked | cross-marker concordance against the reconstructed families, not pedigree truth |
| 9. PCR markers of the structural haplotypes (in-silico design) | yes | DEL (gt:sv.vcf.gz) | benchmarked | Primer3 with the settings recovered from the published primer table of the original study (33 of its 35 pairs returned exactly); screening rules of the published screening table; in-silico design only, no laboratory test |

