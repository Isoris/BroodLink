# Display items: the published figures, drawn by BroodLink

The figures of the original study are drawn by `broodlink run` (and `broodlink report`) from
the tables of the run. Each port keeps the panels, layout, encodings and style of the
original script; family names, counts, thresholds and every printed number come from the
run. `figures/MANIFEST.tsv` lists the figures of a run with a note on any panel left out
because its input is absent; `figures/LEDGER.tsv` lists every number printed in a ported
panel with the table it was computed from.

Status words: **ported** (same panels from the run's tables), **partly** (some panels left
out, reason given), **not ported**.

## Figures of a run

| File | Content | Origin | Status |
|---|---|---|---|
| `fig00_overview` | kinship matrix per marker class (biSNPs, DELs, DUPs, structural haplotypes) and the network of families | the kinship row of the SV kinship figure and the network panel of the shared-parent figure | new composite of two ported panels (front page of the README) |
| `fig01_qc_panels` | call rate, heterozygosity, depth, marker density; heterozygosity against median kinship to the cohort | panel E = the contamination panel of the published cohort figure (kinship of the run instead of KING; analysis set only) | no published counterpart for A-D |
| `fig02_structural_haplotypes_pc1`, `fig03_structural_haplotype_scan`, `fig04_arrangement_tag_variants`, `fig05_barcode_all_parents` | PC1 against heterozygosity per interval; the scan; tag variants; barcode of every link and junction | | no published counterpart |
| `fig06_shared_parent_nomination_audit_network` | A matrix of the auditable families, B nomination score, C Gate 1, D / E Gate 2 per SV class, F both Gate 1 scores, G incidence, H assignment completeness, I opposite homozygotes (biSNP R0, DEL), J network by component | main shared-parent figure A-J (`build_fig3_sharedparent_v13.py`) | ported |
| `fig07_nomination_matrix_all_units` | every family, links to unassigned individuals in the margin | supplementary matrix (`build_suppfig_matrix79_v1.py`) | ported |
| `fig08_sv_genotypes_and_kinship` | A genotypes in family order, B eta squared against allele frequency with matched biSNPs, C eta squared above the null, D-H kinship matrix per marker class | `build_Figure_4.py` | ported |
| `fig09_cross_marker_concordance` | A AUC with bootstrap interval, B kinship within / between families, C SV kinship against masked-biSNP kinship, D pairs recovered; both headlines generated | ED Fig. H4 (`build_ED_Fig_H4.py`) | ported |
| `fig10_relatedness_and_families` | A IBD simplex, B relatedness-class matrix, C family sizes with insets, D within / between kinship, E closure per family | panels A-E of the relatedness figure (`build_hatchery_fig_A_v6.py`) | ported; the ideogram and the population-genetic panels of the combined figure are out of scope (genome displays and statistics outside this pipeline); its SV panel is `fig08` |
| `fig11_barcode_deletion_augmentation` | A tag DELs per interval, B tag strength, C DEL evidence against the barcode state, D sweep, E barcode before and after, F information per link, G junction vote, H evidence per SV class | ED Fig. H5 (`ED_Fig_H5.py`) | partly: panel I (held-out prediction of withheld offspring) is a validation outside the tool |
| `fig12_barcode_gene_level` | A evidence per link and interval, B genes, tag variants and barcode cells of the intervals opened | ED Fig. H6 (`ED_Fig_H6.py`), needs `--genes` | partly: gene-family colours and labels and the candidate centromere mark come from cohort inventories and are not drawn |
| `fig13_funnel_to_unsampled_parent` | A, B closure; C, D nomination funnel; E, F audit; G parental barcode of one link; H, I recurrence across two links and its null | Supp. Fig. S-H4 (`Supp_Fig_S-H4.py`) | ported; two closure classes instead of three; the second (60-draw) null is not drawn |
| `fig14_deletion_genotypes_linked_families_chainN_<family>` | A deletion genotypes of a chain of linked families with the expected cross, B sharing per link and chromosome | Supp. Fig. S-H6 (`build_Supp_Fig_S-H6.py`) | ported, one figure per chain (junction family with its supported partners; at most `figures.max_chain_figures`); the read-support quality bar is not drawn (caller-specific input) |
| `fig15_deletion_genotypes_other_links` | one block per remaining link, per-chromosome summary | Supp. Fig. S-H7 (`Supp_Fig_S-H7.py`) | ported, paged to A4; same omission |

Removed because a port shows the same thing: the old nomination matrix, audit gates,
latent-parent incidence, network, theta against k2, family sizes and closure, family effect,
concordance and composite overview figures of 0.1.0.

## What depends on optional inputs

| Input absent | Left out (the layout closes up; the manifest note says so) |
|---|---|
| SNPs (`--bisnp` / `--snp_indel`) | biSNP and masked-biSNP kinship matrices and distributions, panel C of `fig09` and its headline, matched biSNPs in `fig08` B and C, biSNP R0 in `fig06` I unless the relatedness table has an R0 column |
| SV calls | Gate 2 panels, deletion opposite-homozygote panel, genotype heat maps, `fig11`, `fig14`, `fig15` |
| k0 / k1 / k2 (relatedness from a GRM) | simplex and k2 inset of `fig10` |
| structural-haplotype genotypes or intervals | barcode figures, haplotype rows and matrices |
| `--genes` | `fig12` ("not drawn: no --genes input") |
| audit not run | gate panels; the network is titled "nominated (not audited)" |

## Comparison with the original renders (original cohort)

BroodLink was run on the original cohort from the original tables (the likelihood-based
relatedness table, the deletion and duplication genotype matrices, the karyotypes of the 41
candidate polyINVs as `--haplotypes`, the 41 intervals, the original family labels through
`--families`; seed 1, 20,000 permutations; a second run added the biSNP genotype likelihoods,
every sixth SNP). Each ported figure was opened next to the original render and compared
panel by panel. Nothing was tuned to match. Cohort outputs are not part of this repository.

| Figure | Same as the original | Differs, with the cause |
|---|---|---|
| Shared-parent figure (`fig06`) | A: family order, 253 pairs / 12 nominated / 11 pass Gate 1, support counts 4 / 1 / 1 / 5, the rejected pair. B: 3,081 pairs of 79 units, 90 nominated, null mean 0.023. C: 242 / 11. D, E: 0.779 against 0.701 and 0.663 against 0.532, P = 0.0023. F: 242 / 11 and the rejected pair. G: 14 latent parents over the same 22 families, the same family above two parents, violations 0 / 1 / 4. H: 115 closed, 28 unassigned, 167 core. I: ratios 0.16 / 0.54 / 0.43 / 0.67 (R0) and 0.42 / 0.65 / 0.95 (DEL), with the pair counts of the original tables in every class. J: six components with the same families, sizes and link counts (3 + 7, 3 + 0, 3 + 0, 0 + 3, 1 + 0, 1 + 0), same link colours | H: "not closed 73" instead of unresolved 68 + incompatible 5 (rule not found); attached 17 / unresolved member 4 instead of 16 / 5 (one individual: the membership classes are recomputed, the original read a stored table). G: latent parents are numbered by BroodLink (labels are bookkeeping). J: node positions are computed (the original used stored coordinates), components are numbered by size, family colours are generated. Wording: "unit pairs", "failed Gate 1", "Gate 1 only" |
| Matrix of all units (`fig07`) | family order, every cell, symbol and margin count; 51 + 28 = 79 units, 3,081 comparisons, 90 nominated = 54 + 36, 47 / 4 families, key counts 11 / 4 / 1 / 1 / 5 / 1 / 42 | none found |
| SV genotypes and kinship (`fig08`) | Headline 1 (0.60 / 0.58 / 0.53, null 0.27-0.28). A: 188 fish, 51 families, same family order, 41 intervals, 40 and 30 random SVs. B: medians and deletions per bin (3.4k ... 1.1k), null 0.28. C: excess 0.29 / 0.25 / 0.34. Kinship row: AUC 0.9998 / 0.9971 / 0.9806 with 17,276 / 1,178 / 41 variants; 216 = 188 + 28 | A: the heat maps were compared by eye, not cell by cell (same selection rule and seed). biSNP panels: see the biSNP run below |
| Cross-marker concordance (`fig09`) | A: 0.9998 [0.9996, 1.0000], 0.9971 [0.9952, 0.9986], 0.9806 [0.9739, 0.9859]. B: 363 and 17,215 pairs; medians -0.03 / 0.47, -0.02 / 0.29, -0.02 / 0.47. D: 96% [91, 100], 74% [63, 84], 30% [18, 48] | biSNP rows and panel C: see the biSNP run below |
| Relatedness and families (`fig10`) | A: 20,717 / 1,102 / 1,029 / 13 / 363 pairs. C: 22, 6, 9, 6, 4, 1, 1, 1, 1 families; bands 20,717 / 1,102 / 1,029 / 372 of 23,220; closure curve. D: median 0.238, 4 of 363, 13, 21,819 / 1,038 / 363; 23 families in the inset. E: 38 of 51 closed, 332 of 363, 30 bridged, 1 unbridged, triads 437 / 105 / 16, the same starred families | family colours are generated (the original read a stored colour file); B: the family order comes from the same ordering rule, on BroodLink's table |
| Barcode with deletion evidence (`fig11`) | With the original tag catalogue and links the computation returns the original evidence cell for cell (9 settings x 7 links x 41 intervals), the sweep counts (43 of 45 at k = 2, tag strength >= 0.8) and the junction vote (26 + 9 partial, 4 opposing, 2 without a vote). In the figure: P = 0.002, barcode counts per link, opposing / no vote 4 / 2 | The figure uses the tag variants of step 6 because the original catalogue builder was not found: 1,621 tag DELs instead of 2,080, match 95% (38 of 40) instead of 96% (43 of 45), augmented counts within one interval of the original per link. One 3-fish link that the original entered by hand is absent (links need both families at the minimum size): 6 links instead of 7. The joint profile has 25 called states (the rerun profile) instead of the 24 of an older stored file. Labels: class 0 / class 2 instead of STD / INV, family pairs instead of link names |
| Barcode at gene level (`fig12`) | gene counts of the four intervals (168 / 194 / 105 / 229), all 24 barcode cells of the shared links | obligate counts and tagging DELs per interval are close, not equal, and fewer calls hit coding sequence (headline: 2 of 5 against 12 of 15): same tag-catalogue cause; the original also counted weakly tagging DELs from a cohort table. The four intervals of the original were named in the local config; the automatic rule (most called parent states) opens other intervals |
| Funnel (`fig13`) | Headline 1; A bars; 3,081 / 79 / 90 / 50; D 12 / 11 / 6, 32 / 29 / 0, 36; E thresholds and null point; F 4 / 1 / 1 / 5 / 1, 12 auditable, 11 compatible, P = 0.0023; G 20 / 30 / 8 states; headline 5: 14 of 14, 190 pseudo-chains, P = 0.005; I: 0.320 / 0.094 / 0.181 | B: two closure classes. "supported" 6 instead of "validated" 7 and the 3-fish band 10 / 10 / 0 instead of 10 / 10 / 1 (the hand-entered extension link). H: joint profile 25 called instead of 24 (as above) |
| Deletion genotypes, chain (`fig14`) | every ledger number: 160 of 645 eligible deletions, 22 fish, fish and low-call fish per family, expected-cross counts per family, 157 of 168 squares above the median and 93 above the 95th percentile, row counts; all 168 observed / null values equal the original table | one quality bar instead of two; panel B rows are labelled by family pair and ordered by name |
| Deletion genotypes, other links (`fig15`) | the five original blocks: eligible deletions, "1/1 in both" counts, low-call fish; summary rows 26 / 12, 24 / 15, 28 / 16, 28 / 20 | the port shows every remaining kept link (supported, then Gate 1 only, then provisional) over several pages, so page composition and the headline totals differ |

Summary table of the kinship benchmark (published Supplementary Table S-H3) against
`09_cross_marker_concordance/cross_marker_concordance.tsv`:

| Marker class | Published: variants, AUC [95% CI], recovery % [95% CI] | BroodLink |
|---|---|---|
| DELs | 17,276; 0.99982 [0.99961, 0.99997]; 95.9 [90.6, 99.7] | 17,276; 0.99982 [0.99961, 0.99997]; 95.9 [90.6, 99.7] |
| DUPs | 1,178; 0.99708 [0.99516, 0.99856]; 74.1 [63.2, 84.0] | 1,178; 0.99708 [0.99516, 0.99856]; 74.1 [63.2, 84.0] |
| candidate polyINVs | 41; 0.98058 [0.97395, 0.98590]; 29.8 [17.8, 48.4] | 41; 0.98058 [0.97395, 0.98590]; 29.8 [17.8, 48.4] |
| biSNPs | 1,307,168 (every third SNP); 0.99996 [0.99989, 1.00000]; 99.4 [96.3, 100] | 653,632 (every sixth SNP); 0.99996 [0.99988, 1.00000]; 99.2 [96.3, 100] |
| biSNPs, SV regions masked | 985,007; 0.99990 [0.99975, 0.99998]; 97.2 [90.8, 100] | 492,564; 0.99990 [0.99976, 0.99998]; 97.2 [90.4, 100] |

The three SV rows are equal to the published table to the last printed digit. The biSNP rows come
from the run that read the 28 genotype-likelihood files through `--bisnp` with every sixth SNP
(about 25 minutes; the published matrices used every third), so half as many SNPs enter and the
values agree to three or four digits, not exactly. Given the stored GRMs of the original study, the
same functions return the published AUC, recovery and both intervals to the last digit in all five
rows (`docs/REPRODUCTION.md`).

biSNP panels in that run, against the originals: kinship matrices of biSNPs and masked biSNPs
present beside the DEL / DUP / candidate-polyINV ones (`fig08` D-H, `fig00`); kinship medians
-0.03 / 0.42 for biSNPs as published and for masked biSNPs (published -0.02 / 0.42); Pearson r of SV kinship against masked-biSNP
kinship 0.96 / 0.83 / 0.49 over 23,220 pairs as published, and the same second headline (r = 0.96,
30% recovered); recovery 99% [96, 100] and 97% [90, 100] (published 97% [91, 100]); matched biSNPs
of the family effect: median eta squared 0.61 for the three matched sets (published 0.62 / 0.61 /
0.59) and excess over the null 0.31 / 0.31 / 0.32 (published 0.32 / 0.32 / 0.29): another SNP
subset and, in this local run only, 5,000 deletions and 200 permutations instead of all deletions
and 1,000. The biSNP R0 panel of `fig06` uses the R0 column of the relatedness table in both runs.

## Tables

| Original table | BroodLink |
|---|---|
| Family membership and closure | `T02_families`, `T03_family_membership` |
| Unit-pair nomination and audit edge table | `03_nomination/unit_pair_nomination.tsv`, `04_audit/audit_edges.tsv`, `T04_links`, `T05_audit_gate2` |
| Latent parents, slot limit | `T06_latent_parents`, `T07_slot_limit` |
| Opposite homozygotes by relationship class | `05_network/opposite_homozygotes_by_class.tsv` |
| Interval tiers, karyotype calls | `T08_structural_haplotypes`, `06_structural_haplotypes/*` |
| Arrangement-private deletions | `T09_arrangement_tag_variants` |
| Shared-parent barcode, deletion evidence, recurrence | `T10_barcode`, `07_barcode/barcode_tag_evidence.tsv`, `tag_evidence_sweep.tsv`, `junction_super_barcode.tsv`, `junction_recurrence.tsv` |
| Family effect | `T11_family_effect` |
| Kinship benchmark (S-H3) | `T12_cross_marker_concordance`, `09_cross_marker_concordance/kinship_grm.npz` |
| Sample QC, exclusions | `T01_qc_exclusion_table`, `00_qc/*` |

## Out of scope

Diversity statistics (theta, pi, Tajima's D), runs of homozygosity, heterozygosity
landscapes, population structure and admixture, the ideogram, gene families and functional
annotation of intervals, leave-one-out and held-out validations of the links, PCR marker
design. These belong to the original study, not to this pipeline.
