# Display items: what BroodLink draws and what it does not

Inventory of the figure panels and tables of the subsections of the original study that
BroodLink implements, with the BroodLink output that corresponds to each. Status words:
**implemented** (drawn or written by `broodlink run`), **partly** (the same quantity in a
different form), **planned** (not in 0.1.0), **out of scope**.

BroodLink's figures are generic views drawn from its own tables. They are not ports of the
published figures: layouts, family order, colours and stored network coordinates of the
published figures are specific to one cohort and were not carried over.

## Shared-parent figure of the original study (panels A-J, design v13)

| Panel | Content | BroodLink 0.1.0 | Status |
|---|---|---|---|
| A | Family-pair matrix: nomination below the diagonal, audit above (families of at least 4) | `fig04_nomination_matrix_audit`, `fig14` panel A (all units, or families only) | partly: same two triangles; no size-ordered blocks, no symbol per SV class (the class is in `T04_links`) |
| B | Distribution of the nomination score, all unit pairs and nominated | none (values in `03_nomination/unit_pair_nomination.tsv`) | planned |
| C | Gate 1 score distribution (one shared parent) | `fig05_audit_gates` panel A (scatter of both scores) | partly |
| D, E | Gate 2 sharing of deletions / duplications against the matched null | `fig05_audit_gates` panel B (one column per SV class) | partly: violins of the non-nominated auditable pairs, not the matched-null histogram |
| F | Gate 1, both scores (one parent against both parents) | `fig05_audit_gates` panel A | implemented |
| G | Family x latent-parent incidence with slot-limit violations | `fig06_latent_parent_incidence` | implemented |
| H | Family assignment completeness (block status, membership class) | `fig03_family_sizes_closure` (sizes and closure); counts in `T02`, `T03` | partly |
| I | Opposite-homozygote rates by relationship class | none | planned |
| J | Audited network by component | `fig07_network` (seeded layout) | implemented, different layout |

Supplementary matrix of all units: `fig04_nomination_matrix_audit`. Implemented.

## Families figure of the original study

| Panel | Content | BroodLink 0.1.0 | Status |
|---|---|---|---|
| theta against k2 with thresholds | pairs by relationship class | `fig02_relatedness_theta_k2` | implemented |
| kinship heat map ordered by family | | none | planned |
| family sizes, pairs per theta band | | `fig03_family_sizes_closure`; `01_relatedness/relatedness_class_counts.tsv` | partly |
| closure against the k2 cut | sensitivity inset | `broodlink sweep` (families_k2) | partly: a sweep table, not a panel |
| ideogram of SV density and inversions | | none | out of scope (a genome display, not a relatedness result) |
| family effect per SV class | eta squared against the null | `fig12_family_effect` | implemented |

## Extended Data items

| Item | Content | BroodLink 0.1.0 | Status |
|---|---|---|---|
| SV calls co-segregate with families | eta squared per variant class, share above the null | `fig12_family_effect`, `T11_family_effect` | partly (one panel) |
| Cross-marker kinship concordance | AUC and recovery per marker class; within / between kinship; SV against SNP kinship | `fig13_cross_marker_concordance`, `T12` | partly: no bootstrap intervals, no SV-against-SNP scatter |
| Arrangement-private DEL / DUP and the barcode | tag variants per interval, tag strength, deletion evidence against the barcode, sweep | `fig10_arrangement_tag_variants`, `fig11_barcode`, `T09`, `T10` | partly: the test of deletion evidence against the barcode state is not implemented |
| Haplotype sub-clusters inside arrangements | sub-clusters and deletion support | `06_structural_haplotypes/haplotype_subclusters.tsv` | table only, no figure |

## BroodLink figures without an original counterpart

`fig01_qc_panels`, `fig08_structural_haplotypes_pc1` (PC1 against heterozygosity per
interval), `fig09_structural_haplotype_scan`, `fig14_shared_parent_overview` (a compact
composite of matrix, network, gates, incidence and barcode).

## Tables

| Original table | BroodLink |
|---|---|
| Family membership and closure | `T02_families`, `T03_family_membership` |
| Unit-pair nomination and audit edge table | `03_nomination/unit_pair_nomination.tsv`, `04_audit/audit_edges.tsv`, `T04_links`, `T05_audit_gate2` |
| Latent parents, slot limit | `T06_latent_parents`, `T07_slot_limit` |
| Interval tiers, karyotype calls | `T08_structural_haplotypes`, `06_structural_haplotypes/*` |
| Arrangement-private deletions | `T09_arrangement_tag_variants` |
| Shared-parent barcode | `T10_barcode` |
| Family effect, kinship concordance | `T11_family_effect`, `T12_cross_marker_concordance` |
| Sample QC, exclusions | `T01_qc_exclusion_table`, `00_qc/*` |

## Out of scope

Diversity statistics (theta, pi, Tajima's D), runs of homozygosity, heterozygosity landscapes,
population structure and admixture, gene content of intervals, the recurrence test of one
parent across two links, leave-one-out and held-out validations of the links, PCR marker
design. These belong to the original study, not to this pipeline.
