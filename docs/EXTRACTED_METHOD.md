# Extracted method: where every BroodLink step comes from

Written before any BroodLink code. It records, for each step, the original script and
function, every threshold with its source line, every null model, and every place where the
original code is tied to one cohort. Paths are given relative to two project roots that are
not part of this repository:

- `FAM/` = the family-analysis project folder (`Gariepinus_families_v1.0`)
- `HAT/` = the hatchery manuscript folder (`manuscript_ready/_hatchery`)
- `POP/` = the population-genetics project folder (`fClaHyb_population_v1.0`)

"Methods" = the Methods section of the hatchery manuscript part file, version v40g
(10 November 2026). Line numbers refer to the files as they were on 13 November 2026.

Status words used below: **extracted** (the code was found and is reused), **reimplemented
from Methods** (no code was found; the step is written from the Methods text), **new in
BroodLink** (a rule the original analysis did not have; it is never presented as original).

---

## 0. Shared constants and machinery

Source: `FAM/_project_docs/code/common.py` (toolkit "hatchpg").

| Constant | Value | Line | Meaning |
|---|---|---|---|
| `THETA_DUP` | 0.354 | 11 | duplicate / monozygotic |
| `THETA_FD` | 0.177 | 12 | first degree |
| `THETA_HS` | 0.0884 | 13 | second degree |
| `THETA_3RD` | 0.0442 | 14 | third degree |
| `K2_FS` | 0.15 | 15 | separates full-sib from parent-offspring |
| `GT_ERR` | 0.03 | 16 | per-allele genotype error of the Mendelian models |

- `classify(theta, k2)` (lines 84-89): duplicate/MZ, full-sib (theta >= 0.177 and k2 >= 0.15),
  parent-offspring (theta >= 0.177 and k2 < 0.15), half-sib (theta >= 0.0884), third degree,
  unrelated.
- `transmission_tensor(err)` (lines 98-110): `T[gP, gQ, g_obs]`, Mendelian transmission
  folded with a genotype-error matrix.
- `hwe_prior(af)` (lines 113-115).
- `eta_squared(G, labels, min_called=20)` (lines 118-131): loop version of the family effect.
- `load_gt_matrix` (lines 43-66): SV genotype matrix `CHROM POS END [ID SVTYPE]` + one column
  per sample, values `0/0`, `0/1`, `1/1`; missing = -1; allele-frequency filter is strict
  (`af > lo`, `af < hi`).
- `load_karyotypes` (lines 69-81): individuals x intervals matrix of 0/1/2, interval ids
  `chrom:start-end`.

## 1. Pairwise relatedness

- **Original**: external tool NGSrelate on ANGSD genotype likelihoods (Beagle format) with
  population allele frequencies; output `.res` with columns `a b ida idb nSites J9 J8 J7 ...
  theta ... KING ...`. `J9 = k0`, `J8 = k1`, `J7 = k2` (`get_FS.py` line 32;
  `FAM/_fs_families/code/pipeline_sibship.py` line 53). Methods, "Pairwise relatedness and
  full-sib family reconstruction".
- **Project's own likelihood implementation**: none. `FAM/_sv_inheritance/code/
  masked_relatedness_GL_v01.py` states (lines 9-18) that "ngsRelate's Jacquard EM is not
  reproduced here"; that script computes KING-robust kinship and IBS0 from hard calls and was
  used only as a masked-versus-unmasked contrast, never to call families. BroodLink therefore
  calls NGSrelate when it is installed and otherwise falls to the GRM (section 2). It does
  not contain a likelihood-based k0/k1/k2 estimator. The NGSrelate call (`relatedness.try_ngsrelate`)
  is made by `broodlink run` when the binary is on the PATH, but that branch was NEVER RUN: the
  tool is not installed on the development machine. Without the binary, `--bisnp` alone uses
  the GRM of hard-called genotypes and the run says so.

## 2. GRM-based relatedness and its benchmark

Two variants exist in the project.

**(a) VanRaden GRM with per-pair denominators** — the benchmarked one, used by BroodLink.
`FAM/_manuscript_figures_2026-10-28/code/compute_sv_kinship_v1.py`, function `grm_acc`
(lines 41-48): `Z = (g - 2p) / sqrt(2p(1-p))`, missing -> 0,
`GRM_ij = sum_m Z_im Z_jm / (number of markers called in both)`.
Filters: allele frequency 0.05-0.95 (line 44, strict); called in at least 150 of 216
individuals for SNPs (line 41, `nmin=150`) and at least 100 of 216 for SVs (line 62,
`nmin=100`); SNP hard call = most likely genotype when its likelihood is > 0.90 (line 56);
every third SNP (`STEP`, run with 3).
Benchmark (`.../tables/sv_kinship_validation_summary_v1.tsv`, `.../upgrade_batch_v2_2026-10-29/
tables/sv_kinship_validation_summary_v2.tsv`, code `build_extdata_SV_validation_v2.py` lines
30-58): rank AUC of within-family (363) against between-family (17,215) pairs of assigned
individuals; recovery = share of within-family pairs above the 99.9th percentile of the
between-family pairs of the same class; family-level bootstrap (1,000 resamples).

| Marker class | Markers | AUC [95% CI] | Recovery % [95% CI] |
|---|---|---|---|
| biSNPs | 1,307,168 | 0.99996 [0.99989, 1.00000] | 99.4 [96.3, 100] |
| biSNPs, SV regions masked | 985,007 | 0.99990 [0.99975, 0.99998] | 97.2 [90.8, 100] |
| DELs | 17,276 | 0.99982 [0.99961, 0.99997] | 95.9 [90.6, 99.7] |
| DUPs | 1,178 | 0.99708 [0.99516, 0.99856] | 74.1 [63.2, 84.0] |
| polymorphic INVs | 41 | 0.98058 [0.97395, 0.98590] | 29.8 [17.8, 48.4] |

The families of the benchmark were built from biSNP relatedness, so this is cross-marker
concordance, not pedigree truth (stated in the benchmark caption).

**(b) kinship phi = GRM / 2 with the theta thresholds** —
`FAM/_shared_parent_links/code/grm_verify.py` line 11 (`grm = Z.T@Z / nv / 2.0`) and line 23
(classes at 0.354 / 0.177 / 0.0884 / 0.0442); `relatedness_atlas.py`
`cluster_sibships(relate, samples, fs_lo=0.177, fs_hi=0.45)` (line 678): sibships = connected
components of the graph of pairs with 0.177 <= phi < 0.45. On the DELLY deletion GRM this gave
48 clusters against 52 biSNP families on 226 individuals
(`FAM/_fs_families/from_bundle_MS_Inversions_DEL_vs_biSNP_fullsib_v01/README.md`).

**What the benchmark did not establish.** No GRM-specific full-sib threshold and no minimum
marker count were fixed by any benchmark file. The only threshold applied to a GRM in the
project is variant (b): phi = GRM / 2 >= 0.177. BroodLink uses estimator (a) with the filters
above, converts to phi = GRM / 2 and applies the theta thresholds of section 0, as in (b).
The marker-count check quotes the benchmark rows above (17,276 markers: recovery 95.9%;
1,178 markers: 74.1%); these are observations, not a derived minimum.

## 3. Full-sib families

- **Code**: `FAM/_project_docs/code/get_FS.py` (toolkit) and
  `FAM/_fs_families/code/pipeline_sibship.py` (lines 33-34, 69).
- **Rule** (`get_FS.py` lines 50-70): full-sib edge = theta >= 0.177 and k2 >= 0.15
  (line 50); blocks = connected components of the full-sib-edge graph (lines 56-70);
  core = maximal complete clique of full-sib edges in the block (lines 75-83); attached =
  at least two full-sib edges into the core and no first-degree edge outside the block
  (lines 85-92); a block is closed when the core is the whole block (line 93).
- **Checked on the original relatedness table**: connected components of the graph
  {theta >= 0.177 and k2 >= 0.15} over the 216 retained individuals give exactly the 51
  canonical families (188 individuals; 332 full-sib edges; 28 unassigned). With k2 >= 0.12
  (an older value in `FAM/_fs_families/code/sibling_build.py` line 7 and in two README files)
  the result is 51 components holding 191 individuals, three of which differ. The canonical
  tables therefore correspond to k2 >= 0.15.
- **Wording versus rule.** The Methods say a block is accepted "when its members were
  pairwise compatible on both gates" and report closure separately (332 of 363 within-family
  pairs; 38 closed families, 12 unresolved, 1 incompatible). The 51 families are the blocks
  (components of the full-sib-edge graph), not complete cliques: 13 families have closure
  below 1. BroodLink reproduces the blocks and reports closure and the membership class of
  every individual (core, attached, unresolved member); a config switch offers cliques only.
- Closure classes: closed (closure 1.0), unresolved, incompatible. The rule that separates
  "unresolved" from "incompatible" was **not found** in code; BroodLink reports closed /
  not closed and the closure value only.
- Parent-offspring check: `get_FS.py` lines 43-48 (count of pairs classified
  parent-offspring; warn if any). Methods filter: k1 >= 0.8 and k2 <= 0.05.
- Unrelated subset: maximum independent set of the theta >= 0.177 graph, 78 individuals
  (Methods; exact integer programming; the script was not found). BroodLink uses a
  deterministic greedy independent set (lowest degree first), which can be smaller than the
  exact maximum.

## 4. Half-sib nomination

- **Code (extracted)**: `FAM/_shared_parent_links/rerun_canonical188_2026-10-29/code/
  sp_lib.py`, function `pair_matrix` (lines 31-72). This is itself a reimplementation made
  on 29 October 2026 from `FAM/_shared_parent_links/docs/METHODS_shared_parent.md` and the
  toolkit script `get_HS_audit_get_SP.py`, validated against the frozen tables (the
  original stage scripts were not archived; Methods, "Software and reimplementation").
  An older two-gate script exists at `FAM/_inversions_localPCA_sv_discovery/code/
  sharedparent_two_gate_DELLY.py`.
- **Units**: families plus each unassigned individual as a one-individual unit
  (`run_nomination.py` line 7; unit name `S_<id>`).
- **Statistic**: fraction of cross-unit pairs with theta >= 0.0884 (`sp_lib.py` line 40;
  `THETA_HS`). There is **no upper bound**: the Methods say "second degree or closer". The
  task description of BroodLink gives 0.0884 <= theta < 0.177; the code and the published
  numbers use theta >= 0.0884, so that is the default, and an optional upper bound is a
  config key.
- **Null**: unit labels permuted among individuals, which keeps unit sizes (lines 55-60);
  exceedance counted with a 1e-12 tolerance (line 60); P = (exceedances + 1) / (permutations
  + 1) (line 61); 20,000 permutations and seed 1 in the primary run; Benjamini-Hochberg over
  all unit pairs (line 70, `bh` lines 18-23); nominated at q < 0.05.
- **Order dependence**: individuals are sorted by id (line 36) and units by name (line 42),
  and the permutation stream depends on that order. The total number of nominated pairs
  depends on the number of permutations and on the seed (83-85 at 1,500 permutations, 90 at
  20,000, seed 1); the set between families of at least four individuals did not change in
  20 runs (`nomination_seed_sweep.tsv`).
- Eligibility tiers by the smaller unit: >= 4 primary audit, 3 extension, 2 provisional,
  1 not testable (`run_audit.py` lines 49-50; Methods "Eligibility").

## 5. Audit

- **Code (extracted)**: `sp_lib.py` `oh_llr_vec` (lines 90-112), `af_filter` (lines 84-87),
  `run_audit.py` `gate_stats` (lines 20-36) and `run` (lines 39-89); `make_tables.py`
  `build` (lines 6-38).
- **Markers**: DELLY deletion and duplication genotypes, not used for nomination; allele
  frequency strictly inside 0.02-0.98, computed over the analysed individuals
  (`af_filter`; run with `af_lo = 0.02`).
- **Gate 1** (`oh_llr_vec`): per marker, with Hardy-Weinberg priors and the transmission
  tensor (error 0.03), `L_F(gP) = sum_gQ prior(gQ) prod_i T[gP, gQ, g_i]`;
  one shared parent = `sum_gP prior(gP) L_A(gP) L_B(gP)`; independent parents =
  `(sum prior L_A)(sum prior L_B)`; both shared = `sum_gP sum_gQ prior prior C_A C_B`.
  `LLR_one` and `LLR_both` = sum of the per-marker log ratios divided by the number of
  markers (line 111: markers with a zero likelihood contribute 0 and stay in the
  denominator). Gate 1 passes when `LLR_one > 0` and `LLR_both < 0` (`run_audit.py`
  line 46). Scored for every family pair with both families >= 2 individuals; the audit
  uses pairs with both >= 4.
- **Gate 2** (`gate_stats`, inner function `sh`, lines 30-33): sharing = loci with a carrier in both families /
  loci with a carrier in either, restricted to loci with at least one called genotype in each
  family. Null: for each gate-1-compatible nominated pair, the 40 non-nominated auditable
  pairs nearest in `LLR_one` (`k_match = 40`, line 39; lines 62-65); pooled one-sided
  P = (1 + number of null values >= observed median) / (1 + number of null values)
  (line 71). Run once per class (DEL, DUP).
- **Support rule** (`make_tables.py` lines 12-19): a link is supported when its DEL or DUP
  sharing exceeds the median of the gate-1-compatible nominated pairs with both families
  >= 4 individuals. The bar is relative (about half pass by construction); the alternative
  (above the 99th percentile of the matched null) is reported as a sensitivity column
  (`run_audit.py` lines 67-68, 80-82).
- **Link classes** (`make_tables.py` lines 14-19, 25-30): "DEL + DUP support", "DEL support",
  "DUP support" (supported); "passed OH, no SV support" (gate 1 only); "not testable
  (family < 4 fish)" (provisional); "failed OH gate" (rejected).

## 6. Network, latent parents, slot limit

- **Code (extracted)**: `make_tables.py` lines 35-38 and 57-62; `HAT/figures/main/
  code_sharedparent/build_fig2_halfsib_v8.py` lines 146-164 (function `viol`).
- Post-audit edges: nominated family pairs that did not fail gate 1 (accepted = audited
  pairs that passed gate 1, provisional = not testable; rejected kept in the table).
- Supported network: links with SV support between families >= 4 individuals; families and
  connected components counted on that graph (`make_tables.py` lines 33-37).
- Latent parents: maximal cliques of the graph of non-rejected nominated family pairs
  (`make_tables.py` lines 58-60). They are bookkeeping labels ("H"), not reconstructed
  parents (Methods, "Terminology").
- Slot limit (`viol`): a family has two parental slots, so its link partners must split
  into at most two groups, each of mutually linked families. A family violates the limit
  when the graph "partners that are not linked to each other" is not two-colourable.
  Evaluated at minimum family size 4, 3 and 2.
- Layout file `postaudit_network_layout_v01.tsv`: stored coordinates of one cohort; not
  reused. BroodLink computes a seeded layout.

## 7. Barcode (shared-parent genotype at large structural haplotypes)

- **Code (extracted)**: `FAM/_project_docs/code/get_INV_barcode_SP_get_LP.py` (per link,
  lines 45-82) and `FAM/_shared_parent_links/rerun_canonical188_2026-10-29/code/
  lp01_profile.py` (joint profile of one parent over several families, lines 16-31).
- `Pr(gP | D_A, D_B) ∝ prior(gP) L_A(gP) L_B(gP)`, `L_F` as in gate 1; priors from the
  allele frequency of the arrangement among the individuals in families; error 0.03; a
  state is called at posterior >= 0.90 (`--post`, line 27; `lp01_profile.py` line 31).
  Private parents `Q_A`, `Q_B` conditional on the best `gP` (lines 65-75).
- Junction families (on two or more supported links): `LP_candidate_junctions` (lines
  91-102). The recurrence test of one parent across two links (`refine_LP.py`) exists in
  the toolkit but is **not part of BroodLink 0.1.0**; BroodLink reports the joint profile of
  a junction under the "same parent on every link" assumption and labels it a candidate.

## 8. Family effect (eta squared)

- **Code (extracted)**: `FAM/_manuscript_figures_2026-10-28/code/compute_association_eta2_v1.py`,
  `eta_vec` (lines 21-29) and `assoc` (lines 33-40); older loop version
  `FAM/_sv_inheritance/code/inheritance_DEL_DUP_INV.py`.
- Individuals: only those assigned to a family (line 16-17). Filters: allele frequency
  strictly inside 0.05-0.95, called in at least 20 individuals (line 15).
- Null: family labels permuted 1,000 times, seed 1, genotype matrix fixed (lines 34-36);
  empirical P with +1 (line 37); Benjamini-Hochberg within class (lines 77).
- Reported: median observed eta squared, median of the per-variant null medians, share above
  the null median, above the null 95th percentile, above every permutation, BH count.

## 9. Structural-haplotype (inversion) genotyping by local PCA

- **Genotyping code (extracted)**: `FAM/_inversions_localPCA_sv_discovery/
  unrelated78_rerun_2026-10-16/inv41_unrelated_rerun.py`: drop zero-variance sites,
  standardise each site, SVD, `pc1 = Vt[0] * S[0]` (lines 67-72); one-dimensional k-means,
  k = 3, started at the 0.1 / 0.5 / 0.9 quantiles (`kmeans1d`, lines 28-37); class 0 =
  lowest PC1, 1 = middle, 2 = highest. Input: expected genotype dosage from the genotype
  likelihoods of all SNPs in the interval.
- **Five criteria** (lines 86-93; `README.md` of that folder; Methods): (1) separation
  >= 2 pooled SD between classes 0-1 and 1-2 (`sep`, lines 39-42); (2) the middle class has
  higher within-interval heterozygosity than both others, one-sided Mann-Whitney P < 0.05
  for both (lines 79-81); (3) calls in the unrelated subset agree with the full-cohort calls
  for >= 90% of individuals; (4) >= 5 copies of the minor arrangement in the unrelated
  subset (line 84); (5) family eta squared above its permutation null. All five -> high
  confidence, otherwise candidate. Result on the original cohort: 34 and 7.
- **Order of operations in the original**: the published calls were made on the full cohort
  (216); the unrelated subset (78) was used to recompute the PCA and the calls for the
  tiers. BroodLink keeps both: it fits the clusters on the unrelated subset when one of the
  configured minimum size exists, projects the other individuals onto that axis, and also
  computes the full-cohort calls for criterion 3.
- **Labels**: class 0 / 2 follow the sign of an SVD axis; "STD" and "INV" carry no ancestral
  meaning (`HAT/per_inversion_first_pass/AUDIT_intervals_labels_barcode_2026-10-30.md`,
  section 2). BroodLink orients 0 = the more frequent homokaryotype.
- **Ambiguous individuals**: the original forced a call. The 0.95 posterior flag of the
  fish-level audit (same note, section 5: three-component Gaussian mixture on PC1) is the
  only recorded ambiguity rule; BroodLink sets genotypes below that posterior to missing
  (**new in BroodLink** as a calling rule).
- Criterion 5 column `above_null` of `INV41_eta2_v01.tsv`: the exact comparison (median,
  95th percentile or maximum of the null) was **not found**. BroodLink uses "observed above
  the 95th percentile of the permuted values" and says so.
- Heterozygosity inside the interval (`het_inside` of the karyotype tables): the script that
  computed it was not found. BroodLink uses the share of sites with a heterozygous hard call
  (or dosage between 0.5 and 1.5 when only dosages exist, as in the scan below).

## 10. Automatic scan for large structural haplotypes

- **SNP scan of the published catalogue**: local PCA "in the style of lostruct" (windows of
  100 SNPs, 2 principal components, multidimensional scaling of window distances onto 20
  coordinates, outlier windows at |z(MDS)| >= 3, merging of outlier windows). The scan
  script and its log are **not in the project** (stated in the Methods and in the audit
  note). BroodLink's SNP scan is **reimplemented from Methods**. Not recorded and therefore
  BroodLink defaults: which coordinates are tested (first 2), the merge gap (3 windows),
  the minimum number of outlier windows per region (3).
- **Simpler SNP scans found in code** (`FAM/_inversions_localPCA_sv_discovery/code/
  snp_genome_scan.py`, `scan_large_inv.py`): windows of 2 Mb (step 1 Mb) or 3 Mb, PC1 share
  of variance and the heterozygosity "inverted-U" ratio (heterozygosity of the middle PC1
  tertile over the outer tertiles); candidates at PC1 share > 0.45 and ratio > 2.3 (or 0.40
  and 1.8). BroodLink reports both statistics per candidate interval.
- **Deletion scan (extracted)**: `.../code/step2_inversions.py`: 1-Mb windows, deletions
  with allele frequency 0.1-0.9, at least 5 per window, LD = mean |pairwise correlation| of
  dosages; hit = window at or above the 90th percentile; hits merged across gaps <= 2 Mb;
  a region needs a window at or above the 95th percentile. The sibling-roughness and
  crossover parts of that script are not carried over.
- **Deletion-based arrangement classes (extracted)**: `arrangement_classes` in
  `FAM/_relatedness_ngsrelate/code/relatedness_inversion_pipeline.py` (line 297) and
  `.../code/arrangement_classes.py`: deletions in the interval with allele frequency
  0.1-0.9, at least 8 (5 in the first version); missing imputed to the marker mean;
  standardised; subset co-segregating with the block (mean |r| > 0.3, or the upper half
  when fewer than 5 remain); PC1 by SVD; three classes on PC1 (Gaussian mixture when
  scikit-learn is present, tertiles otherwise). BroodLink uses the same marker selection
  and PC1 and the one-dimensional k-means of section 9 for the classes (scikit-learn is not
  a dependency). It was used in the original project to karyotype candidate regions from
  deletions and to test for missing homozygotes; the published 41 karyotypes are SNP-based.

## 11. Haplotype structure inside arrangements and deletion support

- **Code (extracted)**: `FAM/_manuscript_figures_2026-10-28/code/compute_subclusters_v1.py`
  (`gmm_fit`, `cluster_class`, `merge`; lines 12-55) and
  `build_extdata_haplotype_subclusters_v1.py` (lines 24-33).
- Sub-clusters: per interval and karyotype class with >= 12 individuals, Gaussian mixture on
  PC1 and PC2 (full covariance, k = 1-4, 12 k-means++ restarts, every component >= 5
  individuals, k by BIC); sub-clusters closer than 2.0 within-cluster SD along the line
  joining them are merged; numbered by PC1.
- Deletion support: deletions inside the interval with allele frequency 0.05-0.95 among the
  class; diagnostic = |difference of mean dosage| >= 1.0 between two sub-clusters; null =
  sub-cluster labels permuted inside the class, 200 permutations (the docstring says 500,
  the executed code 200); P with +1. Original result: 32 of 37 STD splits and 21 of 35 INV
  splits supported; median 11 and 7 diagnostic deletions.

## 12. Arrangement-tag deletions

- **Catalogue found** (`FAM/_inversions_localPCA_sv_discovery/tables/
  arrangement_private_DEL_catalogue.tsv`; columns `arr_tagged`, `carrier_STD`,
  `carrier_INV`, `delta`, `n_STD`, `n_INV`); the script that built it was **not found**.
  The rule is read from the columns and from `HAT/LEDGER_sharedparent_rerun_2026-10-29.tsv`
  row SP15 (called in >= 20 STD/STD and >= 10 INV/INV individuals; tag strength = difference
  of carrier fraction; thresholds 0.6 / 0.8 / 0.9 used downstream). BroodLink computes the
  carrier fraction per homokaryotype class, the difference, the tagged arrangement and a
  "fully private" flag; minimum called individuals per class is a config key.

## 13. Quality control, mask and depth

- **Rules found in the Methods** ("Sample quality control", "Callable mask and sequencing
  depth"); **no QC script was found**.
  - Low coverage: call rate below half of the cohort median (0.46 against a median of 0.91).
  - Samples genotyped twice: pairs with k2 >= 0.95 and theta about 0.49.
  - Contamination: cohort-relative outlier in median KING kinship to all other individuals;
    "no numeric threshold was fixed in advance".
  - Mask: uppercase A/C/G/T bases of the repeat-masked reference (N excluded); 576,895,264 bp
    = 59.85% of 963,905,721 bp; the denominator of depth, density and rate statistics.
  - Depth: median depth over the bins of the callable mask (mask-relative; cohort median
    2.54x) and mean depth over the assembly (median 9.64x); both from tables, the depth tool
    was not recorded.
- BroodLink QC is **reimplemented from Methods**. **New in BroodLink** (not original
  rules, each printed with its value): the numeric contamination flag (median relatedness
  to the cohort above median + 5 MAD of the cohort, or heterozygosity above median + 5 MAD);
  the identical-genotype flag at theta >= 0.354 (the duplicate class of section 0) or
  genotype concordance >= 0.99; the sex / batch tests (Kruskal-Wallis on call rate and
  heterozygosity); the SV mask rule (>= 0.5 of the span and both breakpoints callable) and
  the warning when a mask removes more than half of an SV class.

## 14. Places where the original code is tied to one cohort

| Where | What is hard-coded |
|---|---|
| `sp_lib.py` lines 9, 12-15 | absolute paths to the relatedness table and the DELLY matrices |
| `run_nomination.py`, `run_audit.py` lines 8-17 | paths to the two membership tables; layer names `frozen183`, `canonical188`; unit names built from column `Fnum` |
| `make_tables.py` lines 4-5, 52-62 | output folder; comparison with `postaudit_edges_v01.tsv` and `latent_parents_v01.tsv` |
| `run_audit.py` line 56 | classes fixed to `DEL_share`, `DUP_share` |
| `compute_association_eta2_v1.py` lines 12-17, 53 | paths; `assert len(fish) == 188`; `assert len(order) == 226`; 28 linkage-group file names |
| `compute_sv_kinship_v1.py` lines 18-19, 41, 53, 62 | `assert len(fish) == 216`; `nmin = 150` and `100` (absolute counts for 216 individuals); file name pattern `C_gar_LG%02d`; sample order file of 226 |
| `inv41_unrelated_rerun.py` lines 13-26, 63 | paths; unrelated list of 78; interval names parsed from file names; `'C_gar_' + lg` |
| `lp01_profile.py` line 14 | families `F19`, `F1`, `F3` |
| `compute_subclusters_v1.py`, `build_extdata_haplotype_subclusters_v1.py` | paths; `range(41)`; interval labels A1-A41 |
| `build_fig2_halfsib_v8.py` | stored layout, family names, panel positions |
| Sample names | prefix of the cohort's sample identifiers assumed in several figure scripts; Beagle columns `Ind0..Ind225` mapped through a 226-line order file |
| Absolute counts | `min_called = 20` (eta squared), `>= 12` individuals per class and `>= 5` per component (sub-clusters), `k_match = 40`, `>= 5` minor copies: kept as config defaults in BroodLink; the GRM call-rate filters are converted to fractions (150/216 = 0.694, 100/216 = 0.463) |

## 15. Not found (no rule was invented)

1. The original stage scripts of nomination and gates (before the rerun of 29 October 2026).
2. The genome-wide local-PCA scan script and its run log (section 10).
3. The script that made the full-cohort karyotype calls and `het_inside`.
4. The script that built the arrangement-private deletion catalogue.
5. The exact definition of `above_null` for criterion 5.
6. The rule separating "unresolved" from "incompatible" families.
7. The exact-integer-programming script of the unrelated subset.
8. Any QC script (call rate panel, coverage pass, KING contamination screen).
9. A GRM-specific full-sib threshold or minimum marker count.
10. A k2 estimator for hard genotypes: with a GRM, families rest on the kinship threshold alone.

## 16. Display items, tracks, sweeps (added after the first version of this record)

- **Figure scripts read**: the shared-parent figure (`HAT/Fig_sharedparent_nomination_upgrade_2026-11-12/
  code/build_fig3_sharedparent_v13.py`, panels A-J), the matrix of all units
  (`build_suppfig_matrix79_v1.py`), the families figure (`HAT/figures/main/code_relatedness/
  build_hatchery_fig_A_v6.py`) and the captions of the Extended Data figures on kinship
  concordance and arrangement-private deletions. They hard-code family order, colours, panel
  positions in millimetres and stored network layouts of one cohort. Only the slot-limit
  function was extracted (section 6). BroodLink's figures are written new from its own tables;
  `DISPLAY_ITEMS.md` maps each original panel to a BroodLink figure and says which are not drawn.
- **Tracks**: the original project has no track writer for these results. `broodlink/tracks.py`
  is **new in BroodLink** (GFF3 1-based closed, BED 0-based half-open).
- **Coordinates of intervals**: the original interval ids (`chrom:start-end`) are used as
  given; BroodLink reads `--inversions` as BED (0-based start).
- **Sweeps and benchmark**: the original project has a seed sweep of the nomination
  (`nomination_seed_sweep.tsv`, `sweep_nomination.py`) and a k2 sensitivity inset in the
  families figure. `broodlink/sweep.py` and `broodlink/benchmark.py` are new; the provenance of
  every tested value is in `SWEEP.md`.
