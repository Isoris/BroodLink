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
| Absolute counts | `min_called = 20` (eta squared), `>= 12` individuals per class and `>= 5` per component (sub-clusters), `k_match = 40`, `>= 5` minor copies: kept as config defaults in BroodLink; the GRM call-rate filters are converted to fractions (150/216 and 100/216, stored as 0.6944 and 0.4629, rounded down so that a marker called in exactly 150 or 100 individuals is kept, as in the original) |

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

## 17. Ports of the published figures (added with the figure package)

The figures of the original study are now drawn by `broodlink run` from the tables of the run
(`broodlink/figures/`; `docs/DISPLAY_ITEMS.md` has the panel-by-panel status). The house style
(`figures/style.py`: fonts, type scale, palette, placement in millimetres, numbered banners)
is ported from the style module of the original figures (`ros_figstyle.py`, `hat3_style.py`).
Family names, counts, thresholds and every printed number come from the run; the printed
numbers are listed with their source table in `figures/LEDGER.tsv` (the originals wrote one
ledger per figure). Computations added to feed the panels:

| Quantity | Module, step | Original | Status |
|---|---|---|---|
| Opposite homozygotes by relationship class: R0 of the biSNPs, opposite-homozygote rate of the deletions, ratio to unrelated / unlinked pairs | `opposite.py`, step 5 | tables `FAM/_sv_inheritance/tables/opposite_homozygote_test_v01.tsv` (NGSrelate column R0, column `cls`) and `opposite_homozygote_DEL_v01.tsv` (`oh_rate`, `n_shared`, `cls`), described in `FAM/_project_docs/docs/00_FUNNEL_STRATEGY.md`; the script that wrote them was **not found** | reimplemented from the tables and the note. R0 is read from the relatedness table when it has the column; the count from hard-called SNPs (opposite homozygotes / both heterozygous) is **new in BroodLink**. No marker filter is recorded: every deletion of the input is used |
| Frequency-matched SNPs for the family effect | `family_effect.matched_snps`, step 7 | `compute_association_eta2_v1.py`, block "matched SNP pool" (bin width 0.02, seed + 7, null with seed + 1) | extracted |
| GRM per marker class written out; SNPs outside SV regions; bootstrap intervals of AUC and recovery; Pearson correlation with the masked-SNP GRM | `concordance.py`, step 8 | `compute_sv_kinship_v1.py` (`merged`, `in_mask`: intervals + 100 kb, DEL and DUP spans), `build_extdata_SV_validation_v2.py` (`boot_ci`: 1,000 resamples of the families, seed 7) | extracted. The SV spans masked are those of the calls given to the run |
| Closure, status and membership class of imported families (`--families`) | `families.reconstruct(given=...)`, step 2 | same rules as section 3 | extracted (applied to given labels) |
| Deletions that fit a Mendelian cross in every family of a chain, expected-cross strip, thinning, order of individuals | `deletion_paths.py` (`fit_family`, `Context.block`, `order_individuals`, `both_hom`), step 5 | `HAT/code/build_supp_chromosome_DEL_heatmaps_v2.py` (`fit_family`, `block2`, `strips2`, `summarise`), `HAT/code/chromDEL_common.py` (`af_mask`, `order_fish`) | extracted: allele frequency among family individuals inside the audit window; call rate >= 0.90 among well-called individuals (at most 30% blank on the chromosome) of every family; best parental pair of 6 (error 0.02) with simulated P >= 0.05 (300 draws, seed 1, >= 3 called); thinning within 100 kb at >= 90% concordance; even subsample of 160 |
| Sharing per link and chromosome against matched pairs | `deletion_paths.chrom_sharing`, step 5 | `HAT/code/build_supp_chromosome_DEL_heatmaps_v1.py` (`chrom_share_table`), `chromDEL_common.py` (`share`, `matched_null`, `edges`) | extracted. The allele-frequency filter is on the individuals in families (the genome-wide Gate 2 of section 5 uses the analysis set) |
| Chains of linked families, family order in a chain, chromosome shown, paging | `deletion_paths.find_chains`, `choose_chromosome` | typed into the original scripts (one chain, one linkage group "with most DELs") | **new in BroodLink**: a chain = a junction family with its supported partners (`deletion_paths.chain_links: kept` adds families on two or more kept links); chromosome = the one with the most deletions |
| Tag-variant evidence on the shared parent (obligate carriers), sweep over k and tag strength, permutation check, augmented barcode | `tag_evidence.py`, `augmentation.py`, step 7 | `FAM/_manuscript_figures_2026-10-28/code/compute_del_augmentation_v1.py` (`obligate`, `evidence`, `acc`), `HAT/per_inversion_first_pass/relabelled_augmentation/code/compute_augmentation_relabelled_v1.py` (`prep`, `ob`, `evidence`, `acc`) | extracted: both families hold a homozygous carrier, called in >= 2 individuals per family, allele frequency 0.02-0.98 among family individuals; evidence class from >= k obligate tag variants per class at tag strength >= delta (k = 2, delta = 0.8 shown; sweep 1-3 by 0.6 / 0.8 / 0.9); null = evidence classes permuted across intervals inside each link, 500 permutations. Checked cell for cell against the stored original when given the original catalogue and links. The tag catalogue is that of step 6 (section 12: the original catalogue builder was not found) |
| Vote of the links of a junction family ("super barcode") | `augmentation.super_vote` | `ED_Fig_H5.py`, block "G vote logic" (from `HAT/code/build_extdata_SP_DEL_augmentation_v5.py`) | extracted |
| Recurrence of a junction parent across two links, pseudo-chain null | `recurrence.py`, step 7 | `FAM/_project_docs/code/refine_LP.py` | extracted (section 7 listed it as not part of 0.1.0). Partners come from the network; a seeded subset is drawn when the null pairs exceed `barcode.recurrence.n_null_max` |
| Genes of every barcode interval, tag variants on coding sequence | `genes.py`, step 7, needs `--genes` | `ED_Fig_H6.py` (`prepare`, `cell_color`), `HAT/code/del_gene_loci_common.py` (`load_hits`, `calls_table`, `size_label`) | gene membership by overlap extracted; hits on coding sequence computed from the gene file are **new in BroodLink** (the original read a cohort table of hits); the rule choosing the intervals opened is new |

Not ported, with the reason: the held-out prediction of withheld offspring (panel I of the
barcode-augmentation figure, `m_incremental_augmentation.py`: a validation outside the tool);
the second, 60-draw run of the recurrence null; the read-support quality bar of the
deletion-genotype figures (a caller-specific table); gene-family colours and the candidate
centromere mark of the gene-level figure (cohort inventories); the ideogram and the
population-genetic panels of the combined relatedness figure (genome displays and statistics
outside this pipeline); the split of "not closed" families into unresolved and incompatible
(rule not found, section 15 item 6).

## 18. PCR markers of the structural haplotypes (step 10)

`PCR/` = `FAM/_pcr_del_markers` (toolkit `code/`, tables `tables/`, bundle `from_bundle_pcr_del_toolkit_v01/`);
`FUN/` = the functional-divergence manuscript folder (`manuscript_ready/_inv_functional_NLR`: `code/`, `ledger/`, `figures/`).
Published tables: screened deletions (`PCR/tables/Supplementary_Data_S4_...`, 4,771 rows), designed primers
(`...Data_S5_...` = `inversion_primers_ALL.tsv`, 35 pairs, the source of truth according to `FUN/ledger/POLARITY_and_primer_reconciliation_2026-10-29.md`),
markers per interval (`...Table_S6_...`), assay per interval (`...Table_S7_...`, v3 in `FUN/supplementary/tables/`).

### 18.1 What the original design used: Primer3, not the toolkit routine

The script that wrote the published primer table was **not found**. Three scripts exist and none of them is it:
`PCR/code/get_PCR_DEL_markers.py` (screening only, reports primer *windows* of 200 bp, designs nothing),
`PCR/code/deltools.py` with `FAM/_inversions_localPCA_sv_discovery/code/03_design_primers_DELLY.py` (a later port "built from
the MS_Inversions arrangement-PCR pipeline" for family-private and shared-parent deletions, bundle README and
`FAM/_project_docs/tables/audit.txt`), and the check scripts in `FUN/code/`.

Evidence, all re-run for this record:

| Test | Result |
|---|---|
| `deltools.pick_flanking_primers` (default windows) on the split-read consensus of the 35 published deletions, primer3-py absent (its own Tm fallback) | 0 of 35 pairs reproduced (6 designs succeed, 26 `no_primer_in_window`, 3 `breakpoint_too_close_to_end`) |
| the same with primer3-py 2.4.0 installed (Tm from `calc_tm`) | 0 of 35 (25 designs succeed, all different pairs) |
| stored result of `FUN/code/recon_primers.py` (`FUN/ledger/recon_primers_QC.tsv`, column `repo_deltools_reproduces`) | False for 35 of 35; the ledger note says "the repository design tool does not reproduce the primers (0/35)" |
| layout of the published primer table | Primer3's conventions: columns `fpos, flen, rpos, rlen` with `rpos` = last base of the reverse site and product = `rpos - fpos + 1` (Primer3 `PRIMER_LEFT_0`, `PRIMER_RIGHT_0`), and a column `expl` (Primer3 "explain"); `deltools` uses other positions |
| Primer3 (`primer3.design_primers`, primer3-py 2.4.0) on the same consensus, target = the breakpoint, optimum Tm 56, Tm 54-57, minimum GC 30%, product 60-300 bp, everything else default | **33 of 35 pairs identical** (both primers), with identical Tm, band sizes and GC |
| sensitivity of that result | optimum Tm 55.9 gives 21 of 35, 56.1 gives 19, 55.5 gives 3; Primer3 defaults (optimum 60) give 0; without the target 20; minimum GC 20% gives 31; a GC clamp gives 3 |

So the published primers were designed **with Primer3** (the author's recollection is right); the statement in
`recon_primers.py` that the table "reproduces exactly" with the toolkit routine is not borne out by its own output.
The toolkit routine shares only the Tm function with the published table (`primer3.calc_tm` when primer3-py is
installed); its Tm window (54.5-56.6) is the observed range of the published Tm values, not a design setting.

The Primer3 settings are **recovered by reproduction, not read from a script**: optimum Tm 56.0 is sharply
determined; the Tm range (54-57), the minimum GC (30%) and the product range sit on plateaus (several values give
the same 33 pairs). Two published pairs are not returned (one primer of each differs by 2 and by 15 bp of position;
the published primer has a higher Primer3 penalty than the one Primer3 2.4.0 returns): a setting or a Primer3
version that was not recovered. No further setting was tried to force them.

BroodLink therefore uses **Primer3 as the default engine** (`pcr.engine: primer3`, config block `pcr.primer3`) and
the toolkit routine as the fallback when the optional package primer3-py is absent (`pcr.engine: builtin`,
`--primer-engine builtin`; block `pcr.builtin` = the defaults of `pick_flanking_primers`). The fallback is stated
in the run report and does not return the Primer3 pairs.

### 18.2 Marker selection

| Rule | Value | Source | Status |
|---|---|---|---|
| Deletions considered | inside the interval: `POS >= start` and `END <= end` | `get_PCR_DEL_markers.py`, line `inside = ...` | extracted |
| Frequency per class | deletion allele frequency = mean genotype / 2 among the called individuals of the class (homokaryotype 0, heterokaryotype, homokaryotype 2) | same script, function `f` | extracted |
| Called individuals per class | at least 5 in each homokaryotype class | the script's default is `--min-called 10`; the published table holds rows down to 5 and 5 | **from the published table** (smallest value observed); `pcr.selection.min_called_per_class` |
| Diagnostic deletion | frequency difference between the homokaryotype classes >= 0.90 **and** heterokaryotype frequency 0.30-0.70 | published table: column `diagnostic` is true for 307 of 341 rows at difference >= 0.90; the 34 others have a heterokaryotype frequency of 0.05-0.295 or 0.72-0.91 (one at 0.458 with a displayed difference of 0.9, a rounding case), the 307 have 0.300-0.676. The script has `--min-diff 0.80` and a `het_near_midpoint` column (within 0.20 of the midpoint) | **reimplemented from the published table**; the two bounds lie in the gaps 0.295-0.300 and 0.676-0.723 and were set to 0.30 and 0.70 |
| Screening bands | generic pair of 2 x 100 bp flank: 200 bp with the deletion, 200 + length without | published table (`band_deleted_bp` = 200 in every row); the script uses `--flank 200` (400 bp) | from the published table |
| Gel-resolvable | deletion length <= 1,800 bp (long screening band <= 2 kb) | published table: `gel_resolvable` true for 290-1,800 bp, false from 1,803 bp | **reimplemented from the published table**. Lower bound: none can be read (the shortest deletion of that catalogue is 290 bp); BroodLink uses 60 bp, the band-difference rule of `deltools` (`abs(band_undeleted - band_deleted) >= 60`) |
| Ranking | frequency difference, largest first; ties by position | order of the published table (458 of 458 adjacent ties are in position order) | from the published table |
| How many per interval | one | published primer table (35 pairs, one per interval) | `pcr.markers_per_interval` |
| Which candidate | the first diagnostic gel-resolvable deletion, in rank order, for which a pair can be designed | **not recorded**. In 32 of 35 intervals the published deletion is the first such candidate; in 3 it is the next one after a candidate for which the recovered Primer3 settings also return a pair | **reimplemented from Methods** ("one deletion was selected per inversion"); see REPRODUCTION |
| Designed pair, toolkit flag `gel_resolvable` | products differ by >= 60 bp, short product 60-1,200 bp, long product <= 6,000 bp | `deltools.pick_flanking_primers`, last statement | extracted, reported as a check column |

On the original cohort these rules return the published screening table row for row (REPRODUCTION, section (f)).
The arrangement-private catalogue of the same folder (`...Data_S3_...`, carrier frequencies, 2,080 rows, columns
`sizeok`, `fixed`, `pick_score`) is a different, older screen; its builder was not found (section 12) and it is not
what the primers were picked from: step 10 does not use it. Step 10 screens every deletion inside an interval and
marks those that are also step-6 tag variants (`step06_tag`); on the original cohort every diagnostic gel-resolvable deletion is one.

### 18.3 Primer design

| Constraint | Primer3 engine (default) | Built-in routine (fallback) |
|---|---|---|
| Source | recovered by reproduction of the published primer table (18.1) | `deltools.pick_flanking_primers` (signature defaults), unchanged |
| Design sequence | split-read consensus of the junction (DELLY INFO `CONSENSUS`), breakpoint index `CONSBP` | the same |
| Position | pair must contain the two bases around the breakpoint (`SEQUENCE_TARGET = [bp - 1, 2]`) | 3' end of each primer 40-160 nt from the breakpoint (`gap_min`, `gap_max`) |
| Length | 18-27 nt, optimum 20 (Primer3 defaults; published primers are 19-26 nt) | 18-26 nt (`len_min`, `len_max`) |
| Tm | optimum 56.0, range 54.0-57.0 (published primers: 54.5-56.6) | 54.5-56.6, target 55.5 (`tm_lo`, `tm_hi`, `tm_target`) |
| GC | >= 30% (published: 30.8-50%) | 30-60% (`gc_lo`, `gc_hi`) |
| Homopolymer | <= 5 (Primer3 default; a published primer has a run of 5) | <= 4 (`max_homopolymer`) |
| Choice | lowest Primer3 pair penalty (\|Tm - optimum\| + \|length - optimum\| per primer) | per side the window with the smallest \|Tm - target\| + 0.6 when the 3' base is A or T |
| Product on the allele with the deletion | 60-300 bp; = last base of the reverse site - first base of the forward primer + 1 | follows from the windows; `r1 - f0` |
| Product on the allele without the deletion | the above + deletion length (END - POS) | the same (`band_deleted + svlen`) |

Tm model (both engines, and every Tm in the tables): `primer3.calc_tm` with its defaults, i.e. nearest neighbour,
SantaLucia (1998) parameters and salt correction, 50 mM monovalent, 1.5 mM divalent, 0.6 mM dNTP, 50 nM primer.
`pcr.tm_formula` is a pure-Python implementation of the same model, equal to `calc_tm` to 1e-12 on random
oligonucleotides (tested), used when primer3-py is absent, so the Tm scale never depends on what is installed.
The toolkit's own fallback (`deltools.nn_tm`: 250 nM, 16.6 log10[Na+] term, no divalent cations) runs about
10 degrees C lower; it is kept as `pcr.tm_toolkit_fallback` for the record and is **not** used (with it the
built-in routine would pick different primers on machines without primer3-py). The validation script
`FUN/code/primer_tm_gc_check_v1.py` recomputed Tm with another model (Biopython, Owczarzy salt correction, 250 nM:
about 3.4 degrees higher); that is a check of the original, not its design model, and is not ported.

New in BroodLink, not in the original tables: hairpin, self-dimer and cross-dimer Tm and free energy of every pair
(primer3-py `calc_hairpin`, `calc_homodimer`, `calc_heterodimer`; empty without primer3-py); the design sequence
cut from `--reference` (`pcr.reference_flank_bp` = 150 bp on each side of POS / END) for deletions without a
consensus; `--consensus-vcf` and `--pcr-picks`.

### 18.4 Checks

Ported (generic): primers exact substrings of the design sequence, forward before and reverse after the breakpoint
(`validate_primers_v1.py` V5, `primer_in_consensus_check_v1.py`); products recomputed from the primer positions and
the deletion length (V6); frequency difference at or above the threshold (V8); frequencies per class recomputed from
the genotypes (V7: they are computed by the step itself); in-silico PCR against `--reference`
(`primer_in_silico_pcr_v1.py`: <= 2 mismatches per primer found by three exact seeds, matching 3' base, forward on
+ and reverse on - of one sequence, product 60-6,000 bp, "at the locus" = within 2,500 bp of the deletion; V9:
product within 5 bp of the designed long product); Tm difference of the pair, lower GC, 3' base, longest
homopolymer (`primer_tm_gc_check_v1.py`); markers per interval and the reason for an interval without one.

Left out, because they need files of that cohort or a second assembly: V1-V4 (the deletion is in the catalogue VCF,
coordinates and consensus equal the published tables, PRECISE / FILTER flags and read support, which the reader does
not keep); the search in the other parental haplotype of the F1 (`primer_in_silico_mac_*.py`; a second run with
that assembly as `--reference` gives the same numbers); the look-up of cohort SNPs at primer mismatch positions
(`primer_site_cohort_check_v1.py`); the depth and heterozygosity support of the deletion calls
(`del_validation_markers_v1.py`).

### 18.5 Gel simulation and figures

`PCR/code/04_simulate_gels_MPL.py`: three lanes per deletion (no deletion: long band; heterozygote: both; homozygous
deletion: short band), migration on a log size axis, a ladder; "indicative only". The published sheets
(`FUN/code/build_figure_PCR_v3.py`, upgraded in `FUN/figures/upgrades_2026-11-03/round_ed/F-S14`, `F-S15`) draw the
lanes per karyotype class from the designed products with the ladder 100-1,000 bp by 100, 1.5, 2, 3, 4 kb, 18 and 17
gels per sheet. BroodLink: the same gel, `pcr.gels_per_page` = 9 (each published sheet becomes two consecutive
pages), `pcr.ladder_bp`. The assay figure is the port of `round_ed2/F-15/F-15.py` (panels A-H); its columns for the
other parental haplotype are not drawn. The schematic is the port of `PCR/code/05_primer_schematic_MPL.py`. Example
markers (typed by hand in the original): the first marker with the deletion on class 2 (panel A and B), the first
on class 0, the longest and the shortest deletion.

### 18.6 Tied to the cohort in the original

Inversion labels (A1-A41) and the STD / INV names read from stored tables (BroodLink: M01... in interval order,
class 0 / class 2); the four example markers of the figure; chromosome prefix added by hand in the in-silico script;
the path of the catalogue VCF and of the assembly; 35 as an asserted count; the dagger footnotes naming markers.

### 18.7 Not found

The script that designed the published primers and its Primer3 settings (18.1); the rule that chose among candidates
of equal rank (18.2); a lower size bound of "gel-resolvable"; the origin of the 0.30-0.70 heterokaryotype window
(only its effect in the published table). `03_design_primers_DELLY.py` is not in the toolkit's `code/` folder but
in `FAM/_inversions_localPCA_sv_discovery/code/`. A "Supplementary Note S12" cited by a draft does not exist (folder README).
