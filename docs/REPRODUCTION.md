# Reproduction of the original analysis

BroodLink was extracted from analysis scripts that had already produced published numbers
for one cohort (216 retained individuals of 226 sequenced). This page compares BroodLink
0.1.0 with those numbers. Each line is marked **reproduced**, **differs** (with the cause)
or **not run** (with the reason and the command). No threshold was tuned to match.

How it was run: on the original intermediate tables (the likelihood-based relatedness table,
the deletion and duplication genotype matrices, the karyotype table of the 41 intervals, the
stored per-interval SNP dosage, the stored GRMs), through a local, git-ignored configuration
that names the files and the ten excluded samples. No cohort data, sample identifier or
path is part of this repository. Raw reads, alignments and the genotype-likelihood files
were not reprocessed.

Seed 1; 20,000 permutations for the nomination; 1,000 for the family effect.

## (a) Full run

`broodlink run --relatedness <NGSrelate .res> --gt <DEL matrix> --gt <DUP matrix>
--haplotypes <karyotypes of the 41 intervals> --inversions <41 intervals.bed>
--samples <sample order> --config <local config>`

| Quantity | Published | BroodLink 0.1.0 | Status |
|---|---|---|---|
| Individuals analysed | 216 | 216 | reproduced |
| Families / individuals in families / unassigned | 51 / 188 / 28 | 51 / 188 / 28, the same sets of individuals | reproduced |
| Within-family pairs passing both thresholds | 332 of 363 | 332 of 363 | reproduced |
| Closed families | 38 (12 unresolved, 1 incompatible) | 38 closed, 13 not closed | reproduced for "closed"; the rule separating unresolved from incompatible was not found, so BroodLink reports "not closed" only |
| Families of at least 4 | 23 | 23 | reproduced |
| Unrelated subset | 78 (exact maximum set) | 78 (greedy set) | same size; the members were not compared |
| Units / unit pairs | 79 / 3,081 | 79 / 3,081 | reproduced |
| Nominated unit pairs (20,000 permutations) | 90 (54 family-family, 36 with an unassigned individual) | 90 (54 + 36), the same 90 pairs, with the original family labels and with BroodLink's own labels | reproduced |
| Pairs between families of at least 4 | 253 | 253 | reproduced |
| of which nominated / pass Gate 1 / rejected | 12 / 11 / 1 | 12 / 11 / 1 | reproduced |
| Gate 2, deletions: median sharing, matched null, P | 0.779, 0.701, 0.0023 | 0.7794, 0.7011, 0.0023 | reproduced |
| Gate 2, duplications | 0.663, 0.532, 0.0023 | 0.6631, 0.5321, 0.0023 | reproduced |
| Supported links / families / components | 6 / 11 / 5 | 6 / 11 / 5, the same six pairs | reproduced |
| Supported by DEL + DUP / DEL / DUP; Gate 1 only | 4 / 1 / 1; 5 | 4 / 1 / 1; 5 | reproduced |
| Post-audit network, families of at least 3: links, latent parents | 22 edges (11 accepted, 10 provisional, 1 rejected), 14 maximal cliques | 21 kept links + 1 rejected, 14 latent parents | reproduced |
| Slot-limit violations at minimum size 4 / 3 / 2 | 0 / 1 / 4 | 0 / 1 / 4 | reproduced |
| Barcode of the parent shared across three families (41 states) | 25 at posterior >= 0.90, 24 at >= 0.95 | 25 and 24 | reproduced |
| Family effect, median eta squared (null): deletions | 0.585 (0.285), 17,524 tested | 0.5846 (0.2845), 17,524 | reproduced |
| duplications | 0.528 (0.278), 1,182 | 0.5278 (0.2779), 1,182 | reproduced |
| inversions | 0.602 (0.266), 41 | 0.6024 (0.2655), 41 | reproduced |
| Cross-marker concordance: AUC, recovery (deletions 17,276; duplications 1,178; inversions 41) | 0.99982, 95.9%; 0.99708, 74.1%; 0.98058, 29.8% | the same values | reproduced after a correction (below) |

Notes on (a):

- **Monte Carlo.** The count of 90 depends on the number of permutations and on the seed
  (82 at 2,000 permutations in a first check; 84-85 at 1,500 and 88 at 5,000 across three
  seeds). The pairs between families of at least 4 did not change in any run. BroodLink
  names families by size, not by the original numbers; with either naming the same 90 pairs
  were nominated at 20,000 permutations and seed 1.
- **Correction made during this check.** The original GRM filter "called in at least 100 of
  216 individuals" had been stored as the rounded fraction 0.463, which dropped ten
  deletions called in exactly 100 individuals (17,266 instead of 17,276 markers). The
  default fractions were corrected to 0.4629 and 0.6944 and the concordance step rerun.
  This is a fix of a transcription error, not a tuned threshold.
- **Families are connected blocks.** The published Methods describe families as groups of
  mutually consistent pairs; the published 51 families are the connected blocks of full-sib
  pairs, and 13 of them contain a pair that fails a threshold. BroodLink reproduces the
  blocks and reports the closure.
- **No upper bound in the nomination.** The published statistic is "theta >= 0.0884", second
  degree or closer. With an upper bound of 0.177 the count differs (93 against 88 at 5,000
  permutations; see `SWEEP.md`).
- **QC flags.** BroodLink's QC ran on the 226 sequenced individuals and flagged 5 low-coverage
  samples, 3 pairs with identical genotypes and 21 individuals by its numeric contamination
  rule. That rule is new in BroodLink (the original screen fixed no number) and flags far
  more individuals than the original analysis excluded; the ten exclusions of the original
  were taken from the local configuration, not from these flags. Not compared further.
- The total number of barcode states called over all links (230 of 492) has no published
  counterpart: the published count for seven links (137 of 287) includes a link at minimum
  family size 3, which BroodLink's barcode step does not include by default.

## (b) One input at a time

Reference: the 51 families and the 6 supported links of the full run.

| Mode | Families identical to the 51 | Pair precision / recall | Links | Status |
|---|---|---|---|---|
| `--gt` only (deletion and duplication matrices; GRM of deletions, audit on duplications; through the command line) | 42 (49 called) | 0.80 / 0.99 | 6 supported: 4 of the reference six, 2 others; 5 of the six among the kept links | run; differs from the full run because a GRM has no k2 (families sharing a parent merge) and the nomination uses GRM kinship |
| SNPs: stored GRM of every third SNP (1,307,168), nomination and audit on deletions and duplications by the library functions, 5,000 permutations | 46 (49 called) | 0.98 / 0.99 | 7 supported: all six of the reference, 1 other | run outside the command line; not the single-input split mode |
| `--bisnp` only through the command line | | | | **not run**: reading the 28 genotype-likelihood files takes 20-30 minutes; command in `BENCHMARK.md` |
| `--snp_indel` only | | | | **not run**: the small-variant cohort VCF was not among the located inputs. `broodlink run --snp_indel <cohort.vcf.gz> --samples <list> --out out_snp_indel` |

## (c) `--relatedness` alone

`broodlink run --relatedness <NGSrelate .res> --samples <sample order> --config <local config>`

| Quantity | Expected | BroodLink 0.1.0 | Status |
|---|---|---|---|
| Families | the same 51 | 51, identical sets | reproduced |
| Nominated unit pairs | the same 90 | 90, the same pairs | reproduced |
| Audit | none | none; links labelled "nominated, not audited" | as designed |

## (d) Mask and depth

| Quantity | Published | BroodLink 0.1.0 | Status |
|---|---|---|---|
| Callable mask | 576,895,264 bp, 59.85% of 963,905,721 bp | | **not run**: needs the soft-masked reference FASTA. `broodlink run ... --mask <soft-masked.fa> --genome <fa.fai>`; the span and share are printed in every table header |
| Per-sample heterozygosity over the mask | | | **not run**: needs the genotype-likelihood files and the mask |
| Median depth over the mask / over the assembly | 2.54x / 9.64x | | **not run**: the per-sample depth table was not located. `--depth <table or mosdepth folder>` reports both in separate columns |

The mask code (BED and soft-masked FASTA, intersection, exclusion, the structural-variant
rule) is covered by unit tests and by the synthetic run only.

## (e) The 41 intervals

| Quantity | Published | BroodLink 0.1.0 | Status |
|---|---|---|---|
| Karyotypes from SNPs: BroodLink's local-PCA genotyping on the stored per-interval SNP dosage, against the published calls (all 216 individuals, forced calls) | | median 98% of individuals identical per interval; at least 95% in 33 intervals, at least 90% in 35; 5 intervals below 80% (49-69%), all five in the "candidate" tier of the original; identical for all individuals in 5 intervals | **differs**. BroodLink fits the three clusters on the unrelated subset and projects the others; the script that made the published full-cohort calls was not found. Fitting on all individuals gives at least 95% in 38 intervals, minimum 66% |
| Ambiguous genotypes set to missing | none (forced calls) | median 7% of individuals per interval | new rule in BroodLink |
| Tiers | 34 high confidence, 7 candidate | 35 and 6; the same tier for 40 of 41 intervals | **differs** by one interval (criterion 5 is defined as "above the 95th percentile of the null" because the original comparison was not found) |
| Karyotypes from deletions alone, against the published calls | not published | 41 of 41 intervals genotyped; median concordance 0.94; at least 0.90 in 25 intervals, at least 0.80 in 33; minimum 0.45; median 19% of individuals ambiguous | new; deletion-only karyotypes are usable for most intervals, not all |
| Sub-clusters inside arrangements supported by deletions | 32 of 37 splits (one class), 21 of 35 (the other) | 31 of 35 and 19 of 29 | **differs**: BroodLink's run used PC1 / PC2 of the deletion-based local PCA (no SNP input in this run) and sets ambiguous individuals to missing; the published values come from the SNP-based PCA |
| Arrangement-tag deletions | 2,080 catalogued, 875 with tag strength >= 0.8 | 1,621 with difference >= 0.6, 905 with >= 0.8, 70 fully private | **differs**: the catalogue script was not found; BroodLink requires 5 called individuals per class (the catalogue required 20 and 10), keeps only differences >= 0.6 and uses its own deletion-based genotypes in this run |
| Scan: published intervals overlapped by a deletion-scan candidate | not published | 8 of 41; 23 further candidates | weak recovery from deletions alone |
| Scan from SNPs (reimplemented from Methods) | 41 intervals came from the original scan | | **not run** on this cohort (needs the genotype-likelihood files); tested on synthetic data only |

## Visual check of the figures

BroodLink's figures for the full run were rendered (14 of 14 drawn) and compared by eye with
the published shared-parent figure and the matrix of all units for content: the same 12
nominated pairs among the 253, the same rejected pair, six supported links in five
components, and the same family on three latent parents at minimum size 3. Layout, family
order and colours differ by design (`DISPLAY_ITEMS.md`). This was a reading of the figures,
not a pixel or data comparison.

## What this page does not show

It shows that the refactored code returns the numbers of the analysis it was extracted from,
on the same tables. It is not an independent validation of that analysis: the cohort has no
pedigree, and the families, links and barcodes remain inferences.
