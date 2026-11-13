# A quick benchmark

Scope, stated first: this is a small check that every input mode runs and of what each mode
recovers, on **one simulated cohort (80 individuals, one seed)** and on **one real cohort (216
individuals)**. It is not a performance study and not a comparison with other tools. Run
times are wall-clock seconds of a whole run on one laptop, with other jobs running at the
same time; treat them as orders of magnitude.

## 1. Synthetic cohort: recovery of the simulated truth

Reproduce with `broodlink benchmark --out benchmark_out` (about 2 minutes; 2,000
permutations, figures off). Truth: 15 families (72 individuals), 7 pairs of families sharing
one parent, the genotype of each shared parent at 4 structural haplotypes (28 states).

| Mode | Relatedness | Families identical to truth | Pair precision / recall | True links kept (of 7) | False links | Supported, true / false | Barcode: called, correct (of states) | Seconds |
|---|---|---|---|---|---|---|---|---|
| full run (relatedness + bisnp + gt + inversions BED) | table | 15 / 15 | 1.00 / 1.00 | 7 | 0 | 5 / 0 | 17, 17 (28) | 17 |
| relatedness + gt + supplied haplotype genotypes | table | 15 / 15 | 1.00 / 1.00 | 7 | 0 | 5 / 0 | 17, 17 (28) | 4 |
| bisnp + gt, no relatedness table | GRM, 5,451 SNPs | 15 / 15 | 1.00 / 1.00 | 5 | 0 | 3 / 0 | 12, 12 (20) | 31 |
| `--relatedness` only | table | 15 / 15 | 1.00 / 1.00 | 7 | 0 | not audited | no barcode | 6 |
| `--bisnp` only | GRM, 5,451 SNPs | 15 / 15 | 1.00 / 1.00 | 7 | 0 | none supported (no SV class: Gate 2 cannot run) | no barcode | 18 |
| `--gt` only (DEL + DUP in one VCF) | GRM, 605 deletions | 12 / 15 | 0.84 / 0.98 | 3 | 0 | 1 / 0 | no barcode | 7 |
| `--snp_indel` only (hard-call VCF) | GRM, 5,451 SNPs | 15 / 15 | 1.00 / 1.00 | 7 | 0 | none supported (no SV class) | no barcode | 17 |

Reading the table:

- With a relatedness table (theta and k2), families and links are recovered completely. "Supported"
  is 5 of 7 by construction: the Gate-2 bar is the median of the compatible pairs, so about
  half pass. The other two true links stay "Gate-1 only".
- A GRM of hard-called low-coverage SNPs recovers the families but underestimates half-sib
  kinship (mean about 0.08 against 0.125 expected), so some true links are not nominated
  in the two-input mode (5 of 7).
- 605 deletions are below the smallest marker set of the original GRM benchmark (1,178), and
  the run labels the result provisional: three families that share a parent are merged or
  split, and only 3 of 7 links are kept.
- Every state called at posterior >= 0.90 equals the true genotype of the shared parent
  (17 of 17); 11 of 28 states stay unresolved.
- The synthetic data are easy: no population structure, no inbreeding, one generation, clean
  simulated genotypes, and the relatedness table is the realised identity by descent plus small
  noise rather than an estimate from reads.

Table: `benchmark_out/benchmark_synthetic.tsv`.

## 2. Original cohort: each mode against the full run

Reference = the published full run (51 families holding 188 of 216 individuals; 6 supported
links). No pedigree exists, so this measures agreement with the full run, not correctness.
Runs were made through a local, git-ignored configuration; no sample identifier appears here.

| Mode | How it was run | Families identical to the 51 | Pair precision / recall | Supported links: in the 6 / other | Seconds |
|---|---|---|---|---|---|
| full run (relatedness + DEL + DUP + 41 intervals) | `broodlink run`, 20,000 permutations | 51 / 51 | 1.00 / 1.00 | 6 / 0 | not timed as one run (interrupted and resumed step by step) |
| `--relatedness` only | `broodlink run`, 20,000 permutations | 51 / 51 | 1.00 / 1.00 | not audited; all 6 among the 54 nominated family pairs | 183 |
| `--gt` only (DEL and DUP matrices) | `broodlink run`, 20,000 permutations; GRM of 17,276 deletions, audit on duplications | 42 / 51 (49 called) | 0.80 / 0.99 | 4 / 2 | 985 |
| SNPs (stand-in for `--bisnp` with a second input) | stored SNP GRM of the original benchmark (every third SNP, 1,307,168 SNPs), nomination and audit by the library functions, 5,000 permutations; **not** through the command line | 46 / 51 (49 called) | 0.98 / 0.99 | 6 / 1 | not timed |
| `--bisnp` only through the command line | **not run** (see below) | | | | |
| `--snp_indel` only | **not run**: the small-variant cohort VCF was not among the located inputs | | | | |

Families from each stored GRM at kinship >= 0.177 (no k2 exists for a GRM):

| GRM | Families called | Identical to the 51 | Pair precision | Pair recall |
|---|---|---|---|---|
| SNPs, every third (1,307,168) | 49 | 46 | 0.98 | 0.99 |
| SNPs, SV regions masked (985,007) | 48 | 45 | 0.84 | 0.99 |
| Deletions (about 17,300) | 49 | 42 | 0.80 | 0.99 |
| Duplications (about 1,200) | 29 | 10 | 1.00 | 0.32 |
| 41 inversions | 1 | 0 | 0.02 | 0.99 |

Reading the table:

- The relatedness-only run returns the same 51 families and the same 90 nominated unit
  pairs as the full run.
- Without k2, a GRM cannot separate full sibs from other first-degree pairs, and families that
  share a parent merge: precision falls to 0.80-0.98 while recall stays at 0.99. This is the
  main limit of the single-input modes.
- Duplications alone and inversions alone do not reconstruct the families, in line with the
  recovery of 74% and 30% in the original cross-marker benchmark.
- In the deletion-only run the labels of the families differ from the reference, and its "6
  supported" are 4 of the reference links plus 2 others; links from this mode are provisional.

Not run, with the command to run later (about 20-30 minutes on a laptop, mostly reading 28
Beagle files; it was left out to keep every computation short while the data drive could be
unplugged):

```
broodlink run --bisnp <folder of Beagle files> --samples <sample order file> \
    --haplotypes <karyotype table> --inversions <intervals.bed> \
    --config <local config with genotypes.bisnp_step: 3 and the excluded samples> --out out_bisnp_only
```

## 3. The GRM benchmark of the original project (quoted, not rerun)

The VanRaden estimator used by BroodLink was benchmarked in the original project against the
SNP-relatedness families (rank AUC of within-family against between-family pairs; recovery =
share of within-family pairs above the 99.9th percentile of the between-family pairs;
family-level bootstrap):

| Marker class | Markers | AUC [95% CI] | Recovery % [95% CI] |
|---|---|---|---|
| biSNPs | 1,307,168 | 0.99996 [0.99989, 1.00000] | 99.4 [96.3, 100] |
| biSNPs, SV regions masked | 985,007 | 0.99990 [0.99975, 0.99998] | 97.2 [90.8, 100] |
| DELs | 17,276 | 0.99982 [0.99961, 0.99997] | 95.9 [90.6, 99.7] |
| DUPs | 1,178 | 0.99708 [0.99516, 0.99856] | 74.1 [63.2, 84.0] |
| polymorphic INVs | 41 | 0.98058 [0.97395, 0.98590] | 29.8 [17.8, 48.4] |

These are observations of cross-marker concordance. That benchmark set **no GRM-specific
kinship threshold and no minimum marker count**; BroodLink applies kinship >= 0.177 because
an older script of the project did, and quotes the two marker counts above as reference
points in its marker-count check. BroodLink's own cross-marker step (09) on the same cohort
returned the same marker counts (17,276, 1,178 and 41), AUC values (0.99982, 0.99708, 0.98058)
and recovery (95.9%, 74.1%, 29.8%): reproduced. A first run kept 17,266 deletions because the
call-rate filter "100 of 216" had been stored as the rounded fraction 0.463, which excluded
markers called in exactly 100 individuals; the default fractions were corrected to 0.4629 and
0.6944 and the step was rerun. Bootstrap intervals are not computed by BroodLink.
