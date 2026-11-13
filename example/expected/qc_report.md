# BroodLink QC report

**Analysis set: 80 individuals** (0 removed by QC, 3 flagged and kept). Families are relatedness-based groups, not pedigree records. callable span: 47,620,925 bp (79.37% of the 60,000,000-bp assembly).

QC flags samples and continues; it removes a sample only when the sample is confirmed in the config (`samples.exclude_confirmed`) or the run uses `--auto-exclude`.

## Rules

- **low coverage**: call rate below 0.5 x cohort median, per input
- **pairs of samples with identical genotypes**: theta >= 0.354 or genotype concordance >= 0.99 (hard genotypes)
- **possible contamination**: median relatedness or heterozygosity above cohort median + 5.0 MAD (new in BroodLink)
- **removal**: only samples confirmed in the config (samples.exclude_confirmed) or, with --auto-exclude, low-coverage samples and one sample of each identical pair

## Inputs

| input | marker_class | n_samples | n_markers | filters | mean_call_rate | af_median |
|---|---|---|---|---|---|---|
| SNP (bisnp) | SNP | 80 | 5863 | callable mask: 6000 -> 5863 | 0.7697 | 0.5 |
| DEL (gt:sv.vcf.gz) | DEL | 80 | 768 | callable mask: 944 -> 768 | 0.976 | 0.2215 |
| DUP (gt:sv.vcf.gz) | DUP | 80 | 203 | callable mask: 243 -> 203 | 0.9776 | 0.1711 |

## Mask

callable span: 47,620,925 bp (79.37% of the 60,000,000-bp assembly).

| input | marker_class | markers_total | markers_inside_mask | markers_outside_mask | rule | warning |
|---|---|---|---|---|---|---|
| SNP (bisnp) | SNP | 6000 | 5863 | 137 | position callable |  |
| DEL (gt:sv.vcf.gz) | DEL | 944 | 768 | 176 | >= 0.5 of the span and both breakpoints callable |  |
| DUP (gt:sv.vcf.gz) | DUP | 243 | 203 | 40 | >= 0.5 of the span and both breakpoints callable |  |

## Flagged samples

| sample | reason | rule | value | action |
|---|---|---|---|---|
| IND004 | low coverage | call rate below 0.5 x cohort median (0.396) in SNP (bisnp) | 0.098 | flagged, kept |
| IND004 | possible contamination | heterozygosity above cohort median + 5.0 MAD in SNP (bisnp) (new in BroodLink; the original fixed no number) | 0.8115 | flagged, kept |
| IND035 | pair of samples with identical genotypes | theta >= 0.354; with --auto-exclude the sample of the pair with the lower call rate is removed | pair with IND036 | flagged, kept |
| IND036 | pair of samples with identical genotypes | theta >= 0.354; with --auto-exclude the sample of the pair with the lower call rate is removed | pair with IND035 | flagged, kept |

## Pairs of samples with identical genotypes

| sample_1 | sample_2 | theta | genotype_concordance | n_markers_compared | rule |
|---|---|---|---|---|---|
| IND035 | IND036 | 0.5 | nan | nan | theta >= 0.354 |

## Depth

Cohort median depth over the mask: 5.99x; over the assembly: 10.09x (two different denominators, reported in separate columns of `qc_per_sample.tsv`).

## Sex / batch / group effects

| variable | input | statistic | n_groups | test | p_value | flag |
|---|---|---|---|---|---|---|
| sex | SNP (bisnp) | call_rate | 2 | Kruskal-Wallis (new in BroodLink) | 0.6072 |  |
| sex | SNP (bisnp) | heterozygosity | 2 | Kruskal-Wallis (new in BroodLink) | 0.8805 |  |
| sex | DEL (gt:sv.vcf.gz) | call_rate | 2 | Kruskal-Wallis (new in BroodLink) | 0.508 |  |
| sex | DEL (gt:sv.vcf.gz) | heterozygosity | 2 | Kruskal-Wallis (new in BroodLink) | 0.9691 |  |
| sex | DUP (gt:sv.vcf.gz) | call_rate | 2 | Kruskal-Wallis (new in BroodLink) | 0.1684 |  |
| sex | DUP (gt:sv.vcf.gz) | heterozygosity | 2 | Kruskal-Wallis (new in BroodLink) | 0.2926 |  |
| batch | SNP (bisnp) | call_rate | 2 | Kruskal-Wallis (new in BroodLink) | 0.4954 |  |
| batch | SNP (bisnp) | heterozygosity | 2 | Kruskal-Wallis (new in BroodLink) | 0.809 |  |
| batch | DEL (gt:sv.vcf.gz) | call_rate | 2 | Kruskal-Wallis (new in BroodLink) | 0.9574 |  |
| batch | DEL (gt:sv.vcf.gz) | heterozygosity | 2 | Kruskal-Wallis (new in BroodLink) | 0.7828 |  |
| batch | DUP (gt:sv.vcf.gz) | call_rate | 2 | Kruskal-Wallis (new in BroodLink) | 0.8253 |  |
| batch | DUP (gt:sv.vcf.gz) | heterozygosity | 2 | Kruskal-Wallis (new in BroodLink) | 0.8052 |  |
| group | SNP (bisnp) | call_rate | 1 | Kruskal-Wallis (new in BroodLink) | nan |  |
| group | SNP (bisnp) | heterozygosity | 1 | Kruskal-Wallis (new in BroodLink) | nan |  |
| group | DEL (gt:sv.vcf.gz) | call_rate | 1 | Kruskal-Wallis (new in BroodLink) | nan |  |
| group | DEL (gt:sv.vcf.gz) | heterozygosity | 1 | Kruskal-Wallis (new in BroodLink) | nan |  |
| group | DUP (gt:sv.vcf.gz) | call_rate | 1 | Kruskal-Wallis (new in BroodLink) | nan |  |
| group | DUP (gt:sv.vcf.gz) | heterozygosity | 1 | Kruskal-Wallis (new in BroodLink) | nan |  |

## Relatedness

Source: --relatedness, used as is (plain TSV, 80 samples, 3160 pairs).

| class | n_pairs |
|---|---|
| unrelated | 2808 |
| second degree | 190 |
| full-sib | 143 |
| first degree, not full-sib (k2 below gate) | 12 |
| third degree | 6 |
| identical genotypes (theta at or above the identity threshold) | 1 |

Parent-offspring pairs (k1 >= 0.8 and k2 <= 0.05): 0.
