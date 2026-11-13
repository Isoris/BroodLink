# Changelog

## Unreleased

PCR markers (step `10_pcr_markers`, `broodlink pcr-markers`): for every structural haplotype
the deletions inside it are screened, and one flanking primer pair is designed on the sequence
across the junction of the first suitable deletion.

- New module `broodlink/pcr.py`: screening rules of the published screening table (allele
  frequency difference >= 0.90 between the homokaryotype classes, heterokaryotypes at 0.30-0.70,
  deletion <= 1.8 kb), two design engines (Primer3 through the optional package primer3-py, the
  default and the procedure of the original study; the routine of the original deletion toolkit
  as the fallback), Tm as `primer3.calc_tm` (a pure-Python formula equal to it when primer3-py
  is absent), checks per pair, in-silico PCR against `--reference`, hairpin and dimer values
  (new in BroodLink).
- New inputs: `--reference FASTA` (small indexed reader, `.fai` used when present),
  `--consensus-vcf`, `--pcr-picks`, `--primer-engine`. The VCF reader keeps the split-read
  consensus of SV records (DELLY `CONSENSUS`, `CONSBP`).
- New tables `pcr_candidate_deletions.tsv`, `pcr_primers.tsv`, `pcr_markers_per_interval.tsv`
  (marker counts joined with the primers, one row per interval), `pcr_primer_order_list.tsv`;
  also `T13` to `T16` in `tables/` and the XLSX. Run report section and Methods paragraph.
- New figures `fig16_pcr_assay` (port of the assay figure, panels A-H), `fig17_pcr_simulated_gels`
  (paged, `pcr.gels_per_page` = 9) and `fig18_pcr_primer_schematic`.
- Example data: the simulated deletions carry a consensus; `simulate_pedigree.py --reference-only`
  writes the simulated reference (61 MB, not committed).
- `environment.yml` (conda / mamba, with primer3-py) and the extra `broodlink[primers]`.
- Finding recorded in `docs/EXTRACTED_METHOD.md`, section 18: the published primers were
  designed with Primer3, not with the toolkit routine (which returns none of them).

Figures: the published figures of the original study are now drawn by `broodlink run` from
the tables of the run (`docs/DISPLAY_ITEMS.md` has the panel-by-panel status and the
comparison with the original renders).

- `broodlink/figures/` replaces `figures.py`: house style of the originals (`style.py`), the
  views without a published counterpart (`generic.py`) and one module per published figure.
  Sixteen figures (`fig00` to `fig15`), PNG at 450 dpi (default was 300) and PDF, none taller
  than A4. `figures/MANIFEST.tsv` gains a `note` column (panels left out because their input
  is absent); `figures/LEDGER.tsv` lists every number printed in a ported panel with the table
  it came from.
- Ported: the shared-parent figure (matrix, nomination score, Gate 1, Gate 2, incidence,
  assignment completeness, opposite homozygotes, network by component), the matrix of all
  units, SV genotypes and kinship per marker class, cross-marker concordance (four panels and
  both headlines), relatedness and families, the barcode with deletion evidence, the barcode
  at gene level, the funnel to an unsampled parent, deletion genotypes along linked families
  and for the other links. `fig00_overview` (kinship matrices and network) is the front-page
  figure. Removed because a port shows the same thing: the old matrix, gates, incidence,
  network, theta / k2, family sizes, family effect, concordance and composite figures.
- biSNP panels are drawn whenever SNPs are given (kinship matrix, kinship distributions,
  matched biSNPs of the family effect, biSNP R0) and left out cleanly otherwise.
- New input `--genes` (GFF3 or BED): genes of every barcode interval and the gene-level
  barcode figure. Simulated gene file `tests/data/genes.gff3`.
- Step 5 (network) also writes: opposite homozygotes by relationship class
  (`opposite_homozygotes_*.tsv`), deletions that fit a Mendelian cross along chains of linked
  families and sharing per chromosome (`deletion_path_*.tsv`, `deletion_sharing_per_chromosome.tsv`).
- Step 7 (barcode) also writes: tag-variant evidence on the shared parent with its sweep and
  permutation check, the vote of a junction's links, the recurrence test of a junction parent
  across two links (previously listed as not implemented), genes per interval with `--genes`.
- Step 8 (family effect): frequency-matched SNPs as the reference level of every non-SNP class.
- Step 9 (concordance): the GRM of every marker class is written (`kinship_grm.npz` with
  `kinship_grm_individuals.tsv`); the summary table carries bootstrap intervals over families
  for AUC and recovery (previously not implemented), the class "SNPs outside SV regions" and
  the Pearson correlation of each SV GRM with it.
- Families imported with `--families` are described with the closure rules (closure, status,
  membership class) instead of being labelled "imported".
- The relatedness reader keeps the R0 column when the table has one.
- Progress output: step headers with timing, progress lines for long loops (silent for loops
  under 2 s), `--quiet`.
- New config sections `concordance`, `deletion_paths`, `barcode.tag_evidence`,
  `barcode.augmentation`, `barcode.recurrence`, and `figures.individual_noun`,
  `figures.haplotype_label`, `figures.max_chain_figures`, `figures.gene_level_*`.
- Present in the code but not wired and not tested: a likelihood-based relatedness estimator
  (`relatedness.method`, `relatedness.ml`); `--bisnp` without a relatedness table still uses
  NGSrelate when installed, otherwise the GRM.

## 0.1.0 (2026-11-13)

First release. The analysis logic is extracted from the scripts of one hatchery-cohort study
(see `docs/EXTRACTED_METHOD.md`); inputs are generalised and the refactor is checked against
the original results (`docs/REPRODUCTION.md`).

- `broodlink run`: QC, relatedness, full-sib families, half-sib nomination, two-gate audit,
  network with latent parents and slot limit, structural haplotypes (local-PCA genotyping,
  automatic scan, tag variants, sub-clusters), shared-parent barcode, family effect,
  cross-marker concordance; tables (TSV + XLSX), statistics, figures, tracks, reports.
- One subcommand per step, `report`, `sweep`, `benchmark`.
- Synthetic cohort with truth tables, generated by `tests/data/simulate_pedigree.py`.
- Reimplemented from the written Methods (no original code existed): QC rules, the SNP
  local-PCA scan.
- Not implemented: see "Limits" in the README.
