# AlphaGenome assessment TODO for the siTRF1/TRF1 study

Date: 2026-05-24

Purpose: evaluate whether AlphaGenome adds useful, testable evidence to the
existing conclusion that siTRF1 DEGs are not explained by direct TRF1 promoter
binding.

## Current assessment

AlphaGenome is potentially useful as a secondary sequence-to-function annotation
tool, not as primary evidence for TRF1 knockdown causality.

Supported uses for this project:

- annotate hg38 DEG loci with predicted RNA-seq, CAGE/PRO-cap, ATAC/DNase,
  histone ChIP-seq, TF ChIP-seq, and contact-map tracks;
- check whether candidate DEG promoters/enhancers look regulatory in relevant
  biosamples;
- test limited in silico sequence perturbations around candidate regulatory
  motifs, especially stress-response TF motifs, if a concrete variant or edited
  sequence is defined;
- compare candidate loci against local shuffled background sequences as a pilot
  hypothesis-generation step.

Unsupported or weak uses for this project:

- do not use AlphaGenome to infer the effect of TRF1 knockdown directly;
- do not use it to prove TERF1/TRF1 protein binding unless the `CHIP_TF`
  metadata confirms an applicable TERF1/TRF1 track and the biological context is
  relevant;
- do not treat predicted expression tracks as replacement evidence for the
  observed siTRF1 RNA-seq/DEG result;
- do not liftover the legacy hg18 TRF1 peaks for AlphaGenome. For this project,
  use independent newer hg38 resources and gene-centered hg38 loci instead.

## Evidence from official documentation

- AlphaGenome predicts sequence-derived molecular tracks and variant effects
  from DNA sequence, including RNA expression, chromatin accessibility, histone
  marks, TF binding tracks, splicing, and contact maps.
- The model uses human hg38 and mouse mm10 references. This project has TRF1
  ChIP-seq peaks originally identified as hg18, so any AlphaGenome workflow must
  use hg38 coordinates.
- The maximum input interval is 1,048,576 bp. Smaller supported lengths include
  approximately 16 kb, 100 kb, and 500 kb, but official guidance recommends 1 Mb
  where possible.
- The API is non-commercial, rate-limited by demand, and suited to smaller or
  medium-scale analyses rather than very large runs.
- Stated limitations include tissue specificity, long-range interactions,
  species scope, personal genomes, molecular-only predictions, and lack of
  diploid-aware inference.

Sources:

- https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/
- https://www.nature.com/articles/s41586-025-10014-0
- https://www.alphagenomedocs.com/
- https://www.alphagenomedocs.com/faqs.html
- https://www.alphagenomedocs.com/exploring_model_metadata.html

## Required TODOs before any biological claim

- [x] Obtain an AlphaGenome API key under the non-commercial terms.
- [x] Query `output_metadata(HOMO_SAPIENS)` and save the metadata table under
      `results/alphagenome/metadata/`.
- [x] Verify whether any `CHIP_TF` track is TERF1/TRF1-related. If no such
      track exists, explicitly document that AlphaGenome cannot test direct TRF1
      binding for this study.
- [x] Identify the closest available biosample(s) to the original siTRF1 cell
      context. If the original cell line/tissue is absent or weakly matched,
      document the mismatch and treat all outputs as exploratory.
- [x] Use hg38 gene-centered loci and newer hg38 target-gene resources. No
      hg18-to-hg38 conversion was performed.

Completed findings:

- AlphaGenome human metadata has 5,563 rows and 1,617 `CHIP_TF` rows, but
  0 TERF1/TRF1 candidate tracks.
- The closest metadata match to the legacy BJ fibroblast context is BJ
  (`EFO:0002779`), but its available TF ChIP-seq output is CTCF, not TERF1.
- MSigDB/GTRD `TERF1_TARGET_GENES` has 299 target genes and 0 overlap with the
  274 siTRF1 DEGs.
- ChIP-Atlas hg38 TERF1 target-gene tables show 30-32 overlap DEGs at
  TSS +/-1/5/10 kb, but the signal is dominated by one LoVo experiment and
  disappears at the strongest average-score tier (`average_score >= 100`).

## Pilot analysis TODOs

- [x] Run a small pilot on representative DEGs:
      - core stress genes: `DDIT4`, `ASNS`, `CHAC1`, `DDB2`, `PHLDA3`;
      - IFN genes: `IFITM1`, `RNASEL`;
      - mitochondrial/stress-linked genes: `MT-ND1` should be handled
        separately because AlphaGenome's standard human reference workflow is
        nuclear hg38-centric and chrM interpretation needs explicit validation.
- [x] For each gene, request hg38 intervals centered on the TSS and collect
      RNA-seq, DNase, histone ChIP-seq, and TF ChIP-seq outputs where available.
- [x] Run a second AlphaGenome batch on the 32 ChIP-Atlas TERF1-overlap DEG
      candidates.
- [x] Summarize whether predicted regulatory activity is consistent with
      stress-response regulatory architecture, without claiming perturbation
      causality.
- [ ] If a concrete motif hypothesis is selected, perform in silico mutagenesis
      only on defined motifs and compare against matched shuffled controls.
- [x] Add a short section to `docs/FINAL_REPORT.md` only after the pilot results
      are saved and reproducible.

## Interpretation rules

- AlphaGenome results can strengthen the annotation of candidate regulatory
  regions but cannot overturn the existing direct-binding null result by itself.
- A positive regulatory prediction near a DEG is not evidence that TRF1 binds
  there.
- A negative prediction is not definitive proof that the region is inactive in
  the original siTRF1 experimental context.
- Any AlphaGenome-derived claim must name the exact output type, biosample,
  genomic interval, reference build, and scorer/aggregation method.
