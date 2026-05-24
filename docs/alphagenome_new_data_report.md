# AlphaGenome and Newer TERF1 Data Report

Date: 2026-05-24

Scope: evaluate whether newer public TERF1/TRF1 resources and AlphaGenome add
evidence for direct upstream/downstream regulation of the 274 siTRF1 DEGs. No
hg18-to-hg38 liftover was performed.

## Sources Checked

- AlphaGenome official model/API documentation and metadata.
- MSigDB `TERF1_TARGET_GENES`, derived from GTRD v20.06 ChIP-seq harmonization.
- ChIP-Atlas hg38 TERF1 target-gene tables at TSS +/-1 kb, +/-5 kb, and
  +/-10 kb.
- ChIP-Atlas experiment metadata for the TERF1 hg38 target-gene table.

## Reproducible Outputs

- `results/alphagenome/metadata/alphagenome_metadata_probe_summary.md`
- `results/alphagenome/pilot/alphagenome_bj_pilot_summary.md`
- `results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_summary.md`
- `results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_gene_summary.md`
- `results/new_data/terf1_gtrd_deg_overlap_summary.md`
- `results/new_data/chipatlas/chipatlas_terf1_deg_overlap_summary.md`
- `results/new_data/chipatlas/chipatlas_terf1_score_cutoff_summary.tsv`

## Results

### AlphaGenome Metadata

- Human output metadata rows: 5,563.
- Human `CHIP_TF` rows: 1,617.
- Unique `CHIP_TF` transcription factors: 751.
- TERF1/TRF1 candidate tracks: 0.

Conclusion: AlphaGenome cannot currently be used as direct TERF1/TRF1 ChIP-TF
evidence for this study. It can only annotate predicted basal regulatory
context around candidate loci.

### MSigDB / GTRD TERF1 Target Genes

- TERF1 target gene symbols: 299.
- siTRF1 DEG symbols: 274.
- Overlap: 0 genes.
- Hypergeometric p-value over the local annotation universe: 1.

Conclusion: the harmonized GTRD promoter-target gene set gives no support for
direct TERF1 control of the siTRF1 DEGs.

### ChIP-Atlas hg38 TERF1 Target Genes

ChIP-Atlas hg38 target-gene tables contain four TERF1 experiments:

- SRX4957970: HEK293 / HHV6-GFP, TRF1 ChIP-seq.
- SRX359962: LoVo colon adenocarcinoma, TERF1 ChIP-seq.
- SRX472275: lymphoblastoid cell line, TERF1 ChIP-seq.
- SRX4957969: smooth muscle / iciHHV6, TRF1 ChIP-seq.

Overlap with the 274 siTRF1 DEGs:

| TSS window | TERF1 target genes | Overlap DEGs | Direction | p-value | Dominant source | Strong tier (`average_score >= 100`) |
|---|---:|---:|---|---:|---|---:|
| +/-1 kb | 2,566 | 30 | 8 up / 22 down | 1.44e-06 | SRX359962 LoVo, 29 genes | 0 |
| +/-5 kb | 2,792 | 30 | 8 up / 22 down | 7.76e-06 | SRX359962 LoVo, 29 genes | 0 |
| +/-10 kb | 3,211 | 32 | 9 up / 23 down | 1.66e-05 | SRX359962 LoVo, 30 genes | 0 |

Interpretation: this is a real newer-data signal, but it is weak for causal
interpretation because it is dominated by one non-matched LoVo experiment and
the strongest average-score tier has no DEG overlap. The result should be used
to prioritize candidate genes for follow-up, not as proof of direct TRF1
regulation in the siTRF1 system.

### AlphaGenome BJ Candidate Annotation

The 32 ChIP-Atlas overlap DEGs were annotated in AlphaGenome using BJ
(`EFO:0002779`) metadata and 131,072 bp TSS-centered hg38 intervals.

- Candidate genes: 32.
- Successful loci: 32.
- Track rows: 224.
- Failed genes: 0.
- Available BJ TF ChIP-seq track: CTCF only, not TERF1/TRF1.

Top TSS-context candidates:

- DNase at TSS: `CDK19`, `DDIT4`, `ASNS`, `EFCAB14`, `DCBLD2`.
- H3K4me3 at TSS: `ZSCAN5A`, `MLEC`, `ARRDC1`, `ZNF524`, `BLOC1S6`.
- RNA max signal: `COX6A1`, `RPS27L`, `VMAC`, `LMAN2`, `CDIPT`.
- CTCF at TSS: `CDIPT`, `SIRT4`, `RPL22L1`, `PRELID2`, `RPS27L`.

Interpretation: AlphaGenome confirms that several ChIP-Atlas-overlap genes
have plausible regulatory context in BJ-like metadata, especially stress and
housekeeping loci such as `DDIT4`, `ASNS`, `CDK19`, and `COX6A1`. This does not
establish TRF1 binding because AlphaGenome has no TERF1/TRF1 output track and
does not model siTRF1 knockdown.

## Professional Conclusion

The strongest defensible conclusion is still that the siTRF1 DEG program is
mostly downstream/indirect rather than direct TRF1 promoter regulation.

What changed after newer-data analysis:

- GTRD/MSigDB independently supports the null direct-target conclusion
  (`0/274` overlap).
- ChIP-Atlas hg38 adds a candidate list of 30-32 overlapping DEGs, but the
  evidence is context-mismatched and score-sensitive.
- AlphaGenome is useful for candidate prioritization and regulatory annotation,
  not for proving direct TERF1/TRF1 binding.

Priority follow-up candidates from the newer-data plus AlphaGenome pass:

- `DDIT4` and `ASNS`: stress-response genes with ChIP-Atlas overlap and strong
  AlphaGenome BJ DNase/TSS regulatory signal.
- `CDK19`: strongest AlphaGenome BJ DNase-at-TSS candidate among ChIP-Atlas
  overlaps.
- `COX6A1`, `RPS27L`, `LMAN2`, `CDIPT`: strong BJ expression/regulatory signal,
  but direct TERF1 causality remains unproven.

These candidates are appropriate for manual genome-browser inspection and, if
biologically important, experimental validation by TRF1 ChIP-qPCR/CUT&RUN in
the matched siTRF1 cell context.
