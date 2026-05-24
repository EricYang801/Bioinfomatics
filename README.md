# Bioinfomatics

Bioinformatics pipelines centered on **TRF1 ChIP-seq promoter overlap with DEG gene lists**, plus archived utilities and follow-up analyses for TRF1, TERF1, AlphaGenome, ChIP-Atlas, BETA, enrichment, and upstream-sequence extraction.

## File inventory

This inventory covers the reader-facing repository contents and local generated outputs; it excludes `.git/`, macOS `.DS_Store`, and editor metadata.

### Top level

| Path | Type | Use |
|---|---|---|
| **`LICENSE`** | *[doc]* | MIT license for the repository. |
| **`README.md`** | *[doc]* | Main entry point documenting repository contents, setup, inputs, conventions, tests, and legacy workflows. |
| **`.gitignore`** | *[config]* | Ignores large local references, generated result directories, Python caches, virtualenvs, and editor scratch files. |
| **`.gitattributes`** | *[config]* | Stores `github_release_assets/*.tar.gz` through Git LFS. |
| **`environment.yml`** | *[config]* | Conda environment with Python 3.12, bedtools, pybedtools, pyfaidx, pandas, matplotlib, and editable package install. |
| **`Makefile`** | *[config]* | Convenience targets for installing the package, running the TRF1 pipeline, excluding `chrM`, and cleaning outputs. |
| **`pyproject.toml`** | *[config]* | Python package metadata, dependencies, setuptools configuration, and `trf1-promoter-overlap` console entry point. |
| **`Homo_sapiens.GRCh38.108.gtf`** | *[input]* | Local Ensembl release 108 GTF used to regenerate `data/gene_annotation.tsv` when needed; ignored because it is large. |
| **`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`** | *[input]* | Local Ensembl GRCh38 soft-masked primary-assembly FASTA used by archived upstream extraction and local motif scans; ignored because it is large. |
| **`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai`** | *[generated]* | FASTA index for random access by tools such as pyfaidx; ignored through the `*.fai` rule. |
| **`github_release_assets.sha256`** | *[doc]* | SHA256 checksum sidecar for the upstream-sequence release tarball. |
| **`src/`** | *[script]* | Python source root for the importable TRF1 promoter-overlap package. |
| **`scripts/`** | *[script]* | Standalone wrappers and exploratory analysis scripts that use repo inputs and write `results/` outputs. |
| **`tests/`** | *[script]* | Pytest smoke tests for coordinate and chromosome helper logic. |
| **`data/`** | *[input]* | Small version-controlled inputs for the main TRF1 promoter-overlap pipeline. |
| **`docs/`** | *[doc]* | Human-readable reports, format notes, and interpretation documents. |
| **`legacy/`** | *[script]* | Archived upstream-sequence extraction workflow retained for provenance. |
| **`github_release_assets/`** | *[generated]* | Git-LFS release assets for the archived upstream FASTA extraction output. |
| **`results/`** | *[output]* | Gitignored output root for primary pipeline results and follow-up analyses. |

### `src/trf1_overlap/`

| Path | Type | Use |
|---|---|---|
| **`src/trf1_overlap/`** | *[script]* | Importable package implementing the primary promoter-overlap pipeline. |
| **`src/trf1_overlap/__init__.py`** | *[script]* | Package marker and version string. |
| **`src/trf1_overlap/__main__.py`** | *[script]* | Enables `python -m trf1_overlap` by delegating to the CLI. |
| **`src/trf1_overlap/annotation.py`** | *[script]* | Reads BED-style gene annotation TSVs or converts Ensembl GTF gene records into cached annotation TSVs. |
| **`src/trf1_overlap/chroms.py`** | *[script]* | Normalizes chromosome names to UCSC style and strips Ensembl gene ID version suffixes. |
| **`src/trf1_overlap/cli.py`** | *[script]* | Command-line orchestration for reading inputs, building promoters, intersecting peaks, plotting, and reporting. |
| **`src/trf1_overlap/deg.py`** | *[script]* | Reads DEG tables, validates regulation labels, and matches DEGs to annotation by gene ID or gene name. |
| **`src/trf1_overlap/intersect.py`** | *[script]* | Wraps pybedtools and bedtools for promoter-by-peak intersection. |
| **`src/trf1_overlap/peaks.py`** | *[script]* | Parses SISSRs TRF1 peak files and writes cleaned BED coordinates. |
| **`src/trf1_overlap/plots.py`** | *[script]* | Generates matched-gene bar plots and UP/DOWN stacked bar plots with matplotlib. |
| **`src/trf1_overlap/promoter.py`** | *[script]* | Builds 1000, 2000, and 3000 bp TSS-anchored upstream promoter windows and supports chromosome exclusion. |
| **`src/trf1_overlap/report.py`** | *[script]* | Formats overlap TSVs, summarizes matched genes, and writes the text report. |

### `scripts/`

| Path | Type | Use |
|---|---|---|
| **`scripts/trf1_promoter_overlap.py`** | *[script]* | Direct checkout-friendly wrapper for the packaged `trf1_overlap` CLI. |
| **`scripts/nearest_peak_to_deg.py`** | *[script]* | Legacy BETA-run helper that finds the nearest hg18 TRF1 peak center to each DEG TSS. |
| **`scripts/ttaggg_motif_batch.py`** | *[script]* | Batched Ensembl REST scan of DEG promoter windows for TTAGGG and CCCTAA motif counts. |
| **`scripts/scan_ttaggg_promoters.py`** | *[script]* | Local FASTA-backed TTAGGG promoter motif scanner using pyfaidx and `data/gene_annotation.tsv`. |
| **`scripts/trf1_extended_audit.py`** | *[script]* | Extended QC summaries for legacy TRF1 peaks, promoter overlap, BETA distance, and motif-scan outputs. |
| **`scripts/analyze_trf1_binding_error_network.py`** | *[script]* | Checks whether siTRF1 DEGs implicate upstream TRF1-binding regulators or downstream failure markers. |
| **`scripts/analyze_new_terf1_gtrd.py`** | *[script]* | Compares MSigDB/GTRD `TERF1_TARGET_GENES` against the siTRF1 DEG set. |
| **`scripts/analyze_chipatlas_terf1.py`** | *[script]* | Analyzes ChIP-Atlas hg38 TERF1 target-gene tables against siTRF1 DEGs. |
| **`scripts/alphagenome_metadata_probe.py`** | *[script]* | Queries AlphaGenome human output metadata and records TERF1/TRF1 track availability without storing the API key. |
| **`scripts/alphagenome_bj_pilot.py`** | *[script]* | Runs a BJ-context AlphaGenome pilot over selected siTRF1 DEG loci. |
| **`scripts/alphagenome_bj_chipatlas_candidates.py`** | *[script]* | Runs AlphaGenome BJ-context annotation for ChIP-Atlas TERF1-overlap DEG candidates. |
| **`scripts/summarize_alphagenome_chipatlas_candidates.py`** | *[script]* | Collapses AlphaGenome candidate track outputs into per-gene summary rankings. |

### `tests/`

| Path | Type | Use |
|---|---|---|
| **`tests/test_chroms.py`** | *[script]* | Smoke tests for chromosome normalization, comma-separated chromosome parsing, and Ensembl gene ID cleanup. |
| **`tests/test_promoter.py`** | *[script]* | Smoke tests for plus-strand, minus-strand, zero-clamped, and known MT-ND1 promoter-window math. |

### `data/`

| Path | Type | Use |
|---|---|---|
| **`data/README.md`** | *[doc]* | Input provenance, column definitions, coordinate notes, and regeneration guidance. |
| **`data/GSM638201_TRF1-p0.05.bed.gz`** | *[input]* | Compressed SISSRs TRF1 ChIP-seq peak file parsed by `peaks.py`. |
| **`data/siTRF1_DEG_list_gene_name.txt`** | *[input]* | DEG gene list with `gene_id`, `gene_name`, and `regulation` columns used by `deg.py`. |
| **`data/gene_annotation.tsv`** | *[input]* | BED-style gene-level annotation used to compute TSS positions and promoter windows. |

### `docs/`

| Path | Type | Use |
|---|---|---|
| **`docs/trf1_promoter_overlap.md`** | *[doc]* | Detailed usage notes for the TRF1 promoter-overlap pipeline. |
| **`docs/gff_gtf_format.md`** | *[doc]* | Reference notes for GFF/GTF columns, attributes, and coordinate semantics. |
| **`docs/FINAL_REPORT.md`** | *[doc]* | Full siTRF1 and TRF1 ChIP-seq integration report in Markdown. |
| **`docs/FINAL_REPORT.html`** | *[generated]* | HTML rendering of `docs/FINAL_REPORT.md` for browser viewing or sharing. |
| **`docs/alphagenome_new_data_report.md`** | *[doc]* | Report on AlphaGenome, MSigDB/GTRD, and ChIP-Atlas follow-up evidence. |
| **`docs/alphagenome_todo.md`** | *[doc]* | Checklist and evidence notes for using AlphaGenome in this project. |
| **`docs/trf1_binding_error_interpretation.md`** | *[doc]* | Interpretation of upstream TRF1-binding regulator checks and downstream failure markers. |

### `legacy/upstream_extraction/`

| Path | Type | Use |
|---|---|---|
| **`legacy/upstream_extraction/`** | *[script]* | Archived workflow for extracting upstream FASTA windows per gene from GRCh38 Ensembl 108. |
| **`legacy/upstream_extraction/README.md`** | *[doc]* | Workflow overview, run commands, output conventions, and archive location. |
| **`legacy/upstream_extraction/extract_upstream.py`** | *[script]* | Extracts per-gene upstream FASTA windows and writes manifests/checksums. |
| **`legacy/upstream_extraction/find_where.py`** | *[script]* | Searches the Ensembl GTF by text, GTF column, or attribute filters. |
| **`legacy/upstream_extraction/release_notes.md`** | *[doc]* | Release notes for the upstream FASTA archive and validation counts. |
| **`legacy/upstream_extraction/verify_archive.md`** | *[doc]* | Verification procedure and expected SHA256/coordinate checks for the release archive. |

### `github_release_assets/`

| Path | Type | Use |
|---|---|---|
| **`github_release_assets/upstream_by_gene_GRCh38_Ensembl108_reference_orientation.tar.gz`** | *[generated]* | Git-LFS tarball containing per-gene upstream FASTA windows from the archived extraction workflow. |

### `results/` primary pipeline outputs

`results/` is gitignored. These files are expected outputs from the primary `trf1-promoter-overlap` pipeline when using the bundled inputs.

| Path | Type | Use |
|---|---|---|
| **`results/TRF1_cleaned_peaks.bed`** | *[output]* | Cleaned three-column BED coordinates derived from the SISSRs TRF1 peak file. |
| **`results/TRF1_promoter_overlap_1000.tsv`** | *[output]* | Per-overlap table for DEG promoter windows 1000 bp upstream of TSS. |
| **`results/TRF1_promoter_overlap_2000.tsv`** | *[output]* | Per-overlap table for DEG promoter windows 2000 bp upstream of TSS. |
| **`results/TRF1_promoter_overlap_3000.tsv`** | *[output]* | Per-overlap table for DEG promoter windows 3000 bp upstream of TSS. |
| **`results/TRF1_promoter_overlap_summary.tsv`** | *[output]* | Summary counts across the 1000, 2000, and 3000 bp promoter windows. |
| **`results/TRF1_matched_genes_barplot.png`** | *[output]* | Bar plot of matched DEG counts by upstream window size. |
| **`results/TRF1_matched_genes_stacked_barplot.png`** | *[output]* | Stacked UP/DOWN bar plot of matched DEG counts by upstream window size. |
| **`results/report.txt`** | *[output]* | Plain-text run report with command, inputs, counts, notes, and matched genes. |

### `results/integration_analysis/`

`results/integration_analysis/` is also gitignored. It contains local follow-up outputs from BETA, ChEA3, Enrichr, and TTAGGG motif scans.

| Path | Type | Use |
|---|---|---|
| **`results/integration_analysis/`** | *[output]* | Local integration-analysis output root. |
| **`results/integration_analysis/beta/`** | *[output]* | BETA input/output directory containing formatted DEG and peak files, BETA results, and DEG-nearest-peak distances. |
| **`results/integration_analysis/beta/TRF1_peaks.bed`** | *[output]* | BETA-ready BED peak file with peak IDs and fold values. |
| **`results/integration_analysis/beta/TRF1_peaks_full.bed`** | *[output]* | Extended BETA peak table retaining tag count, fold, and p-value fields. |
| **`results/integration_analysis/beta/siTRF1_DEGs.bsf`** | *[output]* | BETA BSF expression input with gene symbol, log2 fold change, and p-value. |
| **`results/integration_analysis/beta/deg_nearest_peak.tsv`** | *[output]* | Per-DEG nearest TRF1 peak-center distance table from `nearest_peak_to_deg.py`. |
| **`results/integration_analysis/beta/out_hg18/`** | *[output]* | BETA basic hg18 output directory. |
| **`results/integration_analysis/beta/out_hg18/siTRF1_hg18_function_prediction.R`** | *[output]* | BETA-generated R plotting script for hg18 function prediction. |
| **`results/integration_analysis/beta/out_hg18/siTRF1_hg18_function_prediction.pdf`** | *[output]* | BETA-generated hg18 function prediction PDF. |
| **`results/integration_analysis/beta/out_minus/`** | *[output]* | BETA minus-direction output directory. |
| **`results/integration_analysis/beta/out_minus/TRF1_minus_targets.txt`** | *[output]* | BETA minus target-gene summary output. |
| **`results/integration_analysis/beta/out_minus/TRF1_minus_targets_associated_peaks.txt`** | *[output]* | BETA minus associated-peak table with target gene, distance, and score fields. |
| **`results/integration_analysis/chea3/`** | *[output]* | ChEA3 transcription-factor enrichment output directory. |
| **`results/integration_analysis/chea3/Integrated_meanRank.tsv`** | *[output]* | ChEA3 integrated mean-rank TF enrichment table. |
| **`results/integration_analysis/enrichr/`** | *[output]* | Enrichr enrichment output directory. |
| **`results/integration_analysis/enrichr/GO_Biological_Process_2026_table.txt`** | *[output]* | Enrichr GO Biological Process 2026 enrichment table. |
| **`results/integration_analysis/enrichr/KEGG_2026_table.txt`** | *[output]* | Enrichr KEGG 2026 enrichment table. |
| **`results/integration_analysis/enrichr/MSigDB_Hallmark_2020_table.txt`** | *[output]* | Enrichr MSigDB Hallmark 2020 enrichment table. |
| **`results/integration_analysis/enrichr/Reactome_Pathways_2024_table.txt`** | *[output]* | Enrichr Reactome Pathways 2024 enrichment table. |
| **`results/integration_analysis/motif_scan/`** | *[output]* | TTAGGG motif-scan output directory. |
| **`results/integration_analysis/motif_scan/ttaggg_motif_scan_full.tsv`** | *[output]* | Per-DEG promoter motif-count table for TTAGGG, tandem TTAGGG, and CCCTAA. |
| **`results/integration_analysis/motif_scan/ttaggg_motif_summary_full.txt`** | *[output]* | Text summary of the full TTAGGG motif scan. |

### Other present `results/` subtrees

These gitignored local outputs are present in this checkout and document follow-up analyses beyond the primary pipeline.

| Path | Type | Use |
|---|---|---|
| **`results/extended_audit/`** | *[output]* | Extended QC output directory for legacy peaks and existing analysis summaries. |
| **`results/extended_audit/trf1_extended_audit_summary.md`** | *[output]* | Markdown summary of the extended TRF1 audit. |
| **`results/extended_audit/legacy_trf1_peaks_parsed.tsv`** | *[output]* | Parsed numeric table of legacy TRF1 peaks. |
| **`results/extended_audit/legacy_trf1_peaks_by_chrom.tsv`** | *[output]* | Chromosome-level peak count and score summary. |
| **`results/extended_audit/legacy_trf1_peak_threshold_grid.tsv`** | *[output]* | Peak-count grid across p-value and fold-change thresholds. |
| **`results/extended_audit/legacy_trf1_top_fold_peaks.tsv`** | *[output]* | Highest-fold TRF1 peak subset for QC review. |
| **`results/extended_audit/existing_promoter_overlap_summary.copy.tsv`** | *[output]* | Copied primary promoter-overlap summary used by the audit. |
| **`results/extended_audit/existing_ttaggg_motif_summary.tsv`** | *[output]* | Condensed motif-scan count summary used by the audit. |
| **`results/extended_audit/existing_beta_nearest_peak_with_bins.tsv`** | *[output]* | BETA nearest-peak distances annotated with distance bins. |
| **`results/extended_audit/existing_beta_nearest_peak_bins.tsv`** | *[output]* | Aggregated nearest-peak distance-bin counts. |
| **`results/new_data/`** | *[output]* | Follow-up output directory for newer TERF1/GTRD/MSigDB and ChIP-Atlas checks. |
| **`results/new_data/TERF1_TARGET_GENES.json`** | *[input]* | MSigDB/GTRD TERF1 target-gene JSON used by `analyze_new_terf1_gtrd.py`. |
| **`results/new_data/terf1_gtrd_deg_overlap.tsv`** | *[output]* | DEG overlap table for the TERF1 MSigDB/GTRD target-gene set. |
| **`results/new_data/terf1_gtrd_deg_overlap_summary.tsv`** | *[output]* | Machine-readable summary metrics for the TERF1 MSigDB/GTRD overlap. |
| **`results/new_data/terf1_gtrd_deg_overlap_summary.md`** | *[output]* | Markdown interpretation summary for the TERF1 MSigDB/GTRD overlap. |
| **`results/new_data/chipatlas/`** | *[output]* | ChIP-Atlas TERF1 target-gene input and overlap-output directory. |
| **`results/new_data/chipatlas/TERF1.1.tsv`** | *[input]* | ChIP-Atlas hg38 TERF1 target-gene table for the 1 kb TSS window. |
| **`results/new_data/chipatlas/TERF1.5.tsv`** | *[input]* | ChIP-Atlas hg38 TERF1 target-gene table for the 5 kb TSS window. |
| **`results/new_data/chipatlas/TERF1.10.tsv`** | *[input]* | ChIP-Atlas hg38 TERF1 target-gene table for the 10 kb TSS window. |
| **`results/new_data/chipatlas/SRX359962.metadata.json`** | *[input]* | ChIP-Atlas experiment metadata for one TERF1 ChIP-seq source. |
| **`results/new_data/chipatlas/SRX472275.metadata.json`** | *[input]* | ChIP-Atlas experiment metadata for one TERF1 ChIP-seq source. |
| **`results/new_data/chipatlas/SRX4957969.metadata.json`** | *[input]* | ChIP-Atlas experiment metadata for one TERF1/TRF1 ChIP-seq source. |
| **`results/new_data/chipatlas/SRX4957970.metadata.json`** | *[input]* | ChIP-Atlas experiment metadata for one TERF1/TRF1 ChIP-seq source. |
| **`results/new_data/chipatlas/chipatlas_terf1_deg_overlap.tsv`** | *[output]* | Per-gene overlap table between ChIP-Atlas TERF1 targets and siTRF1 DEGs. |
| **`results/new_data/chipatlas/chipatlas_terf1_deg_overlap_summary.tsv`** | *[output]* | Summary metrics for ChIP-Atlas TERF1 target overlap. |
| **`results/new_data/chipatlas/chipatlas_terf1_deg_overlap_summary.md`** | *[output]* | Markdown interpretation summary for ChIP-Atlas TERF1 overlap. |
| **`results/new_data/chipatlas/chipatlas_terf1_experiment_metadata.tsv`** | *[output]* | Flattened ChIP-Atlas experiment metadata table. |
| **`results/new_data/chipatlas/chipatlas_terf1_overlap_by_experiment.tsv`** | *[output]* | Overlap summary broken down by ChIP-Atlas experiment source. |
| **`results/new_data/chipatlas/chipatlas_terf1_score_cutoff_summary.tsv`** | *[output]* | ChIP-Atlas overlap summary across score thresholds. |
| **`results/alphagenome/`** | *[output]* | AlphaGenome output root for metadata, pilot, and candidate-locus analyses. |
| **`results/alphagenome/metadata/`** | *[output]* | AlphaGenome output-metadata probe directory. |
| **`results/alphagenome/metadata/alphagenome_human_output_metadata.tsv`** | *[output]* | Full AlphaGenome human output metadata table. |
| **`results/alphagenome/metadata/alphagenome_human_output_type_counts.tsv`** | *[output]* | Counts of AlphaGenome human metadata rows by output type. |
| **`results/alphagenome/metadata/alphagenome_chip_tf_factor_counts.tsv`** | *[output]* | AlphaGenome CHIP_TF transcription-factor count table. |
| **`results/alphagenome/metadata/alphagenome_terf1_trf1_candidate_tracks.tsv`** | *[output]* | AlphaGenome metadata rows matching TERF1/TRF1 search terms. |
| **`results/alphagenome/metadata/alphagenome_biosample_keyword_candidates.tsv`** | *[output]* | AlphaGenome biosample keyword-match summary. |
| **`results/alphagenome/metadata/alphagenome_metadata_probe_summary.md`** | *[output]* | Markdown summary of AlphaGenome metadata availability. |
| **`results/alphagenome/pilot/`** | *[output]* | AlphaGenome BJ-cell pilot output directory. |
| **`results/alphagenome/pilot/alphagenome_bj_pilot_loci.tsv`** | *[output]* | Selected pilot gene loci submitted to AlphaGenome. |
| **`results/alphagenome/pilot/alphagenome_bj_pilot_track_summary.tsv`** | *[output]* | Track-level AlphaGenome pilot signal summaries. |
| **`results/alphagenome/pilot/alphagenome_bj_pilot_compact_summary.tsv`** | *[output]* | Compact per-gene AlphaGenome pilot summary. |
| **`results/alphagenome/pilot/alphagenome_bj_pilot_summary.md`** | *[output]* | Markdown summary of the AlphaGenome BJ pilot. |
| **`results/alphagenome/chipatlas_candidates/`** | *[output]* | AlphaGenome output directory for ChIP-Atlas TERF1-overlap DEG candidates. |
| **`results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_loci.tsv`** | *[output]* | Candidate DEG loci submitted to AlphaGenome. |
| **`results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_track_summary.tsv`** | *[output]* | Track-level AlphaGenome summaries for candidate loci. |
| **`results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_compact_summary.tsv`** | *[output]* | Compact candidate-locus AlphaGenome summary table. |
| **`results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_gene_summary.tsv`** | *[output]* | Per-gene AlphaGenome summary merged with ChIP-Atlas evidence. |
| **`results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_gene_summary.md`** | *[output]* | Markdown ranked candidate-gene AlphaGenome summary. |
| **`results/alphagenome/chipatlas_candidates/alphagenome_bj_chipatlas_candidate_summary.md`** | *[output]* | Markdown summary of the candidate-locus AlphaGenome run. |
| **`results/trf1_binding_error/`** | *[output]* | Output directory for upstream-regulator and downstream-failure interpretation checks. |
| **`results/trf1_binding_error/trf1_binding_error_gene_catalog.tsv`** | *[output]* | Curated TRF1-binding regulator catalog used by the interpretation script. |
| **`results/trf1_binding_error/trf1_binding_error_upstream_deg_hits.tsv`** | *[output]* | DEG hits among candidate upstream TRF1-binding regulators. |
| **`results/trf1_binding_error/trf1_binding_error_downstream_deg_markers.tsv`** | *[output]* | DEG hits among downstream stress, DDR, IFN, apoptosis, and mitochondrial markers. |
| **`results/trf1_binding_error/trf1_binding_error_downstream_class_counts.tsv`** | *[output]* | Downstream marker counts grouped by biological class. |
| **`results/trf1_binding_error/trf1_binding_error_network_summary.md`** | *[output]* | Markdown summary of upstream-regulator and downstream-failure evidence. |

## Pipeline overview

The primary pipeline answers: *for each gene in a DEG list, does the TRF1 ChIP-seq peak set overlap a fixed-size window upstream of its TSS?*

```text
gene annotation -> TSS -> upstream promoter window -> TRF1 peak overlap
```

Three nested upstream window sizes are evaluated together: **1000 / 2000 / 3000 bp**.

## Quick start

### 1. Create the conda environment

`pybedtools` needs the `bedtools` binary, so a conda env is the easiest path:

```bash
conda env create -f environment.yml
conda activate Bioinfomatics
```

`environment.yml` already includes `pip install -e .`, so the
`trf1-promoter-overlap` command and `python -m trf1_overlap` are both
available immediately.

### 2. Run the analysis

```bash
trf1-promoter-overlap \
    --peak data/GSM638201_TRF1-p0.05.bed.gz \
    --deg data/siTRF1_DEG_list_gene_name.txt \
    --annotation data/gene_annotation.tsv \
    --outdir results
```

Or, equivalently:

```bash
make trf1
```

Outputs land in `results/`:

```text
results/
|-- TRF1_cleaned_peaks.bed
|-- TRF1_promoter_overlap_1000.tsv
|-- TRF1_promoter_overlap_2000.tsv
|-- TRF1_promoter_overlap_3000.tsv
|-- TRF1_promoter_overlap_summary.tsv
|-- TRF1_matched_genes_barplot.png
|-- TRF1_matched_genes_stacked_barplot.png
`-- report.txt
```

### 3. Nuclear-only sanity check

Mitochondrial DNA biology differs from nuclear promoters; to exclude `chrM`:

```bash
make trf1-no-chrM
```

## Inputs

| File | Purpose |
|---|---|
| `data/GSM638201_TRF1-p0.05.bed.gz` | TRF1 SISSRs ChIP-seq peaks with header, kept verbatim. |
| `data/siTRF1_DEG_list_gene_name.txt` | siTRF1 DEG list: `gene_id`, `gene_name`, `regulation` (`UP`/`DOWN`). |
| `data/gene_annotation.tsv` | BED-style gene annotation: `chr`, `gene_start`, `gene_end`, `gene_name`, `gene_id`, `strand`. |

See [data/README.md](data/README.md) for column specs and provenance.

If you only have the GTF, pass `--gtf path/to/Homo_sapiens.GRCh38.108.gtf`
and the annotation TSV will be generated automatically and reused on
subsequent runs.

## Coordinate conventions

- **Annotation TSV** uses BED-style **0-based half-open** intervals:
  `gene_start = GTF_start - 1`, `gene_end = GTF_end`.
- **Promoter windows** are also BED 0-based half-open:
  - `+` strand: `[TSS - size, TSS)` where `TSS = gene_start`.
  - `-` strand: `[TSS, TSS + size)` where `TSS = gene_end` (the BED-exclusive end).
- **Chromosome names** are normalized to UCSC style (`chr1`, `chrM`, ...);
  Ensembl `MT` and `chrMT` both collapse to `chrM`.

## Documentation

- [TRF1 promoter overlap usage](docs/trf1_promoter_overlap.md) - detailed pipeline notes.
- [GFF / GTF format reference](docs/gff_gtf_format.md) - column semantics for both formats.
- [Final integration report](docs/FINAL_REPORT.md) - complete siTRF1 and TRF1 ChIP-seq integration report.
- [AlphaGenome and newer data report](docs/alphagenome_new_data_report.md) - AlphaGenome, MSigDB/GTRD, and ChIP-Atlas follow-up.
- [TRF1 binding error interpretation](docs/trf1_binding_error_interpretation.md) - upstream-regulator and downstream-failure interpretation.

## Tests

Pure-logic modules (chromosome normalization, promoter math) ship with smoke tests:

```bash
conda run -n Bioinfomatics python -m pytest tests/
```

## Legacy workflows

The earlier per-gene upstream-FASTA extraction pipeline is preserved under
[`legacy/upstream_extraction/`](legacy/upstream_extraction/README.md), along
with the published archive's SHA256 and verification procedure.

## License

[MIT](LICENSE).
