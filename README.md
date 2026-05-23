# Bioinfomatics

Bioinformatics pipelines centered on **TRF1 ChIP-seq promoter overlap with
DEG gene lists**, plus archived utilities for upstream-sequence extraction.

The primary pipeline answers: *for each gene in a DEG list, does the TRF1
ChIP-seq peak set overlap a fixed-size window upstream of its TSS?*

```text
gene annotation -> TSS -> upstream promoter window -> TRF1 peak overlap
```

Three nested upstream window sizes are evaluated together: **1000 / 2000 / 3000 bp**.

## Repository layout

```text
.
├── src/trf1_overlap/        # importable Python package (the pipeline)
├── scripts/                 # thin CLI wrapper for src/
├── tests/                   # pytest smoke tests for pure-logic modules
├── data/                    # small, version-controlled inputs (peak, DEG, annotation)
├── docs/                    # usage docs and format references
├── legacy/                  # archived workflows (upstream sequence extraction)
├── github_release_assets/   # Git-LFS-tracked tarball + SHA256 sidecar
├── environment.yml          # conda env spec (bedtools + python deps)
├── pyproject.toml           # package metadata, entry point
├── Makefile                 # `make trf1`, `make trf1-no-chrM`, `make clean`
└── results/                 # pipeline output (gitignored)
```

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
├── TRF1_cleaned_peaks.bed
├── TRF1_promoter_overlap_1000.tsv
├── TRF1_promoter_overlap_2000.tsv
├── TRF1_promoter_overlap_3000.tsv
├── TRF1_promoter_overlap_summary.tsv
├── TRF1_matched_genes_barplot.png
├── TRF1_matched_genes_stacked_barplot.png
└── report.txt
```

### 3. Nuclear-only sanity check

Mitochondrial DNA biology differs from nuclear promoters; to exclude `chrM`:

```bash
make trf1-no-chrM
```

## Inputs

| File | Purpose |
|---|---|
| `data/GSM638201_TRF1-p0.05.bed.gz` | TRF1 SISSRs ChIP-seq peaks (with header, kept verbatim) |
| `data/siTRF1_DEG_list_gene_name.txt` | siTRF1 DEG list: `gene_id`, `gene_name`, `regulation` (UP/DOWN) |
| `data/gene_annotation.tsv` | BED-style gene annotation: `chr`, `gene_start`, `gene_end`, `gene_name`, `gene_id`, `strand` |

See [data/README.md](data/README.md) for column specs and provenance.

If you only have the GTF, pass `--gtf path/to/Homo_sapiens.GRCh38.108.gtf`
and the annotation TSV will be generated automatically (and reused on
subsequent runs).

## Coordinate conventions

- **Annotation TSV** uses BED-style **0-based half-open** intervals:
  `gene_start = GTF_start - 1`, `gene_end = GTF_end`.
- **Promoter windows** are also BED 0-based half-open:
  - `+` strand: `[TSS - size, TSS)` where `TSS = gene_start`.
  - `-` strand: `[TSS, TSS + size)` where `TSS = gene_end` (the BED-exclusive end).
- **Chromosome names** are normalized to UCSC style (`chr1`, `chrM`, …);
  Ensembl `MT` and `chrMT` both collapse to `chrM`.

## Documentation

- [TRF1 promoter overlap usage](docs/trf1_promoter_overlap.md) — detailed pipeline notes
- [GFF / GTF format reference](docs/gff_gtf_format.md) — column semantics for both formats

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
