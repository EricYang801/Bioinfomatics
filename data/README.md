# Inputs

Small, version-controlled inputs for the TRF1 promoter overlap pipeline.

## `GSM638201_TRF1-p0.05.bed.gz`

TRF1 ChIP-seq peaks called by SISSRs at p<0.05.

- **Source**: GEO accession [GSM638201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM638201).
- **Format**: SISSRs output with a verbose multi-line header. After the
  literal column header line
  `Chr\tcStart\tcEnd\tNumTags\tFold\tp-value`, each row is a peak.
- **Genome**: GRCh38 (UCSC-style chromosome names: `chr1`, `chrM`, …).
- **Loader**: `src/trf1_overlap/peaks.py` skips the header preamble and
  rejects rows where `peak_end <= peak_start`.

## `siTRF1_DEG_list_gene_name.txt`

DEG list from a siRNA knockdown experiment targeting TRF1 / TERF1.

- **Columns** (tab-separated):

  | Column | Type | Notes |
  |---|---|---|
  | `gene_id` | str | Ensembl gene id (versioned suffixes stripped at parse time) |
  | `gene_name` | str | HGNC symbol; used as fallback when `gene_id` does not match annotation |
  | `regulation` | enum | `UP` or `DOWN`; any other value is a hard error |

- **Loader**: `src/trf1_overlap/deg.py`.

## `gene_annotation.tsv`

BED-style gene-level annotation derived from
`Homo_sapiens.GRCh38.108.gtf` (Ensembl release 108).

- **Columns**: `chr`, `gene_start`, `gene_end`, `gene_name`, `gene_id`, `strand`.
- **Coordinates**: 0-based half-open. `gene_start = GTF_start - 1`, `gene_end = GTF_end`.
- **Chromosome names**: Ensembl names (`1`, `MT`, …) are rewritten to UCSC
  style (`chr1`, `chrM`, …); `chrMT` collapses to `chrM`.
- **Regeneration**: if this file is removed, pass `--gtf path/to/Homo_sapiens.GRCh38.108.gtf`
  to the CLI and it will be rebuilt and cached here.

The GTF itself is **not** committed (it is ~1.4 GB); fetch it from
[Ensembl release 108](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz).
