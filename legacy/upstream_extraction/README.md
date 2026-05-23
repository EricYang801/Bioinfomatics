# Upstream sequence extraction (archived)

Earlier workflow for extracting per-gene upstream FASTA windows from the
GRCh38 Ensembl 108 reference. Kept here for provenance; the active pipeline
is the TRF1 promoter overlap analysis at the repo root.

## Contents

| File | Purpose |
|---|---|
| `extract_upstream.py` | Extract `N` bp upstream of each gene's TSS into per-gene FASTA files, then verify each output against the genome with SHA256. |
| `find_where.py` | Interactive / argv-driven grep over `Homo_sapiens.GRCh38.108.gtf` by text, GTF column, or attribute. |
| `release_notes.md` | Description of the published archive (`github_release_assets/upstream_by_gene_GRCh38_Ensembl108_reference_orientation.tar.gz`). |
| `verify_archive.md` | Procedure for re-verifying the archive against the original FASTA + GTF. |

The Git-LFS-tracked archive itself lives at the repo root in
`github_release_assets/`, alongside its SHA256 sidecar.

## Running `extract_upstream.py`

```bash
# A single gene
python legacy/upstream_extraction/extract_upstream.py \
    --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.release108.full.fa \
    --gtf   Homo_sapiens.GRCh38.108.gtf \
    --outdir upstream_by_gene \
    --gene ATAD3B

# All genes
python legacy/upstream_extraction/extract_upstream.py \
    --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.release108.full.fa \
    --gtf   Homo_sapiens.GRCh38.108.gtf \
    --outdir upstream_by_gene \
    --all
```

Output convention:

- `+` strand: extract `[gene_start - N, gene_start - 1]` (1-based inclusive).
- `-` strand: extract `[gene_end + 1, gene_end + N]` (1-based inclusive).
- Sequences are stored in **reference genome orientation** (not reverse-complemented).

Each run writes `manifest.tsv`, `files.sha256.tsv`, `manifest.tsv.sha256`,
and `summary.tsv` next to the per-gene FASTA files.

## Running `find_where.py`

```bash
# CLI form
python legacy/upstream_extraction/find_where.py --attr gene_name=ATAD3B --where feature=gene
# Interactive form (no args)
python legacy/upstream_extraction/find_where.py
```

The script defaults to looking for `Homo_sapiens.GRCh38.108.gtf` in the
current working directory; pass `--gtf path/to/other.gtf` to override.
