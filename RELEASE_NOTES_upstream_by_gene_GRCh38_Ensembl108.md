# Upstream gene sequences, GRCh38 Ensembl release 108

Archive:
`upstream_by_gene_GRCh38_Ensembl108_reference_orientation.tar.gz`

SHA256:
`fc6c81e87154b3c679ee4b0b7dc1da627d3e1704fdb46df2a86d9ed90c063fb9`

Contents:
- `upstream_by_gene/all_genes/1000bp/{plus,minus}/`
- `upstream_by_gene/all_genes/2000bp/{plus,minus}/`
- `upstream_by_gene/all_genes/3000bp/{plus,minus}/`
- `upstream_by_gene/all_genes/manifest.tsv`
- `upstream_by_gene/all_genes/files.sha256.tsv`
- `upstream_by_gene/all_genes/summary.tsv`

Reference:
- Genome FASTA: `Homo_sapiens.GRCh38.dna_sm.primary_assembly.release108.full.fa`
- Annotation: `Homo_sapiens.GRCh38.108.gtf`
- Feature type: `gene`

Output rules:
- One FASTA file per gene per upstream length.
- Upstream lengths: 1000 bp, 2000 bp, 3000 bp.
- `+` strand genes use `gene_start - N` to `gene_start - 1`.
- `-` strand genes use `gene_end + 1` to `gene_end + N`.
- Sequences are stored in reference genome orientation, not reverse-complemented.
- Results are separated into `plus/` and `minus/` folders.

Validation summary:
- Genes requested: 62,703
- FASTA files written: 188,109
- Verified OK against the original reference FASTA: 188,109
- Verification failed: 0
- Extra ATAD3B FASTA files included: 3
- Archive metadata entries (`._*`): 0
