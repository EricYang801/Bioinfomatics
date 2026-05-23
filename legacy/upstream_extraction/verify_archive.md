# Verify uploaded upstream FASTA archive

This checks the uploaded Git LFS archive against the original Ensembl release 108 FASTA and GTF.

Required local files:
- `github_release_assets/upstream_by_gene_GRCh38_Ensembl108_reference_orientation.tar.gz`
- `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`
- `Homo_sapiens.GRCh38.108.gtf`

Official sources:
- FASTA: `https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`
- GTF: `https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz`

Expected archive SHA256:

```text
fc6c81e87154b3c679ee4b0b7dc1da627d3e1704fdb46df2a86d9ed90c063fb9
```

Verification result from a clean archive:

```text
all_genes_manifest_rows              188109
atad3b_manifest_rows                 3
fasta_members_checked                188112
sequence_matches_original_reference  188112
sequence_mismatches                  0
missing_manifest_rows                0
missing_fasta_members                0
metadata_entries                     0
gtf_coordinate_mismatches            0
region_coordinate_mismatches         0
sha256_mismatches                    0
```
