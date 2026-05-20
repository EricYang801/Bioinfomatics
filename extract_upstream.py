#!/usr/bin/env python3
import argparse
import hashlib
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from pyfaidx import Fasta


DEFAULT_FASTA = Path("Homo_sapiens.GRCh38.dna_sm.primary_assembly.release108.full.fa")
DEFAULT_GTF = Path("Homo_sapiens.GRCh38.108.gtf")
DEFAULT_LENGTHS = (1000, 2000, 3000)
WRAP_WIDTH = 60
CoordinateField = Union[int, str, Path]


@dataclass(frozen=True)
class Gene:
    seqname: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_version: str
    gene_name: str
    gene_biotype: str


def parse_attributes(attributes_text):
    attributes = {}
    for part in attributes_text.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        attributes[key] = value.strip().strip('"')
    return attributes


def read_genes(gtf_path):
    with gtf_path.open() as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) != 9 or fields[2] != "gene":
                continue

            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            yield Gene(
                seqname=fields[0],
                start=int(fields[3]),
                end=int(fields[4]),
                strand=fields[6],
                gene_id=gene_id,
                gene_version=attrs.get("gene_version", ""),
                gene_name=attrs.get("gene_name", gene_id),
                gene_biotype=attrs.get("gene_biotype", ""),
            )


def choose_genes(genes, target_gene):
    if target_gene is None:
        return list(genes)

    selected = [
        gene for gene in genes
        if gene.gene_name == target_gene or gene.gene_id == target_gene
    ]
    if not selected:
        raise SystemExit(f"No gene matched gene_name or gene_id: {target_gene}")
    return selected


def upstream_region(gene, upstream_length, chrom_length):
    if gene.strand == "+":
        start = max(1, gene.start - upstream_length)
        end = gene.start - 1
    elif gene.strand == "-":
        start = gene.end + 1
        end = min(chrom_length, gene.end + upstream_length)
    else:
        return None

    if start > end:
        return None
    return start, end


def fetch_upstream_sequence(genome, gene, region_start, region_end):
    return str(genome[gene.seqname][region_start - 1:region_end]).upper()


def safe_name(text):
    text = text or "unnamed"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text)


def gene_file_path(run_dir, gene, upstream_length):
    gene_label = safe_name(gene.gene_name)
    gene_id = safe_name(gene.gene_id)
    filename = f"{gene_label}__{gene_id}__upstream_{upstream_length}bp.fa"
    return run_dir / f"{upstream_length}bp" / strand_dir(gene.strand) / filename


def strand_dir(strand):
    if strand == "+":
        return "plus"
    if strand == "-":
        return "minus"
    return "unknown_strand"


def fasta_header(gene, upstream_length, region_start, region_end, actual_length):
    return (
        f">{gene.gene_name}|gene_id={gene.gene_id}|gene_version={gene.gene_version}|"
        f"gene_biotype={gene.gene_biotype}|seqname={gene.seqname}|"
        f"gene={gene.start}-{gene.end}|strand={gene.strand}|"
        f"upstream={upstream_length}|genomic_region={region_start}-{region_end}|"
        f"length={actual_length}|sequence_orientation=reference_genome"
    )


def write_fasta(path, header, sequence):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(header + "\n")
        for index in range(0, len(sequence), WRAP_WIDTH):
            handle.write(sequence[index:index + WRAP_WIDTH] + "\n")


def read_fasta_sequence(path):
    sequence_parts = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            sequence_parts.append(line.strip())
    return "".join(sequence_parts).upper()


def sha256_text(text):
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_lengths(text):
    lengths = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        value = int(item)
        if value <= 0:
            raise argparse.ArgumentTypeError("Lengths must be positive integers.")
        lengths.append(value)
    if not lengths:
        raise argparse.ArgumentTypeError("At least one length is required.")
    return tuple(lengths)


def manifest_header():
    return [
        "gene_id",
        "gene_version",
        "gene_name",
        "gene_biotype",
        "seqname",
        "strand",
        "gene_start",
        "gene_end",
        "requested_upstream_bp",
        "upstream_genomic_start",
        "upstream_genomic_end",
        "actual_length",
        "output_file",
        "expected_sequence_sha256",
        "output_sequence_sha256",
        "file_sha256",
        "verify_status",
    ]


def manifest_row(
    gene: Gene,
    upstream_length: int,
    upstream_genomic_start: CoordinateField = ".",
    upstream_genomic_end: CoordinateField = ".",
    actual_length: int = 0,
    output_file: CoordinateField = ".",
    expected_sha: str = ".",
    output_sha: str = ".",
    file_sha: str = ".",
    verify_status: str = ".",
):
    return [
        gene.gene_id,
        gene.gene_version,
        gene.gene_name,
        gene.gene_biotype,
        gene.seqname,
        gene.strand,
        str(gene.start),
        str(gene.end),
        str(upstream_length),
        str(upstream_genomic_start),
        str(upstream_genomic_end),
        str(actual_length),
        str(output_file),
        expected_sha,
        output_sha,
        file_sha,
        verify_status,
    ]


def write_tsv(path, rows):
    with path.open("w") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")


def write_summary(path, stats, manifest_path, manifest_sha):
    lines = [
        f"genes_requested\t{stats['genes_requested']}",
        f"genes_with_reference_sequence\t{stats['genes_with_reference_sequence']}",
        f"genes_missing_reference_sequence\t{stats['genes_missing_reference_sequence']}",
        f"fasta_files_written\t{stats['fasta_files_written']}",
        f"verified_ok\t{stats['verified_ok']}",
        f"verification_failed\t{stats['verification_failed']}",
        f"skipped_intervals\t{stats['skipped_intervals']}",
        f"manifest\t{manifest_path}",
        f"manifest_sha256\t{manifest_sha}",
    ]
    path.write_text("\n".join(lines) + "\n")


def run_name(target_gene):
    return safe_name(target_gene) if target_gene else "all_genes"


def extract_and_verify(fasta_path, gtf_path, outdir, lengths, target_gene, overwrite):
    genes = choose_genes(read_genes(gtf_path), target_gene)
    run_dir = outdir / run_name(target_gene)

    if run_dir.exists():
        if not overwrite:
            raise SystemExit(
                f"Output directory already exists: {run_dir}. "
                "Use --overwrite to replace it."
            )
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    genome = Fasta(str(fasta_path), sequence_always_upper=True)
    reference_names = set(genome.keys())

    manifest_rows = [manifest_header()]
    file_hash_rows = [["file_sha256", "output_file"]]
    stats = {
        "genes_requested": len(genes),
        "genes_with_reference_sequence": 0,
        "genes_missing_reference_sequence": 0,
        "fasta_files_written": 0,
        "verified_ok": 0,
        "verification_failed": 0,
        "skipped_intervals": 0,
    }

    for length in lengths:
        for folder in ("plus", "minus", "unknown_strand"):
            (run_dir / f"{length}bp" / folder).mkdir(
                parents=True,
                exist_ok=True,
            )

    for gene in genes:
        if gene.seqname not in reference_names:
            stats["genes_missing_reference_sequence"] += 1
            for length in lengths:
                manifest_rows.append(
                    manifest_row(
                        gene,
                        length,
                        verify_status="missing_reference_sequence",
                    )
                )
            continue

        stats["genes_with_reference_sequence"] += 1
        chrom_length = len(genome[gene.seqname])

        for length in lengths:
            region = upstream_region(gene, length, chrom_length)
            if region is None:
                stats["skipped_intervals"] += 1
                status = (
                    "unsupported_strand"
                    if gene.strand not in {"+", "-"}
                    else "empty_upstream_interval"
                )
                manifest_rows.append(
                    manifest_row(gene, length, verify_status=status)
                )
                continue

            region_start, region_end = region
            expected_sequence = fetch_upstream_sequence(
                genome,
                gene,
                region_start,
                region_end,
            )
            expected_sha = sha256_text(expected_sequence)
            output_path = gene_file_path(run_dir, gene, length)
            relative_output_path = output_path.relative_to(run_dir)

            write_fasta(
                output_path,
                fasta_header(
                    gene,
                    length,
                    region_start,
                    region_end,
                    len(expected_sequence),
                ),
                expected_sequence,
            )
            stats["fasta_files_written"] += 1

            observed_sequence = read_fasta_sequence(output_path)
            output_sha = sha256_text(observed_sequence)
            file_sha = sha256_file(output_path)
            verify_status = (
                "ok"
                if (
                    observed_sequence == expected_sequence
                    and output_sha == expected_sha
                    and len(observed_sequence) == len(expected_sequence)
                )
                else "mismatch"
            )
            if verify_status == "ok":
                stats["verified_ok"] += 1
            else:
                stats["verification_failed"] += 1

            manifest_rows.append(
                manifest_row(
                    gene=gene,
                    upstream_length=length,
                    upstream_genomic_start=region_start,
                    upstream_genomic_end=region_end,
                    actual_length=len(expected_sequence),
                    output_file=relative_output_path,
                    expected_sha=expected_sha,
                    output_sha=output_sha,
                    file_sha=file_sha,
                    verify_status=verify_status,
                )
            )
            file_hash_rows.append([file_sha, str(relative_output_path)])

    manifest_path = run_dir / "manifest.tsv"
    files_sha_path = run_dir / "files.sha256.tsv"
    write_tsv(manifest_path, manifest_rows)
    write_tsv(files_sha_path, file_hash_rows)

    manifest_sha = sha256_file(manifest_path)
    (run_dir / "manifest.tsv.sha256").write_text(
        f"{manifest_sha}  manifest.tsv\n"
    )
    write_summary(run_dir / "summary.tsv", stats, manifest_path, manifest_sha)

    if stats["verification_failed"]:
        raise SystemExit(
            f"{stats['verification_failed']} output files failed verification. "
            f"See {manifest_path}."
        )

    return run_dir, stats, manifest_path, manifest_sha


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract one FASTA per gene per upstream length from a GTF, "
            "then verify each output against the original genome FASTA."
        )
    )
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--outdir", type=Path, default=Path("upstream_by_gene"))
    parser.add_argument(
        "--lengths",
        type=parse_lengths,
        default=DEFAULT_LENGTHS,
        help="Comma-separated upstream lengths. Default: 1000,2000,3000",
    )
    parser.add_argument(
        "--gene",
        help="Extract only this gene_name or gene_id.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Extract all gene features from the GTF.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace this run's existing output directory.",
    )

    args = parser.parse_args()
    if bool(args.gene) == bool(args.all):
        raise SystemExit("Use exactly one of --gene GENE or --all.")

    run_dir, stats, manifest_path, manifest_sha = extract_and_verify(
        fasta_path=args.fasta,
        gtf_path=args.gtf,
        outdir=args.outdir,
        lengths=args.lengths,
        target_gene=args.gene,
        overwrite=args.overwrite,
    )

    print(f"Output directory: {run_dir}")
    print(f"Genes requested: {stats['genes_requested']}")
    print(f"Genes with reference sequence: {stats['genes_with_reference_sequence']}")
    print(f"Genes missing reference sequence: {stats['genes_missing_reference_sequence']}")
    print(f"FASTA files written: {stats['fasta_files_written']}")
    print(f"Verified OK: {stats['verified_ok']}")
    print(f"Verification failed: {stats['verification_failed']}")
    print(f"Manifest: {manifest_path}")
    print(f"Manifest SHA256: {manifest_sha}")


if __name__ == "__main__":
    main()
