"""GTF parsing and BED-style gene annotation TSV reading."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .chroms import gene_id_base, normalize_chr


ANNOTATION_COLUMNS = [
    "chr",
    "gene_start",
    "gene_end",
    "gene_name",
    "gene_id",
    "strand",
]
REQUIRED_ANNOTATION_COLUMNS = set(ANNOTATION_COLUMNS)


def parse_gtf_attributes(attributes_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in attributes_text.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        attributes[key] = value.strip().strip('"')
    return attributes


def annotation_from_gtf(gtf_path: Path) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    with gtf_path.open() as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) != 9 or fields[2] != "gene":
                continue

            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)

            # GTF is 1-based closed; BED-style interval analysis uses 0-based start.
            gene_start = max(0, int(fields[3]) - 1)
            gene_end = int(fields[4])

            records.append(
                {
                    "chr": normalize_chr(fields[0]),
                    "gene_start": gene_start,
                    "gene_end": gene_end,
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "strand": fields[6],
                }
            )

    if not records:
        raise ValueError(f"No gene records were parsed from GTF: {gtf_path}")
    return pd.DataFrame.from_records(records)


def read_annotation(
    annotation_text: str, gtf_text: str | None
) -> tuple[pd.DataFrame, list[str], str]:
    annotation_path = Path(annotation_text)
    messages: list[str] = []

    if annotation_path.exists():
        annotation = pd.read_csv(annotation_path, sep="\t", dtype=str)
        source = f"TSV: {annotation_path}"
    else:
        if not gtf_text:
            raise FileNotFoundError(
                f"Annotation TSV not found: {annotation_path}; provide --gtf to create it."
            )
        gtf_path = Path(gtf_text)
        if not gtf_path.exists():
            raise FileNotFoundError(f"GTF file not found: {gtf_path}")
        annotation = annotation_from_gtf(gtf_path)
        annotation_path.parent.mkdir(parents=True, exist_ok=True)
        annotation.to_csv(annotation_path, sep="\t", index=False)
        source = f"GTF converted to BED-style TSV: {gtf_path}"
        messages.append(
            f"Annotation TSV {annotation_path} was missing; generated it from {gtf_path}."
        )

    annotation.columns = [str(col).strip() for col in annotation.columns]
    missing = REQUIRED_ANNOTATION_COLUMNS.difference(annotation.columns)
    if missing:
        raise ValueError(
            "Annotation is missing required columns: " + ", ".join(sorted(missing))
        )

    annotation = annotation[ANNOTATION_COLUMNS].copy()
    annotation["chr"] = annotation["chr"].map(normalize_chr)
    annotation["gene_start"] = pd.to_numeric(annotation["gene_start"], errors="raise").astype(int)
    annotation["gene_end"] = pd.to_numeric(annotation["gene_end"], errors="raise").astype(int)
    annotation["gene_name"] = annotation["gene_name"].astype(str)
    annotation["gene_id"] = annotation["gene_id"].astype(str)
    annotation["strand"] = annotation["strand"].astype(str).str.strip()

    invalid_strands = sorted(set(annotation["strand"]) - {"+", "-"})
    if invalid_strands:
        raise ValueError(f"Unsupported strand value(s): {', '.join(invalid_strands)}")

    annotation["_gene_id_base"] = annotation["gene_id"].map(gene_id_base)
    return annotation, messages, source
