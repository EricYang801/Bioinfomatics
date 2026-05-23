"""Chromosome name normalization and Ensembl gene-id helpers."""

from __future__ import annotations

import re


def normalize_chr(chrom: object) -> str:
    text = str(chrom).strip()
    if not text:
        return text
    if text.startswith("chr"):
        return "chrM" if text in {"chrMT", "chrMt"} else text
    if text in {"MT", "Mt", "M"}:
        return "chrM"
    return f"chr{text}"


def parse_chrom_list(chroms_text: str) -> set[str]:
    if not chroms_text:
        return set()
    return {
        normalize_chr(chrom)
        for chrom in chroms_text.split(",")
        if chrom.strip()
    }


def gene_id_base(gene_id: object) -> str:
    text = str(gene_id).strip()
    return re.sub(r"\.\d+$", "", text)
