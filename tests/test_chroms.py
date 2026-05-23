"""Smoke tests for chromosome name normalization and Ensembl gene-id helpers."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from trf1_overlap.chroms import gene_id_base, normalize_chr, parse_chrom_list


def test_normalize_chr_ensembl_to_ucsc():
    assert normalize_chr("1") == "chr1"
    assert normalize_chr("X") == "chrX"
    assert normalize_chr("MT") == "chrM"
    assert normalize_chr("M") == "chrM"


def test_normalize_chr_ucsc_passthrough_and_mt_collapse():
    assert normalize_chr("chr1") == "chr1"
    assert normalize_chr("chrM") == "chrM"
    assert normalize_chr("chrMT") == "chrM"
    assert normalize_chr("chrMt") == "chrM"


def test_parse_chrom_list():
    assert parse_chrom_list("") == set()
    assert parse_chrom_list("chrM") == {"chrM"}
    assert parse_chrom_list("MT,Y") == {"chrM", "chrY"}
    assert parse_chrom_list(" chr1 , 2 ,") == {"chr1", "chr2"}


def test_gene_id_base_strips_version():
    assert gene_id_base("ENSG00000198888") == "ENSG00000198888"
    assert gene_id_base("ENSG00000198888.2") == "ENSG00000198888"
    assert gene_id_base("  ENSG00000198888.10  ") == "ENSG00000198888"
