"""TSS-anchored upstream promoter window construction.

For the plus strand we use ``[TSS - size, TSS)``; for the minus strand we use
``[TSS, TSS + size)``. Both are BED-style 0-based half-open intervals, so
they intersect directly with peaks in the same coordinate system.
"""

from __future__ import annotations

import pandas as pd


PROMOTER_SIZES = (1000, 2000, 3000)
PROMOTER_BED_COLUMNS = [
    "chr",
    "promoter_start",
    "promoter_end",
    "gene_name",
    "gene_id",
    "regulation",
    "TSS",
    "strand",
    "_deg_key",
]


def exclude_chromosomes(
    peaks: pd.DataFrame,
    annotation: pd.DataFrame,
    excluded_chroms: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, str | None]:
    if not excluded_chroms:
        return peaks, annotation, None

    peak_count_before = len(peaks)
    annotation_count_before = len(annotation)
    peaks = peaks[~peaks["peak_chr"].isin(excluded_chroms)].copy()
    annotation = annotation[~annotation["chr"].isin(excluded_chroms)].copy()
    message = (
        "Excluded chromosomes: "
        f"{', '.join(sorted(excluded_chroms))}; "
        f"removed {peak_count_before - len(peaks)} peaks and "
        f"{annotation_count_before - len(annotation)} annotation genes."
    )
    return peaks, annotation, message


def build_promoters(target_genes: pd.DataFrame, size: int) -> pd.DataFrame:
    promoters = target_genes.copy()
    promoters["TSS"] = promoters.apply(
        lambda row: int(row["gene_start"]) if row["strand"] == "+" else int(row["gene_end"]),
        axis=1,
    )
    promoters["promoter_start"] = promoters.apply(
        lambda row: max(0, int(row["TSS"]) - size)
        if row["strand"] == "+"
        else int(row["TSS"]),
        axis=1,
    )
    promoters["promoter_end"] = promoters.apply(
        lambda row: int(row["TSS"])
        if row["strand"] == "+"
        else int(row["TSS"]) + size,
        axis=1,
    )
    promoters = promoters[promoters["promoter_end"] > promoters["promoter_start"]].copy()
    return promoters[PROMOTER_BED_COLUMNS].reset_index(drop=True)
