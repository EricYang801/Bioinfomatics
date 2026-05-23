"""pybedtools wrapper for promoter ∩ peak intersection."""

from __future__ import annotations

import shutil

import pandas as pd

from .peaks import PEAK_BED_COLUMNS
from .promoter import PROMOTER_BED_COLUMNS


def require_pybedtools() -> str:
    if shutil.which("bedtools") is None:
        raise RuntimeError(
            "bedtools executable was not found. Install it in the active environment, "
            "for example: conda install -c bioconda bedtools"
        )
    try:
        import pybedtools  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "pybedtools Python package was not found. Install it in the active "
            "environment, for example: conda install -c bioconda pybedtools"
        ) from None
    return "pybedtools + bedtools"


def build_peak_bedtool(peaks: pd.DataFrame):
    from pybedtools import BedTool

    return BedTool.from_dataframe(peaks[PEAK_BED_COLUMNS])


def intersect_with_pybedtools(promoters: pd.DataFrame, peak_bt) -> pd.DataFrame:
    from pybedtools import BedTool

    promoter_bt = BedTool.from_dataframe(promoters[PROMOTER_BED_COLUMNS])
    result_bt = promoter_bt.intersect(peak_bt, wa=True, wb=True)
    rows = [feature.fields for feature in result_bt]
    names = PROMOTER_BED_COLUMNS + PEAK_BED_COLUMNS
    if not rows:
        return pd.DataFrame(columns=names)
    result = pd.DataFrame(rows, columns=names)
    for col in ["promoter_start", "promoter_end", "TSS", "peak_start", "peak_end"]:
        result[col] = pd.to_numeric(result[col], errors="raise").astype(int)
    return result
