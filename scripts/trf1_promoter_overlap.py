#!/usr/bin/env python3
"""CLI wrapper for the ``trf1_overlap`` package.

This script supports two usage modes:

1. After ``pip install -e .`` — runs the installed package.
2. Direct invocation from a checkout — prepends ``src/`` to ``sys.path``
   so the package resolves without installation.
"""

from __future__ import annotations

import sys
from pathlib import Path


def _bootstrap_src_path() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    src_dir = repo_root / "src"
    if src_dir.is_dir() and str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))


_bootstrap_src_path()

from trf1_overlap.cli import main  # noqa: E402  (path setup must run first)

if __name__ == "__main__":
    raise SystemExit(main())
