"""Command-line entry point for the human proteome/isoform dataset rebuild."""

from __future__ import annotations

import os
import sys


PIPELINE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rbp_pipeline")
if PIPELINE_DIR not in sys.path:
    sys.path.insert(0, PIPELINE_DIR)

from rebuild_from_scratch import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
