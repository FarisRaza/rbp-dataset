"""Shared command-line runner for the small per-family QC entry points."""

from __future__ import annotations

import argparse
import os
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
REPOSITORY_ROOT = os.path.dirname(HERE)
PIPELINE_DIR = os.path.join(REPOSITORY_ROOT, "rbp_pipeline")
if PIPELINE_DIR not in sys.path:
    sys.path.insert(0, PIPELINE_DIR)

import generate_feature_reports
import rebuild_from_scratch


def parser(family):
    command = argparse.ArgumentParser(
        description=(
            f"Generate exploratory figures and run sanity checks for {family}."
        )
    )
    command.add_argument(
        "--work-dir",
        required=True,
        help="build directory containing catalog/ and features/",
    )
    command.add_argument(
        "--output-dir",
        help="report destination (default: WORK_DIR/reports)",
    )
    command.add_argument(
        "--strict",
        action="store_true",
        help="exit with status 1 when any reported sanity check fails",
    )
    return command


def main(family, argv=None):
    args = parser(family).parse_args(argv)
    config = rebuild_from_scratch.layout(args.work_dir)
    catalog = generate_feature_reports._read_table(config["catalog"])
    report_dir = os.path.abspath(
        args.output_dir or os.path.join(args.work_dir, "reports")
    )

    if family == "identifiers":
        report, summary = generate_feature_reports.write_identifier_report(
            catalog, report_dir
        )
    else:
        sidecar = os.path.join(config["features"], family + ".parquet")
        if not os.path.exists(sidecar):
            raise FileNotFoundError(
                f"{sidecar} does not exist; build the {family} feature first"
            )
        report, summary = generate_feature_reports.write_family_report(
            family, catalog, sidecar, report_dir
        )

    failed = int(summary.get("sanity_checks_failed", 0))
    print(f"wrote {report}")
    print(f"sanity checks: {summary.get('sanity_checks_run', 0)} run, {failed} failed")
    return 1 if args.strict and failed else 0

