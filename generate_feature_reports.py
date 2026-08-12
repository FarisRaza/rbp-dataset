"""Generate one Markdown quality-control report per feature family.

The reports are intentionally derived from keyed sidecars rather than from the
very wide final table. This keeps report generation fast and lets a user inspect
one family even when they intentionally requested only a subset of features.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics

import pandas as pd

import rebuild_from_scratch


PRIMARY_COLUMNS = {
    "cider": ["NCPR"],
    "idr": ["IDR_count"],
    "domains": ["Domains_count"],
    "go": ["C_ids", "P_ids", "F_ids"],
    "eclip": [],  # source-specific column names are discovered dynamically
    "interpro": ["InterPro_n_hits"],
    "ptm": [],  # all ptm_<type> indicator columns are summed
    "opentargets": [
        "opentargets_expression_tissue_count",
        "opentargets_disease_count",
    ],
    "cdcode": ["UID"],
    "string": ["string_partners_ensp_by_query"],
    "go_roles": [
        "role_in_transcription",
        "role_in_translation",
        "role_in_mrna_stability",
    ],
    "pslab": ["Delta G [kT]"],
    "rcsb": ["RCSB_PDB_count"],
}

PROVENANCE_WORDS = ("scope", "method", "version", "source", "coordinate")


def _read_table(path):
    if path.lower().endswith(".parquet"):
        return pd.read_parquet(path)
    if path.lower().endswith(".csv"):
        return pd.read_csv(path, low_memory=False)
    if path.lower().endswith(".jsonl"):
        return pd.read_json(path, lines=True)
    raise ValueError(f"unsupported report input: {path}")


def _decode(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        if stripped.startswith(("[", "{")):
            try:
                return json.loads(stripped)
            except json.JSONDecodeError:
                return value
    return value


def _is_present(value):
    value = _decode(value)
    if value is None:
        return False
    if isinstance(value, (list, tuple, set, dict, str)):
        return len(value) > 0
    return True


def _container_size(value):
    value = _decode(value)
    if value is None:
        return 0.0
    if isinstance(value, dict):
        nested = [_container_size(item) for item in value.values()]
        return float(sum(nested)) if nested else 0.0
    if isinstance(value, (list, tuple, set)):
        return float(len(value))
    if isinstance(value, bool):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return 1.0 if str(value).strip() else 0.0


def _metric_columns(family, frame):
    available = [c for c in PRIMARY_COLUMNS[family] if c in frame.columns]
    if family == "ptm":
        available = [
            c for c in frame.columns
            if c.startswith("ptm_")
            and not c.endswith(("_positions", "_residues"))
            and c not in {"ptm_projection_dropped_sites", "ptm_coordinate_system"}
        ]
    elif family == "eclip" and not available:
        preferred = [
            c for c in frame.columns
            if any(token in c.lower() for token in ("has_data", "site_count", "peak_count"))
        ]
        available = preferred or [
            c for c in frame.columns
            if c != "protein_key" and not any(x in c.lower() for x in PROVENANCE_WORDS)
        ]
    return available


def _row_metric(family, frame, columns):
    if not columns:
        return [0.0] * len(frame)
    values = []
    for record in frame[columns].itertuples(index=False, name=None):
        if family == "string":
            total = 0
            for value in record:
                decoded = _decode(value)
                if isinstance(decoded, dict):
                    total += sum(
                        len(partners) if isinstance(partners, dict) else 0
                        for partners in decoded.values()
                    )
            values.append(float(total))
        else:
            values.append(sum(_container_size(value) for value in record))
    return values


def _coverage(frame):
    values = []
    denominator = max(len(frame), 1)
    for column in frame.columns:
        if column == "protein_key" or any(x in column.lower() for x in PROVENANCE_WORDS):
            continue
        present = sum(_is_present(value) for value in frame[column])
        values.append((column, present, 100.0 * present / denominator))
    return sorted(values, key=lambda item: (-item[2], item[0]))


def _plot_reports(family, coverage, metric, report_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = os.path.join(report_dir, "figures")
    os.makedirs(figure_dir, exist_ok=True)

    coverage_path = os.path.join(figure_dir, f"{family}_coverage.png")
    shown = coverage[:12]
    fig, axis = plt.subplots(figsize=(8, max(3.2, 0.36 * len(shown))))
    axis.barh([x[0] for x in shown][::-1], [x[2] for x in shown][::-1], color="#315f72")
    axis.set_xlim(0, 100)
    axis.set_xlabel("Rows with a populated value (%)")
    axis.set_title(f"{family}: column coverage")
    fig.tight_layout()
    fig.savefig(coverage_path, dpi=160)
    plt.close(fig)

    distribution_path = os.path.join(figure_dir, f"{family}_distribution.png")
    finite = [x for x in metric if math.isfinite(x)]
    fig, axis = plt.subplots(figsize=(7, 4.2))
    if finite:
        upper = statistics.quantiles(finite, n=100)[98] if len(finite) >= 100 else max(finite)
        clipped = [min(value, upper) for value in finite]
        axis.hist(clipped, bins=min(40, max(5, int(math.sqrt(len(clipped))))), color="#bd5d38")
        if max(finite) > upper:
            axis.set_xlabel(f"Per-row signal count/value (values above p99 clipped to {upper:g})")
        else:
            axis.set_xlabel("Per-row signal count/value")
    axis.set_ylabel("Rows")
    axis.set_title(f"{family}: primary feature distribution")
    fig.tight_layout()
    fig.savefig(distribution_path, dpi=160)
    plt.close(fig)
    return coverage_path, distribution_path


def write_family_report(family, catalog, sidecar, report_dir):
    frame = _read_table(sidecar)
    if "protein_key" not in frame.columns:
        raise ValueError(f"{sidecar} has no protein_key column")
    catalog_keys = set(catalog["protein_key"].astype(str))
    sidecar_keys = set(frame["protein_key"].astype(str))
    coverage = _coverage(frame)
    metric_columns = _metric_columns(family, frame)
    metric = _row_metric(family, frame, metric_columns)
    coverage_path, distribution_path = _plot_reports(
        family, coverage, metric, report_dir
    )
    positive = sum(value > 0 for value in metric)
    numeric = [value for value in metric if math.isfinite(value)]
    summary = {
        "catalog_rows": len(catalog),
        "sidecar_rows": len(frame),
        "sidecar_columns_excluding_key": len(frame.columns) - 1,
        "duplicate_protein_keys": int(frame["protein_key"].duplicated().sum()),
        "catalog_rows_missing_from_sidecar": len(catalog_keys.difference(sidecar_keys)),
        "sidecar_rows_not_in_catalog": len(sidecar_keys.difference(catalog_keys)),
        "primary_columns": metric_columns,
        "rows_with_nonzero_primary_signal": positive,
        "primary_signal_median": statistics.median(numeric) if numeric else None,
        "primary_signal_mean": statistics.fmean(numeric) if numeric else None,
        "primary_signal_max": max(numeric) if numeric else None,
    }
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(report_dir, f"{family}.md")
    relative_coverage = os.path.relpath(coverage_path, report_dir).replace("\\", "/")
    relative_distribution = os.path.relpath(distribution_path, report_dir).replace("\\", "/")
    lines = [
        f"# {family} feature report",
        "",
        "This file is generated by `generate_feature_reports.py`; rerun it after changing inputs.",
        "",
        "## Sanity-check summary",
        "",
        "| Check | Value |",
        "|---|---:|",
    ]
    for key, value in summary.items():
        shown_value = ", ".join(value) if isinstance(value, list) else value
        if isinstance(shown_value, float):
            shown_value = f"{shown_value:.3f}"
        lines.append(f"| {key.replace('_', ' ')} | {shown_value} |")
    lines.extend([
        "",
        "A valid sidecar should have zero duplicate keys, zero missing catalog keys, and zero keys outside the catalog.",
        "",
        "## Visual checks",
        "",
        f"![Column coverage]({relative_coverage})",
        "",
        f"![Primary distribution]({relative_distribution})",
        "",
        "## Column coverage",
        "",
        "| Column | Populated rows | Coverage |",
        "|---|---:|---:|",
    ])
    for column, present, percentage in coverage:
        lines.append(f"| `{column}` | {present:,} | {percentage:.1f}% |")
    lines.extend([
        "",
        f"Input sidecar: `{os.path.basename(sidecar)}`",
        "",
    ])
    with open(report_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
    return report_path, summary


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work-dir", required=True)
    p.add_argument(
        "--features", default="all",
        help="all, default, or comma-separated feature families",
    )
    p.add_argument("--output-dir", help="default: WORK_DIR/reports")
    p.add_argument("--skip-unavailable", action="store_true")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    config = rebuild_from_scratch.layout(args.work_dir)
    catalog = _read_table(config["catalog"])
    if "protein_key" not in catalog.columns:
        raise ValueError("catalog has no protein_key column")
    families = rebuild_from_scratch.parse_features(args.features)
    report_dir = os.path.abspath(args.output_dir or os.path.join(args.work_dir, "reports"))
    generated = []
    summaries = {}
    for family in families:
        sidecar = os.path.join(config["features"], family + ".parquet")
        if not os.path.exists(sidecar):
            if args.skip_unavailable:
                print(f"SKIP {family}: {sidecar} does not exist")
                continue
            raise FileNotFoundError(sidecar)
        report, summary = write_family_report(family, catalog, sidecar, report_dir)
        generated.append(report)
        summaries[family] = summary
        print(f"wrote {report}")
    with open(os.path.join(report_dir, "report_manifest.json"), "w", encoding="utf-8") as handle:
        json.dump(summaries, handle, indent=2)
    return generated


if __name__ == "__main__":
    main()
