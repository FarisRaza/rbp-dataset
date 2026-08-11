"""Create meeting-ready eCLIP summary statistics and figures.

Consumes ``rbp_master_with_eclip.csv`` and writes source coverage, mean region
fractions, ENCODE-vs-ENCORI agreement statistics, and three PNG figures. Raw
site/peak counts are reported only as provenance; the comparisons use fractions.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd


REGIONS = ["utr3", "utr5", "cds", "ncrna_exon", "intron", "intergenic"]
FAMILIES = {
    "ENCODE published": "encode_published",
    "ENCODE matched": "encode_matched",
    "ENCORI published": "encori_published",
    "ENCORI matched": "encori_matched",
    "POSTAR3": "postar",
    "Skipper": "skipper",
}
COUNT_COLUMNS = {
    "encode_published": "encode_published_n_regions",
    "encode_matched": "encode_matched_n_regions",
    "encori_published": "encori_published_n_regions",
    "encori_matched": "encori_matched_n_regions",
    "postar": "postar_n_sites",
    "skipper": "skipper_n_windows",
}
COLORS = {
    "ENCODE published": "#4C78A8",
    "ENCODE matched": "#72A0C1",
    "ENCORI published": "#F58518",
    "ENCORI matched": "#FFB05A",
    "POSTAR3": "#54A24B",
    "Skipper": "#B279A2",
}


def coverage_table(frame):
    rows = []
    for label, prefix in FAMILIES.items():
        column = COUNT_COLUMNS[prefix]
        measured = frame[column].notna()
        values = pd.to_numeric(frame.loc[measured, column], errors="coerce")
        rows.append({
            "source": label,
            "genes_measured": int(measured.sum()),
            "fraction_of_rbp_table": float(measured.mean()),
            "median_raw_regions_or_windows": float(values.median()),
            "raw_count_warning": "not comparable across sources",
        })
    return pd.DataFrame(rows)


def composition_table(frame):
    rows = []
    for label, prefix in FAMILIES.items():
        for region in REGIONS:
            column = f"{prefix}_frac_{region}"
            if column not in frame:
                continue
            values = pd.to_numeric(frame[column], errors="coerce")
            rows.append({
                "source": label,
                "region": region,
                "n_genes": int(values.notna().sum()),
                "mean_fraction": float(values.mean()),
                "median_fraction": float(values.median()),
            })
    return pd.DataFrame(rows)


def agreement_table(frame):
    rows = []
    for criterion in ("published", "matched"):
        left, right = f"encode_{criterion}", f"encori_{criterion}"
        dominant_left = f"{left}_dominant"
        dominant_right = f"{right}_dominant"
        shared_dom = frame[[dominant_left, dominant_right]].dropna()
        concordance = (
            float((shared_dom[dominant_left] == shared_dom[dominant_right]).mean())
            if len(shared_dom) else None
        )
        for region in REGIONS:
            xcol, ycol = f"{left}_frac_{region}", f"{right}_frac_{region}"
            shared = frame[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna()
            correlation = (
                float(shared[xcol].corr(shared[ycol])) if len(shared) >= 3 else None
            )
            rows.append({
                "criterion": criterion,
                "region": region,
                "shared_genes": len(shared),
                "pearson_r": correlation,
                "encode_mean_fraction": float(shared[xcol].mean()) if len(shared) else None,
                "encori_mean_fraction": float(shared[ycol].mean()) if len(shared) else None,
                "mean_absolute_gene_difference": (
                    float((shared[xcol] - shared[ycol]).abs().mean())
                    if len(shared) else None
                ),
                "dominant_region_concordance": concordance,
            })
    return pd.DataFrame(rows)


def _setup_plotting():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    return plt


def plot_coverage(coverage, path):
    plt = _setup_plotting()
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    ordered = coverage.sort_values("genes_measured")
    bars = ax.barh(
        ordered.source,
        ordered.genes_measured,
        color=[COLORS[name] for name in ordered.source],
    )
    ax.bar_label(bars, padding=4, fmt="%d")
    ax.set_xlabel("RBP genes with measurements")
    ax.set_title("CLIP source coverage in the RBP table")
    ax.set_xlim(0, max(ordered.genes_measured) * 1.16)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_composition(composition, path):
    plt = _setup_plotting()
    selected = ["ENCODE published", "ENCORI published", "POSTAR3", "Skipper"]
    pivot = composition[composition.source.isin(selected)].pivot(
        index="region", columns="source", values="mean_fraction"
    ).reindex(REGIONS)
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    x = np.arange(len(pivot.index))
    width = 0.19
    for offset, source in enumerate(selected):
        if source not in pivot:
            continue
        ax.bar(
            x + (offset - 1.5) * width,
            pivot[source],
            width,
            label=source,
            color=COLORS[source],
        )
    ax.set_xticks(x, [name.replace("_", " ") for name in pivot.index])
    ax.set_ylabel("Mean fraction of assigned binding regions")
    ax.set_title("Binding-region composition by CLIP source")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_encode_encori(frame, path, criterion="published"):
    plt = _setup_plotting()
    left, right = f"encode_{criterion}", f"encori_{criterion}"
    regions = ["utr3", "utr5", "cds", "intron"]
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 8.0), sharex=True, sharey=True)
    for ax, region in zip(axes.flat, regions):
        xcol, ycol = f"{left}_frac_{region}", f"{right}_frac_{region}"
        shared = frame[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna()
        ax.scatter(
            shared[xcol], shared[ycol], s=18, alpha=0.65,
            color="#5B6F8F", edgecolor="none",
        )
        ax.plot([0, 1], [0, 1], linestyle="--", linewidth=1, color="#777777")
        correlation = shared[xcol].corr(shared[ycol]) if len(shared) >= 3 else np.nan
        ax.set_title(f"{region.replace('_', ' ')}: n={len(shared)}, r={correlation:.2f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    fig.supxlabel("ENCODE fraction")
    fig.supylabel("ENCORI fraction")
    fig.suptitle(f"ENCODE versus ENCORI ({criterion} criteria)", y=0.99)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def summarize(input_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    frame = pd.read_csv(input_path, low_memory=False)
    coverage = coverage_table(frame)
    composition = composition_table(frame)
    agreement = agreement_table(frame)
    coverage.to_csv(os.path.join(out_dir, "eclip_source_coverage.csv"), index=False)
    composition.to_csv(os.path.join(out_dir, "eclip_region_composition.csv"), index=False)
    agreement.to_csv(os.path.join(out_dir, "encode_encori_agreement.csv"), index=False)

    plot_coverage(coverage, os.path.join(out_dir, "eclip_source_coverage.png"))
    plot_composition(composition, os.path.join(out_dir, "eclip_region_composition.png"))
    plot_encode_encori(
        frame, os.path.join(out_dir, "encode_vs_encori_published.png"), "published"
    )

    summary = {
        "input": os.path.abspath(input_path),
        "rbp_rows": len(frame),
        "unique_gene_symbols": int(frame["Name"].nunique()),
        "genes_with_any_clip_source": int(
            (pd.to_numeric(frame.get("n_eclip_sources"), errors="coerce") > 0).sum()
        ),
        "coverage": coverage.to_dict("records"),
        "published_encode_encori": agreement[
            agreement.criterion == "published"
        ].to_dict("records"),
        "interpretation": (
            "Use fractions/enrichments for biology; raw counts are source-depth provenance."
        ),
    }
    with open(os.path.join(out_dir, "eclip_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="rbp_master_with_eclip.csv")
    parser.add_argument("--out-dir", default="eclip_summary_results")
    args = parser.parse_args(argv)
    summary = summarize(args.input, args.out_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

