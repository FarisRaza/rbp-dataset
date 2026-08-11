"""ENCORI sites passing BOTH criteria: p <= 1e-3 AND >= 2 independent datasets.

This implements the full two-part criterion:
  1. Significance  -- log10(p) <= -3, matching ENCODE's -log10(p) >= 3.
  2. Reproducibility -- the binding locus is supported by at least N distinct
     source datasets ("two experiments showing the same region").

Reproducibility CANNOT be a row predicate. Sites called from different datasets
almost never share exact coordinates, so requiring identical rows would reject
nearly everything. Instead, overlapping sites for the same RBP are merged into
clusters, and a cluster survives if it draws support from >= N distinct SBDH
dataset ids. (`SBDH<n>-<P|T|M|D>-<idx>`: the leading number is the dataset;
1,968 distinct datasets exist, with per-RBP counts that track profiling depth --
AGO2 226, TARDBP 100, PTBP1 23.)

The reported unit is therefore the MERGED CLUSTER, not the raw site. A cluster is
one binding locus regardless of how many datasets or redundant calls cover it,
which also removes the publication-effort confound that makes raw ENCORI site
counts uninterpretable across RBPs.

Memory stays flat by exploiting the file's sort order (verified: chromosomes are
contiguous blocks, start-sorted within each), so only one chromosome is held at
a time.

Usage:
    python encori_reproducible.py [--max-log10p -3] [--min-datasets 2]
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from annotate_peaks import (INTERGENIC, PRECEDENCE, IntervalIndex,
                            assign_even_fractional)

ENCORI = "sites_annotated.bed"


def merge_with_support(intervals, min_datasets):
    """intervals: list of (start, end, dataset). Must be start-sorted.

    Sweep-merge overlapping intervals, tracking the set of contributing datasets.
    Returns [(start, end, n_datasets), ...] for clusters meeting the threshold.
    """
    out = []
    cs = ce = None
    ds = set()
    for s, e, d in intervals:
        if cs is None:
            cs, ce, ds = s, e, {d}
        elif s <= ce:                      # overlap -> same locus
            ce = max(ce, e)
            ds.add(d)
        else:
            if len(ds) >= min_datasets:
                out.append((cs, ce, len(ds)))
            cs, ce, ds = s, e, {d}
    if cs is not None and len(ds) >= min_datasets:
        out.append((cs, ce, len(ds)))
    return out


def assign_regions(chrom, strand, starts, ends, indexes):
    n = len(starts)
    label = np.full(n, INTERGENIC, dtype=object)
    todo = np.ones(n, dtype=bool)
    for cls in PRECEDENCE:
        if not todo.any():
            break
        hit = indexes[cls].overlaps(chrom, strand, starts, ends) & todo
        if hit.any():
            label[hit] = cls
            todo &= ~hit
    return label


def flush_chrom(chrom, buf, indexes, counts, support, min_datasets, stats):
    """buf: (rbp, strand) -> list of (start, end, dataset)"""
    for (rbp, strand), ivs in buf.items():
        ivs.sort()
        clusters = merge_with_support(ivs, min_datasets)
        stats["sites_in"] += len(ivs)
        if not clusters:
            continue
        s = np.fromiter((c[0] for c in clusters), dtype=np.int64, count=len(clusters))
        e = np.fromiter((c[1] for c in clusters), dtype=np.int64, count=len(clusters))
        nd = np.fromiter((c[2] for c in clusters), dtype=np.int64, count=len(clusters))
        for cls, v in assign_even_fractional(chrom, strand, s, e, indexes).items():
            counts[rbp][cls] += v
        support[rbp].extend(int(x) for x in nd)
        stats["clusters_out"] += len(clusters)
    buf.clear()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-log10p", type=float, default=-3.0)
    ap.add_argument("--min-datasets", type=int, default=2)
    ap.add_argument("--regions", default="gencode_regions")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    out_path = args.out or (
        f"encori_region_counts_p1e{abs(int(args.max_log10p))}_rep{args.min_datasets}.csv")

    print(f"loading region tracks from {args.regions}/ ...", flush=True)
    indexes = {}
    for cls in PRECEDENCE:
        p = os.path.join(args.regions, f"gencode_v46_{cls}.bed")
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        indexes[cls] = IntervalIndex(p)
    print("  ok", flush=True)

    counts = defaultdict(lambda: defaultdict(int))
    support = defaultdict(list)
    buf = defaultdict(list)
    stats = {"sites_in": 0, "clusters_out": 0}

    cur_chrom = None
    prev = None
    n_lines = n_kept = n_filtered = n_bad = 0

    print(f"streaming {ENCORI}  (p <= 1e{int(args.max_log10p)}, "
          f">= {args.min_datasets} datasets) ...", flush=True)
    with open(ENCORI, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            n_lines += 1
            if args.limit and n_lines > args.limit:
                break
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 14:
                n_bad += 1
                continue
            key = tuple(f[:7])
            if key == prev:
                continue
            prev = key

            try:
                score = float(f[4])
                s, e = int(f[1]), int(f[2])
            except ValueError:
                n_bad += 1
                continue
            if score > args.max_log10p:
                n_filtered += 1
                continue
            rbp = f[6].strip()
            if not rbp:
                n_bad += 1
                continue

            chrom = f[0]
            if chrom != cur_chrom:
                if cur_chrom is not None:
                    flush_chrom(cur_chrom, buf, indexes, counts, support,
                                args.min_datasets, stats)
                    print(f"  {cur_chrom}: {stats['sites_in']:,} sites -> "
                          f"{stats['clusters_out']:,} reproducible clusters", flush=True)
                cur_chrom = chrom

            buf[(rbp, f[5])].append((s, e, f[3].split("-")[0]))
            n_kept += 1

    if cur_chrom is not None:
        flush_chrom(cur_chrom, buf, indexes, counts, support, args.min_datasets, stats)
        print(f"  {cur_chrom}: {stats['sites_in']:,} sites -> "
              f"{stats['clusters_out']:,} reproducible clusters", flush=True)

    print(f"\nlines={n_lines:,}  passed p filter={n_kept:,}  "
          f"failed p filter={n_filtered:,}  malformed={n_bad:,}", flush=True)
    print(f"sites in={stats['sites_in']:,} -> reproducible clusters="
          f"{stats['clusters_out']:,} "
          f"({100*stats['clusters_out']/max(stats['sites_in'],1):.1f}%)", flush=True)

    classes = PRECEDENCE + [INTERGENIC]
    rows = []
    for rbp, d in counts.items():
        tot = sum(d.values())
        rec = {"rbp": rbp, "source": "encori", "n_sites": tot,
               "median_dataset_support": float(np.median(support[rbp])) if support[rbp] else 0.0,
               "max_dataset_support": int(max(support[rbp])) if support[rbp] else 0}
        for c in classes:
            rec[f"n_{c}"] = d.get(c, 0)
        for c in classes:
            rec[f"frac_{c}"] = d.get(c, 0) / tot if tot else 0.0
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values("n_sites", ascending=False)
    out["frac_sum"] = out[[f"frac_{c}" for c in classes]].sum(axis=1)
    out.to_csv(out_path, index=False)
    print(f"wrote {out_path}  ({len(out)} RBPs)", flush=True)
    bad = out[(out.n_sites > 0) & (out.frac_sum.sub(1).abs() > 1e-9)]
    print(f"RBPs whose fractions do not sum to 1: {len(bad)}  (must be 0)", flush=True)


if __name__ == "__main__":
    main()


