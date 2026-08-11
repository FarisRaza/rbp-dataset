"""ENCODE peaks passing BOTH criteria: p <= 1e-3 AND >= 2 independent experiments.

The exact parallel of encori_reproducible.py, so the two sources are held to the
same two-part standard:

  criterion          ENCODE                          ENCORI
  ---------          ------                          ------
  significance       -log10(p) >= 3   (col 8)        log10(p) <= -3   (col 5)
  reproducibility    locus called in >= 2            locus called in >= 2
                     biological replicates           distinct SBDH datasets
  merging            overlapping peaks -> cluster    overlapping sites -> cluster

Both use overlap-based clustering rather than exact-coordinate matching, because
independent experiments essentially never call identical boundaries. The reported
unit in both is the merged cluster.

Why not ENCODE's own IDR files instead? They exist (output_type "peaks" with
biological_replicates [1,2]; the name field reads e.g. AARS_K562_IDR) but are far
too sparse for stable region fractions -- AARS has 49 IDR peaks against 62,056 and
70,659 in its two replicates. Replicate intersection at a fixed p threshold keeps
enough peaks to estimate a composition while still requiring reproducibility, and
it is the construction that matches what is possible on the ENCORI side.

NOTE ON WHAT "TWO EXPERIMENTS" MEANS: for ENCODE these are two biological
replicates of ONE experiment in ONE cell line -- a narrower claim than ENCORI's
two independent published datasets, which may differ in lab, protocol, and cell
type. The criterion is structurally parallel but not equally demanding, and the
ENCODE version is the weaker of the two. Do not describe them as equivalent.

Usage:
    python encode_reproducible.py [--min-neglog10p 3] [--min-experiments 2]
"""
import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from annotate_peaks import (INTERGENIC, PRECEDENCE, IntervalIndex,
                            assign_even_fractional)
from encori_reproducible import merge_with_support

NARROWPEAK = ["chrom", "start", "end", "name", "score", "strand",
              "log2fc", "neglog10p", "qvalue", "peak"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", default="encode_peaks_grch38")
    ap.add_argument("--manifest", default="encode_peaks_grch38_manifest.csv")
    ap.add_argument("--regions", default="gencode_regions")
    ap.add_argument("--min-neglog10p", type=float, default=3.0)
    ap.add_argument("--min-log2fc", type=float, default=None,
                    help="optional; omit to match ENCORI, which has no fold-change column")
    ap.add_argument("--min-experiments", type=int, default=2)
    ap.add_argument("--out", default="encode_region_counts_p1e3_rep2.csv")
    args = ap.parse_args()

    print(f"loading region tracks from {args.regions}/ ...", flush=True)
    indexes = {}
    for cls in PRECEDENCE:
        p = os.path.join(args.regions, f"gencode_v46_{cls}.bed")
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        indexes[cls] = IntervalIndex(p)
    print("  ok", flush=True)

    man = pd.read_csv(args.manifest)
    man = man[man.kind == "replicate"]
    print(f"{len(man)} replicate peak files, "
          f"{man.rbp.nunique()} RBPs, "
          f"{man.groupby(['rbp','biosample']).ngroups} RBP x cellline groups\n", flush=True)

    counts = defaultdict(lambda: defaultdict(int))
    support = defaultdict(list)
    n_in = n_out = 0

    groups = list(man.groupby(["rbp", "biosample"]))
    for gi, ((rbp, biosample), grp) in enumerate(groups, 1):
        # pool this RBP+cellline's replicate files, tagging each peak with its
        # source accession -- the "experiment" identity for the >=N test
        pooled = defaultdict(list)          # (chrom, strand) -> [(s, e, acc)]
        n_files = 0
        for _, row in grp.iterrows():
            path = os.path.join(args.peaks, row["file"])
            if not os.path.exists(path):
                continue
            with gzip.open(path, "rt") as fh:
                d = pd.read_csv(fh, sep="\t", header=None, names=NARROWPEAK)
            if not len(d):
                continue
            n_files += 1
            m = d.neglog10p >= args.min_neglog10p
            if args.min_log2fc is not None:
                m &= d.log2fc >= args.min_log2fc
            d = d[m]
            acc = row["accession"]
            for c, s, e, st in zip(d.chrom, d.start, d.end, d.strand):
                pooled[(c, st)].append((s, e, acc))
            n_in += len(d)

        if n_files < args.min_experiments:
            continue                        # cannot satisfy the criterion

        for (chrom, strand), ivs in pooled.items():
            ivs.sort()
            clusters = merge_with_support(ivs, args.min_experiments)
            if not clusters:
                continue
            s = np.fromiter((c[0] for c in clusters), np.int64, len(clusters))
            e = np.fromiter((c[1] for c in clusters), np.int64, len(clusters))
            nd = np.fromiter((c[2] for c in clusters), np.int64, len(clusters))
            for cls, v in assign_even_fractional(chrom, strand, s, e, indexes).items():
                counts[rbp][cls] += v
            support[rbp].extend(int(x) for x in nd)
            n_out += len(clusters)

        if gi % 40 == 0 or gi == len(groups):
            print(f"  {gi}/{len(groups)} groups", flush=True)

    print(f"\npeaks passing p filter={n_in:,} -> reproducible clusters={n_out:,} "
          f"({100*n_out/max(n_in,1):.1f}%)", flush=True)

    classes = PRECEDENCE + [INTERGENIC]
    rows = []
    for rbp, d in counts.items():
        tot = sum(d.values())
        rec = {"rbp": rbp, "source": "encode", "n_sites": tot,
               "median_experiment_support": float(np.median(support[rbp])) if support[rbp] else 0.0}
        for c in classes:
            rec[f"n_{c}"] = d.get(c, 0)
        for c in classes:
            rec[f"frac_{c}"] = d.get(c, 0) / tot if tot else 0.0
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values("n_sites", ascending=False)
    out["frac_sum"] = out[[f"frac_{c}" for c in classes]].sum(axis=1)
    out.to_csv(args.out, index=False)
    print(f"wrote {args.out}  ({len(out)} RBPs)", flush=True)
    bad = out[(out.n_sites > 0) & (out.frac_sum.sub(1).abs() > 1e-9)]
    print(f"RBPs whose fractions do not sum to 1: {len(bad)}  (must be 0)", flush=True)


if __name__ == "__main__":
    main()


