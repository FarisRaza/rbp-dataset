"""ENCORI binding regions built to ENCORI's OWN published methodology.

Replaces our approximations with the definitions from Zhou et al., Nature Methods
2026 ("An encyclopedic regulatory and functional atlas of RNA interactomes"),
Methods, "Construction of comprehensive RBP-binding maps in human":

  1. HIGH-CONFIDENCE CRITERION (verbatim):
       "regions with P < 1.0x10^-5 in cases of only a single experiment, or
        regions observed in at least 50% of available experiments (must be >= 2)
        in the same cell line or tissue. Each RBP was required to have a minimum
        of 300 high-confidence binding regions."
     Implemented as, per merged region:
         keep if  (n_support == 1 and min_log10p <= -5)
               or (n_support >= max(2, ceil(0.5 * n_available_in_that_cellline)))
     The denominator comes from hg38_all_dataset_info_table.txt (1,984 datasets,
     286 RBPs, 108 cell/tissue types), which maps SBDH id -> GeneSymbol +
     CellTissue. Merging happens WITHIN a cell line, as the paper specifies.

  2. MERGE PROTOCOL (verbatim):
       "For single-nucleotide binding sites, loci were extended by +-20 bp before
        merging. For each merged region, the smallest P value from its
        constituent regions was assigned as its P value."
     Site widths confirmed empirically: type M (mutation) and D (deletion) are
     width 1; type P (peak) is ~40 bp. Only width-1 sites are extended.

  3. REGION ASSIGNMENT (verbatim):
       "In cases in which an RBP-binding region overlapped with multiple
        functional regions or REs, the binding event was EVENLY distributed as
        fractional counts across these regions or REs."
     So a region overlapping CDS and 3'UTR contributes 0.5 to each -- an EQUAL
     split across overlapping classes, NOT weighted by overlap length, and NOT
     the precedence rule used elsewhere in this project. This is also the answer
     to the objection that an RBP can bind more than one region type: it can, and
     ENCORI accounts for it fractionally rather than by forcing a single label.

Memory stays flat by processing one chromosome at a time (chromosomes are
contiguous, start-sorted blocks in sites_annotated.bed -- verified).
"""
import argparse
import math
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from annotate_peaks import (INTERGENIC, PRECEDENCE, IntervalIndex,
                            assign_even_fractional)

ENCORI = "sites_annotated.bed"
DATASET_INFO = "encori_dataset_info_hg38.txt"

SINGLE_NT_EXTEND = 20        # paper: +-20 bp for single-nucleotide sites
SINGLE_P_THRESHOLD = -5      # paper: P < 1e-5 for single-experiment regions
MIN_REGIONS_PER_RBP = 300    # paper: >= 300 high-confidence regions per RBP


def load_dataset_info(path):
    """SBDH id -> (rbp, celltissue); and (rbp, celltissue) -> n_available."""
    d = pd.read_csv(path, sep="\t", low_memory=False)
    ds2cell = {}
    for sid, sym, cell in zip(d.DataSetId, d.GeneSymbol, d.CellTissue):
        ds2cell[str(sid)] = (str(sym), str(cell))
    avail = d.groupby(["GeneSymbol", "CellTissue"]).size().to_dict()
    return ds2cell, avail


def merge_regions(ivs):
    """ivs: sorted [(start, end, dataset, log10p)]. -> [(s, e, {datasets}, minp)]"""
    out = []
    cs = ce = None
    ds, mp = set(), 0
    for s, e, d, p in ivs:
        if cs is None:
            cs, ce, ds, mp = s, e, {d}, p
        elif s <= ce:
            ce = max(ce, e)
            ds.add(d)
            mp = min(mp, p)          # smallest P == most negative log10(P)
        else:
            out.append((cs, ce, ds, mp))
            cs, ce, ds, mp = s, e, {d}, p
    if cs is not None:
        out.append((cs, ce, ds, mp))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", default="gencode_regions")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="encori_region_counts_published.csv")
    ap.add_argument("--min-regions", type=int, default=MIN_REGIONS_PER_RBP)
    args = ap.parse_args()

    if not os.path.exists(DATASET_INFO):
        sys.exit(f"missing {DATASET_INFO} (from ENCORI_referenceData.zip)")

    ds2cell, avail = load_dataset_info(DATASET_INFO)
    print(f"dataset info: {len(ds2cell):,} datasets, "
          f"{len(avail):,} (RBP, cell/tissue) groups", flush=True)

    print(f"loading region tracks from {args.regions}/ ...", flush=True)
    indexes = {}
    for cls in PRECEDENCE:
        p = os.path.join(args.regions, f"gencode_v46_{cls}.bed")
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        indexes[cls] = IntervalIndex(p)
    print("  ok", flush=True)

    frac = defaultdict(lambda: defaultdict(float))   # rbp -> class -> fractional count
    n_regions = defaultdict(int)                     # rbp -> high-confidence regions
    support_hist = defaultdict(list)

    buf = defaultdict(list)      # (rbp, cell, strand) -> [(s,e,ds,p)]
    cur_chrom = None
    prev = None
    n_lines = n_sites = n_bad = n_nodataset = 0
    n_merged = n_kept = 0

    def flush(chrom):
        nonlocal n_merged, n_kept
        for (rbp, cell, strand), ivs in buf.items():
            ivs.sort()
            regions = merge_regions(ivs)
            n_merged += len(regions)
            navail = avail.get((rbp, cell), 1)
            need = max(2, math.ceil(0.5 * navail))

            keep_s, keep_e = [], []
            for s, e, ds, mp in regions:
                k = len(ds)
                if (k == 1 and mp <= SINGLE_P_THRESHOLD) or (k >= need):
                    keep_s.append(s)
                    keep_e.append(e)
                    support_hist[rbp].append(k)
            if not keep_s:
                continue
            n_kept += len(keep_s)
            n_regions[rbp] += len(keep_s)
            s = np.asarray(keep_s, dtype=np.int64)
            e = np.asarray(keep_e, dtype=np.int64)
            for cls, v in assign_even_fractional(chrom, strand, s, e, indexes).items():
                frac[rbp][cls] += v
        buf.clear()

    print(f"streaming {ENCORI} ...", flush=True)
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
                s, e = int(f[1]), int(f[2])
                p = int(float(f[4]))
            except ValueError:
                n_bad += 1
                continue

            sid = f[3].split("-")[0]
            meta = ds2cell.get(sid)
            if meta is None:
                n_nodataset += 1
                continue
            _, cell = meta
            rbp = f[6].strip()
            if not rbp:
                n_bad += 1
                continue

            # paper: extend single-nucleotide loci by +-20 bp before merging
            if e - s <= 1:
                s = max(0, s - SINGLE_NT_EXTEND)
                e = e + SINGLE_NT_EXTEND

            chrom = f[0]
            if chrom != cur_chrom:
                if cur_chrom is not None:
                    flush(cur_chrom)
                    print(f"  {cur_chrom}: merged={n_merged:,} kept={n_kept:,}", flush=True)
                cur_chrom = chrom

            buf[(rbp, cell, f[5])].append((s, e, sid, p))
            n_sites += 1

    if cur_chrom is not None:
        flush(cur_chrom)
        print(f"  {cur_chrom}: merged={n_merged:,} kept={n_kept:,}", flush=True)

    print(f"\nlines={n_lines:,} sites={n_sites:,} malformed={n_bad:,} "
          f"unknown-dataset={n_nodataset:,}", flush=True)
    print(f"merged regions={n_merged:,} -> high-confidence={n_kept:,} "
          f"({100*n_kept/max(n_merged,1):.1f}%)", flush=True)

    classes = PRECEDENCE + [INTERGENIC]
    rows = []
    for rbp, d in frac.items():
        tot = sum(d.values())
        if tot <= 0:
            continue
        rec = {"rbp": rbp, "source": "encori",
               "n_sites": n_regions[rbp],
               "n_regions_fractional": tot,
               "median_support": float(np.median(support_hist[rbp])) if support_hist[rbp] else 0.0}
        for c in classes:
            rec[f"n_{c}"] = d.get(c, 0.0)
        for c in classes:
            rec[f"frac_{c}"] = d.get(c, 0.0) / tot
        rows.append(rec)

    out = pd.DataFrame(rows)
    before = len(out)
    out = out[out.n_sites >= args.min_regions]      # paper: >= 300 regions per RBP
    print(f"\nRBPs: {before} -> {len(out)} after >= {args.min_regions} "
          f"high-confidence regions requirement", flush=True)

    out = out.sort_values("n_sites", ascending=False)
    out["frac_sum"] = out[[f"frac_{c}" for c in classes]].sum(axis=1)
    out.to_csv(args.out, index=False)
    print(f"wrote {args.out}  ({len(out)} RBPs)", flush=True)
    bad = out[(out.n_sites > 0) & (out.frac_sum.sub(1).abs() > 1e-9)]
    print(f"RBPs whose fractions do not sum to 1: {len(bad)}  (must be 0)", flush=True)


if __name__ == "__main__":
    main()


