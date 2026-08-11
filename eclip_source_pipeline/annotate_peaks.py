"""Assign eCLIP peaks to exactly one genomic region class, and summarise per RBP.

This replaces the previous annotation, which double-counted peaks because it
intersected against unmerged per-isoform GENCODE intervals (one output line per
overlapping isoform). Here every peak is assigned to EXACTLY ONE region by an
explicit precedence rule, so region fractions sum to 1.0 by construction.

Two assignment modes:
    precedence   -- peak takes the highest-priority class it overlaps at all
    max_overlap  -- peak takes the class it shares the most base pairs with

Peaks overlapping no annotated class are assigned 'intergenic', so nothing is
silently dropped. The previous pipeline lost ~30% of peaks this way (rows summed
to 0.70) without reporting it.

Significance tiers (all emitted so they can be compared):
    idr          -- ENCODE IDR reproducible peaks (biological_replicates [1,2])
    vannostrand  -- log2FC >= 3 AND -log10(p) >= 3   (Van Nostrand et al. 2020, Nature)
    p1e4         -- -log10(p) > 4  (p < 0.0001), no fold-change gate
    all          -- unfiltered, for reference only

Usage:
    python annotate_peaks.py [--peaks DIR] [--regions DIR] [--mode precedence]
"""
import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# Precedence: earlier = higher priority. A peak touching both a 3'UTR and an
# intron is called 3'UTR. Exonic classes outrank intron so that a peak on an
# exon of one isoform and an intron of another is called exonic.
PRECEDENCE = ["utr3", "utr5", "cds", "ncrna_exon", "intron"]
INTERGENIC = "intergenic"

NARROWPEAK = ["chrom", "start", "end", "name", "score", "strand",
              "log2fc", "neglog10p", "x", "y"]


class IntervalIndex:
    """Disjoint sorted intervals per (chrom, strand), with O(log n) overlap tests.

    Relies on the region tracks being MERGED (non-overlapping within a class),
    which build_gencode_regions.py guarantees.
    """

    def __init__(self, bed_path):
        starts, ends = defaultdict(list), defaultdict(list)
        with open(bed_path, "r", encoding="utf-8") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6:
                    continue
                k = (f[0], f[5])
                starts[k].append(int(f[1]))
                ends[k].append(int(f[2]))
        self.idx = {}
        for k in starts:
            s = np.asarray(starts[k], dtype=np.int64)
            e = np.asarray(ends[k], dtype=np.int64)
            o = np.argsort(s, kind="stable")
            s, e = s[o], e[o]
            # prefix[i] = total covered bp in intervals[0..i-1]; enables an
            # O(log n) vectorised "covered bp in [0, x)" query, so overlap_bp
            # needs no per-peak Python loop.
            prefix = np.concatenate(([0], np.cumsum(e - s)))
            self.idx[k] = (s, e, prefix)

    def overlaps(self, chrom, strand, s, e):
        """Boolean array: does each peak overlap any interval?

        Intervals are disjoint and sorted, so only the last interval starting
        before the peak ends can overlap it -- one lookup suffices.
        """
        hit = np.zeros(len(s), dtype=bool)
        arr = self.idx.get((chrom, strand))
        if arr is None:
            return hit
        S, E, _ = arr
        i = np.searchsorted(S, e, side="left") - 1
        ok = i >= 0
        if ok.any():
            hit[ok] = E[i[ok]] > s[ok]
        return hit

    def _covered_to(self, S, E, prefix, x):
        """Vectorised: total covered bp in [0, x) for disjoint sorted intervals."""
        k = np.searchsorted(S, x, side="right") - 1
        ok = k >= 0
        kk = np.where(ok, k, 0)
        partial = np.clip(x - S[kk], 0, E[kk] - S[kk])
        return np.where(ok, prefix[kk] + partial, 0)

    def overlap_bp(self, chrom, strand, s, e):
        """Total overlapping base pairs per peak, fully vectorised."""
        arr = self.idx.get((chrom, strand))
        if arr is None:
            return np.zeros(len(s), dtype=np.int64)
        S, E, prefix = arr
        return self._covered_to(S, E, prefix, e) - self._covered_to(S, E, prefix, s)


def load_peaks(path):
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, names=NARROWPEAK)
    return df


def tier_mask(df, tier):
    if tier == "idr":
        return np.ones(len(df), dtype=bool)          # file is already the IDR set
    if tier == "vannostrand":
        return (df.log2fc >= 3) & (df.neglog10p >= 3)
    if tier == "p1e4":
        return df.neglog10p > 4
    return np.ones(len(df), dtype=bool)


def assign_even_fractional(chrom, strand, starts, ends, indexes):
    """ENCORI's rule: split each region evenly across every class it overlaps.

    From Zhou et al., Nature Methods 2026, Methods: "In cases in which an
    RBP-binding region overlapped with multiple functional regions or REs, the
    binding event was evenly distributed as fractional counts across these
    regions or REs."

    So a region touching CDS and 3'UTR contributes 0.5 to each -- an EQUAL split,
    not weighted by overlap length. Every region contributes exactly 1.0 in
    total, so fractions still sum to 1.

    This is preferred over PRECEDENCE for the same reason precedence was
    questioned: an RBP genuinely can bind more than one region type, and
    precedence silently favours whichever class is ordered first. Ordering utr3
    first inflated 3'UTR-dominant RBP counts by roughly 2x (82 -> 39 on ENCORI).

    Returns dict class -> summed fractional count.
    """
    n = len(starts)
    hits = np.zeros((len(PRECEDENCE), n), dtype=bool)
    for i, cls in enumerate(PRECEDENCE):
        hits[i] = indexes[cls].overlaps(chrom, strand, starts, ends)

    n_classes = hits.sum(axis=0)
    out = {c: 0.0 for c in PRECEDENCE}
    out[INTERGENIC] = float((n_classes == 0).sum())
    weight = np.where(n_classes > 0, 1.0 / np.maximum(n_classes, 1), 0.0)
    for i, cls in enumerate(PRECEDENCE):
        out[cls] += float((hits[i] * weight).sum())
    return out


def tally_multilabel(df, indexes):
    """Count each peak in EVERY region class it touches.

    Answers "what share of this RBP's peaks touch a 3'UTR at all", which is a
    different and equally valid question from the precedence composition. Sums
    to >1 by design, and that is correct -- a peak straddling a stop codon
    really is in both CDS and 3'UTR.

    The numerator is UNIQUE PEAKS touching the class, never intersect output
    lines. Counting output lines is what produced the impossible frac_3UTR=3.50
    in the old pipeline: an isoform-multiplied numerator over a deduplicated
    denominator.
    """
    out = {cls: 0 for cls in PRECEDENCE}
    for (chrom, strand), grp in df.groupby(["chrom", "strand"], sort=False):
        s = grp.start.to_numpy()
        e = grp.end.to_numpy()
        for cls in PRECEDENCE:
            out[cls] += int(indexes[cls].overlaps(chrom, strand, s, e).sum())
    return out


def tally_proportional(df, indexes):
    """Split each peak across classes by the share of its length in each.

    A peak 60% in CDS and 40% in 3'UTR contributes 0.6 and 0.4. Keeps the
    boundary information that precedence discards while still summing to 1.
    """
    out = {cls: 0.0 for cls in PRECEDENCE}
    out[INTERGENIC] = 0.0
    for (chrom, strand), grp in df.groupby(["chrom", "strand"], sort=False):
        s = grp.start.to_numpy()
        e = grp.end.to_numpy()
        width = (e - s).astype(np.float64)
        width[width <= 0] = 1.0

        bp = np.zeros((len(PRECEDENCE), len(s)), dtype=np.float64)
        for i, cls in enumerate(PRECEDENCE):
            bp[i] = indexes[cls].overlap_bp(chrom, strand, s, e)

        # Class tracks can overlap each other (e.g. one isoform's 3'UTR is
        # another's CDS), so covered bp may exceed the peak width. Normalise by
        # the larger of the two: each peak then contributes at most 1.0 in
        # total, and any genuinely uncovered remainder falls to intergenic.
        covered = bp.sum(axis=0)
        denom = np.maximum(covered, width)

        for i, cls in enumerate(PRECEDENCE):
            out[cls] += float((bp[i] / denom).sum())
        out[INTERGENIC] += float((1.0 - covered / denom).sum())
    return out


def assign(df, indexes, mode):
    """Return an array of region labels, one per peak."""
    n = len(df)
    label = np.full(n, INTERGENIC, dtype=object)

    for (chrom, strand), grp in df.groupby(["chrom", "strand"], sort=False):
        pos = grp.index.to_numpy()
        s = grp.start.to_numpy()
        e = grp.end.to_numpy()

        if mode == "precedence":
            unassigned = np.ones(len(pos), dtype=bool)
            for cls in PRECEDENCE:
                if not unassigned.any():
                    break
                hit = indexes[cls].overlaps(chrom, strand, s, e) & unassigned
                if hit.any():
                    label[pos[hit]] = cls
                    unassigned &= ~hit
        else:  # max_overlap
            bp = np.zeros((len(PRECEDENCE), len(pos)), dtype=np.int64)
            for i, cls in enumerate(PRECEDENCE):
                bp[i] = indexes[cls].overlap_bp(chrom, strand, s, e)
            best = bp.argmax(axis=0)
            any_hit = bp.max(axis=0) > 0
            lab = np.array(PRECEDENCE, dtype=object)[best]
            label[pos[any_hit]] = lab[any_hit]

    return label


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", default="encode_peaks_grch38")
    ap.add_argument("--regions", default="gencode_regions")
    ap.add_argument("--manifest", default="encode_peaks_grch38_manifest.csv")
    ap.add_argument("--mode", default="precedence", choices=["precedence", "max_overlap"])
    ap.add_argument("--out-prefix", default="eclip")
    args = ap.parse_args()

    print(f"loading region tracks from {args.regions}/ ...", flush=True)
    indexes = {}
    for cls in PRECEDENCE:
        p = os.path.join(args.regions, f"gencode_v46_{cls}.bed")
        if not os.path.exists(p):
            sys.exit(f"missing region track: {p} (run build_gencode_regions.py first)")
        indexes[cls] = IntervalIndex(p)
        n = sum(len(v[0]) for v in indexes[cls].idx.values())
        print(f"  {cls:11s} {n:>9,} merged intervals", flush=True)

    man = pd.read_csv(args.manifest)
    man = man.set_index("file")

    files = sorted(glob.glob(os.path.join(args.peaks, "*.bed.gz")))
    print(f"\nannotating {len(files)} peak files (mode={args.mode}) ...", flush=True)

    rows = []
    for i, path in enumerate(files, 1):
        base = os.path.basename(path)
        meta = man.loc[base] if base in man.index else None
        if meta is None:
            print(f"  !! {base} not in manifest, skipping", flush=True)
            continue

        df = load_peaks(path)
        if not len(df):
            continue

        kind = meta["kind"]
        tiers = ["idr"] if kind == "idr" else ["all", "vannostrand", "p1e4"]

        for tier in tiers:
            m = tier_mask(df, tier)
            sub = df[m].reset_index(drop=True)
            if not len(sub):
                rows.append({
                    "rbp": meta["rbp"], "biosample": meta["biosample"],
                    "accession": meta["accession"], "kind": kind, "tier": tier,
                    "n_peaks": 0,
                })
                continue

            rec = {
                "rbp": meta["rbp"], "biosample": meta["biosample"],
                "accession": meta["accession"], "kind": kind, "tier": tier,
                "n_peaks": len(sub),
            }

            # scheme 1: precedence -- one class per peak, composition sums to 1
            counts = pd.Series(assign(sub, indexes, args.mode)).value_counts()
            for cls in PRECEDENCE + [INTERGENIC]:
                rec[f"n_{cls}"] = int(counts.get(cls, 0))

            # scheme 2: multi-label -- a peak counts in EVERY class it touches,
            # so these sum to >1. Numerator is unique peaks, never intersect lines.
            for cls, n in tally_multilabel(sub, indexes).items():
                rec[f"ml_n_{cls}"] = n

            # scheme 3: proportional -- peak split by share of its length per class
            for cls, v in tally_proportional(sub, indexes).items():
                rec[f"prop_n_{cls}"] = v

            rows.append(rec)

        if i % 25 == 0 or i == len(files):
            print(f"  {i}/{len(files)}", flush=True)

    out = pd.DataFrame(rows).fillna(0)
    classes = PRECEDENCE + [INTERGENIC]
    denom = out.n_peaks.replace(0, np.nan)

    # precedence composition -- sums to exactly 1
    for cls in classes:
        out[f"frac_{cls}"] = np.where(out.n_peaks > 0, out[f"n_{cls}"] / denom, 0.0)
    out["frac_sum"] = out[[f"frac_{c}" for c in classes]].sum(axis=1)

    # multi-label rates -- "share of peaks touching class X", sums to >1 by design
    for cls in PRECEDENCE:
        out[f"ml_frac_{cls}"] = np.where(out.n_peaks > 0, out[f"ml_n_{cls}"] / denom, 0.0)

    # proportional composition -- also sums to 1
    for cls in classes:
        out[f"prop_frac_{cls}"] = np.where(out.n_peaks > 0, out[f"prop_n_{cls}"] / denom, 0.0)

    # ---- enrichment over a length background -------------------------------
    # Raw fractions are dominated by intron simply because introns are ~93% of
    # annotated genomic space. Enrichment = observed fraction / expected fraction,
    # where expected is that class's share of annotated bp. Values >1 mean the RBP
    # binds the class more than its size alone would predict.
    #
    # CAVEAT: genomic space is a crude background. eCLIP reads come from RNA, where
    # introns are far less abundant than in the genome, so this OVER-corrects for
    # introns. The correct background is the size-matched input (SMInput) library
    # for each experiment, which is what Skipper's beta-binomial uses. Treat these
    # columns as directional, not definitive.
    region_bp = {}
    for cls in PRECEDENCE:
        region_bp[cls] = int(sum(
            (arr[1] - arr[0]).sum() for arr in indexes[cls].idx.values()
        ))
    total_bp = sum(region_bp.values())
    print("\n=== length background (share of annotated genomic space) ===", flush=True)
    for cls in PRECEDENCE:
        exp = region_bp[cls] / total_bp
        print(f"  {cls:11s} {region_bp[cls]/1e6:9.1f} Mb  expected frac = {exp:.4f}", flush=True)
        out[f"enrich_{cls}"] = np.where(exp > 0, out[f"frac_{cls}"] / exp, np.nan)

    path = f"{args.out_prefix}_region_counts_by_file.csv"
    out.to_csv(path, index=False)
    print(f"\nwrote {path}  ({len(out)} rows)", flush=True)

    bad = out[(out.n_peaks > 0) & (out.frac_sum.sub(1).abs() > 1e-9)]
    print(f"rows whose fractions do not sum to 1: {len(bad)}  (must be 0)", flush=True)

    print("\n=== peaks retained by tier ===", flush=True)
    print(out.groupby("tier")["n_peaks"].agg(["count", "sum", "median"]).to_string(), flush=True)


if __name__ == "__main__":
    main()


