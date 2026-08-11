"""Annotate ENCORI and POSTAR3 binding sites to genomic regions, per RBP.

These two sources need different handling from ENCODE eCLIP:

  * Neither retains per-peak statistics. POSTAR3's score column is all zeros;
    ENCORI's is negative integers (min -209,955), not a p-value or fold change.
    So the Van Nostrand "log2FC >= 3 and -log10(p) >= 3" filter CANNOT be applied
    here -- these are curated aggregations of already-called peaks whose
    thresholding happened upstream in each source publication.

  * ENCORI's own region labels are only {CDS, UTR, intron} -- UTR is NOT split
    into 5' and 3'. Since the 5'/3' distinction is the point, sites are
    re-annotated against the merged GENCODE v46 tracks instead of trusting
    column 14.

  * sites_annotated.bed is ~2.85x redundant: bedtools emitted one line per
    overlapping gene annotation, so each site repeats. Duplicates are contiguous
    (verified: adjacent-unique == global-unique), so dedup happens in a single
    streaming pass with flat memory over the 4.5 GB file.

  * postar3_sites.bed has CRLF line endings, putting a trailing \\r on the RBP
    column. Unstripped, every exact name match against the master table fails.

Both files are streamed and annotated in vectorised batches.

Usage:
    python annotate_db_sites.py --source encori
    python annotate_db_sites.py --source postar3
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from annotate_peaks import INTERGENIC, PRECEDENCE, IntervalIndex

SOURCES = {
    "encori": {
        "path": "sites_annotated.bed",
        # chrom start end site_id score strand RBP | gene_chrom gene_start gene_end
        # gene_name . gene_strand region
        "cols": {"chrom": 0, "start": 1, "end": 2, "site_id": 3, "score": 4,
                 "strand": 5, "rbp": 6},
        "dedup_fields": 7,
        "min_fields": 14,
        # Column 5 is log10(p-value) from rbsSeeker -- ENCORI's "CLIP Region
        # P-value". Verified empirically: values are strictly <= -1, unimodal
        # with a mode at -4, decaying smoothly to -209,955, and zero mass above
        # -1. So p < 1e-4 is score <= -4. (The ENCORI bindingSite API returns
        # the identical 6-column BED with this same field.)
        "has_pvalue": True,
    },
    "postar3": {
        "path": "postar3_sites.bed",
        # chrom start end site_id score strand RBP
        "cols": {"chrom": 0, "start": 1, "end": 2, "site_id": 3, "score": 4,
                 "strand": 5, "rbp": 6},
        "dedup_fields": 7,
        "min_fields": 7,
        # POSTAR3's score column is genuinely uninformative -- every value is 0.
        "has_pvalue": False,
    },
}

BATCH = 1_000_000


def annotate_batch(chrom, strand, starts, ends, indexes):
    """Assign each site in one (chrom, strand) block to exactly one region class."""
    n = len(starts)
    label = np.full(n, INTERGENIC, dtype=object)
    unassigned = np.ones(n, dtype=bool)
    for cls in PRECEDENCE:
        if not unassigned.any():
            break
        hit = indexes[cls].overlaps(chrom, strand, starts, ends) & unassigned
        if hit.any():
            label[hit] = cls
            unassigned &= ~hit
    return label


def flush_batch(buf, indexes, counts, source_counts):
    """buf: dict (chrom,strand) -> [starts, ends, rbps, srcs]"""
    for (chrom, strand), (s, e, rbps, srcs) in buf.items():
        s = np.asarray(s, dtype=np.int64)
        e = np.asarray(e, dtype=np.int64)
        lab = annotate_batch(chrom, strand, s, e, indexes)
        for rbp, cls, src in zip(rbps, lab, srcs):
            counts[rbp][cls] += 1
            if src is not None:
                source_counts[rbp][src] += 1
    buf.clear()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, choices=list(SOURCES))
    ap.add_argument("--regions", default="gencode_regions")
    ap.add_argument("--limit", type=int, default=0, help="stop after N lines (testing)")
    ap.add_argument(
        "--max-log10p", type=float, default=None,
        help="keep sites with log10(p) <= this. ENCORI only. Use -4 for p<1e-4, "
             "-3 for p<1e-3. Omit for no significance filter.",
    )
    ap.add_argument("--tag", default="", help="suffix for output filenames")
    # NOTE: a cross-dataset reproducibility filter ("site supported by >= 2
    # independent datasets") is deliberately NOT implemented here. Sites from
    # different SBDH datasets rarely share exact coordinates, so it requires
    # interval merging across datasets -- a genuine second pass, not a row
    # predicate. Shipping it as a flag here would silently do the wrong thing.
    # The per-RBP dataset counts are reported instead.
    args = ap.parse_args()

    if args.max_log10p is not None and not SOURCES[args.source]["has_pvalue"]:
        sys.exit(f"--max-log10p not usable for {args.source}: it has no p-value column")

    cfg = SOURCES[args.source]
    if not os.path.exists(cfg["path"]):
        sys.exit(f"missing {cfg['path']}")

    print(f"loading region tracks from {args.regions}/ ...", flush=True)
    indexes = {}
    for cls in PRECEDENCE:
        p = os.path.join(args.regions, f"gencode_v46_{cls}.bed")
        if not os.path.exists(p):
            sys.exit(f"missing region track: {p} (run build_gencode_regions.py first)")
        indexes[cls] = IntervalIndex(p)
    print("  ok", flush=True)

    c = cfg["cols"]
    nd = cfg["dedup_fields"]

    counts = defaultdict(lambda: defaultdict(int))
    source_counts = defaultdict(lambda: defaultdict(int))
    datasets = defaultdict(set)
    buf = defaultdict(lambda: ([], [], [], []))
    buffered = 0

    prev_key = None
    n_lines = n_sites = n_bad = n_filtered = 0

    print(f"streaming {cfg['path']} ...", flush=True)
    with open(cfg["path"], "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            n_lines += 1
            if args.limit and n_lines > args.limit:
                break

            # CRLF-safe: strip \r as well as \n (POSTAR3 is CRLF)
            f = line.rstrip("\r\n").split("\t")
            if len(f) < cfg["min_fields"]:
                n_bad += 1
                continue

            # streaming dedup -- duplicates are contiguous
            key = tuple(f[:nd])
            if key == prev_key:
                continue
            prev_key = key

            try:
                s = int(f[c["start"]])
                e = int(f[c["end"]])
            except ValueError:
                n_bad += 1
                continue

            # significance filter (ENCORI only; col5 is log10(p))
            if args.max_log10p is not None:
                try:
                    if float(f[c["score"]]) > args.max_log10p:
                        n_filtered += 1
                        continue
                except ValueError:
                    n_bad += 1
                    continue

            rbp = f[c["rbp"]].strip()
            if not rbp:
                n_bad += 1
                continue

            # dataset provenance: ENCORI SBDH<n>, POSTAR3 source family
            if args.source == "encori":
                datasets[rbp].add(f[c["site_id"]].split("-")[0])

            src = None
            if args.source == "postar3":
                # id encodes provenance, e.g. human_RBP_eCLIP_ENCODE_1270957
                parts = f[c["site_id"]].split("_")
                src = "_".join(parts[:3]) if len(parts) >= 3 else f[c["site_id"]]

            k = (f[c["chrom"]], f[c["strand"]])
            b = buf[k]
            b[0].append(s)
            b[1].append(e)
            b[2].append(rbp)
            b[3].append(src)
            buffered += 1
            n_sites += 1

            if buffered >= BATCH:
                flush_batch(buf, indexes, counts, source_counts)
                buffered = 0
                print(f"  {n_lines/1e6:6.1f}M lines -> {n_sites/1e6:5.2f}M unique sites", flush=True)

    flush_batch(buf, indexes, counts, source_counts)

    print(f"\nlines={n_lines:,}  unique sites={n_sites:,}  "
          f"redundancy={n_lines/max(n_sites,1):.2f}x  malformed={n_bad:,}", flush=True)
    if args.max_log10p is not None:
        kept = n_sites / max(n_sites + n_filtered, 1)
        print(f"significance filter log10(p) <= {args.max_log10p} "
              f"(p <= 1e{int(args.max_log10p)}): dropped {n_filtered:,} rows, "
              f"kept {kept:.1%} of unique sites", flush=True)

    classes = PRECEDENCE + [INTERGENIC]
    rows = []
    for rbp, d in counts.items():
        tot = sum(d.values())
        rec = {"rbp": rbp, "source": args.source, "n_sites": tot,
               "n_datasets": len(datasets.get(rbp, ()))}
        for cls in classes:
            rec[f"n_{cls}"] = d.get(cls, 0)
        for cls in classes:
            rec[f"frac_{cls}"] = d.get(cls, 0) / tot if tot else 0.0
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values("n_sites", ascending=False)
    out["frac_sum"] = out[[f"frac_{c}" for c in classes]].sum(axis=1)

    path = f"{args.source}_region_counts{args.tag}.csv"
    out.to_csv(path, index=False)
    print(f"wrote {path}  ({len(out)} RBPs)", flush=True)

    bad = out[(out.n_sites > 0) & (out.frac_sum.sub(1).abs() > 1e-9)]
    print(f"RBPs whose fractions do not sum to 1: {len(bad)}  (must be 0)", flush=True)

    if source_counts:
        srows = []
        for rbp, d in source_counts.items():
            for src, n in d.items():
                srows.append({"rbp": rbp, "site_source": src, "n": n})
        sp = f"{args.source}_site_provenance.csv"
        pd.DataFrame(srows).to_csv(sp, index=False)
        print(f"wrote {sp}", flush=True)


if __name__ == "__main__":
    main()


