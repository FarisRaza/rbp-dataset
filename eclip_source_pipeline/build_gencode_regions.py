"""Build merged, non-overlapping genomic region tracks from GENCODE v46 (GRCh38).

Fixes two defects in the existing gencode_3utr.bed:

  1. NOT MERGED. That file has 214,776 intervals but only 141,864 unique coordinates,
     because every transcript isoform contributed its own record. bedtools intersect
     then emitted one line per overlapping isoform, so a peak over a multi-isoform
     3'UTR was counted 2-4x. That is the source of rbp_3utr_annotation.csv reporting
     frac_3UTR = 3.50 (SUB1: 38,999 "3'UTR peaks" out of 11,143 total).

  2. UTR SIDE WAS ASSUMED, NOT DERIVED. GENCODE's GTF has only a generic `UTR`
     feature (336,514 of them) -- five_prime_UTR / three_prime_UTR appear in the
     GFF3, not the GTF. So 5' vs 3' must be derived strand-aware relative to the
     transcript's CDS span. This script does that explicitly:

         + strand:  UTR entirely before CDS start -> 5'UTR ; after CDS end -> 3'UTR
         - strand:  UTR entirely after  CDS end   -> 5'UTR ; before CDS start -> 3'UTR

Emits one merged BED per region class. Merging is per (chromosome, strand) so that
overlapping isoform annotations collapse to a single interval and each genomic base
is represented once within its class.

The GTF is streamed transcript-by-transcript (GENCODE is ordered gene -> transcript ->
features), so memory stays flat regardless of the 1.6 GB input.

Region classes emitted:
    cds, utr5, utr3, intron, ncrna_exon
Precedence between classes is applied downstream, not here -- these are raw per-class
tracks so the precedence rule can be changed without re-parsing the GTF.
"""
import os
import sys
from collections import defaultdict

GTF = "gencode.gtf"
OUTDIR = "gencode_regions"
PREFIX = "gencode_v46"

os.makedirs(OUTDIR, exist_ok=True)

# gene_types treated as non-coding for the ncrna_exon track
NONCODING_HINTS = (
    "lncRNA", "miRNA", "snRNA", "snoRNA", "scaRNA", "rRNA", "misc_RNA",
    "vault_RNA", "sRNA", "scRNA", "ribozyme", "Mt_rRNA", "Mt_tRNA",
)


def attr(s, key):
    """Pull a GTF attribute value without regex (much faster on 1.6 GB)."""
    i = s.find(key + ' "')
    if i == -1:
        return None
    i += len(key) + 2
    j = s.find('"', i)
    return s[i:j]


def merge(intervals):
    """Merge overlapping/book-ended intervals. intervals: list of (start, end)."""
    if not intervals:
        return []
    intervals.sort()
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= out[-1][1]:          # overlap or adjacent
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


def subtract(span, blocks):
    """span=(s,e) minus a sorted merged list of blocks -> list of gap intervals."""
    s, e = span
    gaps = []
    cur = s
    for bs, be in blocks:
        if be <= cur:
            continue
        if bs >= e:
            break
        if bs > cur:
            gaps.append((cur, min(bs, e)))
        cur = max(cur, be)
        if cur >= e:
            break
    if cur < e:
        gaps.append((cur, e))
    return gaps


# accumulators: class -> (chrom, strand) -> [(start, end), ...]
tracks = {k: defaultdict(list) for k in ("cds", "utr5", "utr3", "intron", "ncrna_exon")}

cur_tid = None
cur = None


def flush(t):
    """Resolve one transcript's features into region classes."""
    if t is None or not t["chrom"]:
        return
    chrom, strand = t["chrom"], t["strand"]
    key = (chrom, strand)

    exons = merge(t["exon"])
    cds = merge(t["cds"])

    if cds:
        for s, e in cds:
            tracks["cds"][key].append((s, e))

        cds_lo = cds[0][0]
        cds_hi = cds[-1][1]

        for s, e in t["utr"]:
            if e <= cds_lo:
                side = "utr5" if strand == "+" else "utr3"
            elif s >= cds_hi:
                side = "utr3" if strand == "+" else "utr5"
            else:
                # UTR overlapping the CDS span -- shouldn't occur in GENCODE; skip
                # rather than guess, and count it for the report.
                t["odd_utr"] += 1
                continue
            tracks[side][key].append((s, e))
    else:
        # no CDS -> non-coding transcript; its exons form the ncRNA track
        if t["is_noncoding"]:
            for s, e in exons:
                tracks["ncrna_exon"][key].append((s, e))

    # introns = transcript span minus exons
    if t["span"] and exons:
        for s, e in subtract(t["span"], exons):
            tracks["intron"][key].append((s, e))


odd_utr_total = 0
n_tx = 0

print(f"streaming {GTF} ...", flush=True)
with open(GTF, "r", encoding="utf-8") as fh:
    for lineno, line in enumerate(fh, 1):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat = f[2]
        if feat not in ("transcript", "exon", "CDS", "UTR"):
            continue

        chrom, start, end, strand, a = f[0], int(f[3]) - 1, int(f[4]), f[6], f[8]
        tid = attr(a, "transcript_id")
        if tid is None:
            continue

        if tid != cur_tid:
            flush(cur)
            odd_utr_total += cur["odd_utr"] if cur else 0
            n_tx += 1 if cur else 0
            gt = attr(a, "gene_type") or ""
            tt = attr(a, "transcript_type") or ""
            cur = {
                "chrom": chrom, "strand": strand, "span": None,
                "exon": [], "cds": [], "utr": [], "odd_utr": 0,
                "is_noncoding": (gt in NONCODING_HINTS) or (tt in NONCODING_HINTS),
            }
            cur_tid = tid

        if feat == "transcript":
            cur["span"] = (start, end)
        elif feat == "exon":
            cur["exon"].append((start, end))
        elif feat == "CDS":
            cur["cds"].append((start, end))
        elif feat == "UTR":
            cur["utr"].append((start, end))

        if lineno % 5_000_000 == 0:
            print(f"  {lineno//1_000_000}M lines", flush=True)

flush(cur)
odd_utr_total += cur["odd_utr"] if cur else 0
n_tx += 1
print(f"parsed {n_tx:,} transcripts", flush=True)
if odd_utr_total:
    print(f"  note: {odd_utr_total:,} UTR features overlapped their CDS span (skipped)", flush=True)

# ---- merge each class across all transcripts and write BED -------------------
print("\n=== merged region tracks ===", flush=True)
summary = []
for cls, bykey in tracks.items():
    path = os.path.join(OUTDIR, f"{PREFIX}_{cls}.bed")
    n_raw = n_merged = bp = 0
    with open(path, "w", encoding="utf-8") as out:
        for (chrom, strand) in sorted(bykey):
            ivs = bykey[(chrom, strand)]
            n_raw += len(ivs)
            for s, e in merge(ivs):
                out.write(f"{chrom}\t{s}\t{e}\t{cls}\t0\t{strand}\n")
                n_merged += 1
                bp += e - s
    summary.append((cls, n_raw, n_merged, bp))
    print(
        f"  {cls:11s} raw={n_raw:>9,}  merged={n_merged:>9,}  "
        f"({100*(1-n_merged/n_raw) if n_raw else 0:5.1f}% collapsed)  {bp/1e6:9.1f} Mb",
        flush=True,
    )

print(f"\nwrote {len(summary)} tracks to {OUTDIR}/", flush=True)


