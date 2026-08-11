"""Calibrate an ENCORI significance threshold against ENCODE's effective depth.

The problem: ENCODE and ENCORI cannot share a filter by construction.
  * ENCODE's criterion is log2FC >= 3 AND -log10(p) >= 3 (Van Nostrand 2020),
    retaining 5.4% of CLIPper's reported clusters (median 4,474 peaks/file).
  * ENCORI has no fold-change column, its log10(p) is INTEGER-valued (one-decade
    resolution, huge ties), and its export has no null mass at all (0 sites with
    p > 0.32), so BH-FDR is not computable on it.
  * The same nominal p<=1e-3 retains 9.1% of ENCODE but 97.3% of ENCORI, because
    each arrives pre-filtered to a different depth.

So "consistent with ENCODE" cannot mean "same p-value". The defensible reading is
SAME EFFECTIVE DEPTH: pick the ENCORI threshold whose surviving site count per RBP
best matches ENCODE's Van Nostrand peak count for that same RBP.

This script makes one streaming pass over sites_annotated.bed, building a per-RBP
histogram of integer log10(p). From that histogram every candidate threshold can
be evaluated without re-reading the file, and both a single global threshold and
per-RBP thresholds are reported.

Outputs:
  encori_pvalue_histogram.csv   per-RBP counts at each integer log10(p)
  encori_calibration.csv        per-RBP: ENCODE target, best global/per-RBP cut
"""
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

ENCORI = "sites_annotated.bed"
ENCODE = "eclip_region_counts_by_file.csv"
MASTER = "rbp_dataframe.csv"

# candidate integer thresholds: keep sites with log10(p) <= T
THRESHOLDS = list(range(-2, -41, -1))


def stream_histogram(path, limit=0):
    """per-RBP Counter over integer log10(p), with adjacent-duplicate dedup."""
    hist = defaultdict(lambda: defaultdict(int))
    prev = None
    n_lines = n_sites = n_bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            n_lines += 1
            if limit and n_lines > limit:
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
                s = int(float(f[4]))
            except ValueError:
                n_bad += 1
                continue
            rbp = f[6].strip()
            if not rbp:
                n_bad += 1
                continue
            hist[rbp][s] += 1
            n_sites += 1
            if n_sites % 2_000_000 == 0:
                print(f"  {n_lines/1e6:5.1f}M lines -> {n_sites/1e6:4.1f}M sites", flush=True)
    print(f"  lines={n_lines:,} sites={n_sites:,} malformed={n_bad:,}", flush=True)
    return hist


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0

    # ---- ENCODE target depth per RBP --------------------------------------
    e = pd.read_csv(ENCODE)
    e = e[e.tier == "vannostrand"]
    # per RBP: median peaks across that RBP's replicate files (not sum -- we want
    # a per-experiment depth, comparable to one ENCORI RBP's site count)
    enc = e.groupby("rbp").n_peaks.median().rename("encode_target")
    print(f"ENCODE Van Nostrand targets: {len(enc)} RBPs, "
          f"median {enc.median():,.0f} peaks\n", flush=True)

    print(f"streaming {ENCORI} ...", flush=True)
    hist = stream_histogram(ENCORI, limit)
    print(f"  {len(hist)} distinct RBPs\n", flush=True)

    # ---- tidy histogram ----------------------------------------------------
    rows = []
    for rbp, h in hist.items():
        for score, n in h.items():
            rows.append({"rbp": rbp, "log10p": score, "n": n})
    H = pd.DataFrame(rows)
    H.to_csv("encori_pvalue_histogram.csv", index=False)

    # counts at each threshold: sites with log10p <= T
    def counts_at(T):
        return H[H.log10p <= T].groupby("rbp").n.sum()

    cum = pd.DataFrame({f"T{T}": counts_at(T) for T in THRESHOLDS}).fillna(0).astype(int)

    # ---- global threshold search ------------------------------------------
    shared = sorted(set(cum.index) & set(enc.index))
    print(f"RBPs in both ENCORI and ENCODE: {len(shared)}\n", flush=True)
    if not shared:
        sys.exit("no shared RBPs -- cannot calibrate")

    tgt = enc.reindex(shared)
    print("=== global threshold search (log-ratio to ENCODE depth) ===", flush=True)
    print(f"  {'thresh':>8} {'median ratio':>13} {'median |log2 ratio|':>21} {'median count':>13}",
          flush=True)
    best, best_score = None, np.inf
    for T in THRESHOLDS:
        c = cum.loc[shared, f"T{T}"].replace(0, np.nan)
        ratio = c / tgt
        score = np.nanmedian(np.abs(np.log2(ratio)))
        if T >= -12 or T % 5 == 0:
            print(f"  {T:>8} {np.nanmedian(ratio):>13.2f} {score:>21.3f} "
                  f"{np.nanmedian(c):>13,.0f}", flush=True)
        if score < best_score:
            best_score, best = score, T
    print(f"\n  BEST GLOBAL THRESHOLD: log10(p) <= {best}  "
          f"(median |log2 ratio| = {best_score:.3f}, "
          f"i.e. typically within {2**best_score:.2f}x of ENCODE depth)", flush=True)

    # ---- per-RBP threshold -------------------------------------------------
    per = {}
    for rbp in shared:
        t = tgt[rbp]
        diffs = {T: abs(np.log2((cum.loc[rbp, f"T{T}"] or np.nan) / t))
                 for T in THRESHOLDS if cum.loc[rbp, f"T{T}"] > 0}
        per[rbp] = min(diffs, key=diffs.get) if diffs else np.nan

    out = pd.DataFrame({
        "encode_target": tgt,
        "encori_at_best_global": cum.loc[shared, f"T{best}"],
        "best_per_rbp_threshold": pd.Series(per),
    })
    out["encori_at_per_rbp"] = [
        cum.loc[r, f"T{int(out.best_per_rbp_threshold[r])}"]
        if pd.notna(out.best_per_rbp_threshold[r]) else np.nan for r in out.index
    ]
    out["ratio_global"] = out.encori_at_best_global / out.encode_target
    out["ratio_per_rbp"] = out.encori_at_per_rbp / out.encode_target
    out.to_csv("encori_calibration.csv")

    print("\n=== how well does each strategy match ENCODE depth? ===", flush=True)
    for col in ("ratio_global", "ratio_per_rbp"):
        r = out[col].replace([np.inf, -np.inf], np.nan).dropna()
        within2 = float(((r >= 0.5) & (r <= 2)).mean())
        print(f"  {col:14s} median {r.median():.2f}x  IQR {r.quantile(.25):.2f}-{r.quantile(.75):.2f}"
              f"  within 2x of ENCODE: {within2:.0%}", flush=True)

    print("\n=== per-RBP thresholds chosen (distribution) ===", flush=True)
    print(out.best_per_rbp_threshold.value_counts().sort_index(ascending=False).head(15).to_string(),
          flush=True)
    print("\nwrote encori_pvalue_histogram.csv, encori_calibration.csv", flush=True)


if __name__ == "__main__":
    main()


