"""Compute Benjamini-Hochberg FDR for ENCODE eCLIP peaks, since ENCODE does not.

ENCODE's eCLIP narrowPeak leaves column 9 (qValue) and column 10 (peak) set to -1
in every released file, so there is no published FDR. Column 8 carries -log10(p)
from CLIPper's input-normalisation test; this script turns that into a q-value.

Two implementation points:

  * -log10(p) is capped at 400, i.e. p = 1e-400, which underflows float64
    (min normal ~1e-308). All arithmetic is therefore done in log10 space:
        log10(q_i) = min over j>=i of [ log10(p_j) + log10(n) - log10(j) ]
    computed as a reverse cumulative minimum over p-sorted values.

  * The implementation is checked for EXACT agreement against
    statsmodels.stats.multitest.multipletests(method='fdr_bh') on the subset
    where p does not underflow, so the log-space version is verified rather
    than assumed.

CAVEAT THAT MATTERS MORE THAN THE ARITHMETIC
--------------------------------------------
BH assumes the supplied p-values are the FULL set of tests, with nulls uniform
on [0,1]. CLIPper reports only candidate clusters it already considered
peak-like, so the null tests have largely been removed before we ever see them.
The observed p distribution is wildly non-uniform (only ~5% of entries have
p > 0.32, where a uniform null predicts ~68%). BH on a pre-filtered set is
ANTI-CONSERVATIVE: it divides by too small an n and reports an FDR that is
optimistic by an unknown factor. The script quantifies the non-uniformity so
the size of that problem is visible rather than hidden.
"""
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd

NARROWPEAK = ["chrom", "start", "end", "name", "score", "strand",
              "log2fc", "neglog10p", "qvalue", "peak"]


def bh_log10(neglog10p):
    """BH q-values in log10 space. Input -log10(p); returns -log10(q)."""
    n = len(neglog10p)
    order = np.argsort(-neglog10p, kind="stable")   # p ascending == -log10p descending
    nl = neglog10p[order]

    rank = np.arange(1, n + 1, dtype=np.float64)
    # log10(q_j) = log10(p_j) + log10(n) - log10(j);  p_j = 10**(-nl_j)
    log10q = -nl + np.log10(n) - np.log10(rank)
    # step-up: q_i = min over j >= i  -> reverse cumulative minimum
    log10q = np.minimum.accumulate(log10q[::-1])[::-1]
    log10q = np.minimum(log10q, 0.0)                # q <= 1

    out = np.empty(n, dtype=np.float64)
    out[order] = -log10q                            # return as -log10(q)
    return out


def verify_against_statsmodels(neglog10p, tol=1e-9):
    """Exact-fit check on the non-underflowing subset."""
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        return None

    keep = neglog10p < 300                          # p >= 1e-300, safe in float64
    if keep.sum() < 100:
        return None
    sub = neglog10p[keep]
    p = np.power(10.0, -sub)

    _, q_sm, _, _ = multipletests(p, method="fdr_bh")
    q_mine = np.power(10.0, -bh_log10(sub))

    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(q_mine - q_sm) / np.where(q_sm > 0, q_sm, 1.0)
    return {
        "n_compared": int(keep.sum()),
        "max_abs_diff": float(np.max(np.abs(q_mine - q_sm))),
        "max_rel_diff": float(np.nanmax(rel)),
        "exact": bool(np.max(np.abs(q_mine - q_sm)) < tol),
    }


def uniformity(neglog10p):
    """How far is the p distribution from the uniform null BH assumes?

    Under a uniform null, P(p > 0.32) = 0.68. Departure quantifies how much
    pre-filtering happened upstream.
    """
    p_gt_032 = float((neglog10p < 0.5).mean())
    p_gt_01 = float((neglog10p < 1.0).mean())
    return {"frac_p_gt_0.32": p_gt_032, "expected_if_uniform": 0.68,
            "frac_p_gt_0.1": p_gt_01, "expected_if_uniform_0.1": 0.90}


def load(path):
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", header=None, names=NARROWPEAK)


def main():
    peaks_dir = sys.argv[1] if len(sys.argv) > 1 else "encode_peaks_grch38"
    files = sorted(glob.glob(os.path.join(peaks_dir, "*_replicate_*.bed.gz")))
    if not files:
        sys.exit(f"no replicate peak files in {peaks_dir}")

    print(f"{len(files)} replicate peak files\n", flush=True)

    # ---- exact-fit verification on the first few files --------------------
    print("=== implementation check vs statsmodels fdr_bh ===", flush=True)
    for path in files[:3]:
        df = load(path)
        v = verify_against_statsmodels(df.neglog10p.to_numpy(dtype=np.float64))
        name = os.path.basename(path)
        if v is None:
            print(f"  {name}: skipped (statsmodels missing or too few rows)", flush=True)
        else:
            print(f"  {name}: n={v['n_compared']:,} max_abs_diff={v['max_abs_diff']:.3e} "
                  f"max_rel_diff={v['max_rel_diff']:.3e} EXACT={v['exact']}", flush=True)

    # ---- null-uniformity diagnostic ---------------------------------------
    print("\n=== is the p-value set BH-eligible? (uniform null check) ===", flush=True)
    for path in files[:3]:
        df = load(path)
        u = uniformity(df.neglog10p.to_numpy(dtype=np.float64))
        print(f"  {os.path.basename(path)}", flush=True)
        print(f"    frac(p>0.32) = {u['frac_p_gt_0.32']:.3f}  (uniform null predicts "
              f"{u['expected_if_uniform']:.2f})", flush=True)
        print(f"    frac(p>0.10) = {u['frac_p_gt_0.1']:.3f}  (uniform null predicts "
              f"{u['expected_if_uniform_0.1']:.2f})", flush=True)

    # ---- FDR across all files ---------------------------------------------
    print("\n=== computing BH q-values for all files ===", flush=True)
    rows = []
    for i, path in enumerate(files, 1):
        df = load(path)
        if not len(df):
            continue
        nl = df.neglog10p.to_numpy(dtype=np.float64)
        nlq = bh_log10(nl)

        base = os.path.basename(path)
        stem = base[: -len(".bed.gz")]
        parts = stem.split("_")
        rbp, biosample = parts[0], parts[1]

        rows.append({
            "file": base, "rbp": rbp, "biosample": biosample, "n_peaks": len(df),
            # how many pass each criterion
            "n_p_lt_1e3": int((nl >= 3).sum()),
            "n_p_lt_1e4": int((nl > 4).sum()),
            "n_vannostrand": int(((df.log2fc >= 3) & (nl >= 3)).sum()),
            "n_q_lt_0.05": int((nlq >= -np.log10(0.05)).sum()),
            "n_q_lt_0.01": int((nlq >= 2).sum()),
            "n_q_lt_0.001": int((nlq >= 3).sum()),
            # what q does the Van Nostrand p cutoff correspond to?
            "q_at_p_1e3": float(np.power(10.0, -nlq[nl >= 3].min())) if (nl >= 3).any() else np.nan,
            "frac_p_gt_032": float((nl < 0.5).mean()),
        })
        if i % 100 == 0 or i == len(files):
            print(f"  {i}/{len(files)}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("encode_fdr_summary.csv", index=False)

    print("\n=== peaks passing each criterion (median per file) ===", flush=True)
    for c in ["n_peaks", "n_p_lt_1e3", "n_p_lt_1e4", "n_vannostrand",
              "n_q_lt_0.05", "n_q_lt_0.01", "n_q_lt_0.001"]:
        print(f"  {c:16s} {out[c].median():>10,.0f}", flush=True)

    print("\n=== q-value implied by the p<=1e-3 (Van Nostrand) cutoff ===", flush=True)
    print(f"  median q at that cutoff: {out.q_at_p_1e3.median():.3e}", flush=True)
    print(f"  IQR: {out.q_at_p_1e3.quantile(.25):.3e} - {out.q_at_p_1e3.quantile(.75):.3e}",
          flush=True)

    print("\n=== null-uniformity across all files ===", flush=True)
    print(f"  median frac(p>0.32) = {out.frac_p_gt_032.median():.3f}  "
          f"(uniform null predicts 0.68)", flush=True)
    print("\nwrote encode_fdr_summary.csv", flush=True)


if __name__ == "__main__":
    main()


