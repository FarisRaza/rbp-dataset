"""Compare the four criterion families, all on even-split region assignment.

    encode_published : log2FC >= 3 AND -log10(p) >= 3          (Van Nostrand 2020)
    encode_matched   : -log10(p) >= 3 AND >= 2 replicates
    encori_published : P < 1e-5 (single) OR >= 50% exps (min 2), >= 300 regions
                                                                (Zhou 2026)
    encori_matched   : log10(p) <= -3 AND >= 2 datasets

Because every family now uses ENCORI's even-split assignment, the frac_* columns
are directly comparable across all four -- only the significance criterion varies.
The key question this answers: does the ENCODE/ENCORI agreement seen under matched
criteria survive, and how much does each source's own published criterion differ
from the matched one?
"""
import os

import numpy as np
import pandas as pd

CLS = ["utr3", "utr5", "cds", "intron", "ncrna_exon", "intergenic"]

FAMILIES = [
    ("encode_published", "encode_region_counts_published.csv"),
    ("encode_matched",   "encode_region_counts_matched.csv"),
    ("encori_published", "encori_region_counts_published.csv"),
    ("encori_matched",   "encori_region_counts_matched.csv"),
]

CONTROLS = ["AGO2", "HNRNPC", "DDX3X", "TARDBP", "U2AF2", "ELAVL1", "PTBP1", "RBFOX2"]


def load(path):
    d = pd.read_csv(path).set_index("rbp")
    f = d[[f"frac_{c}" for c in CLS]]
    f.columns = CLS
    return f, d


def main():
    loaded = {}
    for name, path in FAMILIES:
        if not os.path.exists(path):
            print(f"MISSING {path} -- skipping {name}")
            continue
        loaded[name] = load(path)

    if not loaded:
        return

    print("=== coverage and mean composition (even-split throughout) ===")
    print(f"{'family':18s} {'RBPs':>5} " + " ".join(f"{c:>7s}" for c in CLS))
    for name, (f, _) in loaded.items():
        print(f"{name:18s} {len(f):>5} " + " ".join(f"{f[c].mean():7.3f}" for c in CLS))

    print("\n=== RBPs by dominant region ===")
    print(f"{'family':18s} " + " ".join(f"{c:>10s}" for c in CLS))
    for name, (f, _) in loaded.items():
        d = f.idxmax(axis=1).value_counts()
        print(f"{name:18s} " + " ".join(f"{int(d.get(c,0)):>10d}" for c in CLS))

    # published vs matched, within each source
    for src in ("encode", "encori"):
        a, b = f"{src}_published", f"{src}_matched"
        if a not in loaded or b not in loaded:
            continue
        fa, fb = loaded[a][0], loaded[b][0]
        shared = sorted(set(fa.index) & set(fb.index))
        if not shared:
            continue
        print(f"\n=== {src}: published vs matched  (n={len(shared)} shared RBPs) ===")
        for c in ["utr3", "utr5", "cds", "intron"]:
            x, y = fa.loc[shared, c], fb.loc[shared, c]
            print(f"  {c:11s} published {x.mean():.3f}  matched {y.mean():.3f}  "
                  f"delta {y.mean()-x.mean():+.3f}  r={np.corrcoef(x, y)[0,1]:.3f}")

    # the headline: does ENCODE agree with ENCORI now that both use even-split?
    for tag in ("published", "matched"):
        a, b = f"encode_{tag}", f"encori_{tag}"
        if a not in loaded or b not in loaded:
            continue
        fa, fb = loaded[a][0], loaded[b][0]
        shared = sorted(set(fa.index) & set(fb.index))
        if not shared:
            continue
        print(f"\n=== ENCODE vs ENCORI, both '{tag}'  (n={len(shared)}) ===")
        for c in ["utr3", "utr5", "cds", "intron"]:
            x, y = fa.loc[shared, c], fb.loc[shared, c]
            print(f"  {c:11s} ENCODE {x.mean():.3f}  ENCORI {y.mean():.3f}  "
                  f"gap {abs(x.mean()-y.mean()):.3f}  r={np.corrcoef(x, y)[0,1]:.3f}")
        da = fa.loc[shared].idxmax(axis=1)
        db = fb.loc[shared].idxmax(axis=1)
        print(f"  dominant-region concordance: {(da==db).mean():.1%}")

    print("\n=== positive controls (intron / 3'UTR / 5'UTR) ===")
    hdr = " ".join(f"{n.replace('_','.'):>22s}" for n in loaded)
    print(f"{'RBP':9s} {hdr}")
    for r in CONTROLS:
        cells = []
        for name, (f, _) in loaded.items():
            if r in f.index:
                cells.append(f"{f.loc[r,'intron']:6.3f}/{f.loc[r,'utr3']:5.3f}/{f.loc[r,'utr5']:5.3f}")
            else:
                cells.append(f"{'--':>22s}")
        print(f"{r:9s} " + " ".join(f"{c:>22s}" for c in cells))


if __name__ == "__main__":
    main()


