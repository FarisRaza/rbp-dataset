"""Join all eCLIP sources onto the master RBP table.

Design decisions (set by the user):
  * Per-source column families, NOT a merged consensus. ENCODE eCLIP is a subset
    of both ENCORI and POSTAR3 (POSTAR3 ids literally include
    human_RBP_eCLIP_ENCODE_*), so treating them as independent evidence would be
    pseudo-replication. Keeping them separate also preserves genuine
    disagreements -- e.g. AGO2 is 0.52 3'UTR in ENCORI but 0.25 in POSTAR3, a gap
    that survives significance filtering and is therefore real.
  * Skipper and my own ENCODE pipeline both kept, as parallel families, so the
    two methods can be checked against each other (they correlate r=0.84-0.93).

Join key: master `Name` (1,392 unique HGNC-style symbols, no duplicates).
Aliases are expanded from the master's own `gene_name` column, which is a
UniProt-style space-separated synonym list (e.g. "AARS1 AARS"), plus an optional
HGNC-resolved override file for names that column does not cover.

COUNT COLUMNS ARE NOT COMPARABLE ACROSS SOURCES AND ARE EMITTED AS QC ONLY:
  * encori_n_sites tracks how many datasets exist for that RBP (TARDBP has 100
    ENCORI datasets and 1.24M sites; a once-profiled RBP has a few thousand). It
    is a publication-effort proxy, not a binding measure.
  * encode_n_peaks tracks sequencing depth.
  * The same nominal p<=1e-3 threshold retains 9.1% of ENCODE but 97.3% of
    ENCORI, because each source arrives pre-filtered to a different depth.
Use the frac_* and enrich_* columns for analysis. Counts are for provenance.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

MASTER = "rbp_dataframe.csv"
ALIAS_OVERRIDE = "rbp_alias_overrides.json"   # optional: {"SOURCE_NAME": "MasterName"}
OUT = "rbp_master_with_eclip.csv"
UNMATCHED_REPORT = "eclip_unmatched_rbps.csv"

# ---------------------------------------------------------------------------
# Significance criteria, stated explicitly rather than inferred from which files
# happen to exist on disk. An earlier version picked the ENCORI input by
# os.path.exists() fallback, which would have silently used the UNFILTERED file
# had the filtered one been absent.
#
# ENCODE and ENCORI are held to the SAME significance threshold, p <= 1e-3:
#   ENCODE  -log10(p) >= 3  (plus log2FC >= 3; Van Nostrand et al. 2020, Nature)
#   ENCORI   log10(p) <= -3 (rbsSeeker, column 5 of sites_annotated.bed)
# This matches the one criterion the two sources share. ENCORI has NO
# fold-change column, so ENCODE's effect-size condition has no counterpart and
# is simply absent there -- documented, not worked around.
#
# Deliberately NOT depth-matched: choosing an ENCORI threshold to reproduce
# ENCODE's peak count is post-hoc, and empirically fails anyway (best global
# threshold lands within ~7x of ENCODE depth; only 18% of RBPs within 2x).
# Equal thresholds do NOT mean equal stringency -- p<=1e-3 retains 9.1% of
# ENCODE but 97.6% of ENCORI, because each arrives pre-filtered to a different
# depth upstream. That asymmetry is inherent and is reported, not hidden.
# Four significance families, all sharing ENCORI's EVEN-SPLIT region assignment
# so only the criterion differs and the frac_* columns are directly comparable.
#
#   *_published : each source's own published criterion
#   *_matched   : a common criterion applied to both (p <= 1e-3 AND >= 2 experiments)
#
# Even-split replaced precedence throughout. Precedence ordered utr3 first, which
# inflated 3'UTR-dominant RBP counts roughly 2x (ENCORI: 82 -> 39 once corrected).
SOURCE_FILES = [
    # (prefix, path, count-column name)
    ("encode_published", "encode_region_counts_published.csv", "n_regions"),
    ("encode_matched",   "encode_region_counts_matched.csv",   "n_regions"),
    ("encori_published", "encori_region_counts_published.csv", "n_regions"),
    ("encori_matched",   "encori_region_counts_matched.csv",   "n_regions"),
    ("postar",           "postar3_region_counts.csv",          "n_sites"),
]

# Union families: ENCODE-PRIORITY over the RBP set of both sources.
#
# For an RBP present in both, ENCODE's profile wins. Two independent reasons:
#   * ENCORI pools protocols and cell lines per RBP, which can average away a
#     clean signal. EIF4E is 5'UTR-dominant in ENCODE (0.418) but intron-dominant
#     in ENCORI -- because ENCORI's EIF4E profile is 2 PAR-CLIP datasets in
#     HEK293T with ZERO eCLIP. Same story for EIF3G (2/5 eCLIP) and DDX3X
#     (2/15 eCLIP, 12/15 HEK293).
#   * Pooling raw region counts instead would let ENCORI outvote ENCODE on every
#     shared RBP purely on volume (6.07M vs 2.32M regions), which is a volume
#     artifact, not evidence. Under that scheme DDX3X and EIF3G lose their
#     5'UTR dominance; under ENCODE-priority the 5'UTR set is DDX3X, EIF3G,
#     EIF4E, NCBP2 -- all four cap/5'-end factors, which is the correct biology.
#
# The value of the union is COVERAGE EXTENSION (297/303 RBPs vs ENCODE's 168),
# not arbitration. `union_*_source` records which source supplied each row.
UNION_FAMILIES = [
    ("union_published", "encode_published", "encori_published"),
    ("union_matched",   "encode_matched",   "encori_matched"),
]

CRITERIA_DOC = {
    "encode_published": "log2FC >= 3 AND -log10(p) >= 3 (Van Nostrand et al. 2020)",
    "encode_matched":   "-log10(p) >= 3 AND locus called in >= 2 biological replicates",
    "encori_published": "P < 1e-5 (single exp) OR >= 50% of experiments (min 2); "
                        ">= 300 regions/RBP (Zhou et al. 2026)",
    "encori_matched":   "log10(p) <= -3 AND locus called in >= 2 SBDH datasets",
    "postar":           "UNFILTERED -- POSTAR3's score column is all zeros, no statistic exists",
    "skipper":          "beta-binomial vs SMInput + replicate reproducibility (published output)",
}

CLASSES = ["utr3", "utr5", "cds", "ncrna_exon", "intron", "intergenic"]


def build_alias_map(m):
    """Any known synonym (upper-cased) -> canonical master Name."""
    alias = {}
    for _, r in m.iterrows():
        canon = str(r["Name"]).strip()
        alias.setdefault(canon.upper(), canon)
        for col in ("UNIQUE", "uniprot_name"):
            v = r.get(col)
            if isinstance(v, str) and v.strip():
                alias.setdefault(v.strip().upper(), canon)
        gn = r.get("gene_name")
        if isinstance(gn, str):
            for tok in gn.split():
                alias.setdefault(tok.strip().upper(), canon)

    if os.path.exists(ALIAS_OVERRIDE):
        ov = json.load(open(ALIAS_OVERRIDE))
        names = {str(x).upper(): str(x) for x in m["Name"]}
        applied = 0
        for src, target in ov.items():
            if str(target).upper() in names:
                alias[str(src).upper()] = names[str(target).upper()]
                applied += 1
        print(f"alias overrides applied: {applied}/{len(ov)} from {ALIAS_OVERRIDE}",
              flush=True)
    else:
        print(f"note: {ALIAS_OVERRIDE} not present -- using master-derived aliases only",
              flush=True)
    return alias


def attach(df, alias, name_col="rbp"):
    df = df.copy()
    df["_master"] = df[name_col].astype(str).str.upper().map(alias)
    return df


def agg_counts(df, classes):
    """Per-master-RBP summed region counts (may be fractional under even-split)."""
    ncols = [f"n_{c}" for c in classes if f"n_{c}" in df.columns]
    return df.groupby("_master")[ncols].sum()


def frac_block(df, prefix, classes, count_col):
    """Sum counts per master RBP, recompute fractions so they sum to 1."""
    g = agg_counts(df, classes)
    tot = g.sum(axis=1)
    out = pd.DataFrame(index=g.index)
    out[f"{prefix}_{count_col}"] = tot
    for c in classes:
        col = f"n_{c}"
        if col in g.columns:
            out[f"{prefix}_frac_{c}"] = np.where(tot > 0, g[col] / tot.replace(0, np.nan), np.nan)
    fr = [f"{prefix}_frac_{c}" for c in classes if f"{prefix}_frac_{c}" in out.columns]
    if fr:
        out[f"{prefix}_dominant"] = out[fr].idxmax(axis=1).str.replace(f"{prefix}_frac_", "", regex=False)
    return out


def union_block(raw, prefix, primary, secondary, classes):
    """ENCODE-priority union: primary's profile wins; secondary fills gaps.

    `raw` maps family prefix -> per-master summed counts (from agg_counts).
    Emits a `<prefix>_source` column so every row's provenance is explicit
    rather than implied by the construction.
    """
    if primary not in raw and secondary not in raw:
        return None
    p = raw.get(primary)
    s = raw.get(secondary)

    if p is None:
        combined, src = s.copy(), pd.Series(secondary.split("_")[0], index=s.index)
    elif s is None:
        combined, src = p.copy(), pd.Series(primary.split("_")[0], index=p.index)
    else:
        only_s = s.loc[[i for i in s.index if i not in p.index]]
        combined = pd.concat([p, only_s])
        src = pd.concat([
            pd.Series(primary.split("_")[0], index=p.index),
            pd.Series(secondary.split("_")[0], index=only_s.index),
        ])

    tot = combined.sum(axis=1)
    out = pd.DataFrame(index=combined.index)
    out[f"{prefix}_n_regions"] = tot
    for c in classes:
        col = f"n_{c}"
        if col in combined.columns:
            out[f"{prefix}_frac_{c}"] = np.where(
                tot > 0, combined[col] / tot.replace(0, np.nan), np.nan)
    fr = [f"{prefix}_frac_{c}" for c in classes if f"{prefix}_frac_{c}" in out.columns]
    out[f"{prefix}_dominant"] = out[fr].idxmax(axis=1).str.replace(
        f"{prefix}_frac_", "", regex=False)
    out[f"{prefix}_source"] = src
    return out


def main():
    m = pd.read_csv(MASTER, low_memory=False)
    print(f"master: {m.shape[0]} rows x {m.shape[1]} cols", flush=True)
    alias = build_alias_map(m)
    print(f"alias keys: {len(alias):,}\n", flush=True)

    blocks = []
    unmatched = []

    # ---- the four criterion families + POSTAR3 ----------------------------
    raw = {}
    for prefix, path, count_col in SOURCE_FILES:
        if not os.path.exists(path):
            print(f"  SKIP {prefix:17s} missing {path}", flush=True)
            continue
        d = attach(pd.read_csv(path), alias)
        unmatched += [{"source": prefix, "name": n} for n in
                      sorted(set(d[d._master.isna()].rbp))]
        ok = d[d._master.notna()]
        raw[prefix] = agg_counts(ok, CLASSES)
        blk = frac_block(ok, prefix, CLASSES, count_col)
        for extra in ("median_support", "median_dataset_support",
                      "median_experiment_support"):
            if extra in ok.columns:
                blk[f"{prefix}_{extra}"] = ok.groupby("_master")[extra].median()
                break
        blocks.append(blk)
        print(f"  {prefix:17s} -> {blk.shape[0]:>4} master RBPs   "
              f"({len(set(d[d._master.isna()].rbp))} names unmatched)", flush=True)

    # ---- union families (ENCODE-priority) ---------------------------------
    for prefix, primary, secondary in UNION_FAMILIES:
        blk = union_block(raw, prefix, primary, secondary, CLASSES)
        if blk is None:
            print(f"  SKIP {prefix:17s} (neither source available)", flush=True)
            continue
        blocks.append(blk)
        vc = blk[f"{prefix}_source"].value_counts()
        print(f"  {prefix:17s} -> {blk.shape[0]:>4} master RBPs   "
              f"(from {dict(vc)})", flush=True)

    # ---- Skipper (published) ----------------------------------------------
    if os.path.exists("skipper_region_profile.csv"):
        s = attach(pd.read_csv("skipper_region_profile.csv"), alias)
        unmatched += [{"source": "SKIPPER", "name": n} for n in
                      sorted(set(s[s._master.isna()].rbp))]
        ok = s[s._master.notna()]
        ncols = [c for c in ok.columns if c.startswith("skipper_n_") and c != "skipper_n_windows"]
        g = ok.groupby("_master")[ncols].sum()
        tot = g.sum(axis=1)
        blk = pd.DataFrame(index=g.index)
        blk["skipper_n_windows"] = tot
        for c in ncols:
            cls = c.replace("skipper_n_", "")
            blk[f"skipper_frac_{cls}"] = np.where(tot > 0, g[c] / tot.replace(0, np.nan), np.nan)
        # enrichment: mean across that RBP's datasets (already background-normalised)
        ecols = [c for c in ok.columns if c.startswith("skipper_enrich_")]
        if ecols:
            blk = blk.join(ok.groupby("_master")[ecols].mean())
        blk["skipper_n_celllines"] = ok.groupby("_master").biosample.nunique()
        fr = [c for c in blk.columns if c.startswith("skipper_frac_")]
        blk["skipper_dominant"] = blk[fr].idxmax(axis=1).str.replace("skipper_frac_", "", regex=False)
        blocks.append(blk)
        print(f"SKIPPER -> {blk.shape[0]} master RBPs", flush=True)

    # ---- assemble ----------------------------------------------------------
    out = m.set_index("Name")
    for b in blocks:
        out = out.join(b, how="left")
    out = out.reset_index()

    fams = {p: f"{p}_{c}" for p, _, c in SOURCE_FILES}
    fams["skipper"] = "skipper_n_windows"
    for fam, col in fams.items():
        if col in out.columns:
            out[f"has_{fam}"] = out[col].notna()

    have = [f"has_{f}" for f in fams if f"has_{f}" in out.columns]
    out["n_eclip_sources"] = out[have].sum(axis=1) if have else 0

    # union families are DERIVED from encode_*/encori_*, so they get their own
    # has_ flags but are deliberately excluded from n_eclip_sources -- counting
    # them would double-count the same underlying evidence.
    for prefix, _, _ in UNION_FAMILIES:
        col = f"{prefix}_n_regions"
        if col in out.columns:
            out[f"has_{prefix}"] = out[col].notna()

    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}  {out.shape}", flush=True)

    print("\n=== master-table coverage ===", flush=True)
    for f in fams:
        c = f"has_{f}"
        if c in out.columns:
            print(f"  {f:8s} {int(out[c].sum()):>5} / {len(out)} master RBPs", flush=True)
    print(f"  {'ANY':8s} {int((out.n_eclip_sources > 0).sum()):>5} / {len(out)}", flush=True)
    print("\n  RBPs by number of sources:", flush=True)
    print(out.n_eclip_sources.value_counts().sort_index().to_string(), flush=True)

    if unmatched:
        u = pd.DataFrame(unmatched).drop_duplicates()
        u.to_csv(UNMATCHED_REPORT, index=False)
        print(f"\nwrote {UNMATCHED_REPORT}  "
              f"({u.name.nunique()} distinct names not in master)", flush=True)


if __name__ == "__main__":
    main()


