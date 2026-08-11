"""Turn the published Skipper output into per-RBP region profiles.

Source: Figshare DOI 10.6084/m9.figshare.21206009 (Boyle et al., CC BY 4.0),
the Skipper reprocessing of ENCODE eCLIP -- 219 datasets, 148 distinct RBPs,
GRCh38 / GENCODE v38.

Two products:

1. Per-RBP region composition and enrichment, from
   `encode3_eclip_enrichment.reference.tsv`. This is a dense 219 x 103 matrix of
   (dataset, feature) with `clip_fraction` and `global_fraction`. Enrichment is
   clip_fraction / global_fraction against Skipper's own transcriptome-wide
   background -- which is the right background, unlike the genomic-length
   approximation used elsewhere in this project (introns are 93% of genomic
   space but far less of the transcriptome, so length over-corrects them).

2. An RBP -> target gene side table, from the per-dataset
   `reproducible_enriched_windows` TSVs inside the two tars. These carry
   `gene_name` directly, so no GTF intersection is needed.

Skipper's 103-feature vocabulary is finer than the 5-class scheme used for the
other sources -- it splits by transcript biotype AND includes splice-site
proximity classes (SS5_ADJ, SS3_PROX, SSB_*) that have no equivalent elsewhere.
Both the native profile and a collapsed 5-class mapping are emitted so Skipper
can be compared against ENCODE/ENCORI/POSTAR3 without discarding its extra
resolution.

CAVEAT: the Figshare per-window files are NOT blacklist-filtered. Rows
overlapping encode3_eclip_blacklist.bed must be removed to match the published
filtering (up to 17.7% of rows for some datasets). The blacklist ships in the
Skipper repo, not in the Figshare archive -- if absent, this script says so
rather than silently skipping the step.
"""
import gzip
import io
import os
import re
import sys
import tarfile
from collections import defaultdict

import numpy as np
import pandas as pd

SKDIR = "skipper_published"
ENRICH = os.path.join(SKDIR, "encode3_eclip_enrichment.reference.tsv")
TARS = {
    "K562": os.path.join(SKDIR, "ENCODE_reproducible_enriched_windows.K562.tar"),
    "HepG2": os.path.join(SKDIR, "ENCODE_reproducible_enriched_windows.HepG2.tar"),
}
BLACKLIST = os.path.join(SKDIR, "encode3_eclip_blacklist.bed")

# Collapse Skipper's biotype:region vocabulary onto the 5 classes used for the
# other sources. Splice-site classes are kept separate -- folding them into
# "intron" would overstate intronic binding for splicing factors, which is
# exactly the signal of interest.
SMALL_RNA = {"snoRNA", "snRNA", "scaRNA", "tRNA", "Y_RNA", "7SK", "Mt_RNA",
             "rRNA", "misc_RNA", "vault_RNA", "srpRNA", "RNase_P_RNA",
             "RNase_MRP_RNA", "telomerase_RNA"}
REPEAT_HINTS = ("Alu", "L1", "L2", "LTR", "MIR", "ERV", "SVA", "Satellite",
                "Simple_repeat", "Low_complexity", "DNA", "Antisense_")


def collapse(feature):
    # Features with no biotype prefix are RNA-class or repeat-family labels.
    # These matter enormously: together they are only ~0.3% of the transcriptome
    # background but take ~8% of eCLIP binding -- ~28x enrichment, the strongest
    # signal in the dataset. Folding them into a catch-all "other" buries the
    # single most interesting result, so they get their own classes.
    if ":" not in feature:
        if feature in SMALL_RNA or feature.endswith("_RNA"):
            return "small_rna"
        if any(h in feature for h in REPEAT_HINTS):
            return "repeat"
        return "other"
    biotype, region = feature.split(":", 1)
    if region in ("PRIMIRNA", "EXON_SMALL"):
        return "small_rna"
    if region in ("UTR3",):
        return "utr3"
    if region in ("UTR5",):
        return "utr5"
    if region in ("CDS", "CDS_START", "CDS_STOP"):
        return "cds"
    if region.startswith("SS"):
        return "splice_site"
    if region == "INTRON":
        return "intron"
    if region in ("EXON_LNCRNA", "EXON_PSEUDO"):
        return "ncrna_exon"
    if region == "EXON_MRNA":
        return "cds" if biotype == "protein_coding" else "ncrna_exon"
    return "other"


def build_enrichment():
    df = pd.read_csv(ENRICH, sep="\t")
    # id is "<RBP>:<H|K>"
    parts = df["id"].str.rsplit(":", n=1, expand=True)
    df["rbp"] = parts[0]
    df["biosample"] = parts[1].map({"H": "HepG2", "K": "K562"})
    df["class"] = df["feature"].map(collapse)

    print(f"enrichment reference: {len(df):,} rows, "
          f"{df.rbp.nunique()} RBPs, {df.feature.nunique()} features", flush=True)

    # sanity: clip_fraction should sum to 1 per dataset
    s = df.groupby("id").clip_fraction.sum()
    print(f"  clip_fraction sums per dataset: min={s.min():.6f} max={s.max():.6f}", flush=True)

    # ---- native 103-feature wide matrix -----------------------------------
    native = df.pivot_table(index=["rbp", "biosample"], columns="feature",
                            values="clip_fraction", fill_value=0.0)
    native.columns = [f"skipper_frac_{c}" for c in native.columns]
    native.to_csv("skipper_native_feature_fractions.csv")
    print(f"  wrote skipper_native_feature_fractions.csv {native.shape}", flush=True)

    # ---- collapsed 5+1 class composition and enrichment --------------------
    g = df.groupby(["rbp", "biosample", "class"]).agg(
        clip_fraction=("clip_fraction", "sum"),
        global_fraction=("global_fraction", "sum"),
        clip_count=("clip_count", "sum"),
    ).reset_index()
    g["enrichment"] = np.where(g.global_fraction > 0,
                               g.clip_fraction / g.global_fraction, np.nan)

    frac = g.pivot_table(index=["rbp", "biosample"], columns="class",
                         values="clip_fraction", fill_value=0.0)
    enr = g.pivot_table(index=["rbp", "biosample"], columns="class",
                        values="enrichment")
    cnt = g.pivot_table(index=["rbp", "biosample"], columns="class",
                        values="clip_count", fill_value=0)

    frac.columns = [f"skipper_frac_{c}" for c in frac.columns]
    enr.columns = [f"skipper_enrich_{c}" for c in enr.columns]
    cnt.columns = [f"skipper_n_{c}" for c in cnt.columns]

    out = frac.join(enr).join(cnt).reset_index()
    out["skipper_n_windows"] = cnt.sum(axis=1).values
    out.to_csv("skipper_region_profile.csv", index=False)
    print(f"  wrote skipper_region_profile.csv {out.shape}", flush=True)
    return out


def build_targets():
    """RBP -> target gene side table from the per-dataset window files."""
    bl = None
    if os.path.exists(BLACKLIST):
        b = pd.read_csv(BLACKLIST, sep="\t", header=None,
                        names=["chrom", "start", "end", "name", "score", "strand"])
        bl = set(zip(b.chrom, b.start, b.end, b.strand))
        print(f"blacklist: {len(bl):,} intervals", flush=True)
    else:
        print(f"WARNING: {BLACKLIST} not found -- windows are NOT blacklist-filtered.\n"
              f"         Published Skipper results ARE filtered; up to 17.7% of rows in\n"
              f"         some datasets are blacklisted. Get it from the Skipper repo:\n"
              f"         git clone https://github.com/YeoLab/skipper (annotations/)",
              flush=True)

    rows = []
    per_dataset = []
    for biosample, tarpath in TARS.items():
        if not os.path.exists(tarpath):
            print(f"missing {tarpath}", flush=True)
            continue
        with tarfile.open(tarpath) as tf:
            members = [m for m in tf.getmembers() if m.name.endswith(".tsv.gz")]
            print(f"{biosample}: {len(members)} datasets", flush=True)
            for i, m in enumerate(members, 1):
                base = os.path.basename(m.name)
                rbp = base.split(".")[0].rsplit("_", 1)[0]
                fh = tf.extractfile(m)
                if fh is None:
                    continue
                try:
                    with gzip.open(io.BytesIO(fh.read()), "rt") as gz:
                        d = pd.read_csv(gz, sep="\t", low_memory=False)
                except Exception as exc:
                    print(f"  !! {base}: {exc}", flush=True)
                    continue

                n_raw = len(d)
                n_bl = 0
                if bl is not None and {"chrom", "start", "end", "strand"} <= set(d.columns):
                    keep = ~pd.Series(
                        list(zip(d.chrom, d.start, d.end, d.strand)), index=d.index
                    ).isin(bl)
                    n_bl = int((~keep).sum())
                    d = d[keep]

                per_dataset.append({"rbp": rbp, "biosample": biosample,
                                    "n_windows_raw": n_raw,
                                    "n_windows_blacklisted": n_bl,
                                    "n_windows": len(d)})

                if "gene_name" in d.columns:
                    gcol = "gene_name"
                elif "gene" in d.columns:
                    gcol = "gene"
                else:
                    continue
                agg = d.groupby(gcol).size().reset_index(name="n_windows")
                agg["rbp"] = rbp
                agg["biosample"] = biosample
                agg = agg.rename(columns={gcol: "target_gene"})
                rows.append(agg)

                if i % 25 == 0 or i == len(members):
                    print(f"  {i}/{len(members)}", flush=True)

    if rows:
        targets = pd.concat(rows, ignore_index=True)
        targets = targets[["rbp", "biosample", "target_gene", "n_windows"]]
        targets.to_csv("skipper_rbp_targets.csv.gz", index=False, compression="gzip")
        print(f"wrote skipper_rbp_targets.csv.gz  "
              f"({len(targets):,} RBP-target pairs, "
              f"{targets.rbp.nunique()} RBPs, {targets.target_gene.nunique():,} genes)",
              flush=True)

    if per_dataset:
        pd_df = pd.DataFrame(per_dataset)
        pd_df.to_csv("skipper_dataset_summary.csv", index=False)
        print("\n=== window counts per dataset ===", flush=True)
        print(pd_df.n_windows.describe()[["count", "mean", "50%", "min", "max"]].to_string(),
              flush=True)
        if bl is not None:
            print(f"blacklisted rows removed: {pd_df.n_windows_blacklisted.sum():,} "
                  f"({100*pd_df.n_windows_blacklisted.sum()/pd_df.n_windows_raw.sum():.1f}%)",
                  flush=True)


if __name__ == "__main__":
    build_enrichment()
    print()
    build_targets()


