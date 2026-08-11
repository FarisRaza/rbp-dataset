"""Re-download ENCODE eCLIP peaks with correct file selection.

Replaces download_encode_peaks.py, which had two defects:

  1. It filtered files on `preferred_default is True`. ENCODE sets that flag on the
     SINGLE-REPLICATE rep2 file, and sets it once per assembly. The `break` therefore
     picked hg19 or GRCh38 depending on list order (225 hg19 / 27 GRCh38 in practice)
     AND discarded replicate 1 entirely.
  2. It recorded neither assembly nor biosample, so the resulting peak set could not
     be audited.

This version selects on (assembly == GRCh38, output_type == peaks, file_format == bed)
and keeps every such file, classifying each by its biological_replicates field:

     len(biological_replicates) == 1  -> per-replicate peaks  (for log2FC/p cutoffs)
     len(biological_replicates) == 2  -> IDR reproducible set (name field ends _IDR)

Both are needed: IDR is the most rigorous reproducibility filter but is extremely
sparse for eCLIP (AARS: 49 peaks vs 62k-70k per replicate), so per-replicate peaks
under the Van Nostrand et al. 2020 cutoff (log2FC >= 3 and -log10(p) >= 3) are the
practical basis for region-distribution analysis.

Writes a full manifest so every downstream number is traceable to a file accession.
Safe to re-run: existing files are skipped.
"""
import csv
import os
import sys
import time

import requests

OUTDIR = "encode_peaks_grch38"
MANIFEST = "encode_peaks_grch38_manifest.csv"
ASSEMBLY = "GRCh38"

os.makedirs(OUTDIR, exist_ok=True)

session = requests.Session()
session.headers.update({"accept": "application/json"})


def get_json(url, params=None, tries=4):
    for attempt in range(tries):
        try:
            r = session.get(url, params=params, timeout=90)
            r.raise_for_status()
            return r.json()
        except Exception as exc:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))


# ---- Step 1: enumerate ENCORE eCLIP experiments -----------------------------
search = get_json(
    "https://www.encodeproject.org/search/",
    params={
        "type": "Experiment",
        "assay_title": "eCLIP",
        "assembly": ASSEMBLY,
        "internal_tags": "ENCORE",
        "format": "json",
        "limit": "all",
    },
)
experiments = search["@graph"]
print(f"Found {len(experiments)} eCLIP experiments", flush=True)

rows = []
downloaded = skipped = failed = 0

for i, exp in enumerate(experiments, 1):
    exp_acc = exp["@id"].strip("/").split("/")[-1]
    target = exp.get("target") or {}
    rbp = target.get("label", "unknown") if isinstance(target, dict) else "unknown"

    try:
        exp_data = get_json(f"https://www.encodeproject.org{exp['@id']}?format=json")
    except Exception as exc:
        print(f"  !! {exp_acc} metadata failed: {exc}", flush=True)
        failed += 1
        continue

    biosample = (exp_data.get("biosample_ontology") or {}).get("term_name", "")
    controls = [c.get("accession", "") for c in exp_data.get("possible_controls", [])]

    # ---- Step 2: keep EVERY GRCh38 peaks bed, classified by replicate set ----
    for f in exp_data.get("files", []):
        if f.get("output_type") != "peaks":
            continue
        if f.get("file_format") != "bed":
            continue
        if f.get("assembly") != ASSEMBLY:
            continue

        breps = f.get("biological_replicates") or []
        kind = "idr" if len(breps) > 1 else "replicate"
        acc = f["accession"]

        fname = os.path.join(OUTDIR, f"{rbp}_{biosample}_{kind}_{acc}.bed.gz")
        rows.append(
            {
                "file": os.path.basename(fname),
                "rbp": rbp,
                "biosample": biosample,
                "kind": kind,
                "accession": acc,
                "experiment": exp_acc,
                "assembly": f.get("assembly", ""),
                "biological_replicates": ";".join(map(str, breps)),
                "technical_replicates": ";".join(f.get("technical_replicates") or []),
                "preferred_default": f.get("preferred_default", ""),
                "controls": ";".join(controls),
                "file_size": f.get("file_size", ""),
            }
        )

        if os.path.exists(fname) and os.path.getsize(fname) > 0:
            skipped += 1
            continue

        url = "https://www.encodeproject.org" + f["href"]
        try:
            resp = session.get(url, timeout=300)
            resp.raise_for_status()
            with open(fname, "wb") as out:
                out.write(resp.content)
            downloaded += 1
        except Exception as exc:
            print(f"  !! {acc} download failed: {exc}", flush=True)
            failed += 1

    if i % 20 == 0 or i == len(experiments):
        print(
            f"  {i}/{len(experiments)} experiments | "
            f"{downloaded} new, {skipped} cached, {failed} failed",
            flush=True,
        )

# ---- Step 3: manifest --------------------------------------------------------
with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

from collections import Counter

print("\n=== summary ===", flush=True)
print(f"  peak files recorded : {len(rows)}", flush=True)
print(f"  distinct RBPs       : {len({r['rbp'] for r in rows})}", flush=True)
print(f"  by kind             : {dict(Counter(r['kind'] for r in rows))}", flush=True)
print(f"  by biosample        : {dict(Counter(r['biosample'] for r in rows))}", flush=True)
print(f"  by assembly         : {dict(Counter(r['assembly'] for r in rows))}", flush=True)
print(f"\n  downloaded={downloaded} skipped={skipped} failed={failed}", flush=True)
print(f"  manifest -> {MANIFEST}", flush=True)


