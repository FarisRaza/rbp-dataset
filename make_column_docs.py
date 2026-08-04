"""Emit a documentation row for every column this pipeline adds beyond the 157.

Writes ``New_Columns_Documentation.csv`` with the same shape as the project's
documentation spreadsheet: Family, Column Name, Information, Source, Granularity,
Example. Regenerate it after changing any column set so the docs cannot drift
from the code.

    python make_column_docs.py
"""

import csv
import os

import eclip
import go_roles
import interpro
import paths

REGION_LABEL = {
    "utr3": "3' UTR",
    "utr5": "5' UTR",
    "cds": "coding sequence",
    "ncrna_exon": "non-coding exon",
    "intron": "intron",
    "intergenic": "intergenic",
    "other": "other region",
    "repeat": "repeat element",
    "small_rna": "small RNA",
    "splice_site": "splice site",
}

SOURCE_LABEL = {
    "encode_published": "ENCODE eCLIP, published significance criterion "
                        "(-log10 p >= 3 and log2FC >= 3; Van Nostrand 2020)",
    "encode_matched": "ENCODE eCLIP, matched criterion (p <= 1e-3 only, so it "
                      "is directly comparable to ENCORI)",
    "encori_published": "ENCORI/starBase, published significance criterion",
    "encori_matched": "ENCORI, matched criterion (p <= 1e-3), comparable to ENCODE",
    "postar": "POSTAR3 aggregated CLIP sites",
    "skipper": "Skipper reprocessing of ENCODE eCLIP",
}

COUNT_WARNING = ("PROVENANCE ONLY -- not comparable across sources: it tracks "
                 "how many datasets exist for that RBP and how deeply it was "
                 "sequenced, not how much it binds.")


def eclip_rows(columns):
    rows = []
    for column in columns:
        info = None

        if column == "n_eclip_sources":
            info = ("How many of the five CLIP sources have data for this gene "
                    "(0-5). 0 means never profiled, not 'does not bind'.")
        elif column.startswith("has_union_"):
            which = column.replace("has_union_", "")
            info = (f"1 if any source has data for this gene under the {which} "
                    f"significance criterion.")
        elif column.startswith("has_"):
            src = column[4:]
            info = (f"1 if {SOURCE_LABEL.get(src, src)} has data for this gene. "
                    f"Empty means the gene was never profiled, which is distinct "
                    f"from a measured zero.")
        else:
            for src, label in SOURCE_LABEL.items():
                if not column.startswith(src + "_"):
                    continue
                rest = column[len(src) + 1:]
                if rest.startswith("frac_"):
                    region = rest[5:]
                    info = (f"Fraction of this RBP's binding sites falling in the "
                            f"{REGION_LABEL.get(region, region)}. Source: {label}. "
                            f"Use these, not the counts, for cross-source analysis.")
                elif rest.startswith("enrich_"):
                    region = rest[7:]
                    info = (f"Enrichment of this RBP's binding in the "
                            f"{REGION_LABEL.get(region, region)} over background. "
                            f"Source: {label}. Enrichment, not raw fraction, is what "
                            f"distinguishes real preference from genomic real estate "
                            f"-- almost every RBP overlaps some 3' UTR by chance.")
                elif rest == "dominant":
                    info = (f"The region holding the largest fraction of this RBP's "
                            f"sites. Source: {label}.")
                elif rest in ("n_regions", "n_sites", "n_windows"):
                    info = f"Number of binding regions. Source: {label}. {COUNT_WARNING}"
                elif rest == "n_celllines":
                    info = f"Number of cell lines this RBP was profiled in. Source: {label}."
                elif rest.startswith("median_"):
                    info = (f"Median number of experiments/datasets supporting a site "
                            f"for this RBP. Source: {label}. A reproducibility proxy.")
                break

        rows.append({
            "Family": "eCLIP",
            "Column Name": column,
            "Information": info or f"eCLIP-derived column from {column}.",
            "Source": "rbp_master_with_eclip.csv (ENCODE / ENCORI / POSTAR3 / Skipper)",
            "Granularity": "Gene (broadcast to every isoform row of that gene)",
        })
    return rows


INTERPRO_INFO = {
    "InterPro_domains": "List of distinct domain names found by InterProScan, "
                        "restricted to domain-level member databases (Pfam, SMART, "
                        "ProSite, Gene3D, SUPERFAMILY, CDD, PIRSF, NCBIfam, Hamap, SFLD).",
    "InterPro_count": "Dict of domain name -> number of occurrences of that domain.",
    "InterPro_range": "Dict of domain name -> list of (start, end) positions, "
                      "0-based half-open, matching the Domains_range convention.",
    "InterPro_accessions": "Dict of domain name -> InterPro entry accessions (IPR…) "
                           "the signature maps to.",
    "InterPro_databases": "Dict of domain name -> which member databases called it. "
                          "A domain called by several is more confident.",
    "InterPro_go_terms": "Dict of domain name -> GO terms implied by the signature "
                         "(from the -goterms flag). Domain-level, so more specific "
                         "than the protein-level C_/P_/F_ columns.",
    "InterPro_n_hits": "Total number of significant domain hits on this protein.",
}

GO_ROLE_INFO = {
    "role_in_transcription": "1 if any biological-process or molecular-function GO "
                             "term for this protein mentions transcription. Broad by "
                             "design -- also catches 'reverse transcription' and "
                             "'transcription factor binding'.",
    "role_in_translation": "1 if any P or F term indicates regulation of translation "
                           "(regulation of translation, translational repressor/"
                           "activator activity, and their initiation/elongation forms).",
    "role_in_mrna_stability": "1 if any P or F term indicates mRNA stabilisation or "
                              "destabilisation (including 3'-UTR-mediated and "
                              "CRD-mediated stabilisation, miRNA-mediated destabilisation).",
    "role_in_translation_stability": "1 if either role_in_translation or "
                                     "role_in_mrna_stability. Retained because the "
                                     "original analysis used this combined definition.",
}


def build():
    rows = []

    _, eclip_columns = eclip.load(kappel_dir=paths.KAPPEL)
    rows.extend(eclip_rows(eclip_columns))

    for column in interpro.COLUMNS:
        rows.append({
            "Family": "InterPro",
            "Column Name": column,
            "Information": INTERPRO_INFO[column],
            "Source": "InterProScan 5.77-108.0 -> df_np_unique.interpro.tsv",
            "Granularity": "Protein (joined via the RefSeq accession in ProteinHGVS)",
        })

    for column in go_roles.ROLE_COLUMNS:
        rows.append({
            "Family": "GO functional roles",
            "Column Name": column,
            "Information": GO_ROLE_INFO[column],
            "Source": "Derived from the existing P_descriptions / F_descriptions columns",
            "Granularity": "Protein",
        })

    return rows


def main():
    rows = build()
    out_path = os.path.join(paths.KAPPEL, "New_Columns_Documentation.csv")
    fields = ["Family", "Column Name", "Information", "Source", "Granularity"]
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} column descriptions -> {out_path}")
    for family in ("eCLIP", "InterPro", "GO functional roles"):
        n = sum(1 for r in rows if r["Family"] == family)
        print(f"  {family:22s} {n:3d} columns")


if __name__ == "__main__":
    main()
