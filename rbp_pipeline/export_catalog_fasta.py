"""Export canonical, NCBI-isoform, or all clean-catalog sequences as FASTA."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os

import catalog_selection
from catalog_io import read_rows


ROW_KIND = {
    "canonical": "swissprot_canonical",
    "isoform": "ncbi_isoform",
}


def _identifier(row, style):
    if style == "protein-key":
        return row["protein_key"]
    if style == "uniprot":
        return row.get("uniprot_id")
    if style == "refseq":
        values = row.get("refseq_protein_ids") or []
        return values[0] if values else None
    if row.get("row_kind") == "swissprot_canonical" and row.get("uniprot_id"):
        return row["uniprot_id"]
    values = row.get("refseq_protein_ids") or []
    return values[0] if values else row["protein_key"]


def export(catalog, output, row_kind="all", identifier="auto", wrap=80):
    rows = read_rows(catalog)
    if row_kind != "all":
        rows = [row for row in rows if row.get("row_kind") == ROW_KIND[row_kind]]
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    seen = set()
    count = 0
    with open(output, "w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            accession = _identifier(row, identifier)
            if not accession:
                raise ValueError(
                    f"{row['protein_key']} has no identifier for --identifier {identifier}"
                )
            if accession in seen:
                raise ValueError(f"duplicate FASTA identifier: {accession}")
            seen.add(accession)
            sequence = str(row["sequence"]).upper().rstrip("*")
            handle.write(f">{accession} protein_key={row['protein_key']}\n")
            for start in range(0, len(sequence), wrap):
                handle.write(sequence[start:start + wrap] + "\n")
            count += 1
    manifest = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "catalog": os.path.abspath(catalog),
        "catalog_sha256": catalog_selection.file_sha256(catalog),
        "output": os.path.abspath(output),
        "output_sha256": catalog_selection.file_sha256(output),
        "row_kind": row_kind,
        "identifier": identifier,
        "sequence_count": count,
    }
    with open(output + ".manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--row-kind", choices=["canonical", "isoform", "all"], default="all"
    )
    parser.add_argument(
        "--identifier",
        choices=["auto", "protein-key", "uniprot", "refseq"],
        default="auto",
    )
    parser.add_argument("--wrap", type=int, default=80)
    args = parser.parse_args(argv)
    if args.wrap < 1:
        parser.error("--wrap must be positive")
    manifest = export(
        args.catalog, args.output, args.row_kind, args.identifier, args.wrap
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
