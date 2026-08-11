"""Append compact RCSB PDB columns to the already extended 242-column table.

The element-level secondary structures remain normalized in
``Human_Proteome_RCSB_Secondary_Structure.csv.gz`` and join on ``uniprot_id``.
Only compact per-protein lists and counts are broadcast to every isoform row.

    python append_rcsb_to_extended.py
"""

import argparse
import csv
import os
import time

import paths
import rcsb

csv.field_size_limit(1 << 30)


def append(input_path, summary_path, output_path, progress_every=5000):
    by_accession = rcsb.load_summary(summary_path)
    empty = rcsb.empty_summary()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    part = output_path + ".part"
    written = matched = 0
    started = time.time()

    with open(input_path, newline="", encoding="utf-8", errors="replace") as src, \
            open(part, "w", newline="", encoding="utf-8") as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst)
        header = next(reader)
        duplicates = set(header).intersection(rcsb.SUMMARY_COLUMNS)
        if duplicates:
            raise ValueError(f"input already contains RCSB columns: {sorted(duplicates)}")
        try:
            accession_index = header.index("uniprot_id")
        except ValueError as exc:
            raise ValueError(f"{input_path} has no uniprot_id column") from exc
        writer.writerow(header + rcsb.SUMMARY_COLUMNS)

        for row in reader:
            if len(row) != len(header):
                raise ValueError(
                    f"row {written + 2} has {len(row)} fields; expected {len(header)}"
                )
            accession = row[accession_index]
            values = by_accession.get(accession)
            if values is None:
                values = empty
            else:
                matched += 1
            writer.writerow(row + [values[column] for column in rcsb.SUMMARY_COLUMNS])
            written += 1
            if progress_every and written % progress_every == 0:
                rate = written / max(time.time() - started, 0.001)
                print(f"  {written:,} rows ({rate:,.0f}/s)", flush=True)

    os.replace(part, output_path)
    return written, matched, len(header) + len(rcsb.SUMMARY_COLUMNS)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", default=paths.EXTENDED_TABLE)
    p.add_argument("--summary", default=paths.RCSB_SUMMARY)
    p.add_argument("--output", default=paths.RCSB_ENRICHED_TABLE)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    rows, matched, columns = append(args.input, args.summary, args.output)
    print(f"wrote {os.path.abspath(args.output)}")
    print(f"  {rows:,} rows x {columns:,} columns")
    print(f"  {matched:,} rows matched a reviewed human UniProt accession")


if __name__ == "__main__":
    main()

