"""Produce the expanded master table.

Two jobs, in one streaming pass so that a 7.4 GB file is never held in memory:

  1. copy every existing row across, updating its two
     ``PPI_*_in_Dataframe`` columns to account for the newly-added proteins;
  2. append the new rows.

Why existing rows have to be touched at all
-------------------------------------------
Three STRING columns are defined relative to the set of proteins present in the
table -- ``PPI_UniProt_Partners`` as well as the two ``_in_Dataframe`` ones,
since translating a partner ENSP to UniProt requires that partner to be in the
table. Adding proteins makes previously-ineligible partners count, so all three
go stale on existing rows unless recomputed. The set only grows, so the update
is purely additive: a row's dicts gain entries and never lose any.

Everything else on an existing row is copied through byte-for-byte.
"""

import ast
import csv
import os
import pickle
import sys

csv.field_size_limit(1 << 30)

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import schema
import string_ppi
from build_rows import format_value


class Unparseable(Exception):
    """A stored cell could not be read back, so it must not be rewritten."""


def parse_dict(text):
    """Read a stored dict cell.

    Empty means an empty dict. Anything present but unreadable raises, rather
    than degrading to ``{}`` -- these cells hold up to ~1,400 interaction
    partners each, and quietly substituting an empty dict would delete them all
    while looking like a successful run. The caller leaves such rows untouched.

    (Every one of the 20,000 rows sampled parses cleanly today. This guard is
    for the day a source refresh writes numpy reprs into these columns, as it
    already has for the CD-CODE ones.)
    """
    if text is None or not text.strip():
        return {}
    try:
        value = ast.literal_eval(text)
    except (ValueError, SyntaxError) as exc:
        raise Unparseable(text[:120]) from exc
    if not isinstance(value, dict):
        raise Unparseable(f"expected a dict, got {type(value).__name__}")
    return value


def run(master_path, new_rows_path, out_path, reverse_path, string_path,
        progress_every=5000):
    with open(reverse_path, "rb") as fh:
        reverse = pickle.load(fh)
    with open(string_path, "rb") as fh:
        string_data = pickle.load(fh)
    ensp_to_uniprot = string_data["ensp_to_uniprot"]

    schema.validate_against(master_path)
    schema.validate_against(new_rows_path)

    i_ensp = schema.COLUMNS.index("ENSP_clean")
    i_uni = schema.COLUMNS.index("PPI_UniProt_Partners")
    i_ensp_in = schema.COLUMNS.index("PPI_ENSP_Partners_in_Dataframe")
    i_uni_in = schema.COLUMNS.index("PPI_UniProt_Partners_in_Dataframe")

    copied = updated = appended = skipped = 0

    with open(out_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(schema.COLUMNS)

        with open(master_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh)
            next(reader)
            for row in reader:
                gained = reverse.get(row[i_ensp])
                if gained:
                    try:
                        uni, ensp_in, uni_in = string_ppi.extend_table_relative(
                            parse_dict(row[i_uni]),
                            parse_dict(row[i_ensp_in]),
                            parse_dict(row[i_uni_in]),
                            gained,
                            ensp_to_uniprot,
                        )
                    except Unparseable as exc:
                        # Copy the row through untouched rather than replace
                        # real partner data with a partial dict.
                        skipped += 1
                        if skipped <= 5:
                            print(f"  !! row {copied}: {exc}; left unmodified",
                                  flush=True)
                    else:
                        row[i_uni] = format_value(uni)
                        row[i_ensp_in] = format_value(ensp_in)
                        row[i_uni_in] = format_value(uni_in)
                        updated += 1
                writer.writerow(row)
                copied += 1
                if progress_every and copied % progress_every == 0:
                    print(f"  copied {copied} existing rows ({updated} PPI-updated)",
                          flush=True)

        with open(new_rows_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh)
            next(reader)
            for row in reader:
                writer.writerow(row)
                appended += 1

    return copied, updated, appended, skipped


if __name__ == "__main__":
    import paths

    scratch = sys.argv[1] if len(sys.argv) > 1 else paths.SCRATCH
    master = paths.MASTER_TABLE
    new_rows = paths.NEW_ROWS
    out = paths.EXPANDED_TABLE

    copied, updated, appended, skipped = run(
        master, new_rows, out,
        os.path.join(scratch, "string_reverse.pkl"),
        os.path.join(scratch, "string_new.pkl"),
    )
    print(f"\nwrote {out}")
    print(f"  {copied} existing rows copied, of which {updated} had their "
          f"PPI_*_in_Dataframe columns extended")
    if skipped:
        print(f"  {skipped} rows left unmodified because their stored PPI "
              f"columns could not be parsed")
    print(f"  {appended} new rows appended")
    print(f"  {copied + appended} rows total")
