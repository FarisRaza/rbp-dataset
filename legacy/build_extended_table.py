"""Build the final table: 157 original columns + 85 new ones, in one pass.

Streams the master, applies everything that needs applying, appends the new
rows, and writes ``..._Extended.csv``. Doing it in a single pass matters: each
pass over these files reads and writes ~15 GB.

What happens to each row:

  * the three table-relative STRING columns are brought up to date, since adding
    proteins makes previously-ineligible partners count (see
    ``string_ppi.extend_table_relative``);
  * 74 eCLIP columns are joined on gene symbol;
  * 7 InterPro columns are joined on the RefSeq accession in ``ProteinHGVS``;
  * 4 GO functional-role flags are derived from the GO columns already present.

Everything else is copied through unchanged.

    python build_extended_table.py
"""

import csv
import os
import pickle
import sys
import time

csv.field_size_limit(1 << 30)

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import eclip
import go_roles
import interpro
import paths
import schema
import string_ppi
from append_to_master import Unparseable, parse_dict
from build_rows import format_value

OUTPUT = os.path.join(
    paths.KAPPEL, "Isoform_Post_Merge_PSLab_OpenTargets_Extended.csv"
)


def build(master_path, new_rows_path, out_path, scratch, progress_every=5000):
    print("loading annotation indexes")
    by_gene, eclip_columns = eclip.load(kappel_dir=paths.KAPPEL)
    print(f"  eCLIP: {len(by_gene)} genes, {len(eclip_columns)} columns")

    interpro_path = os.path.join(paths.KAPPEL, interpro.DEFAULT_TSV)
    if os.path.exists(interpro_path):
        by_accession = interpro.index_by_accession(interpro_path)
        print(f"  InterPro: {len(by_accession)} proteins")
    else:
        by_accession = {}
        print(f"  InterPro: {interpro_path} absent, columns will be empty")

    with open(os.path.join(scratch, "string_reverse.pkl"), "rb") as fh:
        reverse = pickle.load(fh)
    with open(os.path.join(scratch, "string_new.pkl"), "rb") as fh:
        ensp_to_uniprot = pickle.load(fh)["ensp_to_uniprot"]
    print(f"  STRING: {len(reverse)} existing ENSPs gain partners")

    extension_columns = eclip_columns + interpro.COLUMNS + go_roles.ROLE_COLUMNS
    header = schema.COLUMNS + extension_columns
    print(f"\noutput: {len(header)} columns "
          f"({len(schema.COLUMNS)} original + {len(extension_columns)} new)")

    index = {c: i for i, c in enumerate(schema.COLUMNS)}
    i_ensp = index["ENSP_clean"]
    i_uni = index["PPI_UniProt_Partners"]
    i_ensp_in = index["PPI_ENSP_Partners_in_Dataframe"]
    i_uni_in = index["PPI_UniProt_Partners_in_Dataframe"]
    i_name = index["Name"]
    i_hgvs = index["ProteinHGVS"]

    stats = {"eclip": 0, "interpro": 0, "ppi_updated": 0, "ppi_skipped": 0}
    written = 0
    start = time.time()

    def extend(row):
        """Return the 85 extension values for one row, in header order."""
        values = []

        record = by_gene.get(row[i_name])
        if record is not None:
            stats["eclip"] += 1
        values.extend(
            eclip.columns_for(row[i_name], by_gene, eclip_columns)[c]
            for c in eclip_columns
        )

        ip = interpro.columns_for(row[i_hgvs], by_accession)
        if ip["InterPro_n_hits"]:
            stats["interpro"] += 1
        values.extend(format_value(ip[c]) for c in interpro.COLUMNS)

        flags = go_roles.flags_for({
            "P_descriptions": row[index["P_descriptions"]],
            "F_descriptions": row[index["F_descriptions"]],
        })
        values.extend(str(flags[c]) for c in go_roles.ROLE_COLUMNS)
        return values

    with open(out_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(header)

        with open(master_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh)
            next(reader)
            for row in reader:
                gained = reverse.get(row[i_ensp])
                if gained:
                    try:
                        uni, ensp_in, uni_in = string_ppi.extend_table_relative(
                            parse_dict(row[i_uni]), parse_dict(row[i_ensp_in]),
                            parse_dict(row[i_uni_in]), gained, ensp_to_uniprot,
                        )
                    except Unparseable:
                        stats["ppi_skipped"] += 1
                    else:
                        row[i_uni] = format_value(uni)
                        row[i_ensp_in] = format_value(ensp_in)
                        row[i_uni_in] = format_value(uni_in)
                        stats["ppi_updated"] += 1
                writer.writerow(row + extend(row))
                written += 1
                if progress_every and written % progress_every == 0:
                    rate = written / (time.time() - start)
                    print(f"  {written} rows  ({rate:.0f}/s)", flush=True)

        with open(new_rows_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh)
            next(reader)
            for row in reader:
                writer.writerow(row + extend(row))
                written += 1

    return written, stats, len(header)


if __name__ == "__main__":
    scratch = sys.argv[1] if len(sys.argv) > 1 else paths.SCRATCH
    schema.validate_against(paths.MASTER_TABLE)
    schema.validate_against(paths.NEW_ROWS)

    written, stats, ncols = build(
        paths.MASTER_TABLE, paths.NEW_ROWS, OUTPUT, scratch
    )
    print(f"\nwrote {OUTPUT}")
    print(f"  {written} rows x {ncols} columns")
    print(f"  {stats['eclip']} rows got eCLIP data")
    print(f"  {stats['interpro']} rows got InterPro domains")
    print(f"  {stats['ppi_updated']} rows had their STRING columns extended")
    if stats["ppi_skipped"]:
        print(f"  {stats['ppi_skipped']} rows left unmodified (unparseable PPI cells)")
