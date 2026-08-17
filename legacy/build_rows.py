"""Assemble the 157-column rows for the proteins missing from the master table.

Reads the per-family intermediates produced by the ``run_*.py`` scripts and
emits ``Missing_Proteins_Rows.csv``: same header, same column order and same
value encoding as ``Isoform_Post_Merge_PSLab_OpenTargets.csv``, so the two can
be concatenated directly.

Value encoding
--------------
The master table was written by ``DataFrame.to_csv``, so every non-scalar cell
is a Python ``repr``: lists as ``['a', 'b']``, dicts as ``{'RRM': 3}``, and
missing values as an empty field. `format_value` reproduces that, with one
deliberate improvement -- numpy scalars are cast to Python types first, so the
CD-CODE columns come out as ``[9606]`` rather than the ``[np.int64(9606)]``
that ``ast.literal_eval`` chokes on in the existing rows.
"""

import csv
import json
import os
import pickle
import sys

csv.field_size_limit(1 << 30)

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import schema
import cider
import cdcode
import identity
import opentargets as ot
import paths
import string_ppi

KAPPEL = paths.KAPPEL
MASTER = paths.MASTER_TABLE


def format_value(value):
    """Render one cell the way pandas' to_csv would have.

    Scalars print bare; containers print as Python reprs; None and NaN become an
    empty field.
    """
    if value is None:
        return ""
    if isinstance(value, float) and value != value:  # NaN
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple, dict)):
        return repr(_plain(value))
    return str(_plain(value))


def _plain(value):
    """Recursively replace numpy scalars with Python equivalents."""
    if isinstance(value, dict):
        return {_plain(k): _plain(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_plain(v) for v in value]
    if isinstance(value, tuple):
        return tuple(_plain(v) for v in value)
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            return value.item()
        except (ValueError, AttributeError):
            return value
    return value


def load_intermediates(scratch):
    """Read everything the run_* scripts produced."""
    with open(os.path.join(scratch, "missing_proteins.json")) as fh:
        proteins = json.load(fh)
    with open(os.path.join(scratch, "missing_idmap.json")) as fh:
        idmap = json.load(fh)
    with open(os.path.join(scratch, "missing_ensg.json")) as fh:
        ids = json.load(fh)
    with open(os.path.join(scratch, "missing_idr.json")) as fh:
        idrs = json.load(fh)
    with open(os.path.join(scratch, "missing_domains.json")) as fh:
        domains = json.load(fh)
    with open(os.path.join(scratch, "missing_go.json")) as fh:
        go_terms = json.load(fh)
    with open(os.path.join(scratch, "missing_pslab.json")) as fh:
        pslab_values = json.load(fh)
    with open(os.path.join(scratch, "missing_opentargets.pkl"), "rb") as fh:
        open_targets = pickle.load(fh)
    with open(os.path.join(scratch, "string_new.pkl"), "rb") as fh:
        string_data = pickle.load(fh)
    return dict(
        proteins=proteins, idmap=idmap, ids=ids, idrs=idrs, domains=domains,
        go_terms=go_terms, pslab_values=pslab_values, open_targets=open_targets,
        string_data=string_data,
    )


def build(scratch, out_path):
    data = load_intermediates(scratch)
    proteins = data["proteins"]
    accessions = proteins["missing"]
    seqs = proteins["seqs"]
    headers = proteins["headers"]

    census = identity.load_rbp_census(os.path.join(KAPPEL, "RBPs.xlsx"))

    # CD-CODE condensate identity rests on file ordering, so check it before
    # using it rather than trusting it. A mismatch here means every condensate
    # assignment past the first divergence is wrong, which is not something to
    # discover after the run.
    misaligned = cdcode.verify_alignment(KAPPEL)
    if misaligned:
        raise SystemExit(
            f"CD-CODE member tables are misaligned ({len(misaligned)} problems, "
            f"first: {misaligned[0]}). Refusing to assign condensates. See "
            f"cdcode.verify_alignment()."
        )
    protein_to_condensates, condensate_attributes = cdcode.build_index(KAPPEL)

    otd = data["open_targets"]
    target_columns = otd["target_columns"]
    sd = data["string_data"]

    stats = {
        "with_idr": 0, "with_domain": 0, "with_go": 0, "with_condensate": 0,
        "with_ppi": 0, "with_opentargets": 0, "with_ensg": 0,
    }

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(schema.COLUMNS)

        for accession in accessions:
            sequence = seqs[accession]
            row = {}

            # --- identity + bookkeeping ---------------------------------
            row.update(identity.build_row(
                accession, sequence, headers[accession],
                data["idmap"].get(accession, {}),
                data["ids"]["ensg"].get(accession),
                census,
            ))
            if row["ID_list"]:
                stats["with_ensg"] += 1

            # --- CIDER on the whole sequence ----------------------------
            row.update(cider.whole_sequence(sequence))

            # --- IDR / FOLD geometry ------------------------------------
            geometry = data["idrs"][accession]
            row.update(geometry)
            idr_seqs = geometry["IDR_discrete_seq"]
            if idr_seqs:
                stats["with_idr"] += 1

            # --- CIDER per IDR ------------------------------------------
            per_idr = cider.per_idr(idr_seqs)
            for metric, values in per_idr.items():
                row[f"IDR_{metric}"] = values

            # --- Domains + CIDER per domain -----------------------------
            domain_row = data["domains"].get(accession, dict(empty_domains()))
            row.update(domain_row)
            if domain_row.get("Domains"):
                stats["with_domain"] += 1
            per_domain = cider.per_domain(domain_row.get("Domains_discrete_seq") or {})
            for metric, values in per_domain.items():
                row[f"Domains_{metric}"] = values

            # --- STRING -------------------------------------------------
            ppi = string_ppi.columns_for(
                sd["new_ensp"].get(accession),
                sd["adjacency"],
                sd["ensp_to_uniprot"],
                sd["all_ensps"],
                sd["all_uniprots"],
            )
            row.update(ppi)
            if ppi["PPI_ENSP_Partners"]:
                stats["with_ppi"] += 1

            # --- CD-CODE ------------------------------------------------
            condensates = cdcode.lookup(accession, protein_to_condensates,
                                        condensate_attributes)
            row.update(condensates)
            if condensates["Condensate Name"]:
                stats["with_condensate"] += 1

            # --- Gene Ontology ------------------------------------------
            gterms = data["go_terms"].get(accession, dict(empty_go()))
            row.update(gterms)
            if gterms["C_ids"] or gterms["P_ids"] or gterms["F_ids"]:
                stats["with_go"] += 1

            # --- Open Targets -------------------------------------------
            open_target_columns = ot.columns_for(
                row["ID_list"], otd["associations"], otd["expression"],
                otd["targets"], target_columns,
            )
            row.update(open_target_columns)
            if any(v for v in open_target_columns["approvedSymbol"].values()):
                stats["with_opentargets"] += 1

            # --- PSLab, one entry per IDR -------------------------------
            row.update(collect_pslab(accession, len(idr_seqs), data["pslab_values"]))

            writer.writerow([format_value(row.get(c)) for c in schema.COLUMNS])

    return stats, len(accessions)


def collect_pslab(accession, n_idrs, pslab_values):
    """Gather the 11 PSLab columns as lists ordered to match IDR_discrete_seq."""
    out = {c: [] for c in schema.PSLAB_COLUMNS}
    for i in range(1, n_idrs + 1):
        entry = pslab_values.get(f"{accession}_1_{i}")
        for column in schema.PSLAB_COLUMNS:
            out[column].append(entry.get(column) if entry else None)
    return out


def empty_domains():
    import domains as domains_module
    return domains_module.EMPTY


def empty_go():
    import go as go_module
    return go_module.EMPTY


if __name__ == "__main__":
    scratch = sys.argv[1] if len(sys.argv) > 1 else paths.SCRATCH
    out_path = sys.argv[2] if len(sys.argv) > 2 else paths.NEW_ROWS

    schema.validate_against(MASTER)
    print("schema validated against master header")

    stats, total = build(scratch, out_path)
    print(f"\nwrote {total} rows -> {out_path}")
    for key, value in stats.items():
        print(f"  {key:20s} {value:5d}  ({value/total*100:5.1f}%)")
