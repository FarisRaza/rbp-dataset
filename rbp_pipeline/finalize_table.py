"""Write the sorted, historical-layout-compatible final protein table.

The clean catalog and feature sidecars deliberately retain typed list columns
and explicit provenance.  This final step adds the familiar identity aliases
from the historical table, restores its feature-family block order, and keeps
all clean columns after those blocks so no information is discarded.

``dominant_isoform`` is computed from amino-acid sequence equality against the
canonical Swiss-Prot sequence of a row's mapped UniProt parent(s). It is not
inferred from an NCBI isoform label. The historical table used the capitalization
``Dominant_Isoform``; the lowercase spelling is used here because Parquet/DuckDB
column lookup is case-insensitive and cannot safely carry both spellings.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os

import go_roles
import interpro
import opentargets
import ptm
import rcsb
import schema


IDENTITY_COLUMNS = [
    "uniprot_id",
    "dominant_isoform",
    "sequence",
    "UNIQUE",
    "ProteinHGVS",
    "HGVSDescription",
    "ENSG",
    "ENST",
    "ENSP",
    "ncbi_gene_id",
    "protein_key",
    "ID",
    "Name",
    "Description",
]

STRING_COMPATIBILITY_COLUMNS = list(schema.STRING_COLUMNS)


def _quote_identifier(value):
    return '"' + str(value).replace('"', '""') + '"'


def _sql_string(value):
    return "'" + os.path.abspath(value).replace("'", "''") + "'"


def _sha256(path, chunk=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _output_columns(input_columns):
    """Return the stable family-block order followed by clean/provenance data."""
    ordered = list(IDENTITY_COLUMNS)
    preferred = (
        schema.CIDER_COLUMNS
        + schema.IDR_GEOMETRY_COLUMNS
        + schema.IDR_CIDER_COLUMNS
        + schema.DOMAIN_GEOMETRY_COLUMNS
        + schema.DOMAIN_CIDER_COLUMNS
        + STRING_COMPATIBILITY_COLUMNS
        + schema.CDCODE_COLUMNS
        + schema.GO_COLUMNS
        + ["isoform_number", "ID_list"]
        + opentargets.COLUMNS
        + schema.PSLAB_COLUMNS
    )
    # eCLIP columns vary with the supplied compilation, so preserve their
    # source-table order by taking them from the assembled input schema.
    eclip_columns = [
        column for column in input_columns
        if column.startswith(("encode_", "encori_", "postar_", "skipper_", "has_"))
        or column == "n_eclip_sources"
    ]
    preferred += eclip_columns
    preferred += interpro.COLUMNS + go_roles.ROLE_COLUMNS
    preferred += ptm.PTM_COLUMNS + ptm.PROVENANCE_COLUMNS
    preferred += rcsb.SUMMARY_COLUMNS

    available = set(input_columns)
    generated = set(IDENTITY_COLUMNS + STRING_COMPATIBILITY_COLUMNS + [
        "isoform_number", "ID_list"
    ])
    for column in preferred:
        if column in available or column in generated:
            if column not in ordered:
                ordered.append(column)

    # Keep every clean identity, raw feature and scope/provenance column.  Only
    # source columns replaced by a same-named compatibility alias are skipped.
    replaced = {"uniprot_id", "sequence", "ncbi_gene_id"}
    for column in input_columns:
        if column not in replaced and column not in ordered:
            ordered.append(column)
    return ordered


def _selected_columns(all_columns, requested):
    if not requested:
        return all_columns
    if isinstance(requested, str):
        requested = [item.strip() for item in requested.split(",") if item.strip()]
    unknown = sorted(set(requested).difference(all_columns))
    if unknown:
        raise ValueError(f"unknown final output columns: {unknown}")
    # protein_key is the stable join key and is always retained in a narrowed
    # table, even if it was not typed explicitly.
    return ["protein_key"] + [
        column for column in requested if column != "protein_key"
    ]


def finalize(input_path, output_path, manifest_path=None, columns=None):
    """Finalize one assembled Parquet table and return a verification manifest."""
    import duckdb

    input_path = os.path.abspath(input_path)
    output_path = os.path.abspath(output_path)
    if not input_path.lower().endswith(".parquet"):
        raise ValueError("finalization input must be Parquet")
    if not output_path.lower().endswith((".parquet", ".csv")):
        raise ValueError("final output must end in .parquet or .csv")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    connection = duckdb.connect()
    source = f"read_parquet({_sql_string(input_path)})"
    input_columns = [
        row[0]
        for row in connection.execute(f"DESCRIBE SELECT * FROM {source}").fetchall()
    ]
    required = {
        "protein_key", "row_kind", "sequence", "sequence_sha256",
        "uniprot_id", "uniprot_parent_ids", "refseq_protein_ids",
        "ncbi_isoform_names", "ensembl_gene_ids", "ensembl_transcript_ids",
        "ensembl_protein_ids", "ncbi_gene_id", "gene_symbol", "gene_description",
    }
    missing = sorted(required.difference(input_columns))
    if missing:
        raise ValueError(f"assembled input is missing required columns: {missing}")

    canonical_only_representatives = {
        "UniProt domains": "Domains",
        "domain CIDER": "Domains_FCR",
        "UniProt GO": "C_ids",
        "GO-derived roles": "role_in_transcription",
        "PTM": "ptm_phosphorylation",
        "CD-CODE": "UID",
        "RCSB/PDB": "RCSB_PDB_IDs",
    }
    scope_validation = {}
    for family, column in canonical_only_representatives.items():
        if column not in input_columns:
            continue
        leaked = connection.execute(
            f"SELECT count(*) FROM {source} "
            f"WHERE row_kind <> 'swissprot_canonical' "
            f"AND {_quote_identifier(column)} IS NOT NULL"
        ).fetchone()[0]
        scope_validation[family] = {
            "representative_column": column,
            "noncanonical_nonnull_rows": leaked,
        }
        if leaked:
            raise ValueError(
                f"{family} leaked into {leaked} sequence-distinct rows via {column}"
            )

    sequence_derived_representatives = {
        "whole-sequence CIDER": "FCR",
        "metapredict IDR": "IDR_count",
        "IDR CIDER": "IDR_FCR",
        "PSLab": "mean_lambda",
    }
    sequence_validation = {}
    for family, column in sequence_derived_representatives.items():
        if column not in input_columns:
            continue
        missing_values = connection.execute(
            f"SELECT count(*) FROM {source} "
            f"WHERE {_quote_identifier(column)} IS NULL"
        ).fetchone()[0]
        sequence_validation[family] = {
            "representative_column": column,
            "null_rows": missing_values,
        }
        if missing_values:
            raise ValueError(
                f"{family} is missing on {missing_values} sequence rows via {column}"
            )

    string_aliases = []
    compatibility_sources = {
        "ENSP_clean": "string_query_ensp_ids",
        "PPI_ENSP_Partners": "string_partners_ensp_by_query",
        "PPI_UniProt_Partners": "string_partners_uniprot_by_query",
        "PPI_ENSP_Partners_in_Dataframe": (
            "string_partners_ensp_in_catalog_by_query"
        ),
        "PPI_UniProt_Partners_in_Dataframe": (
            "string_partners_uniprot_by_query"
        ),
    }
    for alias, source_column in compatibility_sources.items():
        if source_column in input_columns:
            string_aliases.append(
                f"b.{_quote_identifier(source_column)} AS {_quote_identifier(alias)}"
            )
        else:
            string_aliases.append(f"NULL::VARCHAR AS {_quote_identifier(alias)}")

    unique_expression = (
        "b.rbp_census_unique"
        if "rbp_census_unique" in input_columns else "NULL::VARCHAR"
    )

    exclusions = ", ".join(_quote_identifier(column) for column in (
        "uniprot_id", "sequence", "ncbi_gene_id"
    ))
    prepared_query = f"""
        WITH source_table AS (
            SELECT * FROM {source}
        ),
        canonical AS (
            SELECT uniprot_id, sequence_sha256
            FROM source_table
            WHERE row_kind = 'swissprot_canonical' AND uniprot_id IS NOT NULL
        ),
        dominance AS (
            SELECT
                b.protein_key,
                max(CASE WHEN c.sequence_sha256 = b.sequence_sha256 THEN 1 ELSE 0 END)
                    ::INTEGER AS dominant_isoform
            FROM source_table b
            LEFT JOIN UNNEST(b.uniprot_parent_ids) p(parent_id) ON TRUE
            LEFT JOIN canonical c ON c.uniprot_id = p.parent_id
            GROUP BY b.protein_key
        ),
        prepared AS (
            SELECT
                COALESCE(
                    b.uniprot_id,
                    NULLIF(array_to_string(b.uniprot_parent_ids, ';'), '')
                ) AS uniprot_id,
                d.dominant_isoform,
                b.sequence,
                {unique_expression} AS "UNIQUE",
                array_to_string(b.refseq_protein_ids, ',') AS "ProteinHGVS",
                array_to_string(b.ncbi_isoform_names, '; ') AS "HGVSDescription",
                array_to_string(b.ensembl_gene_ids, ';') AS "ENSG",
                array_to_string(b.ensembl_transcript_ids, ';') AS "ENST",
                array_to_string(b.ensembl_protein_ids, ';') AS "ENSP",
                b.ncbi_gene_id,
                array_to_string(b.ensembl_gene_ids, ';') AS "ID",
                b.gene_symbol AS "Name",
                b.gene_description AS "Description",
                row_number() OVER (
                    PARTITION BY COALESCE(
                        b.uniprot_id,
                        list_extract(b.uniprot_parent_ids, 1),
                        b.protein_key
                    )
                    ORDER BY d.dominant_isoform DESC, b.protein_key
                ) AS isoform_number,
                CAST(to_json(b.ensembl_gene_ids) AS VARCHAR) AS "ID_list",
                b.* EXCLUDE ({exclusions}),
                {', '.join(string_aliases)},
                COALESCE(
                    b.uniprot_id,
                    list_extract(b.uniprot_parent_ids, 1),
                    b.protein_key
                ) AS __sort_uniprot,
                d.dominant_isoform AS __sort_dominant,
                COALESCE(array_to_string(b.refseq_protein_ids, ';'), b.protein_key)
                    AS __sort_refseq
            FROM source_table b
            JOIN dominance d USING (protein_key)
        )
    """
    output_columns = _selected_columns(_output_columns(input_columns), columns)
    projection = ", ".join(_quote_identifier(column) for column in output_columns)
    final_query = (
        prepared_query
        + f" SELECT {projection} FROM prepared"
        + " ORDER BY __sort_uniprot, __sort_dominant DESC, __sort_refseq, protein_key"
    )

    temporary = output_path + ".part"
    if os.path.exists(temporary):
        os.remove(temporary)
    if output_path.lower().endswith(".parquet"):
        options = "FORMAT PARQUET, COMPRESSION ZSTD"
        output_reader = f"read_parquet({_sql_string(temporary)})"
    else:
        options = "FORMAT CSV, HEADER true"
        output_reader = (
            f"read_csv_auto({_sql_string(temporary)}, header=true, sample_size=-1)"
        )
    connection.execute(f"COPY ({final_query}) TO {_sql_string(temporary)} ({options})")

    input_rows = connection.execute(f"SELECT count(*) FROM {source}").fetchone()[0]
    checks = connection.execute(prepared_query + """
        SELECT
            count(*) AS output_rows,
            sum(CASE WHEN dominant_isoform = 1 THEN 1 ELSE 0 END) AS dominant_rows,
            sum(CASE WHEN row_kind = 'swissprot_canonical' THEN 1 ELSE 0 END)
                AS canonical_rows,
            sum(CASE WHEN dominant_isoform NOT IN (0, 1) THEN 1 ELSE 0 END)
                AS invalid_dominant_rows,
            sum(CASE WHEN dominant_isoform = 1
                          AND row_kind <> 'swissprot_canonical' THEN 1 ELSE 0 END)
                AS noncanonical_dominant_rows,
            count(DISTINCT protein_key) AS unique_protein_keys
        FROM prepared
    """).fetchone()
    output_counts = connection.execute(
        f"SELECT count(*), count(DISTINCT protein_key) FROM {output_reader}"
    ).fetchone()
    connection.close()
    if checks[0] != input_rows or checks[5] != input_rows \
            or output_counts[0] != input_rows or output_counts[1] != input_rows:
        raise ValueError("finalization changed row count or duplicated protein_key")
    if checks[1] != checks[2] or checks[3] or checks[4]:
        raise ValueError(
            "dominant_isoform validation failed: "
            f"dominant={checks[1]}, canonical={checks[2]}, "
            f"invalid={checks[3]}, noncanonical_dominant={checks[4]}"
        )
    os.replace(temporary, output_path)

    manifest = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "input": input_path,
        "output": output_path,
        "rows": input_rows,
        "columns": len(output_columns),
        "canonical_rows": checks[2],
        "dominant_isoform_rows": checks[1],
        "dominant_isoform_definition": (
            "1 iff the row sequence SHA-256 equals the canonical Swiss-Prot "
            "sequence SHA-256 of at least one mapped UniProt parent; otherwise 0"
        ),
        "sort": (
            "first mapped UniProt accession; dominant sequence first; then "
            "RefSeq accession/protein_key"
        ),
        "canonical_only_scope_validation": scope_validation,
        "sequence_derived_coverage_validation": sequence_validation,
        "output_sha256": _sha256(output_path),
    }
    manifest_path = manifest_path or output_path + ".manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="assembled clean Parquet table")
    parser.add_argument("--output", required=True, help="final Parquet or CSV table")
    parser.add_argument("--manifest")
    parser.add_argument(
        "--columns",
        help="comma-separated final columns; protein_key is always retained",
    )
    args = parser.parse_args(argv)
    manifest = finalize(args.input, args.output, args.manifest, args.columns)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
