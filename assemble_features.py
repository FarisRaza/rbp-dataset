"""Assemble a clean catalog and keyed feature sidecars into the final table.

DuckDB performs streaming left joins, so the full wide table need not fit in
Python memory.  Every sidecar must contain exactly one row per ``protein_key``;
duplicate keys are rejected before writing output.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os


def _sql_string(value):
    return "'" + os.path.abspath(value).replace("'", "''") + "'"


def _table_expression(path):
    lower = path.lower()
    quoted = _sql_string(path)
    if lower.endswith(".parquet"):
        return f"read_parquet({quoted})"
    if lower.endswith(".csv"):
        return f"read_csv_auto({quoted}, header=true, sample_size=-1)"
    if lower.endswith(".jsonl"):
        return f"read_json_auto({quoted}, format='newline_delimited')"
    raise ValueError(f"unsupported table type: {path}")


def _sha256(path, chunk=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def assemble(catalog, sidecars, output, manifest_path=None):
    """Left-join ``sidecars`` to ``catalog`` in the supplied family order."""
    import duckdb

    catalog = os.path.abspath(catalog)
    sidecars = [os.path.abspath(path) for path in sidecars]
    output = os.path.abspath(output)
    os.makedirs(os.path.dirname(output), exist_ok=True)

    connection = duckdb.connect()
    base = _table_expression(catalog)
    base_counts = connection.execute(
        f"SELECT count(*) n, count(DISTINCT protein_key) u FROM {base}"
    ).fetchone()
    if base_counts[0] != base_counts[1]:
        raise ValueError("catalog contains duplicate protein_key values")

    select = ["b.*"]
    joins = []
    seen_columns = {
        row[0]
        for row in connection.execute(f"DESCRIBE SELECT * FROM {base}").fetchall()
    }
    family_stats = []
    for number, path in enumerate(sidecars):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        expression = _table_expression(path)
        alias = f"f{number}"
        counts = connection.execute(
            f"SELECT count(*) n, count(DISTINCT protein_key) u FROM {expression}"
        ).fetchone()
        if counts[0] != counts[1]:
            raise ValueError(f"feature sidecar has duplicate protein_key values: {path}")
        columns = [
            row[0]
            for row in connection.execute(
                f"DESCRIBE SELECT * FROM {expression}"
            ).fetchall()
            if row[0] != "protein_key"
        ]
        overlap = sorted(seen_columns.intersection(columns))
        if overlap:
            raise ValueError(
                f"feature sidecar {path} repeats existing columns: {overlap}"
            )
        seen_columns.update(columns)
        select.append(f"{alias}.* EXCLUDE (protein_key)")
        joins.append(f"LEFT JOIN {expression} {alias} USING (protein_key)")
        family_stats.append({"path": path, "rows": counts[0], "columns": len(columns)})

    query = f"SELECT {', '.join(select)} FROM {base} b {' '.join(joins)}"
    temporary = output + ".part"
    if os.path.exists(temporary):
        os.remove(temporary)
    if output.lower().endswith(".parquet"):
        copy_options = "FORMAT PARQUET, COMPRESSION ZSTD"
        temporary_expression = f"read_parquet({_sql_string(temporary)})"
    elif output.lower().endswith(".csv"):
        copy_options = "FORMAT CSV, HEADER true"
        temporary_expression = (
            f"read_csv_auto({_sql_string(temporary)}, header=true, sample_size=-1)"
        )
    else:
        raise ValueError("final output must end in .parquet or .csv")
    connection.execute(
        f"COPY ({query}) TO {_sql_string(temporary)} ({copy_options})"
    )
    final_counts = connection.execute(
        f"SELECT count(*) n FROM {temporary_expression}"
    ).fetchone()[0]
    if final_counts != base_counts[0]:
        raise ValueError(
            f"assembled row count changed from {base_counts[0]} to {final_counts}"
        )
    connection.close()
    os.replace(temporary, output)

    manifest = {
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "catalog": catalog,
        "catalog_rows": base_counts[0],
        "sidecars": family_stats,
        "output": output,
        "output_rows": final_counts,
        "output_columns": len(seen_columns),
        "output_sha256": _sha256(output),
    }
    manifest_path = manifest_path or output + ".manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--feature", action="append", default=[], dest="features")
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest")
    args = parser.parse_args(argv)
    manifest = assemble(args.catalog, args.features, args.output, args.manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
