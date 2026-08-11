"""Shared table I/O for the clean rebuild and independent feature scripts."""

from __future__ import annotations

import csv
import json
import math
import os

from rebuild_schema import BASE_COLUMNS, LIST_COLUMNS


def _json_default(value):
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _json_ready(value):
    """Convert nested numpy values and non-finite floats to strict JSON."""
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        return _json_ready(value.item())
    return value


def write_rows(rows, path, columns=None):
    """Write rows to CSV, JSONL or Parquet based on ``path`` suffix."""
    rows = list(rows)
    columns = list(columns or (rows[0].keys() if rows else BASE_COLUMNS))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lower = path.lower()
    if lower.endswith(".jsonl"):
        with open(path, "w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, default=_json_default, sort_keys=True) + "\n")
        return path
    if lower.endswith(".parquet"):
        import pandas as pd

        pd.DataFrame(rows, columns=columns).to_parquet(path, index=False)
        return path
    if not lower.endswith(".csv"):
        raise ValueError("output must end in .csv, .jsonl or .parquet")
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            encoded = {}
            for column in columns:
                value = row.get(column)
                if isinstance(value, (list, dict, tuple)):
                    value = json.dumps(value, default=_json_default, sort_keys=True)
                encoded[column] = value
            writer.writerow(encoded)
    return path


def read_rows(path):
    """Read CSV, JSONL or Parquet into a list of dictionaries."""
    lower = path.lower()
    if lower.endswith(".jsonl"):
        with open(path, encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]
    if lower.endswith(".parquet"):
        import pandas as pd

        rows = pd.read_parquet(path).to_dict("records")
    else:
        with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
            rows = list(csv.DictReader(handle))
    for row in rows:
        for column in LIST_COLUMNS:
            if column not in row:
                continue
            value = row.get(column)
            if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
                value = value.tolist()
                row[column] = value
            if value is None or (isinstance(value, float) and math.isnan(value)):
                row[column] = []
                continue
            if isinstance(value, tuple):
                row[column] = list(value)
                continue
            if isinstance(value, list):
                continue
            if not isinstance(value, str):
                continue
            try:
                row[column] = json.loads(value) if value else []
            except json.JSONDecodeError:
                row[column] = [value] if value else []
        for column in ("length_aa", "tax_id"):
            value = row.get(column)
            if value in (None, "") or (
                isinstance(value, float) and math.isnan(value)
            ):
                continue
            row[column] = int(value)
        if "is_swissprot_canonical" in row:
            row["is_swissprot_canonical"] = str(
                row["is_swissprot_canonical"]
            ).lower() in {"1", "true", "yes"}
    return rows


def write_feature_rows(rows, path, columns):
    """Write a feature sidecar with all nested values encoded as strict JSON.

    Encoding feature containers as strings avoids a Parquet pitfall where a
    dict keyed by ENSG is inferred as a struct with one field per human gene.
    It also gives CSV and Parquet identical, language-neutral cell contents.
    """
    def encode(row):
        item = {}
        for column in columns:
            value = row.get(column)
            if isinstance(value, (list, dict, tuple)):
                value = json.dumps(_json_ready(value), sort_keys=True)
            elif isinstance(value, float) and not math.isfinite(value):
                value = None
            item[column] = value
        return item

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lower = path.lower()
    encoded_rows = (encode(row) for row in rows)
    if lower.endswith(".csv"):
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(encoded_rows)
        return path
    if lower.endswith(".jsonl"):
        with open(path, "w", encoding="utf-8") as handle:
            for row in encoded_rows:
                handle.write(json.dumps(row, default=_json_default) + "\n")
        return path
    if not lower.endswith(".parquet"):
        raise ValueError("feature output must end in .csv, .jsonl or .parquet")

    # Stream Parquet in bounded batches. A prebuffer establishes scalar types;
    # nested values are already JSON strings, so Open Targets dict keys do not
    # become tens of thousands of Parquet struct fields.
    import itertools
    import pyarrow as pa
    import pyarrow.parquet as pq

    prebuffer = list(itertools.islice(encoded_rows, 10_000))
    type_by_column = {}
    for column in columns:
        kinds = {
            type(row[column])
            for row in prebuffer
            if row.get(column) is not None
        }
        if column == "protein_key" or str in kinds:
            arrow_type = pa.string()
        elif float in kinds:
            arrow_type = pa.float64()
        elif int in kinds:
            arrow_type = pa.int64()
        elif bool in kinds:
            arrow_type = pa.bool_()
        else:
            arrow_type = pa.string()
        type_by_column[column] = arrow_type
    schema = pa.schema([(column, type_by_column[column]) for column in columns])

    def normalize(row):
        item = {}
        for column in columns:
            value = row.get(column)
            if hasattr(value, "item") and not isinstance(value, (str, bytes)):
                value = value.item()
            arrow_type = type_by_column[column]
            if value is not None and pa.types.is_string(arrow_type):
                value = str(value)
            item[column] = value
        return item

    temporary = path + ".part"
    if os.path.exists(temporary):
        os.remove(temporary)
    writer = pq.ParquetWriter(temporary, schema, compression="zstd")
    try:
        source = itertools.chain(prebuffer, encoded_rows)
        while True:
            batch = list(itertools.islice(source, 2_000))
            if not batch:
                break
            table = pa.Table.from_pylist([normalize(row) for row in batch], schema=schema)
            writer.write_table(table)
    finally:
        writer.close()
    os.replace(temporary, path)
    return path


def read_feature_rows(path):
    """Read a feature sidecar and decode cells that contain JSON containers."""
    rows = read_rows(path)
    for row in rows:
        for column, value in list(row.items()):
            if not isinstance(value, str):
                continue
            stripped = value.lstrip()
            if not stripped.startswith(("[", "{")):
                continue
            try:
                row[column] = json.loads(value)
            except json.JSONDecodeError:
                pass
    return rows
