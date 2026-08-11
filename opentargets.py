"""Clean Open Targets tissue-expression and disease annotations.

The clean rebuild deliberately exposes only the two Open Targets products that
are useful in this table:

* tissue-level RNA/protein expression; and
* target--disease/condition associations with readable disease metadata.

Open Targets annotations are gene-level.  They are joined through every row's
``ensembl_gene_ids`` and may therefore be shared by a Swiss-Prot canonical row
and its NCBI isoforms.  An absent record means "not present in this Open Targets
release", not evidence that expression or association is zero.

The preferred inputs are the typed Parquet downloads published by Open Targets.
For compatibility, the loaders also understand this project's historical
``expression_all.csv`` and ``association_by_datatype_direct_full.csv`` exports.
The expression CSV contains multiline numpy representations, so it is parsed
into ordinary Python containers instead of being copied into cells verbatim.

The old 33-column target dump is retained only through the legacy helper
functions at the bottom of this file.  :mod:`feature_runner` uses the clean
interface and never reads ``target_full.csv``.
"""

from __future__ import annotations

import ast
import csv
import datetime as _datetime
import json
import os
import re
import shutil
import sqlite3
import tempfile
import urllib.error
import urllib.request
from collections import defaultdict
from functools import lru_cache


csv.field_size_limit(1 << 30)

DEFAULT_RELEASE = "25.12"
FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"

ASSOCIATION_FILE = "association_by_datatype_direct_full.csv"
EXPRESSION_FILE = "expression_all.csv"
TARGET_FILE = "target_full.csv"
ASSOCIATION_COLUMNS = ["diseaseId", "datatypeId", "score", "evidenceCount"]

COLUMNS = [
    "opentargets_tissue_expression",
    "opentargets_expression_tissue_count",
    "opentargets_disease_associations",
    "opentargets_disease_count",
    "opentargets_disease_names",
    "opentargets_therapeutic_areas",
    "opentargets_release",
    "opentargets_annotation_scope",
]

EXPRESSION_LONG_COLUMNS = [
    "ensembl_gene_id",
    "tissue_id",
    "tissue_name",
    "organs",
    "anatomical_systems",
    "rna_value",
    "rna_unit",
    "rna_zscore",
    "rna_level",
    "protein_level",
    "protein_reliability",
    "protein_cell_types",
    "opentargets_release",
]

DISEASE_LONG_COLUMNS = [
    "ensembl_gene_id",
    "disease_id",
    "disease_name",
    "disease_description",
    "therapeutic_area_ids",
    "therapeutic_area_names",
    "datatype_id",
    "score",
    "evidence_count",
    "opentargets_release",
]


def _parquet_files(path):
    """Return Parquet files represented by a file or Spark-style directory."""
    if os.path.isfile(path):
        return [path]
    if not os.path.isdir(path):
        raise FileNotFoundError(path)
    files = []
    for root, _, names in os.walk(path):
        files.extend(
            os.path.join(root, name)
            for name in names
            if name.lower().endswith((".parquet", ".snappy.parquet"))
        )
    if not files:
        raise FileNotFoundError(f"no Parquet files found under {path}")
    return sorted(files)


def _iter_parquet(path, columns):
    import pyarrow.dataset as ds

    dataset = ds.dataset(_parquet_files(path), format="parquet")
    for batch in dataset.scanner(columns=list(columns), batch_size=4096).to_batches():
        yield from batch.to_pylist()


def _iter_rows(path, columns):
    if os.path.isdir(path) or path.lower().endswith((".parquet", ".snappy.parquet")):
        yield from _iter_parquet(path, columns)
        return
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        missing = set(columns).difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        yield from reader


def _insert_missing_object_commas(text):
    """Repair numpy object-array reprs such as ``[{...}\n {...}]``.

    Numpy separates object-array elements with whitespace rather than commas.
    That display is readable but is not a Python literal.  This scanner inserts
    a comma only between an unquoted closing brace and the next unquoted opening
    brace, leaving quoted tissue labels untouched.
    """
    out = []
    quote = None
    escaped = False
    index = 0
    while index < len(text):
        char = text[index]
        out.append(char)
        if quote:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == quote:
                quote = None
            index += 1
            continue
        if char in {"'", '"'}:
            quote = char
            index += 1
            continue
        if char == "}":
            look = index + 1
            while look < len(text) and text[look].isspace():
                look += 1
            if look < len(text) and text[look] == "{":
                out.append(",")
        index += 1
    return "".join(out)


class _NumpyLiteralCleaner(ast.NodeTransformer):
    """Drop numpy ``array(..., dtype=object)`` wrappers from an expression."""

    def visit_Call(self, node):  # noqa: N802 - ast visitor API
        is_array = (
            isinstance(node.func, ast.Name) and node.func.id == "array"
        ) or (
            isinstance(node.func, ast.Attribute) and node.func.attr == "array"
        )
        if not is_array or not node.args:
            raise ValueError("unsupported call in Open Targets expression literal")
        return self.visit(node.args[0])

    def visit_Name(self, node):  # noqa: N802 - ast visitor API
        if node.id in {"nan", "NaN"}:
            return ast.copy_location(ast.Constant(None), node)
        if node.id in {"True", "False", "None"}:
            return node
        raise ValueError(f"unsupported name in expression literal: {node.id}")


def parse_tissues(value):
    """Parse one historical CSV ``tissues`` cell into strict containers."""
    if value is None or value == "":
        return []
    if isinstance(value, list):
        return value
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return value.tolist()
    text = _insert_missing_object_commas(str(value).strip())
    tree = ast.parse(text, mode="eval")
    tree = _NumpyLiteralCleaner().visit(tree)
    ast.fix_missing_locations(tree)
    result = ast.literal_eval(tree)
    if isinstance(result, tuple):
        result = list(result)
    if not isinstance(result, list):
        raise ValueError("Open Targets tissues value is not a list")
    return result


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return value.tolist()
    return [value]


def _number(value, cast):
    if value in (None, ""):
        return None
    try:
        return cast(value)
    except (TypeError, ValueError):
        return None


def _clean_tissue(ensg, tissue):
    rna = tissue.get("rna") or {}
    protein = tissue.get("protein") or {}
    cell_types = []
    for cell in _as_list(protein.get("cell_type")):
        if not isinstance(cell, dict):
            continue
        cell_types.append({
            "name": cell.get("name"),
            "level": _number(cell.get("level"), int),
            "reliability": (
                bool(cell.get("reliability"))
                if cell.get("reliability") is not None else None
            ),
        })
    return {
        "ensembl_gene_id": ensg,
        "tissue_id": tissue.get("efo_code"),
        "tissue_name": tissue.get("label"),
        "organs": [str(item) for item in _as_list(tissue.get("organs"))],
        "anatomical_systems": [
            str(item) for item in _as_list(tissue.get("anatomical_systems"))
        ],
        "rna_value": _number(rna.get("value"), float),
        "rna_unit": rna.get("unit") or None,
        "rna_zscore": _number(rna.get("zscore"), int),
        "rna_level": _number(rna.get("level"), int),
        "protein_level": _number(protein.get("level"), int),
        "protein_reliability": (
            bool(protein.get("reliability"))
            if protein.get("reliability") is not None else None
        ),
        "protein_cell_types": cell_types,
    }


def load_tissue_expression(path, wanted_ensgs):
    """Return ``{ENSG: [flat tissue-expression records]}`` for wanted genes."""
    wanted = set(wanted_ensgs)
    out = {}
    for row in _iter_rows(path, ("id", "tissues")):
        ensg = row.get("id")
        if ensg not in wanted:
            continue
        tissues = parse_tissues(row.get("tissues"))
        records = [
            _clean_tissue(ensg, tissue)
            for tissue in tissues
            if isinstance(tissue, dict)
        ]
        records.sort(key=lambda item: (
            item.get("tissue_name") or "", item.get("tissue_id") or ""
        ))
        out[ensg] = records
    return out


def load_association_records(path, wanted_ensgs):
    """Return raw association records grouped by Ensembl gene identifier."""
    wanted = set(wanted_ensgs)
    out = defaultdict(list)
    fields = ("diseaseId", "targetId", "datatypeId", "score", "evidenceCount")
    for row in _iter_rows(path, fields):
        ensg = row.get("targetId")
        if ensg not in wanted:
            continue
        out[ensg].append({
            "disease_id": row.get("diseaseId"),
            "datatype_id": row.get("datatypeId"),
            "score": _number(row.get("score"), float),
            "evidence_count": _number(row.get("evidenceCount"), int),
        })
    return dict(out)


def load_disease_metadata(path):
    """Load names, descriptions and named therapeutic areas by disease ID."""
    fields = ("id", "name", "description", "therapeuticAreas")
    raw = {}
    for row in _iter_rows(path, fields):
        disease_id = row.get("id")
        if not disease_id:
            continue
        areas = row.get("therapeuticAreas")
        if isinstance(areas, str):
            try:
                areas = json.loads(areas)
            except json.JSONDecodeError:
                areas = [item for item in re.split(r"[;,]", areas) if item]
        raw[disease_id] = {
            "name": row.get("name"),
            "description": row.get("description"),
            "therapeutic_area_ids": [str(item) for item in _as_list(areas)],
        }
    for record in raw.values():
        record["therapeutic_areas"] = [
            {"id": area_id, "name": raw.get(area_id, {}).get("name")}
            for area_id in record.pop("therapeutic_area_ids")
        ]
    return raw


class OpenTargetsStore:
    """Bounded-memory SQLite index used by the full-proteome feature runner."""

    def __init__(self, path, delete_on_close=False):
        self.path = path
        self.delete_on_close = delete_on_close
        self.connection = sqlite3.connect(path)

    @lru_cache(maxsize=256)
    def expression_for(self, ensg):
        row = self.connection.execute(
            "SELECT records_json FROM expression WHERE ensg = ?", (ensg,)
        ).fetchone()
        return json.loads(row[0]) if row else []

    @lru_cache(maxsize=256)
    def associations_for(self, ensg):
        return self.connection.execute(
            "SELECT disease_id, datatype_id, score, evidence_count "
            "FROM association WHERE ensg = ?", (ensg,)
        ).fetchall()

    def ensgs_with_expression(self):
        for (ensg,) in self.connection.execute("SELECT ensg FROM expression ORDER BY ensg"):
            yield ensg

    def ensgs_with_associations(self):
        query = "SELECT DISTINCT ensg FROM association ORDER BY ensg"
        for (ensg,) in self.connection.execute(query):
            yield ensg

    def close(self):
        self.expression_for.cache_clear()
        self.associations_for.cache_clear()
        self.connection.close()
        if self.delete_on_close and os.path.exists(self.path):
            os.remove(self.path)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def build_store(
    expression_path, association_path, wanted_ensgs, destination=None
):
    """Index large Open Targets sources on disk instead of retaining them in RAM."""
    wanted = set(wanted_ensgs)
    delete_on_close = destination is None
    if destination is None:
        descriptor, temporary = tempfile.mkstemp(
            prefix="opentargets_", suffix=".sqlite"
        )
        os.close(descriptor)
    else:
        destination = os.path.abspath(destination)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        temporary = destination + ".part"
        if os.path.exists(temporary):
            os.remove(temporary)

    connection = sqlite3.connect(temporary)
    try:
        connection.executescript(
            "PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; "
            "CREATE TABLE expression (ensg TEXT PRIMARY KEY, records_json TEXT NOT NULL);"
            "CREATE TABLE association ("
            "ensg TEXT NOT NULL, disease_id TEXT, datatype_id TEXT, "
            "score REAL, evidence_count INTEGER);"
        )
        expression_batch = []
        for row in _iter_rows(expression_path, ("id", "tissues")):
            ensg = row.get("id")
            if ensg not in wanted:
                continue
            records = [
                _clean_tissue(ensg, tissue)
                for tissue in parse_tissues(row.get("tissues"))
                if isinstance(tissue, dict)
            ]
            records.sort(key=lambda item: (
                item.get("tissue_name") or "", item.get("tissue_id") or ""
            ))
            expression_batch.append((ensg, json.dumps(records, separators=(",", ":"))))
            if len(expression_batch) >= 250:
                connection.executemany(
                    "INSERT OR REPLACE INTO expression VALUES (?, ?)", expression_batch
                )
                expression_batch.clear()
        if expression_batch:
            connection.executemany(
                "INSERT OR REPLACE INTO expression VALUES (?, ?)", expression_batch
            )

        association_batch = []
        fields = ("diseaseId", "targetId", "datatypeId", "score", "evidenceCount")
        for row in _iter_rows(association_path, fields):
            ensg = row.get("targetId")
            if ensg not in wanted:
                continue
            association_batch.append((
                ensg,
                row.get("diseaseId"),
                row.get("datatypeId"),
                _number(row.get("score"), float),
                _number(row.get("evidenceCount"), int),
            ))
            if len(association_batch) >= 10_000:
                connection.executemany(
                    "INSERT INTO association VALUES (?, ?, ?, ?, ?)",
                    association_batch,
                )
                association_batch.clear()
        if association_batch:
            connection.executemany(
                "INSERT INTO association VALUES (?, ?, ?, ?, ?)", association_batch
            )
        connection.execute("CREATE INDEX association_ensg ON association (ensg)")
        connection.commit()
    finally:
        connection.close()

    if destination is not None:
        os.replace(temporary, destination)
        temporary = destination
    return OpenTargetsStore(temporary, delete_on_close=delete_on_close)


def _disease_records(ensg, associations, disease_metadata):
    by_disease = defaultdict(list)
    for record in associations.get(ensg, []):
        disease_id = record.get("disease_id")
        if disease_id:
            by_disease[disease_id].append(record)

    result = []
    for disease_id, rows in by_disease.items():
        metadata = disease_metadata.get(disease_id, {})
        by_datatype = {}
        for row in rows:
            datatype = row.get("datatype_id") or "unknown"
            item = by_datatype.setdefault(datatype, {
                "datatype_id": datatype, "score": None, "evidence_count": 0
            })
            score = row.get("score")
            if score is not None and (item["score"] is None or score > item["score"]):
                item["score"] = score
            if row.get("evidence_count") is not None:
                item["evidence_count"] += row["evidence_count"]
        evidence = sorted(by_datatype.values(), key=lambda item: item["datatype_id"])
        scores = [item["score"] for item in evidence if item["score"] is not None]
        result.append({
            "ensembl_gene_id": ensg,
            "disease_id": disease_id,
            "disease_name": metadata.get("name"),
            "disease_description": metadata.get("description"),
            "therapeutic_areas": metadata.get("therapeutic_areas", []),
            # This is intentionally named max_datatype_score: it is not the
            # Open Targets overall association score.
            "max_datatype_score": max(scores) if scores else None,
            "evidence_count_total": sum(
                item["evidence_count"] for item in evidence
            ),
            "evidence_by_datatype": evidence,
        })
    result.sort(key=lambda item: (
        -(item["max_datatype_score"] if item["max_datatype_score"] is not None else -1),
        item.get("disease_name") or "",
        item["disease_id"],
    ))
    return result


def _records_from_compact(ensg, rows):
    """Expose compact SQLite tuples through the mapping used by aggregation."""
    return {ensg: [
        {
            "disease_id": disease_id,
            "datatype_id": datatype_id,
            "score": score,
            "evidence_count": evidence_count,
        }
        for disease_id, datatype_id, score, evidence_count in rows
    ]}


def clean_columns_for(
    ensg_list, expression, associations, disease_metadata, release=DEFAULT_RELEASE
):
    """Build the clean Open Targets columns for one protein/isoform row."""
    ensgs = list(dict.fromkeys(ensg_list or []))
    tissue_records = [record for ensg in ensgs for record in expression.get(ensg, [])]
    disease_records = [
        record
        for ensg in ensgs
        for record in _disease_records(ensg, associations, disease_metadata)
    ]
    disease_names = []
    seen_disease = set()
    therapeutic_areas = []
    seen_area = set()
    for record in disease_records:
        disease_id = record["disease_id"]
        if disease_id not in seen_disease:
            seen_disease.add(disease_id)
            disease_names.append({"id": disease_id, "name": record.get("disease_name")})
        for area in record.get("therapeutic_areas", []):
            area_id = area.get("id")
            if area_id not in seen_area:
                seen_area.add(area_id)
                therapeutic_areas.append(area)
    therapeutic_areas.sort(key=lambda item: (item.get("name") or "", item.get("id") or ""))
    return {
        "opentargets_tissue_expression": tissue_records,
        "opentargets_expression_tissue_count": len({
            (item["ensembl_gene_id"], item.get("tissue_id"), item.get("tissue_name"))
            for item in tissue_records
        }),
        "opentargets_disease_associations": disease_records,
        "opentargets_disease_count": len(seen_disease),
        "opentargets_disease_names": disease_names,
        "opentargets_therapeutic_areas": therapeutic_areas,
        "opentargets_release": release,
        "opentargets_annotation_scope": (
            "gene-level; keyed by ENSG and broadcast to rows sharing that gene; "
            "absence is not a negative result"
        ),
    }


def clean_columns_from_store(
    ensg_list, store, disease_metadata, release=DEFAULT_RELEASE
):
    """Bounded-memory counterpart to :func:`clean_columns_for`."""
    ensgs = list(dict.fromkeys(ensg_list or []))
    expression = {ensg: store.expression_for(ensg) for ensg in ensgs}
    associations = {}
    for ensg in ensgs:
        associations.update(_records_from_compact(ensg, store.associations_for(ensg)))
    return clean_columns_for(
        ensgs, expression, associations, disease_metadata, release
    )


def expression_long_rows(expression, release=DEFAULT_RELEASE):
    """Yield one analysis-ready row per Ensembl gene and tissue."""
    for ensg in sorted(expression):
        for record in expression[ensg]:
            yield {**record, "opentargets_release": release}


def disease_long_rows(associations, disease_metadata, release=DEFAULT_RELEASE):
    """Yield one analysis-ready row per gene, disease and evidence datatype."""
    for ensg in sorted(associations):
        for record in _disease_records(ensg, associations, disease_metadata):
            areas = record.get("therapeutic_areas", [])
            for evidence in record["evidence_by_datatype"]:
                yield {
                    "ensembl_gene_id": ensg,
                    "disease_id": record["disease_id"],
                    "disease_name": record.get("disease_name"),
                    "disease_description": record.get("disease_description"),
                    "therapeutic_area_ids": [item.get("id") for item in areas],
                    "therapeutic_area_names": [item.get("name") for item in areas],
                    "datatype_id": evidence["datatype_id"],
                    "score": evidence["score"],
                    "evidence_count": evidence["evidence_count"],
                    "opentargets_release": release,
                }


def expression_long_rows_from_store(store, release=DEFAULT_RELEASE):
    for ensg in store.ensgs_with_expression():
        for record in store.expression_for(ensg):
            yield {**record, "opentargets_release": release}


def disease_long_rows_from_store(store, disease_metadata, release=DEFAULT_RELEASE):
    for ensg in store.ensgs_with_associations():
        associations = _records_from_compact(ensg, store.associations_for(ensg))
        yield from disease_long_rows(associations, disease_metadata, release)


def _download(url, destination, force=False):
    if os.path.exists(destination) and not force:
        return destination
    os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)
    temporary = destination + ".part"
    with urllib.request.urlopen(url, timeout=120) as source, open(temporary, "wb") as sink:
        shutil.copyfileobj(source, sink, length=1024 * 1024)
    os.replace(temporary, destination)
    return destination


def _dataset_filenames(url):
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            listing = response.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as error:
        if error.code == 404 and url.rstrip("/").endswith("/expression"):
            raise ValueError(
                "this Open Targets release does not publish the compatible "
                "'expression' dataset; release 26.06 introduced the different "
                "'baseline_expression' schema, so use the pinned 25.12 default "
                "or supply a compatible expression export explicitly"
            ) from error
        raise
    return sorted(set(re.findall(r'href="([^"/]+\.parquet)"', listing, re.I)))


def download_release(destination, release=DEFAULT_RELEASE, force=False):
    """Download the three pinned Open Targets datasets needed by this family."""
    root = os.path.abspath(destination)
    base = f"{FTP_ROOT}/{release}/output"
    paths = {
        "associations": os.path.join(root, "association_by_datatype_direct"),
        "expression": os.path.join(root, "expression"),
        "diseases": os.path.join(root, "disease.parquet"),
    }
    urls = []
    for name, remote_dir in (
        ("associations", "association_by_datatype_direct"),
        ("expression", "expression"),
    ):
        directory_url = f"{base}/{remote_dir}/"
        filenames = _dataset_filenames(directory_url)
        if not filenames:
            raise FileNotFoundError(f"no Parquet files listed at {directory_url}")
        for filename in filenames:
            url = directory_url + filename
            _download(url, os.path.join(paths[name], filename), force=force)
            urls.append(url)
    disease_url = f"{base}/disease/disease.parquet"
    _download(disease_url, paths["diseases"], force=force)
    urls.append(disease_url)
    manifest = {
        "opentargets_release": release,
        "downloaded_at_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "urls": urls,
    }
    os.makedirs(root, exist_ok=True)
    with open(os.path.join(root, "download_manifest.json"), "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return paths


# ---------------------------------------------------------------------------
# Legacy compatibility for the historical expansion scripts.  The clean
# feature runner above does not call any function below this line.

_NUMERIC = {"score": float, "evidenceCount": int}


def load_associations(path, wanted_ensgs):
    """Historical parallel-list association representation."""
    wanted = set(wanted_ensgs)
    out = {}
    fields = ("diseaseId", "targetId", "datatypeId", "score", "evidenceCount")
    for row in _iter_rows(path, fields):
        ensg = row.get("targetId")
        if ensg not in wanted:
            continue
        bucket = out.setdefault(ensg, {column: [] for column in ASSOCIATION_COLUMNS})
        for column in ASSOCIATION_COLUMNS:
            value = row.get(column)
            cast = _NUMERIC.get(column)
            if cast is not None:
                value = _number(value, cast)
            bucket[column].append(value)
    return out


def load_expression(path, wanted_ensgs):
    """Historical raw ``tissues`` representation."""
    wanted = set(wanted_ensgs)
    out = {}
    for row in _iter_rows(path, ("id", "tissues")):
        if row.get("id") in wanted:
            out[row["id"]] = row.get("tissues")
    return out


def load_targets(path, wanted_ensgs):
    """Historical target-wide representation used only by old scripts."""
    if os.path.isdir(path) or path.lower().endswith((".parquet", ".snappy.parquet")):
        import pyarrow.dataset as ds

        dataset = ds.dataset(_parquet_files(path), format="parquet")
        columns = [name for name in dataset.schema.names if name != "id"]
        iterator = _iter_parquet(path, ["id"] + columns)
    else:
        handle = open(path, newline="", encoding="utf-8", errors="replace")
        reader = csv.DictReader(handle)
        columns = [name for name in reader.fieldnames if name != "id"]
        iterator = reader
    wanted = set(wanted_ensgs)
    out = {}
    try:
        for row in iterator:
            if row.get("id") in wanted:
                out[row["id"]] = {column: row.get(column) for column in columns}
    finally:
        if "handle" in locals():
            handle.close()
    return out, columns


def columns_for(ensg_list, associations, expression, targets, target_columns):
    """Historical 33-column output retained for old expansion scripts."""
    out = {}
    for column in ASSOCIATION_COLUMNS:
        out[column] = {
            ensg: associations.get(ensg, {}).get(column, []) for ensg in ensg_list
        }
    out["tissues"] = {
        ensg: expression.get(ensg) if expression.get(ensg) is not None else []
        for ensg in ensg_list
    }
    for column in target_columns:
        out[column] = {
            ensg: _target_value(targets, ensg, column) for ensg in ensg_list
        }
    return out


def _target_value(targets, ensg, column):
    record = targets.get(ensg)
    if record is None:
        return []
    value = record.get(column)
    return float("nan") if value is None or value == "" else value


def empty_columns(target_columns):
    out = {column: {} for column in ASSOCIATION_COLUMNS}
    out["tissues"] = {}
    out.update({column: {} for column in target_columns})
    return out
