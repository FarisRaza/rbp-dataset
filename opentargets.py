"""Open Targets feature family -- disease associations, expression, target annotation.

Sources
-------
Three CSV exports of the Open Targets Platform, all keyed on Ensembl gene id:

  association_by_datatype_direct_full.csv  one row per (target, disease, datatype)
  expression_all.csv                       one row per gene, column ``tissues``
  target_full.csv                          one row per gene, 28 annotation columns

Produces schema.OPENTARGETS_COLUMNS (33 columns).

Shape
-----
Every one of the 33 columns is a **dict keyed by ENSG**, not a bare value::

    approvedSymbol = {'ENSG00000148584': 'A1CF'}
    diseaseId      = {'ENSG00000148584': ['DOID_0050890', 'EFO_0000094', ...]}

That is because the table's ``ID`` column can hold several semicolon-separated
ENSGs for one protein (39 of the original genes do). ``ID_list`` is ``ID`` split
on ";", and each Open Targets column maps every element of ``ID_list`` to its
value.

Nulls are encoded two different ways, and the distinction is real -- see
`_target_value`. A gene missing from the source entirely gives ``[]``; a gene
present with that particular field blank gives ``nan``.

The four association columns are lists gathered by grouping the association
table on ``targetId``; ``tissues`` is a single value per gene; the remaining 28
are every column of ``target_full.csv`` except its ``id`` key.

Memory
------
These files total roughly 2.9 GB and ``expression_all.csv``'s ``tissues`` field
alone can exceed 50 KB per gene, so nothing is loaded wholesale. Each loader
streams its file and keeps only rows whose ENSG is actually wanted.

A caveat inherited from the source files: many values are stringified numpy
reprs (``"['ENST00000373993' 'ENST00000373997']"``) rather than clean CSV. They
are carried through verbatim so the new rows match the existing ones.
"""

import csv

csv.field_size_limit(1 << 30)

ASSOCIATION_FILE = "association_by_datatype_direct_full.csv"
EXPRESSION_FILE = "expression_all.csv"
TARGET_FILE = "target_full.csv"

# Columns gathered into per-target lists from the association table.
ASSOCIATION_COLUMNS = ["diseaseId", "datatypeId", "score", "evidenceCount"]

# Numeric association fields, cast on the way in so the stored lists hold
# numbers rather than strings.
_NUMERIC = {"score": float, "evidenceCount": int}


def load_associations(path, wanted_ensgs):
    """Group the association table on ``targetId``, keeping only wanted genes.

    Returns ``{ensg: {column: [values...]}}``. The four lists for a gene are
    positionally parallel: index *i* of ``diseaseId``, ``datatypeId``, ``score``
    and ``evidenceCount`` all describe the same association record.
    """
    wanted = set(wanted_ensgs)
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ensg = row.get("targetId")
            if ensg not in wanted:
                continue
            bucket = out.setdefault(ensg, {c: [] for c in ASSOCIATION_COLUMNS})
            for column in ASSOCIATION_COLUMNS:
                value = row.get(column)
                cast = _NUMERIC.get(column)
                if cast is not None:
                    try:
                        value = cast(value)
                    except (TypeError, ValueError):
                        value = None
                bucket[column].append(value)
    return out


def load_expression(path, wanted_ensgs):
    """Return ``{ensg: tissues_string}`` for wanted genes only."""
    wanted = set(wanted_ensgs)
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ensg = row.get("id")
            if ensg in wanted:
                out[ensg] = row.get("tissues")
    return out


def load_targets(path, wanted_ensgs):
    """Return ``({ensg: {column: value}}, column_names)`` for wanted genes only.

    `column_names` is every column of the file except ``id``, in file order --
    which is the order those 28 columns appear in the table schema.
    """
    wanted = set(wanted_ensgs)
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        columns = [c for c in reader.fieldnames if c != "id"]
        for row in reader:
            ensg = row.get("id")
            if ensg in wanted:
                out[ensg] = {c: row.get(c) for c in columns}
    return out, columns


def columns_for(ensg_list, associations, expression, targets, target_columns):
    """Build all 33 Open Targets columns for one protein.

    `ensg_list` is the protein's ``ID_list``. An empty list yields empty dicts
    throughout, which is what the table stores for proteins with no ENSG.
    """
    out = {}

    for column in ASSOCIATION_COLUMNS:
        out[column] = {
            ensg: associations.get(ensg, {}).get(column, []) for ensg in ensg_list
        }

    # Historical quirk preserved: a gene absent from the expression table stores
    # an empty list here rather than None, because a late pass replaced every
    # None in the dict-valued columns with [].
    out["tissues"] = {
        ensg: (expression.get(ensg) if expression.get(ensg) is not None else [])
        for ensg in ensg_list
    }

    for column in target_columns:
        out[column] = {
            ensg: _target_value(targets, ensg, column) for ensg in ensg_list
        }

    return out


def _target_value(targets, ensg, column):
    """One target-annotation cell, using the master's two distinct null encodings.

    The table distinguishes them, and downstream code that tests ``pd.isna``
    would read them differently, so they are reproduced rather than collapsed:

      gene absent from target_full.csv entirely -> ``[]``
      gene present but this field blank         -> ``nan``

    The second case dominates for sparsely-populated columns: across a 4,000-row
    sample of the master, ``tep`` holds ``nan`` in 3,816 cells and ``[]`` in 193
    -- and that 193 is identical across every target column, i.e. exactly the
    genes missing from the source.
    """
    record = targets.get(ensg)
    if record is None:
        return []
    value = record.get(column)
    if value is None or value == "":
        return float("nan")
    return value


def empty_columns(target_columns):
    """The 33 columns for a protein with no ENSG at all."""
    out = {c: {} for c in ASSOCIATION_COLUMNS}
    out["tissues"] = {}
    for column in target_columns:
        out[column] = {}
    return out
