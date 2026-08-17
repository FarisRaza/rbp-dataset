# Open Targets expression and disease associations

`annotate_opentargets.py` replaces the historical mixed target dump with clean,
normalized tissue-expression and disease/condition summaries from a pinned Open
Targets Platform release.

Column definitions: [Open Targets expression and disease associations](../COLUMN_REFERENCE.md#open-targets-expression-and-disease-associations).

Rows are keyed by version-normalized ENSG. The sidecar records structured per-tissue RNA/protein
values, normalized disease records, therapeutic areas, counts, and release.
Separate normalized long Parquet tables are also written for efficient analysis.

Sanity checks: release is constant; disease counts equal unique disease IDs;
tissue counts equal unique tissue IDs; scores retain their explicit datatype
meaning; catalog rows without an ENSG remain present with empty annotations.

The report plots expression/disease coverage and the combined per-row tissue
plus disease count. Analyze the long tables when comparing individual tissues
or conditions.

Run `python qc/check_opentargets.py --work-dir WORK_DIR`.
