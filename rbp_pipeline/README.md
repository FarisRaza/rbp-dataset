# Pipeline source

`rebuild_from_scratch.py` is the master implementation called by the repository
root's `build_dataset.py` command. Each `annotate_<family>.py` file is an
independently runnable feature-family entry point; its scientific logic and
source notes live in the matching module (`cider.py`, `ptm.py`, and so on).

The catalog is built first, every requested family writes one Parquet sidecar
keyed by `protein_key`, `assemble_features.py` performs checked left joins, and
`finalize_table.py` restores the historical family-block layout and sorts rows
by UniProt accession.
