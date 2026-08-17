# Identifiers and catalog rows

The catalog contains one guaranteed canonical sequence row per reviewed human
UniProtKB/Swiss-Prot accession plus sequence-unique curated NCBI RefSeq `NP_`
products mapped through NCBI Gene ID. Exact amino-acid-sequence matches recover
additional Ensembl identifiers when direct cross-references are incomplete.

The final table adds compact compatibility aliases, while the clean catalog
retains explicit lists, mapping methods, ambiguity flags, source releases, and
the stable `protein_key`. See the [column reference](../COLUMN_REFERENCE.md#identifiers-and-row-provenance).

Sanity checks verify unique row keys, sequence lengths and hashes, a UniProt
accession on every canonical row, and the expected `NP_` prefix on all retained
RefSeq protein accessions.

Run `python qc/check_identifiers.py --work-dir WORK_DIR` to produce catalog
composition, identifier-coverage figures, and a machine-readable check summary.

