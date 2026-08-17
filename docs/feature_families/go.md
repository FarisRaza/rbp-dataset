# Gene Ontology

`annotate_go.py` parses cellular-component, biological-process, and
molecular-function GO cross-references from reviewed UniProt records, including
GO IDs, names, and evidence codes.

Canonical rows receive direct UniProt annotations. Sequence-distinct NCBI
isoforms remain null because the source records are attached to the canonical
UniProt protein rather than to those RefSeq sequences.

Sanity checks: ID/name/evidence lists have equal lengths within each aspect;
GO IDs start with `GO:`; noncanonical rows remain null; sidecar keys are unique
and cover the catalog.

The report compares coverage across C, P, and F annotations and plots the total
GO-term count per row.
