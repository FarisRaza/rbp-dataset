# Gene Ontology

`annotate_go.py` parses cellular-component, biological-process, and
molecular-function GO cross-references from reviewed UniProt records, including
GO IDs, names, and evidence codes.

Canonical rows receive direct UniProt annotations. NCBI isoforms inherit the
union of mapped Swiss-Prot parent annotations, labeled as gene/protein-level
rather than isoform-specific evidence.

Sanity checks: ID/name/evidence lists have equal lengths within each aspect;
GO IDs start with `GO:`; source parent accessions are present; sidecar keys are
unique and cover the catalog.

The report compares coverage across C, P, and F annotations and plots the total
GO-term count per row.
