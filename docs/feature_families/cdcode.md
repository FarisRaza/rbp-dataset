# CD-CODE condensates

`annotate_cdcode.py` joins proteins to biomolecular condensates using canonical
UniProt accessions. Ten parallel lists describe condensate name, UID, type,
species, molecular contents, condensatopathy, and confidence.

Sequence-distinct NCBI isoforms remain null because membership is keyed to the
canonical UniProt protein and is not isoform-specific.

Sanity checks: the ten lists have identical lengths; UID values are unique per
row; the source/member file ordering passes `verify_alignment`; sidecar keys
cover the catalog.

The generated report plots field coverage and condensate memberships per row.
Exact source reproduction requires the validated external snapshot described in
`DATA_SOURCES.md`.
