# Feature-family guides

Each feature family has independent source code, input provenance, annotation
scope, and sanity checks. These pages explain what its generated report means.
The consolidated [column reference](../COLUMN_REFERENCE.md) integrates the
original spreadsheet definitions with the current schema. Run an individual
`qc/check_<family>.py` command after building a sidecar, or use
`python rbp_pipeline/generate_feature_reports.py`, to create data-specific
Markdown and figures under `WORK_DIR/reports/`.

- [Identifiers and catalog](identifiers.md)
- [CIDER](cider.md)
- [IDRs](idr.md)
- [UniProt domains](domains.md)
- [Gene Ontology](go.md)
- [eCLIP](eclip.md)
- [InterPro](interpro.md)
- [PTMs](ptm.md)
- [Open Targets](opentargets.md)
- [CD-CODE](cdcode.md)
- [STRING](string.md)
- [GO-derived regulatory roles](go_roles.md)
- [PSLab](pslab.md)
- [RCSB/PDB](rcsb.md)

See [`qc/README.md`](../../qc/README.md) for every runnable QC command and
examples of strict release validation.
