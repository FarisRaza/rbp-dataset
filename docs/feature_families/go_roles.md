# GO-derived regulatory roles

`annotate_go_roles.py` derives flags for transcription regulation, translation
regulation, and mRNA-stability regulation by applying documented text criteria
to GO biological-process and molecular-function names. A combined legacy flag
is retained.

Column definitions: [GO-derived regulatory roles](../COLUMN_REFERENCE.md#go-derived-regulatory-roles).

This family has no additional external source; it consumes the GO sidecar and
inherits its annotation scope.

Sanity checks: every flag is Boolean; the combined flag equals the logical OR
of the component flags; sidecar keys exactly cover the catalog.

The generated report compares role prevalences and plots the number of positive
role flags per row. Treat these as transparent rule-based labels, not an
independent experimental assay.

Run `python qc/check_go_roles.py --work-dir WORK_DIR`.
