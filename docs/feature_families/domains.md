# Curated UniProt domains

`annotate_domains.py` extracts reviewed UniProt `DOMAIN` and `ZN_FING` feature
coordinates and their sequences from the same Swiss-Prot flat file used to
build canonical rows.

These coordinates are canonical-sequence-specific and are not broadcast to
sequence-distinct NCBI isoforms. `uniprot_domain_annotation_scope` makes this
distinction explicit.

Sanity checks: direct annotations occur only on canonical rows; ranges remain
within sequence length; the domain count dictionary agrees with ranges and
discrete sequences; all catalog keys are retained.

The generated report plots domain-field coverage and the distribution of domain
occurrences per row. Canonical-only coverage is expected, not missing data.
