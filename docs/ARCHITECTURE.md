# Rebuild architecture

`build_dataset.py` runs four resumable stages: download, catalog, features, and
assembly/finalization.

## Catalog rules

The default catalog contains current curated human RefSeq `NP_` proteins and
reviewed human UniProtKB/Swiss-Prot canonical proteins.

1. Every reviewed Swiss-Prot accession receives one guaranteed canonical row.
2. NCBI products are mapped to reviewed UniProt parents by explicit UniProt
   RefSeq cross-reference, then shared NCBI Gene ID.
3. NCBI products with the same `(GeneID, amino-acid sequence)` are collapsed
   into one row while all RefSeq and transcript identifiers are retained.
4. An NCBI sequence exactly equal to its mapped Swiss-Prot canonical sequence
   enriches that canonical row instead of creating a duplicate.
5. A reviewed UniProt protein with no mapped NCBI product remains represented
   by its Swiss-Prot sequence.
6. NCBI products without a reviewed parent are reported separately and are not
   silently assigned.

The catalog reports include unmapped NPs, UniProt canonical fallbacks, and
ambiguous parent mappings. The stable key is `sp:<accession>` for canonical
rows and `refseq:<GeneID>:<sequence-hash>` for sequence-distinct NCBI rows.

## Identifier mapping

Each row retains list-valued mappings for UniProt parents/isoforms, RefSeq
proteins/transcripts, ENSG, ENST, ENSP, HGNC, and NCBI Gene IDs. Direct source
cross-references are preferred. Missing isoform-level Ensembl identifiers may
be filled only by exact full-protein sequence matching against the current
Ensembl peptide FASTA, restricted by ENSG where possible.

`identifier_mapping_methods` and `identifier_ambiguity` make these decisions
auditable. The pipeline never selects one ENSP and broadcasts it to every
protein isoform.

## Feature sidecars

Each requested family writes exactly one row per `protein_key` to
`features/<family>.parquet`. Feature files can be regenerated independently and
are joined only after duplicate-key checks.

Sequence-derived features are computed for all rows. UniProt-coordinate or
canonical-protein annotations remain null for sequence-distinct rows. Gene-level
Open Targets and eCLIP annotations are broadcast to mapped isoforms with an
explicit scope field.

Open Targets emits both compact nested master-table columns and normalized long
tables for tissue expression and disease/datatype associations. Its release is
pinned in the run manifest.

## Assembly and historical compatibility

`assemble_features.py` performs bounded-memory DuckDB left joins and verifies
that row count and `protein_key` uniqueness do not change. The clean assembly
retains native catalog identifiers and all provenance fields.

`finalize_table.py` then:

- adds the familiar identity aliases and feature-family block order from
  `Isoform_Post_Merge_PSLab_OpenTargets.csv`;
- replaces the historical Open Targets dump with the revised expression and
  disease columns;
- computes `dominant_isoform` by exact sequence equality to the canonical
  sequence of a mapped UniProt parent;
- places the canonical row first and sorts rows by UniProt accession;
- optionally retains only exact columns requested with `--columns`.

Parquet is the primary output. CSV is supported for compatibility but is much
larger because sequences and nested annotations are serialized as text.

## Reproducibility

Every run preserves a catalog audit, run manifest, per-family sidecars, and
output checksum manifest. Exact historical-source snapshots that cannot be
downloaded automatically are documented in [DATA_SOURCES.md](../DATA_SOURCES.md)
and should be archived with checksums outside Git.
