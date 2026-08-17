"""Schema for the clean-room human protein/isoform rebuild.

The historical 157-column table grew organically and overloaded several
columns (for example ``ProteinHGVS`` can contain more than one accession).
The clean rebuild keeps identifiers explicit and uses one row per unique
protein sequence per NCBI gene, plus one guaranteed canonical row per reviewed
Swiss-Prot accession.

Feature-family columns are added by :mod:`feature_families`; this module owns
only the stable identity/provenance columns that every row must have.
"""

from __future__ import annotations


BASE_COLUMNS = [
    # Stable row identity and sequence
    "protein_key",
    "row_kind",  # swissprot_canonical | ncbi_isoform
    "sequence",
    "length_aa",
    "sequence_sha256",
    "sequence_source",
    # Gene identity
    "tax_id",
    "gene_symbol",
    "gene_synonyms",
    "gene_description",
    "ncbi_gene_id",
    "ncbi_gene_ids",
    "hgnc_ids",
    # UniProt identity
    "uniprot_id",
    "uniprot_secondary_accessions",
    "uniprot_entry_name",
    "uniprot_parent_ids",
    "uniprot_isoform_ids",
    "swissprot_canonical_accessions",
    "is_swissprot_canonical",
    # RefSeq identity. Lists are JSON-encoded in CSV and native lists in Parquet.
    "refseq_protein_ids",
    "refseq_transcript_ids",
    "ncbi_isoform_names",
    # Ensembl identity (versioned IDs where the source supplies a version)
    "ensembl_gene_ids",
    "ensembl_transcript_ids",
    "ensembl_protein_ids",
    # Mapping audit
    "identifier_mapping_methods",
    "identifier_ambiguity",
    "canonical_match_method",
    # Reproducibility
    "swissprot_release",
    "ncbi_annotation_release",
    "ensembl_release",
    "build_timestamp_utc",
]


LIST_COLUMNS = {
    "hgnc_ids",
    "gene_synonyms",
    "ncbi_gene_ids",
    "uniprot_secondary_accessions",
    "uniprot_parent_ids",
    "uniprot_isoform_ids",
    "swissprot_canonical_accessions",
    "refseq_protein_ids",
    "refseq_transcript_ids",
    "ncbi_isoform_names",
    "ensembl_gene_ids",
    "ensembl_transcript_ids",
    "ensembl_protein_ids",
    "identifier_mapping_methods",
    "identifier_ambiguity",
}


def validate_base_row(row):
    """Raise ``ValueError`` when a catalog row violates a hard invariant."""
    missing = [c for c in BASE_COLUMNS if c not in row]
    if missing:
        raise ValueError(f"base row is missing columns: {missing}")
    if not row["protein_key"]:
        raise ValueError("protein_key may not be blank")
    if row["row_kind"] not in {"swissprot_canonical", "ncbi_isoform"}:
        raise ValueError(f"invalid row_kind: {row['row_kind']!r}")
    seq = row["sequence"]
    if not isinstance(seq, str) or not seq:
        raise ValueError(f"{row['protein_key']}: sequence is empty")
    if int(row["length_aa"]) != len(seq):
        raise ValueError(f"{row['protein_key']}: length_aa disagrees with sequence")
    if row["row_kind"] == "swissprot_canonical" and not row["uniprot_id"]:
        raise ValueError("canonical Swiss-Prot rows require uniprot_id")
    for column in LIST_COLUMNS:
        if not isinstance(row[column], list):
            raise ValueError(f"{row['protein_key']}: {column} must be a list")
