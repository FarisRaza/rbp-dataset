"""Select catalog rows by identifiers or exact sequences.

Identifier selectors are deliberately broad: a canonical UniProt accession is
also an alias of its mapped NCBI isoform rows, so selecting ``P12345`` keeps the
protein family (canonical plus mapped isoforms).  FASTA selectors are narrow:
they retain rows whose amino-acid sequence is an exact match.
"""

from __future__ import annotations

import hashlib
import os
import re


IDENTIFIER_COLUMNS = [
    "protein_key",
    "gene_symbol",
    "ncbi_gene_id",
    "ncbi_gene_ids",
    "hgnc_ids",
    "uniprot_id",
    "uniprot_secondary_accessions",
    "uniprot_parent_ids",
    "uniprot_isoform_ids",
    "swissprot_canonical_accessions",
    "refseq_protein_ids",
    "refseq_transcript_ids",
    "ensembl_gene_ids",
    "ensembl_transcript_ids",
    "ensembl_protein_ids",
]


def _normalize_identifier(value):
    return str(value or "").strip().upper()


def _versionless(value):
    return re.sub(r"\.\d+$", "", value)


def _identifier_variants(value):
    value = _normalize_identifier(value)
    if not value:
        return set()
    variants = {value, _versionless(value)}
    if ":" in value:
        variants.add(value.split(":", 1)[1])
    return {item for item in variants if item}


def row_aliases(row):
    """Return normalized identifier aliases for one catalog row."""
    aliases = set()
    for column in IDENTIFIER_COLUMNS:
        value = row.get(column)
        values = value if isinstance(value, (list, tuple, set)) else [value]
        for item in values:
            aliases.update(_identifier_variants(item))
    return aliases


def parse_inline_identifiers(value):
    """Parse comma/semicolon/whitespace-separated identifiers."""
    if not value:
        return []
    return [item for item in re.split(r"[,;\s]+", value) if item]


def read_identifier_file(path):
    """Read identifiers from a text file; blank lines and comments are ignored."""
    identifiers = []
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            identifiers.extend(parse_inline_identifiers(line))
    return identifiers


def read_fasta(path):
    """Yield ``(header, normalized_sequence)`` records without extra packages."""
    header = None
    sequence = []
    with open(path, encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper().rstrip("*")
                header = line[1:].strip()
                sequence = []
            elif header is None:
                raise ValueError(f"FASTA sequence encountered before a header in {path}")
            else:
                sequence.append("".join(line.split()))
    if header is not None:
        yield header, "".join(sequence).upper().rstrip("*")


def file_sha256(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_spec(proteins=None, protein_list=None, protein_fasta=None):
    """Stable selector description stored in the catalog/run manifests."""
    return {
        "inline_identifiers": sorted(
            {_normalize_identifier(x) for x in parse_inline_identifiers(proteins)}
        ),
        "identifier_file": (
            {"path": os.path.abspath(protein_list), "sha256": file_sha256(protein_list)}
            if protein_list else None
        ),
        "fasta_file": (
            {"path": os.path.abspath(protein_fasta), "sha256": file_sha256(protein_fasta)}
            if protein_fasta else None
        ),
    }


def has_selection(spec):
    return bool(
        spec.get("inline_identifiers")
        or spec.get("identifier_file")
        or spec.get("fasta_file")
    )


def select_rows(
    rows,
    proteins=None,
    protein_list=None,
    protein_fasta=None,
    strict=False,
):
    """Return selected rows and an audit dictionary.

    The three selectors are combined by union. Identifier matching is
    case-insensitive and version-insensitive for RefSeq/Ensembl accessions.
    FASTA records match normalized amino-acid sequences exactly.
    """
    rows = list(rows)
    requested = parse_inline_identifiers(proteins)
    if protein_list:
        requested.extend(read_identifier_file(protein_list))
    normalized = sorted(
        {_normalize_identifier(item) for item in requested if _normalize_identifier(item)}
    )

    matched_identifiers = set()
    selected_keys = set()
    if normalized:
        requested_variants = {
            identifier: _identifier_variants(identifier) for identifier in normalized
        }
        for row in rows:
            aliases = row_aliases(row)
            for identifier, variants in requested_variants.items():
                if aliases.intersection(variants):
                    selected_keys.add(row["protein_key"])
                    matched_identifiers.add(identifier)

    fasta_records = list(read_fasta(protein_fasta)) if protein_fasta else []
    fasta_hashes = {
        hashlib.sha256(sequence.encode("ascii")).hexdigest(): header
        for header, sequence in fasta_records
        if sequence
    }
    matched_fasta_hashes = set()
    for row in rows:
        digest = row.get("sequence_sha256")
        if digest in fasta_hashes:
            selected_keys.add(row["protein_key"])
            matched_fasta_hashes.add(digest)

    selector_used = bool(normalized or fasta_records)
    selected = (
        [row for row in rows if row["protein_key"] in selected_keys]
        if selector_used else rows
    )
    unmatched_identifiers = sorted(set(normalized).difference(matched_identifiers))
    unmatched_fasta = [
        header for digest, header in fasta_hashes.items()
        if digest not in matched_fasta_hashes
    ]
    audit = {
        "selector_used": selector_used,
        "input_catalog_rows": len(rows),
        "selected_rows": len(selected),
        "requested_identifier_count": len(normalized),
        "matched_identifier_count": len(matched_identifiers),
        "unmatched_identifiers": unmatched_identifiers,
        "requested_fasta_records": len(fasta_records),
        "matched_fasta_sequences": len(matched_fasta_hashes),
        "unmatched_fasta_headers": unmatched_fasta,
    }
    if strict and (unmatched_identifiers or unmatched_fasta):
        raise ValueError(
            "selection did not match every request: "
            f"identifiers={unmatched_identifiers}, FASTA={unmatched_fasta}"
        )
    if selector_used and not selected:
        raise ValueError("protein selection matched zero catalog rows")
    return selected, audit
