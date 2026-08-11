"""Build the canonical Swiss-Prot + sequence-unique NCBI isoform catalog.

Hard invariants
---------------
* Every reviewed human Swiss-Prot accession creates one canonical row.
* NCBI rows are unique by ``(GeneID, amino-acid sequence)``.
* RefSeq transcripts/proteins that encode the same sequence are aggregated.
* An NCBI sequence identical to a canonical protein enriches that canonical
  row instead of creating a duplicate isoform row.
* ENSG/ENST/ENSP values are attached only with stated evidence; one identifier
  is never silently broadcast across unrelated isoforms.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
from collections import defaultdict

import ncbi_isoforms
import swissprot_source
from rebuild_schema import BASE_COLUMNS, LIST_COLUMNS, validate_base_row


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def versionless(accession):
    return re.sub(r"\.\d+$", "", str(accession or ""))


def _unique(values):
    return sorted({value for value in values if value not in (None, "")})


def _base_row(**values):
    row = {column: [] if column in LIST_COLUMNS else None for column in BASE_COLUMNS}
    row.update(values)
    return row


def _merge_lists(row, **updates):
    for column, values in updates.items():
        row[column] = _unique(list(row.get(column) or []) + list(values or []))


def _canonical_xrefs(meta):
    accepted_tags = {None, meta.accession, f"{meta.accession}-1"}
    refseq = [x for x in meta.refseq if x.get("uniprot_isoform") in accepted_tags]
    ensembl = [x for x in meta.ensembl if x.get("uniprot_isoform") in accepted_tags]
    return refseq, ensembl


def _gene_release(gene_reports):
    names = _unique(
        name for report in gene_reports.values() for name in report.annotation_names
    )
    return ";".join(names) if names else None


def build_catalog(
    swissprot_path,
    ncbi_gene_report,
    ncbi_product_report,
    ncbi_protein_fasta,
    include_predicted=False,
    include_orphan_refseq=False,
    ensembl_index=None,
    ensembl_release=None,
):
    """Return clean catalog rows from current Swiss-Prot and NCBI artifacts."""
    build_time = dt.datetime.now(dt.timezone.utc).isoformat()
    swiss_manifest = swissprot_source.load_manifest(swissprot_path)
    swiss_release = swiss_manifest.get("uniprot_release")
    gene_reports = ncbi_isoforms.load_gene_reports(ncbi_gene_report)
    products = ncbi_isoforms.load_product_reports(
        ncbi_product_report, include_predicted=include_predicted
    )
    sequences = ncbi_isoforms.load_protein_sequences(
        ncbi_protein_fasta, include_predicted=include_predicted
    )
    ncbi_release = _gene_release(gene_reports)

    metas = [meta for meta, _record in swissprot_source.iter_records(swissprot_path)]
    meta_by_accession = {meta.accession: meta for meta in metas}
    rows_by_accession = {}
    accessions_by_gene = defaultdict(set)
    refseq_to_uniprot = defaultdict(list)
    ensembl_by_isoform = defaultdict(list)

    for meta in metas:
        for gene_id in meta.gene_ids:
            accessions_by_gene[gene_id].add(meta.accession)
        for xref in meta.refseq:
            if xref.get("protein"):
                refseq_to_uniprot[versionless(xref["protein"])].append(
                    (meta.accession, xref.get("uniprot_isoform"))
                )
        for xref in meta.ensembl:
            tag = xref.get("uniprot_isoform")
            if tag:
                ensembl_by_isoform[tag].append(xref)

    # NCBI's gene report provides a second, independent GeneID->Swiss-Prot map.
    for gene_id, report in gene_reports.items():
        for accession in report.swiss_prot_accessions:
            if accession in meta_by_accession:
                accessions_by_gene[gene_id].add(accession)

    # Guarantee every reviewed canonical protein before considering NCBI.
    for meta in metas:
        matching_gene_reports = [
            gene_reports[gene_id] for gene_id in meta.gene_ids if gene_id in gene_reports
        ]
        canonical_refseq, canonical_ensembl = _canonical_xrefs(meta)
        ensg = [x.get("gene") for x in canonical_ensembl]
        ensg.extend(
            value
            for report in matching_gene_reports
            for value in report.ensembl_gene_ids
        )
        row = _base_row(
            protein_key=f"sp:{meta.accession}",
            row_kind="swissprot_canonical",
            sequence=meta.sequence,
            length_aa=len(meta.sequence),
            sequence_sha256=sequence_sha256(meta.sequence),
            sequence_source="UniProtKB/Swiss-Prot canonical",
            tax_id=9606,
            gene_symbol=meta.gene_symbol,
            gene_synonyms=_unique(
                meta.gene_synonyms
                + [value for report in matching_gene_reports for value in report.synonyms]
            ),
            gene_description=meta.description,
            ncbi_gene_id=meta.gene_ids[0] if len(meta.gene_ids) == 1 else None,
            ncbi_gene_ids=_unique(meta.gene_ids),
            hgnc_ids=_unique(
                meta.hgnc_ids
                + [value for report in matching_gene_reports for value in report.hgnc_ids]
            ),
            uniprot_id=meta.accession,
            uniprot_secondary_accessions=_unique(meta.secondary_accessions),
            uniprot_entry_name=meta.entry_name,
            uniprot_parent_ids=[meta.accession],
            uniprot_isoform_ids=[f"{meta.accession}-1"],
            swissprot_canonical_accessions=[meta.accession],
            is_swissprot_canonical=True,
            refseq_protein_ids=_unique(x.get("protein") for x in canonical_refseq),
            refseq_transcript_ids=_unique(x.get("transcript") for x in canonical_refseq),
            ncbi_isoform_names=[],
            ensembl_gene_ids=_unique(ensg),
            ensembl_transcript_ids=_unique(x.get("transcript") for x in canonical_ensembl),
            ensembl_protein_ids=_unique(x.get("protein") for x in canonical_ensembl),
            identifier_mapping_methods=["UniProtKB_cross_reference"],
            identifier_ambiguity=(
                ["multiple_NCBI_GeneIDs"] if len(meta.gene_ids) > 1 else []
            ),
            canonical_match_method="Swiss-Prot canonical source",
            swissprot_release=swiss_release,
            ncbi_annotation_release=ncbi_release,
            ensembl_release=ensembl_release,
            build_timestamp_utc=build_time,
        )
        rows_by_accession[meta.accession] = row

    # Multiple transcript records can point to one product accession, and
    # multiple product accessions can encode exactly the same sequence.
    groups = defaultdict(list)
    missing_fasta = set()
    for product in products:
        sequence = sequences.get(product.refseq_protein)
        if sequence is None:
            missing_fasta.add(product.refseq_protein)
            continue
        groups[(product.gene_id, sequence_sha256(sequence))].append(product)

    isoform_rows = []
    skipped_orphan_groups = 0
    for (gene_id, digest), group in sorted(groups.items()):
        sequence = sequences[group[0].refseq_protein]
        gene_report = gene_reports.get(gene_id)
        gene_candidates = sorted(accessions_by_gene.get(gene_id, set()))
        exact = [
            accession for accession in gene_candidates
            if meta_by_accession[accession].sequence == sequence
        ]
        protein_ids = _unique(product.refseq_protein for product in group)
        transcript_ids = _unique(product.refseq_transcript for product in group)
        names = _unique(product.isoform_name for product in group)
        product_enst = _unique(product.ensembl_transcript for product in group)
        product_ensp = _unique(product.ensembl_protein for product in group)

        explicit_pairs = []
        for protein_id in protein_ids:
            explicit_pairs.extend(refseq_to_uniprot.get(versionless(protein_id), []))
        explicit_parents = _unique(pair[0] for pair in explicit_pairs)
        isoform_ids = _unique(pair[1] for pair in explicit_pairs)

        # The requested universe is the reviewed Swiss-Prot proteome plus its
        # RefSeq isoforms, not every curated NCBI protein lacking a reviewed
        # Swiss-Prot parent. Keep the latter opt-in and count what was excluded.
        if not include_orphan_refseq and not (gene_candidates or explicit_parents):
            skipped_orphan_groups += 1
            continue

        ensembl_from_uniprot = [
            xref for isoform_id in isoform_ids
            for xref in ensembl_by_isoform.get(isoform_id, [])
        ]
        ensg = _unique(
            ([value for value in (gene_report.ensembl_gene_ids if gene_report else [])]
             + [xref.get("gene") for xref in ensembl_from_uniprot])
        )
        enst = _unique(product_enst + [x.get("transcript") for x in ensembl_from_uniprot])
        ensp = _unique(product_ensp + [x.get("protein") for x in ensembl_from_uniprot])

        if exact:
            for accession in exact:
                row = rows_by_accession[accession]
                _merge_lists(
                    row,
                    gene_synonyms=(gene_report.synonyms if gene_report else []),
                    ncbi_gene_ids=[gene_id],
                    hgnc_ids=(gene_report.hgnc_ids if gene_report else []),
                    refseq_protein_ids=protein_ids,
                    refseq_transcript_ids=transcript_ids,
                    ncbi_isoform_names=names,
                    ensembl_gene_ids=ensg,
                    ensembl_transcript_ids=enst,
                    ensembl_protein_ids=ensp,
                    identifier_mapping_methods=[
                        "NCBI_Gene_product_report", "exact_sequence_within_NCBI_Gene"
                    ],
                )
                row["canonical_match_method"] = "exact NCBI/Swiss-Prot sequence"
                if not row.get("ncbi_gene_id"):
                    row["ncbi_gene_id"] = gene_id
                if not row.get("gene_symbol") and gene_report:
                    row["gene_symbol"] = gene_report.symbol
                if not row.get("gene_description") and gene_report:
                    row["gene_description"] = gene_report.description
            continue

        parents = explicit_parents or gene_candidates
        ambiguities = []
        if not parents:
            ambiguities.append("no_Swiss-Prot_parent_for_NCBI_gene")
        elif len(parents) > 1:
            ambiguities.append("multiple_Swiss-Prot_parent_candidates")
        if explicit_parents and set(explicit_parents) != set(gene_candidates) and gene_candidates:
            ambiguities.append("UniProt_RefSeq_and_NCBI_gene_parent_maps_disagree")

        symbol = gene_report.symbol if gene_report else group[0].gene_symbol
        description = (
            gene_report.description if gene_report else group[0].gene_description
        )
        row = _base_row(
            protein_key=f"refseq:{gene_id}:{digest[:16]}",
            row_kind="ncbi_isoform",
            sequence=sequence,
            length_aa=len(sequence),
            sequence_sha256=digest,
            sequence_source="NCBI RefSeq protein.faa",
            tax_id=9606,
            gene_symbol=symbol,
            gene_synonyms=_unique(gene_report.synonyms if gene_report else []),
            gene_description=description,
            ncbi_gene_id=gene_id,
            ncbi_gene_ids=[gene_id],
            hgnc_ids=_unique(gene_report.hgnc_ids if gene_report else []),
            uniprot_id=None,
            uniprot_secondary_accessions=[],
            uniprot_entry_name=None,
            uniprot_parent_ids=parents,
            uniprot_isoform_ids=isoform_ids,
            swissprot_canonical_accessions=gene_candidates,
            is_swissprot_canonical=False,
            refseq_protein_ids=protein_ids,
            refseq_transcript_ids=transcript_ids,
            ncbi_isoform_names=names,
            ensembl_gene_ids=ensg,
            ensembl_transcript_ids=enst,
            ensembl_protein_ids=ensp,
            identifier_mapping_methods=_unique(
                ["NCBI_Gene_product_report"]
                + (["UniProtKB_RefSeq_isoform_cross_reference"] if explicit_pairs else [])
            ),
            identifier_ambiguity=ambiguities,
            canonical_match_method="not identical to Swiss-Prot canonical at same NCBI Gene",
            swissprot_release=swiss_release,
            ncbi_annotation_release=ncbi_release,
            ensembl_release=ensembl_release,
            build_timestamp_utc=build_time,
        )
        isoform_rows.append(row)

    rows = list(rows_by_accession.values()) + isoform_rows

    if ensembl_index is not None:
        import ensembl_ids

        for row in rows:
            mappings = ensembl_ids.exact_mappings(
                row["sequence"], ensembl_index, row["ensembl_gene_ids"]
            )
            if not mappings:
                continue
            _merge_lists(
                row,
                ensembl_gene_ids=[x.get("gene") for x in mappings],
                ensembl_transcript_ids=[x.get("transcript") for x in mappings],
                ensembl_protein_ids=[x.get("protein") for x in mappings],
                identifier_mapping_methods=["exact_sequence_to_Ensembl_peptide_FASTA"],
            )
            if len(mappings) > 1:
                _merge_lists(row, identifier_ambiguity=["multiple_exact_Ensembl_peptides"])

    rows.sort(key=lambda row: (row["row_kind"] != "swissprot_canonical", row["protein_key"]))
    seen = set()
    for row in rows:
        for column in LIST_COLUMNS:
            row[column] = _unique(row[column])
        validate_base_row(row)
        if row["protein_key"] in seen:
            raise ValueError(f"duplicate protein_key: {row['protein_key']}")
        seen.add(row["protein_key"])

    audit = {
        "rows": len(rows),
        "swissprot_canonical_rows": len(rows_by_accession),
        "ncbi_sequence_unique_isoform_rows": len(isoform_rows),
        "ncbi_product_records": len(products),
        "ncbi_unique_gene_sequence_groups": len(groups),
        "ncbi_orphan_gene_sequence_groups_skipped": skipped_orphan_groups,
        "ncbi_accessions_missing_from_fasta": sorted(missing_fasta),
    }
    return rows, audit


def main(argv=None):
    import argparse
    from catalog_io import write_rows

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--swissprot", required=True)
    parser.add_argument("--ncbi-gene-report", required=True)
    parser.add_argument("--ncbi-product-report", required=True)
    parser.add_argument("--ncbi-protein-fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit-json")
    parser.add_argument("--include-predicted", action="store_true")
    parser.add_argument("--include-orphan-refseq", action="store_true")
    args = parser.parse_args(argv)
    rows, audit = build_catalog(
        args.swissprot,
        args.ncbi_gene_report,
        args.ncbi_product_report,
        args.ncbi_protein_fasta,
        include_predicted=args.include_predicted,
        include_orphan_refseq=args.include_orphan_refseq,
    )
    write_rows(rows, args.output, BASE_COLUMNS)
    audit_path = args.audit_json or args.output + ".audit.json"
    with open(audit_path, "w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
