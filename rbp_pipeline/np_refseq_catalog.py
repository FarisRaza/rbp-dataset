"""Build an NP_-only human RefSeq catalog with Swiss-Prot fallbacks.

Every reviewed human Swiss-Prot accession is guaranteed a canonical row.
Current NCBI ``NP_`` products are assigned to reviewed UniProt parents using,
in order, explicit UniProt RefSeq cross-references and shared NCBI GeneIDs.
Sequence-identical products within a gene are aggregated.  Every NP accession
that cannot be assigned to a reviewed UniProt parent is returned in a detailed
audit report rather than silently discarded.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import os
from collections import defaultdict

import isoform_catalog
import ncbi_refseq_ftp
import swissprot_source
from rebuild_schema import BASE_COLUMNS, LIST_COLUMNS, validate_base_row


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
    accepted = {None, meta.accession, f"{meta.accession}-1"}
    refseq = [xref for xref in meta.refseq if xref.get("uniprot_isoform") in accepted]
    ensembl = [xref for xref in meta.ensembl if xref.get("uniprot_isoform") in accepted]
    return refseq, ensembl


def _write_csv(path, rows, columns):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                column: json.dumps(row.get(column), separators=(",", ":"))
                if isinstance(row.get(column), (list, dict)) else row.get(column)
                for column in columns
            })


def build_catalog(
    swissprot_path,
    protein_dir,
    gff_paths,
    accession_product_report=None,
    ensembl_index=None,
    ensembl_release=None,
):
    build_time = dt.datetime.now(dt.timezone.utc).isoformat()
    swiss_manifest = swissprot_source.load_manifest(swissprot_path)
    swiss_release = swiss_manifest.get("uniprot_release")
    metas = [meta for meta, _record in swissprot_source.iter_records(swissprot_path)]
    meta_by_accession = {meta.accession: meta for meta in metas}

    uniprot_by_gene = defaultdict(set)
    uniprot_by_refseq = defaultdict(list)
    ensembl_by_isoform = defaultdict(list)
    for meta in metas:
        for gene_id in meta.gene_ids:
            uniprot_by_gene[gene_id].add(meta.accession)
        for xref in meta.refseq:
            protein = xref.get("protein")
            if protein and protein.startswith("NP_"):
                uniprot_by_refseq[ncbi_refseq_ftp.versionless(protein)].append(
                    (meta.accession, xref.get("uniprot_isoform"))
                )
        for xref in meta.ensembl:
            if xref.get("uniprot_isoform"):
                ensembl_by_isoform[xref["uniprot_isoform"]].append(xref)

    gff_products, release_lines = ncbi_refseq_ftp.load_gff_products(gff_paths)
    gff_only_products = set(gff_products)
    accession_products = ncbi_refseq_ftp.load_accession_product_reports(
        accession_product_report
    )
    for accession, incoming in accession_products.items():
        product = gff_products.setdefault(accession, ncbi_refseq_ftp.GffProduct())
        product.gene_ids.update(incoming.gene_ids)
        product.gene_symbols.update(incoming.gene_symbols)
        product.transcript_ids.update(incoming.transcript_ids)
        product.isoform_names.update(incoming.isoform_names)
        product.assemblies.update(incoming.assemblies)
        product.ensembl_gene_ids.update(incoming.ensembl_gene_ids)
        product.ensembl_transcript_ids.update(incoming.ensembl_transcript_ids)
        product.ensembl_protein_ids.update(incoming.ensembl_protein_ids)
    proteins = list(ncbi_refseq_ftp.iter_np_proteins(protein_dir))

    rows_by_accession = {}
    for meta in metas:
        canonical_refseq, canonical_ensembl = _canonical_xrefs(meta)
        row = _base_row(
            protein_key=f"sp:{meta.accession}",
            row_kind="swissprot_canonical",
            sequence=meta.sequence,
            length_aa=len(meta.sequence),
            sequence_sha256=isoform_catalog.sequence_sha256(meta.sequence),
            sequence_source="UniProtKB/Swiss-Prot canonical",
            tax_id=9606,
            gene_symbol=meta.gene_symbol,
            gene_synonyms=_unique(meta.gene_synonyms),
            gene_description=meta.description,
            ncbi_gene_id=meta.gene_ids[0] if len(meta.gene_ids) == 1 else None,
            ncbi_gene_ids=_unique(meta.gene_ids),
            hgnc_ids=_unique(meta.hgnc_ids),
            uniprot_id=meta.accession,
            uniprot_secondary_accessions=_unique(meta.secondary_accessions),
            uniprot_entry_name=meta.entry_name,
            uniprot_parent_ids=[meta.accession],
            uniprot_isoform_ids=[f"{meta.accession}-1"],
            swissprot_canonical_accessions=[meta.accession],
            is_swissprot_canonical=True,
            refseq_protein_ids=[],
            refseq_transcript_ids=[],
            ncbi_isoform_names=[],
            ensembl_gene_ids=_unique(x.get("gene") for x in canonical_ensembl),
            ensembl_transcript_ids=_unique(
                x.get("transcript") for x in canonical_ensembl
            ),
            ensembl_protein_ids=_unique(x.get("protein") for x in canonical_ensembl),
            identifier_mapping_methods=["UniProtKB_cross_reference"],
            identifier_ambiguity=(
                ["multiple_NCBI_GeneIDs"] if len(meta.gene_ids) > 1 else []
            ),
            canonical_match_method="Swiss-Prot canonical source",
            swissprot_release=swiss_release,
            ncbi_annotation_release=";".join(release_lines) or None,
            ensembl_release=None,
            build_timestamp_utc=build_time,
        )
        rows_by_accession[meta.accession] = row

    mapped_parent_accessions = set()
    orphan_rows = []
    ambiguous_rows = []
    mapped_products = []
    direct_count = 0
    gene_count = 0
    for protein in proteins:
        base = ncbi_refseq_ftp.versionless(protein.accession)
        annotation = gff_products.get(base)
        gene_ids = _unique(annotation.gene_ids if annotation else [])
        explicit_pairs = uniprot_by_refseq.get(base, [])
        explicit_parents = _unique(parent for parent, _tag in explicit_pairs)
        gene_parents = _unique(
            parent for gene_id in gene_ids for parent in uniprot_by_gene.get(gene_id, [])
        )
        parents = explicit_parents or gene_parents
        methods = []
        if explicit_parents:
            methods.append("UniProtKB_RefSeq_cross_reference")
            direct_count += 1
        if gene_parents:
            methods.append("shared_NCBI_GeneID")
            if not explicit_parents:
                gene_count += 1
        isoform_ids = _unique(tag for _parent, tag in explicit_pairs)

        if not parents:
            reason = (
                "NCBI_Gene_has_no_reviewed_UniProt_parent"
                if gene_ids else "no_GeneID_or_UniProt_RefSeq_mapping"
            )
            orphan_rows.append({
                "refseq_protein": protein.accession,
                "length_aa": len(protein.sequence),
                "sequence_sha256": isoform_catalog.sequence_sha256(protein.sequence),
                "description": protein.description,
                "ncbi_gene_ids": gene_ids,
                "gene_symbols": _unique(annotation.gene_symbols if annotation else []),
                "refseq_transcript_ids": _unique(
                    annotation.transcript_ids if annotation else []
                ),
                "ncbi_isoform_names": _unique(
                    annotation.isoform_names if annotation else []
                ),
                "assemblies": _unique(annotation.assemblies if annotation else []),
                "reason": reason,
            })
            continue

        mapped_parent_accessions.update(parents)
        record = {
            "protein": protein,
            "gene_ids": gene_ids,
            "gene_symbols": _unique(annotation.gene_symbols if annotation else []),
            "transcripts": _unique(annotation.transcript_ids if annotation else []),
            "names": _unique(annotation.isoform_names if annotation else []),
            "assemblies": _unique(annotation.assemblies if annotation else []),
            "ensembl_gene_ids": _unique(
                annotation.ensembl_gene_ids if annotation else []
            ),
            "ensembl_transcript_ids": _unique(
                annotation.ensembl_transcript_ids if annotation else []
            ),
            "ensembl_protein_ids": _unique(
                annotation.ensembl_protein_ids if annotation else []
            ),
            "parents": parents,
            "isoform_ids": isoform_ids,
            "methods": methods,
        }
        mapped_products.append(record)
        if len(parents) > 1:
            ambiguous_rows.append({
                "refseq_protein": protein.accession,
                "ncbi_gene_ids": gene_ids,
                "candidate_uniprot_ids": parents,
                "mapping_methods": methods,
                "reason": "multiple_reviewed_UniProt_parent_candidates",
            })

    groups = defaultdict(list)
    for record in mapped_products:
        protein = record["protein"]
        digest = isoform_catalog.sequence_sha256(protein.sequence)
        if len(record["gene_ids"]) == 1:
            identity = "gene:" + record["gene_ids"][0]
        else:
            identity = "parents:" + ",".join(record["parents"])
        groups[(identity, digest)].append(record)

    isoform_rows = []
    exact_canonical_accessions = set()
    for (identity, digest), group in sorted(groups.items()):
        sequence = group[0]["protein"].sequence
        parents = _unique(parent for record in group for parent in record["parents"])
        gene_ids = _unique(value for record in group for value in record["gene_ids"])
        protein_ids = _unique(record["protein"].accession for record in group)
        transcripts = _unique(value for record in group for value in record["transcripts"])
        names = _unique(value for record in group for value in record["names"])
        isoform_ids = _unique(value for record in group for value in record["isoform_ids"])
        methods = _unique(value for record in group for value in record["methods"])
        symbols = _unique(value for record in group for value in record["gene_symbols"])
        product_ensg = _unique(
            value for record in group for value in record["ensembl_gene_ids"]
        )
        product_enst = _unique(
            value for record in group for value in record["ensembl_transcript_ids"]
        )
        product_ensp = _unique(
            value for record in group for value in record["ensembl_protein_ids"]
        )
        exact = [
            parent for parent in parents
            if meta_by_accession[parent].sequence == sequence
        ]
        ensembl = [
            xref for isoform_id in isoform_ids
            for xref in ensembl_by_isoform.get(isoform_id, [])
        ]
        if exact:
            for parent in exact:
                row = rows_by_accession[parent]
                _merge_lists(
                    row,
                    ncbi_gene_ids=gene_ids,
                    refseq_protein_ids=protein_ids,
                    refseq_transcript_ids=transcripts,
                    ncbi_isoform_names=names,
                    uniprot_isoform_ids=isoform_ids,
                    ensembl_gene_ids=(
                        product_ensg + [x.get("gene") for x in ensembl]
                    ),
                    ensembl_transcript_ids=(
                        product_enst + [x.get("transcript") for x in ensembl]
                    ),
                    ensembl_protein_ids=(
                        product_ensp + [x.get("protein") for x in ensembl]
                    ),
                    identifier_mapping_methods=methods + ["exact_sequence_to_Swiss-Prot"],
                )
                row["canonical_match_method"] = "exact NP_/Swiss-Prot sequence"
                exact_canonical_accessions.add(parent)
            continue

        parent_metas = [meta_by_accession[parent] for parent in parents]
        ambiguities = []
        if len(parents) > 1:
            ambiguities.append("multiple_Swiss-Prot_parent_candidates")
        if len(gene_ids) > 1:
            ambiguities.append("multiple_NCBI_GeneIDs")
        symbol = symbols[0] if len(symbols) == 1 else None
        if not symbol:
            parent_symbols = _unique(meta.gene_symbol for meta in parent_metas)
            symbol = parent_symbols[0] if len(parent_symbols) == 1 else None
        description = (
            parent_metas[0].description if len(parent_metas) == 1
            else group[0]["protein"].description
        )
        parent_ensg = _unique(
            xref.get("gene")
            for meta in parent_metas
            for xref in meta.ensembl
        )
        key_identity = gene_ids[0] if len(gene_ids) == 1 else identity.replace(":", "_")
        row = _base_row(
            protein_key=f"refseq:{key_identity}:{digest[:16]}",
            row_kind="ncbi_isoform",
            sequence=sequence,
            length_aa=len(sequence),
            sequence_sha256=digest,
            sequence_source="NCBI RefSeq human NP_ protein release",
            tax_id=9606,
            gene_symbol=symbol,
            gene_synonyms=_unique(value for meta in parent_metas for value in meta.gene_synonyms),
            gene_description=description,
            ncbi_gene_id=gene_ids[0] if len(gene_ids) == 1 else None,
            ncbi_gene_ids=gene_ids,
            hgnc_ids=_unique(value for meta in parent_metas for value in meta.hgnc_ids),
            uniprot_id=parents[0] if len(parents) == 1 else None,
            uniprot_secondary_accessions=[],
            uniprot_entry_name=(parent_metas[0].entry_name if len(parents) == 1 else None),
            uniprot_parent_ids=parents,
            uniprot_isoform_ids=isoform_ids,
            swissprot_canonical_accessions=parents,
            is_swissprot_canonical=False,
            refseq_protein_ids=protein_ids,
            refseq_transcript_ids=transcripts,
            ncbi_isoform_names=names,
            ensembl_gene_ids=_unique(
                product_ensg + parent_ensg + [x.get("gene") for x in ensembl]
            ),
            ensembl_transcript_ids=_unique(
                product_enst + [x.get("transcript") for x in ensembl]
            ),
            ensembl_protein_ids=_unique(
                product_ensp + [x.get("protein") for x in ensembl]
            ),
            identifier_mapping_methods=methods,
            identifier_ambiguity=ambiguities,
            canonical_match_method="not identical to mapped Swiss-Prot canonical",
            swissprot_release=swiss_release,
            ncbi_annotation_release=";".join(release_lines) or None,
            ensembl_release=ensembl_release,
            build_timestamp_utc=build_time,
        )
        isoform_rows.append(row)

    fallback_rows = []
    for meta in metas:
        if meta.accession in mapped_parent_accessions:
            continue
        fallback_rows.append({
            "uniprot_id": meta.accession,
            "uniprot_entry_name": meta.entry_name,
            "gene_symbol": meta.gene_symbol,
            "ncbi_gene_ids": meta.gene_ids,
            "length_aa": len(meta.sequence),
            "sequence_sha256": isoform_catalog.sequence_sha256(meta.sequence),
            "reason": "no_current_NP_product_mapped_to_reviewed_UniProt_entry",
        })

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
                ensembl_gene_ids=[mapping.get("gene") for mapping in mappings],
                ensembl_transcript_ids=[
                    mapping.get("transcript") for mapping in mappings
                ],
                ensembl_protein_ids=[mapping.get("protein") for mapping in mappings],
                identifier_mapping_methods=["exact_sequence_to_Ensembl_peptide_FASTA"],
            )
            row["ensembl_release"] = ensembl_release
            if len(mappings) > 1:
                _merge_lists(
                    row, identifier_ambiguity=["multiple_exact_Ensembl_peptides"]
                )
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
        "scope": "current human curated RefSeq NP_ products + reviewed Swiss-Prot canonical fallbacks",
        "reviewed_human_uniprot_accessions": len(metas),
        "current_ncbi_np_accessions": len(proteins),
        "current_ncbi_np_unique_sequences": len({
            isoform_catalog.sequence_sha256(protein.sequence) for protein in proteins
        }),
        "np_accessions_mapped_to_reviewed_uniprot": len(mapped_products),
        "np_accessions_without_reviewed_uniprot": len(orphan_rows),
        "np_accessions_with_direct_uniprot_refseq_mapping": direct_count,
        "np_accessions_mapped_only_by_shared_geneid": gene_count,
        "np_accessions_with_ambiguous_uniprot_parent": len(ambiguous_rows),
        "reviewed_uniprot_accessions_with_at_least_one_mapped_np": len(mapped_parent_accessions),
        "reviewed_uniprot_accessions_without_mapped_np_canonical_fallback": len(fallback_rows),
        "canonical_rows": len(rows_by_accession),
        "canonical_rows_exactly_matched_to_np": len(exact_canonical_accessions),
        "sequence_unique_noncanonical_np_rows": len(isoform_rows),
        "catalog_rows": len(rows),
        "catalog_rows_with_ensembl_gene": sum(
            bool(row["ensembl_gene_ids"]) for row in rows
        ),
        "catalog_rows_with_ensembl_transcript": sum(
            bool(row["ensembl_transcript_ids"]) for row in rows
        ),
        "catalog_rows_with_ensembl_protein": sum(
            bool(row["ensembl_protein_ids"]) for row in rows
        ),
        "gff_mapped_np_accessions": sum(
            1 for protein in proteins
            if ncbi_refseq_ftp.versionless(protein.accession) in gff_only_products
        ),
        "accession_lookup_mapped_np_accessions": sum(
            1 for protein in proteins
            if ncbi_refseq_ftp.versionless(protein.accession) in accession_products
        ),
        "ncbi_geneid_mapped_np_accessions": sum(
            1 for protein in proteins
            if ncbi_refseq_ftp.versionless(protein.accession) in gff_products
        ),
        "gff_annotation_release": release_lines,
    }
    reports = {
        "np_without_reviewed_uniprot": orphan_rows,
        "uniprot_without_mapped_np": fallback_rows,
        "ambiguous_np_uniprot_mappings": ambiguous_rows,
    }
    return rows, audit, reports


def write_reports(report_dir, reports):
    os.makedirs(report_dir, exist_ok=True)
    paths = {}
    definitions = {
        "np_without_reviewed_uniprot": [
            "refseq_protein", "length_aa", "sequence_sha256", "description",
            "ncbi_gene_ids", "gene_symbols", "refseq_transcript_ids",
            "ncbi_isoform_names", "assemblies", "reason",
        ],
        "uniprot_without_mapped_np": [
            "uniprot_id", "uniprot_entry_name", "gene_symbol", "ncbi_gene_ids",
            "length_aa", "sequence_sha256", "reason",
        ],
        "ambiguous_np_uniprot_mappings": [
            "refseq_protein", "ncbi_gene_ids", "candidate_uniprot_ids",
            "mapping_methods", "reason",
        ],
    }
    for name, columns in definitions.items():
        path = os.path.join(report_dir, name + ".csv")
        _write_csv(path, reports[name], columns)
        paths[name] = path
    return paths
