"""Standard runners for every independently reproducible feature family.

Each runner consumes the clean catalog and writes a sidecar keyed by
``protein_key``.  Sidecars never mutate the catalog, which makes a failed or
updated family safe to rerun.  :mod:`assemble_features` performs the final
left joins.

The small ``annotate_<family>.py`` entry points in this repository call this
module with family-specific help.  The scientific calculations remain in the
corresponding family modules (``cider.py``, ``go.py``, ``interpro.py``, etc.).
"""

from __future__ import annotations

import copy
import os

from catalog_io import read_feature_rows, read_rows, write_feature_rows


def _by_key(path):
    if not path or not os.path.exists(path):
        return {}
    return {row["protein_key"]: row for row in read_feature_rows(path)}


def _require(path, label):
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"{label} is required, received: {path!r}")
    return path


def run_cider(rows, idr_features=None, domain_features=None, **_):
    """Whole-sequence, per-IDR and per-UniProt-domain localCIDER metrics."""
    import cider
    import schema

    idrs = _by_key(idr_features)
    domains = _by_key(domain_features)
    for row in rows:
        out = {"protein_key": row["protein_key"]}
        out.update(cider.whole_sequence(row["sequence"]))

        idr_sequences = idrs.get(row["protein_key"], {}).get("IDR_discrete_seq", [])
        for metric, values in cider.per_idr(idr_sequences or []).items():
            out[f"IDR_{metric}"] = values

        domain_record = domains.get(row["protein_key"], {})
        domain_sequences = domain_record.get("Domains_discrete_seq")
        if domain_sequences is None:
            # These metrics depend on UniProt's canonical DOMAIN/ZN_FING
            # coordinates, so they must remain null on sequence-distinct rows.
            out.update({column: None for column in schema.DOMAIN_CIDER_COLUMNS})
        else:
            for metric, values in cider.per_domain(domain_sequences or {}).items():
                out[f"Domains_{metric}"] = values
        out["cider_sequence_sanitization"] = (
            "U->C; other noncanonical residues removed only while scoring"
        )
        yield out

    # references retained to make the returned schema obvious to readers
    _ = schema.CIDER_COLUMNS, schema.IDR_CIDER_COLUMNS, schema.DOMAIN_CIDER_COLUMNS


def columns_cider(**_):
    import schema

    return (
        ["protein_key"]
        + schema.CIDER_COLUMNS
        + schema.IDR_CIDER_COLUMNS
        + schema.DOMAIN_CIDER_COLUMNS
        + ["cider_sequence_sanitization"]
    )


def run_idr(rows, **_):
    """Metapredict V3 IDR/fold geometry using the documented 0.5 threshold."""
    import idr

    idr.install_python_fallback()
    predictions = idr.predict(
        [row["sequence"] for row in rows], batch_size=200, show_progress=False
    )
    for row, prediction in zip(rows, predictions):
        yield {
            "protein_key": row["protein_key"],
            **prediction,
            "idr_method": "metapredict V3; disorder_threshold=0.5",
        }


def columns_idr(**_):
    import schema

    return ["protein_key"] + schema.IDR_GEOMETRY_COLUMNS + ["idr_method"]


def run_domains(rows, swissprot=None, **_):
    """Curated Swiss-Prot DOMAIN/ZN_FING coordinates on canonical rows."""
    import domains
    import schema
    import swissprot_source

    _require(swissprot, "--swissprot")
    wanted = {
        row["uniprot_id"]: row
        for row in rows
        if row.get("row_kind") == "swissprot_canonical" and row.get("uniprot_id")
    }
    annotations = {}
    for meta, record in swissprot_source.iter_records(swissprot):
        if meta.accession not in wanted:
            continue
        geometry = domains.domains_for_record(record)
        annotations[meta.accession] = domains.attach_sequences(geometry, meta.sequence)
    for row in rows:
        if row.get("row_kind") == "swissprot_canonical":
            value = annotations.get(row.get("uniprot_id"))
            scope = "direct Swiss-Prot canonical feature coordinates"
        else:
            value = None
            scope = "not applicable: UniProt coordinates are canonical-specific"
        yield {
            "protein_key": row["protein_key"],
            **(
                copy.deepcopy(value)
                if value is not None
                else (
                    copy.deepcopy(domains.EMPTY)
                    if row.get("row_kind") == "swissprot_canonical"
                    else {column: None for column in schema.DOMAIN_GEOMETRY_COLUMNS}
                )
            ),
            "uniprot_domain_annotation_scope": scope,
        }
    _ = schema.DOMAIN_GEOMETRY_COLUMNS


def columns_domains(**_):
    import schema

    return ["protein_key"] + schema.DOMAIN_GEOMETRY_COLUMNS + [
        "uniprot_domain_annotation_scope"
    ]


def _merge_go(records):
    import go

    out = copy.deepcopy(go.EMPTY)
    for aspect in "CPF":
        triples = set()
        for record in records:
            triples.update(
                zip(
                    record[f"{aspect}_ids"],
                    record[f"{aspect}_descriptions"],
                    record[f"{aspect}_evidence"],
                )
            )
        ordered = sorted(triples, key=lambda item: tuple(str(x) for x in item))
        out[f"{aspect}_ids"] = [item[0] for item in ordered]
        out[f"{aspect}_descriptions"] = [item[1] for item in ordered]
        out[f"{aspect}_evidence"] = [item[2] for item in ordered]
    return out


def run_go(rows, swissprot=None, **_):
    """UniProt GO cross-references on canonical Swiss-Prot rows only."""
    import go
    import swissprot_source

    _require(swissprot, "--swissprot")
    by_accession = {
        meta.accession: go.go_for_record(record)
        for meta, record in swissprot_source.iter_records(swissprot)
    }
    for row in rows:
        canonical = row.get("row_kind") == "swissprot_canonical"
        parents = [row.get("uniprot_id")] if canonical else []
        scope = (
            "direct UniProtKB canonical protein annotation"
            if canonical else
            "not applicable: UniProt GO annotation is canonical-specific"
        )
        records = [by_accession[parent] for parent in parents if parent in by_accession]
        yield {
            "protein_key": row["protein_key"],
            **(
                _merge_go(records)
                if records else (
                    copy.deepcopy(go.EMPTY)
                    if canonical else {column: None for column in go.EMPTY}
                )
            ),
            "go_annotation_scope": scope,
            "go_source_uniprot_ids": [p for p in parents if p in by_accession],
        }


def columns_go(**_):
    import schema

    return ["protein_key"] + schema.GO_COLUMNS + [
        "go_annotation_scope", "go_source_uniprot_ids"
    ]


def run_eclip(rows, eclip_table=None, **_):
    """ENCODE/ENCORI/POSTAR/Skipper CLIP summaries joined by gene symbol."""
    import eclip

    _require(eclip_table, "--eclip-table")
    by_gene, columns = eclip.load(eclip_table)
    for row in rows:
        symbol = row.get("gene_symbol")
        record = by_gene.get(symbol)
        measured = False
        if record is not None:
            for column, value in record.items():
                if not column.startswith("has_") and column != "n_eclip_sources":
                    continue
                if str(value).strip().lower() not in {"", "0", "0.0", "false", "none", "nan"}:
                    measured = True
                    break
        yield {
            "protein_key": row["protein_key"],
            **eclip.columns_for(symbol, by_gene, columns),
            "eclip_annotation_scope": (
                "gene-level CLIP measurement broadcast to protein isoforms"
                if measured else (
                    "RBP census member not measured in the supplied CLIP compilation"
                    if record is not None else
                    "gene absent from the supplied RBP/CLIP compilation"
                )
            ),
        }


def columns_eclip(eclip_table=None, **_):
    import eclip

    _require(eclip_table, "--eclip-table")
    _records, columns = eclip.load(eclip_table)
    return ["protein_key"] + columns + ["eclip_annotation_scope"]


def run_interpro(rows, interpro_tsv=None, **_):
    """InterProScan domains joined to every matching RefSeq protein accession."""
    import interpro

    _require(interpro_tsv, "--interpro-tsv")
    index = interpro.index_by_accession(interpro_tsv)
    versionless = {}
    for accession, record in index.items():
        versionless.setdefault(accession.split(".")[0], record)
    for row in rows:
        record = None
        source = None
        for accession in row.get("refseq_protein_ids") or []:
            record = index.get(accession) or versionless.get(accession.split(".")[0])
            if record is not None:
                source = accession
                break
        yield {
            "protein_key": row["protein_key"],
            **(copy.deepcopy(record) if record is not None else copy.deepcopy(interpro.EMPTY)),
            "interpro_source_refseq_protein": source,
            "interpro_version": interpro.VERSION,
            "interpro_annotation_scope": (
                "direct InterProScan hits for a matching RefSeq protein accession"
                if record is not None else
                "no hit or accession absent from supplied InterProScan TSV; "
                "the hit-only TSV cannot distinguish these cases"
            ),
        }


def columns_interpro(**_):
    import interpro

    return ["protein_key"] + interpro.COLUMNS + [
        "interpro_source_refseq_protein", "interpro_version",
        "interpro_annotation_scope",
    ]


def run_ptm(rows, ptm_csv=None, **_):
    """Canonical UniProt PTM sites; sequence-distinct isoforms remain null."""
    import ptm

    _require(ptm_csv, "--ptm-csv")
    source = ptm.load_wide_csv(ptm_csv)
    yield from ptm.annotate_rows(rows, source, project_isoforms=False)


def columns_ptm(**_):
    import ptm

    return ["protein_key"] + ptm.PTM_COLUMNS + ptm.PROVENANCE_COLUMNS + [
        "ptm_annotation_scope"
    ]


def run_opentargets(
    rows,
    opentargets_associations=None,
    opentargets_expression=None,
    opentargets_diseases=None,
    opentargets_release=None,
    opentargets_expression_long_output=None,
    opentargets_disease_long_output=None,
    opentargets_cache=None,
    **_,
):
    """Clean Open Targets tissue expression and disease records keyed by ENSG."""
    import opentargets

    _require(opentargets_associations, "--opentargets-associations")
    _require(opentargets_expression, "--opentargets-expression")
    _require(opentargets_diseases, "--opentargets-diseases")
    release = opentargets_release or opentargets.DEFAULT_RELEASE
    wanted = {
        normalized
        for row in rows
        for identifier in row.get("ensembl_gene_ids") or []
        if (normalized := opentargets.normalize_ensg(identifier))
    }
    diseases = opentargets.load_disease_metadata(opentargets_diseases)
    store = opentargets.build_store(
        opentargets_expression,
        opentargets_associations,
        wanted,
        destination=opentargets_cache,
    )
    try:
        if opentargets_expression_long_output:
            write_feature_rows(
                opentargets.expression_long_rows_from_store(store, release),
                opentargets_expression_long_output,
                opentargets.EXPRESSION_LONG_COLUMNS,
            )
        if opentargets_disease_long_output:
            write_feature_rows(
                opentargets.disease_long_rows_from_store(store, diseases, release),
                opentargets_disease_long_output,
                opentargets.DISEASE_LONG_COLUMNS,
            )
        for row in rows:
            ensg = [
                normalized
                for identifier in row.get("ensembl_gene_ids") or []
                if (normalized := opentargets.normalize_ensg(identifier))
            ]
            yield {
                "protein_key": row["protein_key"],
                **opentargets.clean_columns_from_store(
                    ensg, store, diseases, release
                ),
            }
    finally:
        store.close()


def columns_opentargets(**_):
    import opentargets

    return ["protein_key"] + opentargets.COLUMNS


def run_cdcode(rows, cdcode_root=None, **_):
    """CD-CODE condensate membership via canonical UniProt parent accessions."""
    import cdcode

    _require(cdcode_root, "--cdcode-root")
    problems = cdcode.verify_alignment(cdcode_root)
    if problems:
        raise ValueError(f"CD-CODE source/member ordering failed validation: {problems[:3]}")
    index, attributes = cdcode.build_index(cdcode_root)
    for row in rows:
        canonical = row.get("row_kind") == "swissprot_canonical"
        parents = [row.get("uniprot_id")] if canonical else []
        records = [cdcode.lookup(parent, index, attributes) for parent in parents if parent]
        merged = copy.deepcopy(cdcode.EMPTY)
        seen_uid = set()
        for record in records:
            for item_index, uid in enumerate(record["UID"]):
                if uid in seen_uid:
                    continue
                seen_uid.add(uid)
                for column in merged:
                    merged[column].append(record[column][item_index])
        yield {
            "protein_key": row["protein_key"],
            **(
                merged if canonical
                else {column: None for column in cdcode.EMPTY}
            ),
            "cdcode_annotation_scope": (
                "direct canonical UniProt membership"
                if canonical
                else "not applicable: CD-CODE membership is canonical-specific"
            ),
        }


def columns_cdcode(**_):
    import schema

    return ["protein_key"] + schema.CDCODE_COLUMNS + ["cdcode_annotation_scope"]


def run_string(rows, string_links=None, **_):
    """STRING v12 interaction partners for every ENSP represented in the catalog."""
    import string_ppi

    _require(string_links, "--string-links")
    query_by_key = {
        row["protein_key"]: sorted(
            {string_ppi.clean_ensp(x) for x in row.get("ensembl_protein_ids") or [] if x}
        )
        for row in rows
    }
    focus = {ensp for values in query_by_key.values() for ensp in values}
    adjacency, _reverse = string_ppi.scan(string_links, focus)
    ensp_to_uniprot = {}
    for row in rows:
        parents = row.get("uniprot_parent_ids") or []
        if len(parents) != 1:
            continue
        for ensp in query_by_key[row["protein_key"]]:
            ensp_to_uniprot.setdefault(ensp, parents[0])
    for row in rows:
        queries = query_by_key[row["protein_key"]]
        by_query = {query: adjacency.get(query, {}) for query in queries}
        translated = {
            query: {
                ensp_to_uniprot[partner]: score
                for partner, score in partners.items()
                if partner in ensp_to_uniprot
            }
            for query, partners in by_query.items()
        }
        in_catalog = {
            query: {
                partner: score
                for partner, score in partners.items()
                if partner in ensp_to_uniprot
            }
            for query, partners in by_query.items()
        }
        yield {
            "protein_key": row["protein_key"],
            "string_query_ensp_ids": queries,
            "string_partners_ensp_by_query": by_query,
            "string_partners_uniprot_by_query": translated,
            "string_partners_ensp_in_catalog_by_query": in_catalog,
            "string_version": "12.0",
        }


def columns_string(**_):
    return [
        "protein_key",
        "string_query_ensp_ids",
        "string_partners_ensp_by_query",
        "string_partners_uniprot_by_query",
        "string_partners_ensp_in_catalog_by_query",
        "string_version",
    ]


def run_pslab(rows, idr_features=None, pspred_repo=None, **_):
    """PSLab phase-separation features and predictions for each predicted IDR."""
    import pslab

    _require(idr_features, "--idr-features")
    _require(pspred_repo, "--pspred-repo")
    idrs = _by_key(idr_features)
    models, residues, nu_file = pslab.load_models(pspred_repo)
    for row in rows:
        sequences = idrs.get(row["protein_key"], {}).get("IDR_discrete_seq", []) or []
        predictions = list(pslab.predict(sequences, models, residues, nu_file))
        output = {column: [] for column in columns_pslab()[1:-1]}
        for prediction in predictions:
            for column in output:
                output[column].append(prediction.get(column))
        yield {
            "protein_key": row["protein_key"],
            **output,
            "pslab_annotation_scope": "one prediction per metapredict IDR, same order",
        }


def columns_pslab(**_):
    import schema

    return ["protein_key"] + schema.PSLAB_COLUMNS + ["pslab_annotation_scope"]


def run_go_roles(rows, go_features=None, **_):
    """Derived transcription/translation/mRNA-stability flags from GO names."""
    import go_roles

    _require(go_features, "--go-features")
    by_key = _by_key(go_features)
    for row in rows:
        canonical = row.get("row_kind") == "swissprot_canonical"
        yield {
            "protein_key": row["protein_key"],
            **(
                go_roles.flags_for(by_key.get(row["protein_key"], {}))
                if canonical
                else {column: None for column in go_roles.ROLE_COLUMNS}
            ),
            "go_role_annotation_scope": (
                "derived from canonical UniProt GO terms"
                if canonical else
                "not applicable: source GO annotation is canonical-specific"
            ),
        }


def columns_go_roles(**_):
    import go_roles

    return ["protein_key"] + go_roles.ROLE_COLUMNS + ["go_role_annotation_scope"]


def run_rcsb(rows, rcsb_summary=None, **_):
    """Join compact experimental PDB/secondary-structure summaries to canonicals."""
    import rcsb

    _require(rcsb_summary, "--rcsb-summary")
    index = rcsb.load_summary(rcsb_summary)
    for row in rows:
        accession = (
            row.get("uniprot_id")
            if row.get("row_kind") == "swissprot_canonical"
            else None
        )
        record = index.get(accession) if accession else None
        yield {
            "protein_key": row["protein_key"],
            **(
                copy.deepcopy(record)
                if record is not None
                else {column: None for column in rcsb.SUMMARY_COLUMNS}
            ),
            "rcsb_annotation_scope": (
                "direct SIFTS mapping for Swiss-Prot canonical sequence"
                if record is not None else (
                    "canonical accession absent from supplied RCSB/SIFTS summary"
                    if accession else
                    "not inherited to sequence-distinct NCBI isoforms"
                )
            ),
        }


def columns_rcsb(**_):
    import rcsb

    return ["protein_key"] + rcsb.SUMMARY_COLUMNS + ["rcsb_annotation_scope"]


RUNNERS = {
    "cider": run_cider,
    "idr": run_idr,
    "domains": run_domains,
    "go": run_go,
    "eclip": run_eclip,
    "interpro": run_interpro,
    "ptm": run_ptm,
    "opentargets": run_opentargets,
    "cdcode": run_cdcode,
    "string": run_string,
    "pslab": run_pslab,
    "go_roles": run_go_roles,
    "rcsb": run_rcsb,
}

COLUMN_FUNCTIONS = {
    name: globals()[f"columns_{name}"] for name in RUNNERS
}


def run_family(family, input_path, output_path, **kwargs):
    """Run one named family and write its keyed sidecar."""
    if family not in RUNNERS:
        raise ValueError(f"unknown family {family!r}; choose from {sorted(RUNNERS)}")
    rows = read_rows(input_path)
    columns = COLUMN_FUNCTIONS[family](**kwargs)
    write_feature_rows(RUNNERS[family](rows, **kwargs), output_path, columns)
    return output_path


def family_main(family, argv=None):
    import argparse

    descriptions = {name: runner.__doc__ for name, runner in RUNNERS.items()}
    parser = argparse.ArgumentParser(description=descriptions[family])
    parser.add_argument("--input", required=True, help="clean catalog")
    parser.add_argument("--output", required=True, help="feature sidecar CSV/Parquet")
    parser.add_argument("--swissprot")
    parser.add_argument("--idr-features")
    parser.add_argument("--domain-features")
    parser.add_argument("--eclip-table")
    parser.add_argument("--interpro-tsv")
    parser.add_argument("--ptm-csv")
    parser.add_argument("--opentargets-associations")
    parser.add_argument("--opentargets-expression")
    parser.add_argument("--opentargets-diseases")
    parser.add_argument("--opentargets-release")
    parser.add_argument("--opentargets-expression-long-output")
    parser.add_argument("--opentargets-disease-long-output")
    parser.add_argument("--opentargets-cache")
    parser.add_argument("--cdcode-root")
    parser.add_argument("--string-links")
    parser.add_argument("--pspred-repo")
    parser.add_argument("--go-features")
    parser.add_argument("--rcsb-summary")
    args = parser.parse_args(argv)
    options = vars(args)
    input_path = options.pop("input")
    output_path = options.pop("output")
    run_family(family, input_path, output_path, **options)
    print(f"wrote {family} features -> {output_path}")


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=sorted(RUNNERS))
    args, remainder = parser.parse_known_args(argv)
    family_main(args.family, remainder)


if __name__ == "__main__":
    main()
