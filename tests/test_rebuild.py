"""Fast offline tests for the clean catalog, PTM projection and assembly."""

import csv
import gzip
import json
import os
import sys
import tempfile
import unittest

PIPELINE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "rbp_pipeline"
)
if PIPELINE_DIR not in sys.path:
    sys.path.insert(0, PIPELINE_DIR)

import assemble_features
import catalog_selection
import export_catalog_fasta
import finalize_table
import generate_feature_reports
import isoform_catalog
import np_refseq_catalog
import opentargets
import ptm
import rebuild_from_scratch
from catalog_io import read_rows, write_feature_rows, write_rows
from rebuild_schema import BASE_COLUMNS, LIST_COLUMNS


SWISSPROT_FIXTURE = """ID   TEST1_HUMAN             Reviewed;           5 AA.
AC   P00001;
DT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.
DT   01-JAN-2000, sequence version 1.
DT   01-JAN-2000, entry version 1.
DE   RecName: Full=Test protein one;
GN   Name=TEST1;
OS   Homo sapiens (Human).
OC   Eukaryota.
OX   NCBI_TaxID=9606;
DR   GeneID; 1;
DR   HGNC; HGNC:1; TEST1.
DR   RefSeq; NP_000001.1; NM_000001.1. [P00001-1]
DR   Ensembl; ENST000001.1; ENSP000001.1; ENSG000001.1. [P00001-1]
SQ   SEQUENCE   5 AA;  500 MW;  0000000000000000 CRC64;
     MAAAA
//
ID   TEST2_HUMAN             Reviewed;           5 AA.
AC   P00002;
DT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.
DT   01-JAN-2000, sequence version 1.
DT   01-JAN-2000, entry version 1.
DE   RecName: Full=Test protein two;
GN   Name=TEST2;
OS   Homo sapiens (Human).
OC   Eukaryota.
OX   NCBI_TaxID=9606;
DR   GeneID; 2;
DR   HGNC; HGNC:2; TEST2.
SQ   SEQUENCE   5 AA;  500 MW;  0000000000000000 CRC64;
     MCCCC
//
"""


class CatalogTests(unittest.TestCase):
    def test_np_ftp_catalog_maps_geneids_and_reports_orphans(self):
        with tempfile.TemporaryDirectory() as root:
            swiss = os.path.join(root, "human.dat")
            protein_dir = os.path.join(root, "protein")
            gff = os.path.join(root, "human.gff.gz")
            os.makedirs(protein_dir)
            with open(swiss, "w", encoding="utf-8") as handle:
                handle.write(SWISSPROT_FIXTURE)
            fasta = os.path.join(protein_dir, "human.1.protein.faa.gz")
            with gzip.open(fasta, "wt", encoding="ascii") as handle:
                handle.write(
                    ">NP_000001.1 canonical product [Homo sapiens]\nMAAAA\n"
                    ">NP_000002.1 alternative product [Homo sapiens]\nMAAAT\n"
                    ">NP_000003.1 orphan product [Homo sapiens]\nMGGGG\n"
                )
            with gzip.open(gff, "wt", encoding="utf-8") as handle:
                handle.write("#!annotation-source fixture-release\n")
                handle.write(
                    "NC_1\tRefSeq\tCDS\t1\t15\t.\t+\t0\t"
                    "ID=cds-NP_000001.1;Parent=rna-NM_000001.1;"
                    "Dbxref=GeneID:1;gene=TEST1;product=canonical product;"
                    "protein_id=NP_000001.1\n"
                )
                handle.write(
                    "NC_1\tRefSeq\tCDS\t20\t34\t.\t+\t0\t"
                    "ID=cds-NP_000002.1;Parent=rna-NM_000002.1;"
                    "Dbxref=GeneID:1;gene=TEST1;product=alternative product;"
                    "protein_id=NP_000002.1\n"
                )
                handle.write(
                    "NC_1\tRefSeq\tCDS\t40\t54\t.\t+\t0\t"
                    "ID=cds-NP_000003.1;Parent=rna-NM_000003.1;"
                    "Dbxref=GeneID:3;gene=ORPHAN;product=orphan product;"
                    "protein_id=NP_000003.1\n"
                )

            rows, audit, reports = np_refseq_catalog.build_catalog(
                swiss, protein_dir, [gff]
            )
            self.assertEqual(audit["current_ncbi_np_accessions"], 3)
            self.assertEqual(audit["np_accessions_mapped_to_reviewed_uniprot"], 2)
            self.assertEqual(audit["np_accessions_without_reviewed_uniprot"], 1)
            self.assertEqual(
                audit["reviewed_uniprot_accessions_without_mapped_np_canonical_fallback"],
                1,
            )
            self.assertEqual(audit["catalog_rows"], 3)
            self.assertEqual(
                reports["np_without_reviewed_uniprot"][0]["refseq_protein"],
                "NP_000003.1",
            )
            fallback = reports["uniprot_without_mapped_np"]
            self.assertEqual([row["uniprot_id"] for row in fallback], ["P00002"])
            isoform = next(row for row in rows if row["row_kind"] == "ncbi_isoform")
            self.assertEqual(isoform["uniprot_id"], "P00001")
            self.assertEqual(isoform["refseq_protein_ids"], ["NP_000002.1"])

    def test_canonical_guarantee_and_sequence_deduplication(self):
        with tempfile.TemporaryDirectory() as root:
            swiss = os.path.join(root, "human.dat")
            gene = os.path.join(root, "data_report.jsonl")
            product = os.path.join(root, "product_report.jsonl")
            fasta = os.path.join(root, "protein.faa")
            with open(swiss, "w", encoding="utf-8") as handle:
                handle.write(SWISSPROT_FIXTURE)
            with open(gene, "w", encoding="utf-8") as handle:
                handle.write(json.dumps({
                    "gene_id": "1", "tax_id": "9606", "symbol": "TEST1",
                    "description": "test", "ensembl_gene_ids": ["ENSG000001"],
                    "swiss_prot_accessions": ["P00001"],
                    "nomenclature_authority": {"identifier": "HGNC:1"},
                    "annotations": [{"annotation_name": "fixture-release"}],
                }) + "\n")
            with open(product, "w", encoding="utf-8") as handle:
                handle.write(json.dumps({
                    "gene_id": "1", "symbol": "TEST1", "description": "test",
                    "transcripts": [
                        {"accession_version": "NM_000001.1", "protein": {
                            "accession_version": "NP_000001.1", "isoform_name": "isoform 1"
                        }},
                        {"accession_version": "NM_000002.1", "protein": {
                            "accession_version": "NP_000002.1", "isoform_name": "isoform 2"
                        }},
                        {"accession_version": "NM_000003.1", "protein": {
                            "accession_version": "NP_000003.1", "isoform_name": "isoform 2"
                        }},
                    ],
                }) + "\n")
            with open(fasta, "w", encoding="utf-8") as handle:
                handle.write(
                    ">NP_000001.1 canonical\nMAAAA\n"
                    ">NP_000002.1 alternative\nMAATA\n"
                    ">NP_000003.1 same alternative\nMAATA\n"
                )

            rows, audit = isoform_catalog.build_catalog(swiss, gene, product, fasta)
            self.assertEqual(audit["swissprot_canonical_rows"], 2)
            self.assertEqual(audit["ncbi_sequence_unique_isoform_rows"], 1)
            self.assertEqual(len(rows), 3)
            canonical = next(row for row in rows if row["protein_key"] == "sp:P00001")
            self.assertIn("NP_000001.1", canonical["refseq_protein_ids"])
            isoform = next(row for row in rows if row["row_kind"] == "ncbi_isoform")
            self.assertEqual(
                isoform["refseq_protein_ids"], ["NP_000002.1", "NP_000003.1"]
            )

    def test_uniprot_selector_retains_canonical_and_mapped_isoforms(self):
        rows = [
            {
                "protein_key": "sp:P00001", "uniprot_id": "P00001",
                "uniprot_parent_ids": ["P00001"], "sequence": "MAAAA",
                "sequence_sha256": isoform_catalog.sequence_sha256("MAAAA"),
            },
            {
                "protein_key": "ncbi:1:abc", "uniprot_id": None,
                "uniprot_parent_ids": ["P00001"], "sequence": "MAATA",
                "sequence_sha256": isoform_catalog.sequence_sha256("MAATA"),
            },
            {
                "protein_key": "sp:P00002", "uniprot_id": "P00002",
                "uniprot_parent_ids": ["P00002"], "sequence": "MCCCC",
                "sequence_sha256": isoform_catalog.sequence_sha256("MCCCC"),
            },
        ]
        selected, audit = catalog_selection.select_rows(rows, proteins="p00001")
        self.assertEqual(
            {row["protein_key"] for row in selected},
            {"sp:P00001", "ncbi:1:abc"},
        )
        self.assertEqual(audit["unmatched_identifiers"], [])

    def test_fasta_selector_is_exact_sequence_specific(self):
        with tempfile.TemporaryDirectory() as root:
            fasta = os.path.join(root, "selected.fasta")
            with open(fasta, "w", encoding="utf-8") as handle:
                handle.write(">wanted\nMAATA\n")
            rows = [
                {
                    "protein_key": "sp:P00001", "sequence": "MAAAA",
                    "sequence_sha256": isoform_catalog.sequence_sha256("MAAAA"),
                },
                {
                    "protein_key": "ncbi:1:abc", "sequence": "MAATA",
                    "sequence_sha256": isoform_catalog.sequence_sha256("MAATA"),
                },
            ]
            selected, audit = catalog_selection.select_rows(
                rows, protein_fasta=fasta, strict=True
            )
            self.assertEqual([row["protein_key"] for row in selected], ["ncbi:1:abc"])
            self.assertEqual(audit["matched_fasta_sequences"], 1)

    def test_feature_dependencies_are_computed_but_not_requested(self):
        requested = rebuild_from_scratch.parse_features("pslab,go_roles")
        self.assertEqual(requested, ["go_roles", "pslab"])
        self.assertEqual(
            rebuild_from_scratch.execution_features(requested),
            ["idr", "go", "go_roles", "pslab"],
        )

    def test_feature_report_writes_markdown_and_figures(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as root:
            sidecar = os.path.join(root, "idr.parquet")
            report_dir = os.path.join(root, "reports")
            write_feature_rows(
                [
                    {"protein_key": "sp:P1", "IDR_count": 0},
                    {"protein_key": "sp:P2", "IDR_count": 2},
                ],
                sidecar,
                ["protein_key", "IDR_count"],
            )
            catalog = pd.DataFrame({"protein_key": ["sp:P1", "sp:P2"]})
            report, summary = generate_feature_reports.write_family_report(
                "idr", catalog, sidecar, report_dir
            )
            self.assertTrue(os.path.exists(report))
            self.assertTrue(
                os.path.exists(os.path.join(report_dir, "figures", "idr_coverage.png"))
            )
            self.assertEqual(summary["duplicate_protein_keys"], 0)
            self.assertEqual(summary["rows_with_nonzero_primary_signal"], 1)

    def test_catalog_fasta_export_uses_refseq_identifier(self):
        with tempfile.TemporaryDirectory() as root:
            catalog = os.path.join(root, "catalog.jsonl")
            output = os.path.join(root, "isoforms.fasta")
            row = {column: [] if column in LIST_COLUMNS else None for column in BASE_COLUMNS}
            row.update(
                protein_key="ncbi:1:abc", row_kind="ncbi_isoform",
                sequence="MAATA", length_aa=5,
                sequence_sha256=isoform_catalog.sequence_sha256("MAATA"),
                refseq_protein_ids=["NP_000002.1"],
                is_swissprot_canonical=False,
            )
            write_rows([row], catalog, BASE_COLUMNS)
            manifest = export_catalog_fasta.export(
                catalog, output, row_kind="isoform", identifier="refseq"
            )
            with open(output, encoding="ascii") as handle:
                text = handle.read()
            self.assertTrue(text.startswith(">NP_000002.1 "))
            self.assertIn("MAATA", text)
            self.assertEqual(manifest["sequence_count"], 1)

    def test_ptm_projection_drops_deleted_site(self):
        source = "MAKST"
        target = "MAST"  # canonical K2 is deleted; S3 is retained at target index 2
        annotation = ptm.empty_annotation()
        annotation["ptm_acetylation"] = 1
        annotation["ptm_acetylation_positions"] = [2]
        annotation["ptm_acetylation_residues"] = ["K"]
        annotation["ptm_phosphorylation"] = 1
        annotation["ptm_phosphorylation_positions"] = [3]
        annotation["ptm_phosphorylation_residues"] = ["S"]
        projected = ptm.project_annotation(annotation, source, target)
        self.assertEqual(projected["ptm_acetylation"], 0)
        self.assertEqual(projected["ptm_phosphorylation_positions"], [2])
        self.assertEqual(projected["ptm_projection_dropped_sites"], 1)

    def test_duckdb_assembly_preserves_catalog_rows(self):
        with tempfile.TemporaryDirectory() as root:
            catalog = os.path.join(root, "catalog.parquet")
            sidecar = os.path.join(root, "feature.parquet")
            output = os.path.join(root, "final.parquet")
            base = {
                column: ([] if column in LIST_COLUMNS else None)
                for column in BASE_COLUMNS
            }
            rows = []
            for number in (1, 2):
                row = dict(base)
                row.update(
                    protein_key=f"sp:P{number}", row_kind="swissprot_canonical",
                    sequence="MA", length_aa=2, sequence_sha256=str(number),
                    uniprot_id=f"P{number}", is_swissprot_canonical=True,
                )
                rows.append(row)
            write_rows(rows, catalog, BASE_COLUMNS)
            write_feature_rows(
                [{"protein_key": "sp:P1", "score": 7}],
                sidecar, ["protein_key", "score"],
            )
            manifest = assemble_features.assemble(catalog, [sidecar], output)
            self.assertEqual(manifest["output_rows"], 2)
            final = read_rows(output)
            self.assertEqual(len(final), 2)

    def test_finalizer_computes_dominance_and_groups_uniprot(self):
        with tempfile.TemporaryDirectory() as root:
            assembled = os.path.join(root, "assembled.parquet")
            output = os.path.join(root, "final.parquet")
            base = {
                column: ([] if column in LIST_COLUMNS else None)
                for column in BASE_COLUMNS
            }
            canonical = dict(base)
            canonical.update(
                protein_key="sp:P00001", row_kind="swissprot_canonical",
                sequence="MAAAA", length_aa=5,
                sequence_sha256=isoform_catalog.sequence_sha256("MAAAA"),
                uniprot_id="P00001", uniprot_parent_ids=["P00001"],
                is_swissprot_canonical=True, gene_symbol="TEST1",
                ensembl_gene_ids=["ENSG1"], ncbi_gene_id="1",
            )
            isoform = dict(base)
            isoform.update(
                protein_key="refseq:1:test", row_kind="ncbi_isoform",
                sequence="MAATA", length_aa=5,
                sequence_sha256=isoform_catalog.sequence_sha256("MAATA"),
                uniprot_id="P00001", uniprot_parent_ids=["P00001"],
                is_swissprot_canonical=False, gene_symbol="TEST1",
                refseq_protein_ids=["NP_000002.1"],
                ensembl_gene_ids=["ENSG1"], ncbi_gene_id="1",
            )
            write_rows([isoform, canonical], assembled, BASE_COLUMNS)
            manifest = finalize_table.finalize(assembled, output)
            rows = read_rows(output)
            self.assertEqual(manifest["dominant_isoform_rows"], 1)
            self.assertEqual([row["dominant_isoform"] for row in rows], [1, 0])
            self.assertEqual([row["uniprot_id"] for row in rows], ["P00001", "P00001"])
            self.assertEqual(rows[1]["ProteinHGVS"], "NP_000002.1")
            self.assertEqual(rows[0]["ENSG"], "ENSG1")

            narrow = os.path.join(root, "narrow.parquet")
            finalize_table.finalize(
                assembled, narrow, columns="uniprot_id,dominant_isoform,ENSG"
            )
            narrow_rows = read_rows(narrow)
            self.assertEqual(
                list(narrow_rows[0]),
                ["protein_key", "uniprot_id", "dominant_isoform", "ENSG"],
            )


class OpenTargetsTests(unittest.TestCase):
    def test_numpy_expression_export_is_normalized(self):
        raw = """[{'efo_code': 'UBERON_1', 'label': 'brain',
        'organs': array(['brain'], dtype=object),
        'anatomical_systems': array(['nervous system'], dtype=object),
        'rna': {'value': 12.5, 'zscore': 2, 'level': 4, 'unit': 'TPM'},
        'protein': {'reliability': True, 'level': 2,
        'cell_type': array([{'name': 'neurons', 'reliability': True, 'level': 2}], dtype=object)}}
        {'efo_code': 'UBERON_2', 'label': 'liver',
        'organs': array(['liver'], dtype=object),
        'anatomical_systems': array(['digestive system'], dtype=object),
        'rna': {'value': 1.0, 'zscore': -1, 'level': 1, 'unit': 'TPM'},
        'protein': {'reliability': False, 'level': -1,
        'cell_type': array([], dtype=object)}}]"""
        tissues = opentargets.parse_tissues(raw)
        self.assertEqual([item["label"] for item in tissues], ["brain", "liver"])
        clean = opentargets._clean_tissue("ENSG1", tissues[0])
        self.assertEqual(clean["rna_value"], 12.5)
        self.assertEqual(clean["protein_cell_types"][0]["name"], "neurons")

    def test_clean_disease_columns_name_scores_honestly(self):
        expression = {"ENSG1": [{
            "ensembl_gene_id": "ENSG1", "tissue_id": "UBERON_1",
            "tissue_name": "brain",
        }]}
        associations = {"ENSG1": [
            {"disease_id": "EFO_1", "datatype_id": "literature",
             "score": 0.3, "evidence_count": 2},
            {"disease_id": "EFO_1", "datatype_id": "genetic_association",
             "score": 0.8, "evidence_count": 1},
        ]}
        diseases = {"EFO_1": {
            "name": "test condition", "description": "fixture",
            "therapeutic_areas": [{"id": "TA_1", "name": "test area"}],
        }}
        result = opentargets.clean_columns_for(
            ["ENSG1"], expression, associations, diseases, release="test"
        )
        record = result["opentargets_disease_associations"][0]
        self.assertEqual(record["disease_name"], "test condition")
        self.assertEqual(record["max_datatype_score"], 0.8)
        self.assertNotIn("score", record)
        self.assertEqual(result["opentargets_disease_count"], 1)
        self.assertEqual(result["opentargets_release"], "test")
        self.assertEqual(set(result), set(opentargets.COLUMNS))

    def test_disk_store_keeps_clean_family_runnable_in_bounded_memory(self):
        with tempfile.TemporaryDirectory() as root:
            expression_path = os.path.join(root, "expression.csv")
            association_path = os.path.join(root, "association.csv")
            cache_path = os.path.join(root, "opentargets.sqlite")
            with open(expression_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["id", "tissues"])
                writer.writeheader()
                writer.writerow({
                    "id": "ENSG1",
                    "tissues": "[{'efo_code':'UBERON_1','label':'brain',"
                    "'organs':array(['brain'],dtype=object),"
                    "'anatomical_systems':array(['nervous'],dtype=object),"
                    "'rna':{'value':2.0,'zscore':1,'level':3,'unit':'TPM'},"
                    "'protein':{'reliability':True,'level':2,"
                    "'cell_type':array([],dtype=object)}}]",
                })
            with open(association_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "diseaseId", "targetId", "datatypeId", "score",
                    "evidenceCount",
                ])
                writer.writeheader()
                writer.writerow({
                    "diseaseId": "EFO_1", "targetId": "ENSG1",
                    "datatypeId": "literature", "score": "0.5",
                    "evidenceCount": "3",
                })
            store = opentargets.build_store(
                expression_path, association_path, {"ENSG1"}, cache_path
            )
            try:
                result = opentargets.clean_columns_from_store(
                    ["ENSG1"], store,
                    {"EFO_1": {"name": "condition", "description": None,
                                "therapeutic_areas": []}},
                )
            finally:
                store.close()
            self.assertTrue(os.path.exists(cache_path))
            self.assertEqual(result["opentargets_expression_tissue_count"], 1)
            self.assertEqual(result["opentargets_disease_count"], 1)


if __name__ == "__main__":
    unittest.main()
