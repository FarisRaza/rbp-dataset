"""Fast offline tests for the clean catalog, PTM projection and assembly."""

import csv
import json
import os
import tempfile
import unittest

import assemble_features
import isoform_catalog
import opentargets
import ptm
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
