"""Validate RCSB/SIFTS code and, optionally, every generated output row.

    python validate_rcsb.py
    python validate_rcsb.py --full
"""

import argparse
import ast
from collections import Counter
import csv
import gzip
import json
import os
import tempfile

import append_rcsb_to_extended
import rcsb


def test_coordinate_mapping():
    segments = (
        rcsb.MappingSegment(10, 19, 101, 110),
        rcsb.MappingSegment(25, 30, 201, 206),
    )
    assert rcsb.map_element_to_uniprot(12, 17, segments) == (
        [[103, 108]], 6, "complete"
    )
    assert rcsb.map_element_to_uniprot(17, 27, segments) == (
        [[108, 110], [201, 203]], 6, "partial"
    )
    assert rcsb.map_element_to_uniprot(1, 5, segments) == (
        [], 0, "outside_uniprot_segment"
    )


def test_length_mismatch_is_not_fabricated():
    segments = (rcsb.MappingSegment(1, 10, 50, 58),)
    assert rcsb.map_element_to_uniprot(2, 5, segments) == (
        [], 0, "unmapped_length_mismatch"
    )


def test_release_and_known_mapping(sifts_path):
    index = rcsb.load_sifts(sifts_path, ("P04637",))
    # TP53 has many experimental structures; 1TUP is the canonical DNA-binding
    # domain structure and is stable enough for a mapping smoke test.
    assert "1TUP" in index.by_accession["P04637"]["pdb_ids"]
    assert index.by_chain[("1TUP", "A")]["P04637"]


def test_summary_appender():
    with tempfile.TemporaryDirectory() as directory:
        input_path = os.path.join(directory, "input.csv")
        summary_path = os.path.join(directory, "summary.csv")
        output_path = os.path.join(directory, "output.csv")
        with open(input_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["uniprot_id", "value"])
            writer.writerow(["P04637", "known"])
            writer.writerow(["NOT_IN_SUMMARY", "unknown"])
        with open(summary_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["uniprot_id"] + rcsb.SUMMARY_COLUMNS
            )
            writer.writeheader()
            row = {"uniprot_id": "P04637", **rcsb.empty_summary()}
            row["RCSB_PDB_IDs"] = "['1TUP']"
            row["RCSB_PDB_count"] = "1"
            writer.writerow(row)

        rows, matched, columns = append_rcsb_to_extended.append(
            input_path, summary_path, output_path, progress_every=0
        )
        assert (rows, matched, columns) == (2, 1, 2 + len(rcsb.SUMMARY_COLUMNS))
        with open(output_path, newline="", encoding="utf-8") as fh:
            records = list(csv.DictReader(fh))
        assert records[0]["RCSB_PDB_IDs"] == "['1TUP']"
        assert records[1]["RCSB_PDB_IDs"] == "[]"


def validate_outputs(summary_path, elements_path, metadata_path):
    """Stream both full outputs and reconcile every count and coordinate."""
    with open(metadata_path, encoding="utf-8") as fh:
        metadata = json.load(fh)

    summary_rows = 0
    accessions = set()
    observation_sum = 0
    proteins_with_pdb = 0
    tp53_pdbs = None
    with open(summary_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        expected = ["uniprot_id"] + rcsb.SUMMARY_COLUMNS
        assert reader.fieldnames == expected, (reader.fieldnames, expected)
        for row in reader:
            summary_rows += 1
            accession = row["uniprot_id"]
            assert accession and accession not in accessions
            accessions.add(accession)
            pdb_ids = ast.literal_eval(row["RCSB_PDB_IDs"])
            chains = ast.literal_eval(row["RCSB_PDB_chains"])
            assert isinstance(pdb_ids, list) and isinstance(chains, dict)
            assert int(row["RCSB_PDB_count"]) == len(pdb_ids)
            assert sorted(pdb_ids) == pdb_ids
            if pdb_ids:
                proteins_with_pdb += 1
            all_count = int(row["RCSB_secondary_structure_observation_count"])
            helix = int(row["RCSB_helix_observation_count"])
            strand = int(row["RCSB_beta_strand_observation_count"])
            complete = int(row["RCSB_secondary_structure_complete_mapping_count"])
            partial = int(row["RCSB_secondary_structure_partial_mapping_count"])
            assert all_count == helix + strand
            assert all_count == complete + partial
            observation_sum += all_count
            if accession == "P04637":
                tp53_pdbs = pdb_ids

    assert tp53_pdbs is not None and "1TUP" in tp53_pdbs

    element_rows = 0
    types = Counter()
    statuses = Counter()
    saw_tp53_1tup = False
    with gzip.open(elements_path, "rt", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        assert reader.fieldnames == rcsb.ELEMENT_COLUMNS
        for row in reader:
            element_rows += 1
            start = int(row["pdb_seq_start"])
            end = int(row["pdb_seq_end"])
            length = int(row["element_length"])
            mapped = int(row["uniprot_mapped_residue_count"])
            assert start > 0 and end >= start and length == end - start + 1
            assert 0 <= mapped <= length
            ranges = json.loads(row["uniprot_ranges"])
            assert all(a > 0 and b >= a for a, b in ranges)
            element_type = row["element_type"]
            status = row["uniprot_mapping_status"]
            assert element_type in ("helix", "beta_strand")
            assert status in ("complete", "partial", "unmapped_length_mismatch")
            if status == "complete":
                assert mapped == length
            types[element_type] += 1
            statuses[status] += 1
            if row["uniprot_id"] == "P04637" and row["pdb_id"] == "1TUP":
                saw_tp53_1tup = True

    expected_counts = metadata["counts"]
    assert summary_rows == expected_counts["human_swissprot_accessions"]
    assert proteins_with_pdb == expected_counts["proteins_with_pdb"]
    assert element_rows == expected_counts["secondary_structure_element_rows"]
    assert observation_sum == element_rows
    assert saw_tp53_1tup
    return {
        "summary_rows": summary_rows,
        "unique_accessions": len(accessions),
        "proteins_with_pdb": proteins_with_pdb,
        "element_rows": element_rows,
        "element_types": dict(types),
        "mapping_statuses": dict(statuses),
        "tp53_1tup_present": saw_tp53_1tup,
    }


def validate_joined_table(joined_path, summary_path):
    """Verify every appended cell in the 252-column extended table."""
    by_accession = rcsb.load_summary(summary_path)
    empty = rcsb.empty_summary()
    rows = matched = 0
    with open(joined_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        assert header[-len(rcsb.SUMMARY_COLUMNS):] == rcsb.SUMMARY_COLUMNS
        assert len(header) == len(set(header)), "joined table has duplicate columns"
        accession_index = header.index("uniprot_id")
        extension_start = len(header) - len(rcsb.SUMMARY_COLUMNS)
        for row in reader:
            rows += 1
            assert len(row) == len(header), (rows + 1, len(row), len(header))
            expected = by_accession.get(row[accession_index])
            if expected is None:
                expected = empty
            else:
                matched += 1
            actual = row[extension_start:]
            wanted = [expected[column] for column in rcsb.SUMMARY_COLUMNS]
            assert actual == wanted, (rows + 1, row[accession_index])
    return {"rows": rows, "columns": len(header), "matched_rows": matched}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--full", action="store_true", help="stream and verify all outputs")
    p.add_argument("--summary")
    p.add_argument("--elements")
    p.add_argument("--metadata")
    p.add_argument("--joined", action="store_true", help="verify every joined table row")
    p.add_argument("--joined-path")
    return p


def main(argv=None):
    import paths

    args = parser().parse_args(argv)
    test_coordinate_mapping()
    test_length_mismatch_is_not_fabricated()
    test_release_and_known_mapping(paths.RCSB_SIFTS)
    test_summary_appender()
    print("RCSB code validation passed")
    if args.full:
        report = validate_outputs(
            args.summary or paths.RCSB_SUMMARY,
            args.elements or paths.RCSB_ELEMENTS,
            args.metadata or paths.RCSB_METADATA,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
    if args.joined:
        report = validate_joined_table(
            args.joined_path or paths.RCSB_ENRICHED_TABLE,
            args.summary or paths.RCSB_SUMMARY,
        )
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
