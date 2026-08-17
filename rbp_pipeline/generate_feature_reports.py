"""Generate one Markdown quality-control report per feature family.

The reports are intentionally derived from keyed sidecars rather than from the
very wide final table. This keeps report generation fast and lets a user inspect
one family even when they intentionally requested only a subset of features.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import statistics

import pandas as pd

import rebuild_from_scratch


PRIMARY_COLUMNS = {
    "cider": ["NCPR"],
    "idr": ["IDR_count"],
    "domains": ["Domains_count"],
    "go": ["C_ids", "P_ids", "F_ids"],
    "eclip": [],  # source-specific column names are discovered dynamically
    "interpro": ["InterPro_n_hits"],
    "ptm": [],  # all ptm_<type> indicator columns are summed
    "opentargets": [
        "opentargets_expression_tissue_count",
        "opentargets_disease_count",
    ],
    "cdcode": ["UID"],
    "string": ["string_partners_ensp_by_query"],
    "go_roles": [
        "role_in_transcription",
        "role_in_translation",
        "role_in_mrna_stability",
    ],
    "pslab": ["Delta G [kT]"],
    "rcsb": ["RCSB_PDB_count"],
}

PROVENANCE_WORDS = ("scope", "method", "version", "source", "coordinate")


def _read_table(path):
    if path.lower().endswith(".parquet"):
        return pd.read_parquet(path)
    if path.lower().endswith(".csv"):
        return pd.read_csv(path, low_memory=False)
    if path.lower().endswith(".jsonl"):
        return pd.read_json(path, lines=True)
    raise ValueError(f"unsupported report input: {path}")


def _decode(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        if stripped.startswith(("[", "{")):
            try:
                return json.loads(stripped)
            except json.JSONDecodeError:
                return value
    return value


def _as_list(value):
    value = _decode(value)
    if value is None:
        return []
    if isinstance(value, tuple):
        return list(value)
    return value if isinstance(value, list) else []


def _as_dict(value):
    value = _decode(value)
    return value if isinstance(value, dict) else {}


def _number(value):
    value = _decode(value)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _nonnull(value):
    return _decode(value) is not None


def _add_check(checks, name, failures, checked, expectation):
    checks.append({
        "check": name,
        "status": "PASS" if failures == 0 else "FAIL",
        "failed_rows": int(failures),
        "checked_rows": int(checked),
        "expectation": expectation,
    })


def _parallel_list_failures(frame, columns):
    failures = 0
    checked = 0
    available = [column for column in columns if column in frame]
    if len(available) < 2:
        return failures, checked
    for record in frame[available].itertuples(index=False, name=None):
        decoded = [_decode(value) for value in record]
        present = [value for value in decoded if value is not None]
        if not present:
            continue
        checked += 1
        if any(not isinstance(value, list) for value in present):
            failures += 1
            continue
        if len({len(value) for value in present}) != 1:
            failures += 1
    return failures, checked


def _canonical_leakage(frame, catalog_by_key, columns):
    failures = 0
    checked = 0
    available = [column for column in columns if column in frame]
    for record in frame[["protein_key"] + available].to_dict("records"):
        meta = catalog_by_key.get(str(record["protein_key"]), {})
        if meta.get("row_kind") == "swissprot_canonical":
            continue
        checked += 1
        if any(_nonnull(record.get(column)) for column in available):
            failures += 1
    return failures, checked


def _family_sanity_checks(family, catalog, frame):
    """Return machine-readable hard checks for one feature sidecar."""
    catalog_keys = set(catalog["protein_key"].astype(str))
    sidecar_keys = set(frame["protein_key"].astype(str))
    catalog_by_key = {
        str(record["protein_key"]): record
        for record in catalog.to_dict("records")
    }
    checks = []
    _add_check(
        checks, "protein_key is unique", int(frame["protein_key"].duplicated().sum()),
        len(frame), "zero duplicate sidecar keys",
    )
    _add_check(
        checks, "all catalog rows are represented",
        len(catalog_keys.difference(sidecar_keys)), len(catalog_keys),
        "zero catalog keys missing from the sidecar",
    )
    _add_check(
        checks, "sidecar contains no foreign rows",
        len(sidecar_keys.difference(catalog_keys)), len(sidecar_keys),
        "zero sidecar keys outside the catalog",
    )

    if family == "idr":
        geometry_failures = 0
        coordinate_failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            meta = catalog_by_key.get(str(record["protein_key"]), {})
            length = _number(meta.get("length_aa"))
            if length is None:
                continue
            checked += 1
            idr_count = int(_number(record.get("IDR_count")) or 0)
            fold_count = int(_number(record.get("FOLD_count")) or 0)
            idr_ranges = _as_list(record.get("IDR_range"))
            fold_ranges = _as_list(record.get("FOLD_range"))
            idr_sequences = _as_list(record.get("IDR_discrete_seq"))
            fold_sequences = _as_list(record.get("FOLD_discrete_seq"))
            idr_total = _number(record.get("IDR_total_size"))
            fold_total = _number(record.get("FOLD_total_size"))
            if (
                len(idr_ranges) != idr_count
                or len(idr_sequences) != idr_count
                or len(fold_ranges) != fold_count
                or len(fold_sequences) != fold_count
                or idr_total is None or fold_total is None
                or int(idr_total + fold_total) != int(length)
            ):
                geometry_failures += 1
            ranges = idr_ranges + fold_ranges
            valid = all(
                isinstance(item, list) and len(item) == 2
                and _number(item[0]) is not None and _number(item[1]) is not None
                and 0 <= _number(item[0]) <= _number(item[1]) <= length
                for item in ranges
            )
            if not valid:
                coordinate_failures += 1
        _add_check(
            checks, "IDR/fold counts and lengths agree", geometry_failures, checked,
            "counts match list lengths and IDR_total_size + FOLD_total_size equals sequence length",
        )
        _add_check(
            checks, "IDR/fold coordinates are in bounds", coordinate_failures, checked,
            "all zero-based half-open ranges lie within the sequence",
        )

    elif family == "cider":
        bounded = {
            "FCR": (0.0, 1.0),
            "NCPR": (-1.0, 1.0),
            "fraction_negative": (0.0, 1.0),
            "fraction_positive": (0.0, 1.0),
            "fraction_expanding": (0.0, 1.0),
            "fraction_disorder_promoting": (0.0, 1.0),
        }
        failures = 0
        checked = 0
        for column, (lower, upper) in bounded.items():
            if column not in frame:
                continue
            for value in frame[column]:
                number = _number(value)
                if number is None:
                    continue
                checked += 1
                failures += int(not lower <= number <= upper)
        _add_check(
            checks, "whole-sequence CIDER fractions are bounded", failures, checked,
            "fractions are in [0,1] and NCPR is in [-1,1]",
        )
        idr_columns = [column for column in frame if column.startswith("IDR_")]
        failures, checked = _parallel_list_failures(frame, idr_columns)
        _add_check(
            checks, "per-IDR CIDER lists are parallel", failures, checked,
            "all populated IDR metric lists have equal length within a row",
        )
        domain_columns = [column for column in frame if column.startswith("Domains_")]
        failures = 0
        checked = 0
        for record in frame[domain_columns].to_dict("records"):
            decoded = [_decode(record[column]) for column in domain_columns]
            present = [value for value in decoded if value is not None]
            if not present:
                continue
            checked += 1
            if any(not isinstance(value, dict) for value in present):
                failures += 1
                continue
            key_sets = [set(value) for value in present]
            if len({tuple(sorted(keys)) for keys in key_sets}) != 1:
                failures += 1
                continue
            for domain in key_sets[0]:
                values = [value[domain] for value in present]
                if any(not isinstance(value, list) for value in values) \
                        or len({len(value) for value in values}) != 1:
                    failures += 1
                    break
        _add_check(
            checks, "per-domain CIDER dictionaries are parallel", failures, checked,
            "domain keys and per-region metric-list lengths agree within a row",
        )

    elif family == "domains":
        geometry = [
            "Domains", "Domains_count", "Domains_avg_size", "Domains_total_size",
            "Domains_range", "Domains_discrete_seq", "Domains_concat_seq",
        ]
        failures, checked = _canonical_leakage(frame, catalog_by_key, geometry)
        _add_check(
            checks, "UniProt domain coordinates stay canonical-only", failures, checked,
            "all sequence-distinct isoform domain fields are null",
        )
        failures = 0
        coordinate_failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            meta = catalog_by_key.get(str(record["protein_key"]), {})
            if meta.get("row_kind") != "swissprot_canonical":
                continue
            counts = _as_dict(record.get("Domains_count"))
            ranges = _as_dict(record.get("Domains_range"))
            sequences = _as_dict(record.get("Domains_discrete_seq"))
            length = _number(meta.get("length_aa"))
            checked += 1
            keys = set(counts) | set(ranges) | set(sequences)
            for name in keys:
                named_ranges = _as_list(ranges.get(name))
                named_sequences = _as_list(sequences.get(name))
                if int(_number(counts.get(name)) or 0) != len(named_ranges) \
                        or len(named_ranges) != len(named_sequences):
                    failures += 1
                    break
            if length is not None:
                valid = all(
                    isinstance(item, list) and len(item) == 2
                    and _number(item[0]) is not None and _number(item[1]) is not None
                    and 0 <= _number(item[0]) <= _number(item[1]) <= length
                    and len(str(sequence)) == int(_number(item[1]) - _number(item[0]))
                    for name in keys
                    for item, sequence in zip(
                        _as_list(ranges.get(name)), _as_list(sequences.get(name))
                    )
                )
                coordinate_failures += int(not valid)
        _add_check(
            checks, "domain counts match ranges and sequences", failures, checked,
            "per-name counts equal range and discrete-sequence list lengths",
        )
        _add_check(
            checks, "domain coordinates and sequence slices agree", coordinate_failures, checked,
            "ranges are in bounds and their lengths equal stored domain-sequence lengths",
        )

    elif family == "go":
        failures = 0
        id_failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            for aspect in "CPF":
                values = [
                    _as_list(record.get(f"{aspect}_ids")),
                    _as_list(record.get(f"{aspect}_descriptions")),
                    _as_list(record.get(f"{aspect}_evidence")),
                ]
                checked += 1
                failures += int(len({len(value) for value in values}) != 1)
                id_failures += sum(
                    not re.fullmatch(r"GO:\d{7}", str(identifier))
                    for identifier in values[0]
                )
        _add_check(
            checks, "GO ID/name/evidence lists are parallel", failures, checked,
            "equal list lengths within each GO aspect",
        )
        _add_check(
            checks, "GO identifiers are well formed", id_failures, checked,
            "every populated identifier matches GO: followed by seven digits",
        )
        failures, checked = _canonical_leakage(
            frame, catalog_by_key,
            [f"{aspect}_{suffix}" for aspect in "CPF" for suffix in ("ids", "descriptions", "evidence")],
        )
        _add_check(
            checks, "UniProt GO annotations stay canonical-only", failures, checked,
            "all sequence-distinct isoform GO fields are null",
        )

    elif family == "eclip":
        failures = 0
        checked = 0
        for column in [name for name in frame if "_frac_" in name.lower()]:
            for value in frame[column]:
                number = _number(value)
                if number is None:
                    continue
                checked += 1
                failures += int(not 0 <= number <= 1)
        _add_check(
            checks, "CLIP regional fractions are bounded", failures, checked,
            "every populated *_frac_* value is in [0,1]",
        )
        failures = 0
        checked = 0
        allowed = {"", "0", "1", "false", "true", "none", "nan"}
        for column in [name for name in frame if name.startswith("has_")]:
            for value in frame[column]:
                checked += 1
                failures += int(str(value).strip().lower() not in allowed)
        _add_check(
            checks, "CLIP availability flags are Boolean-like", failures, checked,
            "has_* values are empty/unknown or Boolean",
        )

    elif family == "interpro":
        failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            ranges = _as_dict(record.get("InterPro_range"))
            observed = sum(len(_as_list(value)) for value in ranges.values())
            expected = _number(record.get("InterPro_n_hits"))
            if expected is None:
                continue
            checked += 1
            failures += int(observed != int(expected))
        _add_check(
            checks, "InterPro hit counts match retained ranges", failures, checked,
            "InterPro_n_hits equals the total number of stored range records",
        )
        failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            meta = catalog_by_key.get(str(record["protein_key"]), {})
            length = _number(meta.get("length_aa"))
            if length is None:
                continue
            checked += 1
            ranges = [
                item for values in _as_dict(record.get("InterPro_range")).values()
                for item in _as_list(values)
            ]
            valid = all(
                isinstance(item, list) and len(item) == 2
                and _number(item[0]) is not None and _number(item[1]) is not None
                and 0 <= _number(item[0]) <= _number(item[1]) <= length
                for item in ranges
            )
            failures += int(not valid)
        _add_check(
            checks, "InterPro coordinates are in bounds", failures, checked,
            "all retained zero-based half-open ranges lie within the row sequence",
        )

    elif family == "ptm":
        ptm_types = sorted({
            column[len("ptm_"):-len("_positions")]
            for column in frame if column.startswith("ptm_") and column.endswith("_positions")
        })
        failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            for ptm_type in ptm_types:
                indicator = _number(record.get(f"ptm_{ptm_type}"))
                if indicator is None:
                    continue
                checked += 1
                positions = _as_list(record.get(f"ptm_{ptm_type}_positions"))
                residues = _as_list(record.get(f"ptm_{ptm_type}_residues"))
                if indicator not in (0, 1) or len(positions) != len(residues) \
                        or bool(indicator) != bool(positions):
                    failures += 1
        _add_check(
            checks, "PTM indicators and site lists agree", failures, checked,
            "binary flag matches nonempty, equal-length position/residue lists",
        )
        failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            meta = catalog_by_key.get(str(record["protein_key"]), {})
            sequence = str(meta.get("sequence") or "")
            if not sequence:
                continue
            for ptm_type in ptm_types:
                positions = _as_list(record.get(f"ptm_{ptm_type}_positions"))
                residues = _as_list(record.get(f"ptm_{ptm_type}_residues"))
                for position, residue in zip(positions, residues):
                    checked += 1
                    number = _number(position)
                    if number is None or int(number) != number \
                            or not 0 <= int(number) < len(sequence) \
                            or sequence[int(number)] != str(residue):
                        failures += 1
        _add_check(
            checks, "PTM coordinates match canonical residues", failures, checked,
            "every zero-based position is in bounds and stores the matching residue",
        )
        ptm_columns = [
            column for column in frame
            if column.startswith("ptm_") and column != "ptm_annotation_scope"
        ]
        failures, checked = _canonical_leakage(frame, catalog_by_key, ptm_columns)
        _add_check(
            checks, "PTM annotations stay canonical-only", failures, checked,
            "all sequence-distinct isoform PTM fields are null",
        )

    elif family == "opentargets":
        expression_failures = 0
        disease_failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            tissues = _as_list(record.get("opentargets_tissue_expression"))
            tissue_ids = {
                (item.get("ensembl_gene_id"), item.get("tissue_id"), item.get("tissue_name"))
                for item in tissues if isinstance(item, dict)
            }
            diseases = _as_list(record.get("opentargets_disease_associations"))
            disease_ids = {
                item.get("disease_id") for item in diseases
                if isinstance(item, dict) and item.get("disease_id")
            }
            checked += 1
            expression_failures += int(
                int(_number(record.get("opentargets_expression_tissue_count")) or 0)
                != len(tissue_ids)
            )
            disease_failures += int(
                int(_number(record.get("opentargets_disease_count")) or 0)
                != len(disease_ids)
            )
        _add_check(
            checks, "Open Targets tissue counts are normalized", expression_failures, checked,
            "tissue count equals unique ENSG/tissue records",
        )
        _add_check(
            checks, "Open Targets disease counts are normalized", disease_failures, checked,
            "disease count equals unique disease identifiers",
        )
        releases = {
            str(value) for value in frame.get("opentargets_release", [])
            if _nonnull(value) and str(value).strip()
        }
        _add_check(
            checks, "Open Targets release is pinned", max(0, len(releases) - 1),
            len(frame), "exactly one nonempty release value is used",
        )
        eligible = sum(
            bool(_as_list(record.get("ensembl_gene_ids")))
            for record in catalog.to_dict("records")
        )
        positive = sum(
            (_number(record.get("opentargets_expression_tissue_count")) or 0) > 0
            or (_number(record.get("opentargets_disease_count")) or 0) > 0
            for record in frame.to_dict("records")
        )
        catastrophic_zero = int(eligible >= 100 and positive == 0)
        _add_check(
            checks, "Open Targets has biological coverage", catastrophic_zero,
            eligible, "a catalog with at least 100 ENSG-mapped rows must not yield zero annotated rows",
        )

    elif family == "cdcode":
        columns = [
            "Condensate Name", "UID", "Condensate Type", "Species Tax Id",
            "Proteins", "DNA", "RNA", "C-mods", "Condensatopathy",
            "Confidence Score",
        ]
        failures, checked = _parallel_list_failures(frame, columns)
        _add_check(
            checks, "CD-CODE fields are positionally parallel", failures, checked,
            "all ten populated condensate lists have equal length",
        )
        uid_failures = 0
        checked = 0
        if "UID" in frame:
            for value in frame["UID"]:
                uids = _as_list(value)
                if not uids:
                    continue
                checked += 1
                uid_failures += int(len(uids) != len({str(uid) for uid in uids}))
        _add_check(
            checks, "CD-CODE UIDs are unique within rows", uid_failures, checked,
            "no condensate UID is duplicated for one protein",
        )
        failures, checked = _canonical_leakage(frame, catalog_by_key, columns)
        _add_check(
            checks, "CD-CODE membership stays canonical-only", failures, checked,
            "all sequence-distinct isoform CD-CODE fields are null",
        )

    elif family == "string":
        failures = 0
        checked = 0
        column = "string_partners_ensp_by_query"
        if column in frame:
            for cell in frame[column]:
                for partners in _as_dict(cell).values():
                    if not isinstance(partners, dict):
                        failures += 1
                        checked += 1
                        continue
                    for score in partners.values():
                        checked += 1
                        number = _number(score)
                        failures += int(number is None or not 0 <= number <= 1000)
        _add_check(
            checks, "STRING confidence scores are bounded", failures, checked,
            "every retained score is numeric and lies in [0,1000]",
        )

    elif family == "go_roles":
        columns = [
            "role_in_transcription", "role_in_translation",
            "role_in_mrna_stability", "role_in_translation_stability",
        ]
        failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            values = [_number(record.get(column)) for column in columns]
            if all(value is None for value in values):
                continue
            checked += 1
            if any(value not in (0, 1) for value in values) or values[3] != max(values[1], values[2]):
                failures += 1
        _add_check(
            checks, "GO-role flags are internally consistent", failures, checked,
            "all flags are binary and the combined flag is translation OR mRNA stability",
        )
        failures, checked = _canonical_leakage(frame, catalog_by_key, columns)
        _add_check(
            checks, "GO-role labels stay canonical-only", failures, checked,
            "all sequence-distinct isoform role fields are null",
        )

    elif family == "pslab":
        columns = [
            column for column in frame
            if column not in {"protein_key", "pslab_annotation_scope"}
        ]
        failures, checked = _parallel_list_failures(frame, columns)
        _add_check(
            checks, "PSLab per-IDR outputs are parallel", failures, checked,
            "all populated PSLab lists have equal length within a row",
        )
        failures = 0
        checked = 0
        for column in (
            "Saturation concentration [mg/mL]",
            "Saturation concentration [uM]",
        ):
            if column not in frame:
                continue
            for cell in frame[column]:
                for value in _as_list(cell):
                    number = _number(value)
                    if number is None:
                        continue
                    checked += 1
                    failures += int(number < 0)
        _add_check(
            checks, "PSLab concentrations are nonnegative", failures, checked,
            "all finite predicted saturation concentrations are at least zero",
        )

    elif family == "rcsb":
        failures = 0
        checked = 0
        for record in frame.to_dict("records"):
            count = _number(record.get("RCSB_PDB_count"))
            if count is None:
                continue
            checked += 1
            failures += int(count < 0 or int(count) != len(_as_list(record.get("RCSB_PDB_IDs"))))
        _add_check(
            checks, "RCSB entry counts match PDB ID lists", failures, checked,
            "RCSB_PDB_count is nonnegative and equals the number of unique PDB IDs",
        )
        failures = 0
        checked = 0
        count_columns = [column for column in frame if column.startswith("RCSB_") and column.endswith("_count")]
        for column in count_columns:
            for value in frame[column]:
                number = _number(value)
                if number is None:
                    continue
                checked += 1
                failures += int(number < 0 or int(number) != number)
        _add_check(
            checks, "RCSB counts are nonnegative integers", failures, checked,
            "all populated *_count values are whole numbers at least zero",
        )
        columns = [column for column in frame if column.startswith("RCSB_")]
        failures, checked = _canonical_leakage(frame, catalog_by_key, columns)
        _add_check(
            checks, "RCSB mappings stay canonical-only", failures, checked,
            "all sequence-distinct isoform RCSB fields are null",
        )

    return checks


def _identifier_sanity_checks(catalog):
    checks = []
    _add_check(
        checks, "protein_key is unique", int(catalog["protein_key"].duplicated().sum()),
        len(catalog), "zero duplicate catalog keys",
    )
    length_failures = 0
    hash_failures = 0
    canonical_failures = 0
    np_failures = 0
    for record in catalog.to_dict("records"):
        sequence = str(record.get("sequence") or "")
        length = _number(record.get("length_aa"))
        length_failures += int(length is None or int(length) != len(sequence))
        expected_hash = record.get("sequence_sha256")
        if expected_hash:
            hash_failures += int(
                hashlib.sha256(sequence.encode("ascii")).hexdigest() != str(expected_hash)
            )
        if record.get("row_kind") == "swissprot_canonical":
            canonical_failures += int(not record.get("uniprot_id"))
        np_failures += sum(
            not str(accession).startswith("NP_")
            for accession in _as_list(record.get("refseq_protein_ids"))
        )
    _add_check(
        checks, "sequence lengths agree", length_failures, len(catalog),
        "length_aa equals the amino-acid string length",
    )
    _add_check(
        checks, "sequence hashes agree", hash_failures, len(catalog),
        "sequence_sha256 equals SHA-256 of the amino-acid sequence",
    )
    _add_check(
        checks, "canonical rows have UniProt accessions", canonical_failures, len(catalog),
        "every swissprot_canonical row has uniprot_id",
    )
    _add_check(
        checks, "RefSeq protein identifiers are curated NP_ accessions", np_failures, len(catalog),
        "every populated RefSeq protein accession starts with NP_",
    )
    canonical_ids = [
        str(record.get("uniprot_id"))
        for record in catalog.to_dict("records")
        if record.get("row_kind") == "swissprot_canonical" and record.get("uniprot_id")
    ]
    _add_check(
        checks, "canonical UniProt accessions are unique",
        len(canonical_ids) - len(set(canonical_ids)), len(canonical_ids),
        "exactly one canonical catalog row per reviewed UniProt accession",
    )
    return checks


def _is_present(value):
    value = _decode(value)
    if value is None:
        return False
    if isinstance(value, (list, tuple, set, dict, str)):
        return len(value) > 0
    return True


def _container_size(value):
    value = _decode(value)
    if value is None:
        return 0.0
    if isinstance(value, dict):
        nested = [_container_size(item) for item in value.values()]
        return float(sum(nested)) if nested else 0.0
    if isinstance(value, (list, tuple, set)):
        return float(len(value))
    if isinstance(value, bool):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return 1.0 if str(value).strip() else 0.0


def _metric_columns(family, frame):
    available = [c for c in PRIMARY_COLUMNS[family] if c in frame.columns]
    if family == "ptm":
        available = [
            c for c in frame.columns
            if c.startswith("ptm_")
            and not c.endswith(("_positions", "_residues"))
            and c not in {
                "ptm_projection_dropped_sites", "ptm_coordinate_system",
                "ptm_annotation_scope",
            }
            and not c.startswith("ptm_projection_")
        ]
    elif family == "eclip" and not available:
        preferred = [
            c for c in frame.columns
            if any(token in c.lower() for token in ("has_data", "site_count", "peak_count"))
        ]
        available = preferred or [
            c for c in frame.columns
            if c != "protein_key" and not any(x in c.lower() for x in PROVENANCE_WORDS)
        ]
    return available


def _row_metric(family, frame, columns):
    if not columns:
        return [0.0] * len(frame)
    values = []
    for record in frame[columns].itertuples(index=False, name=None):
        if family == "string":
            total = 0
            for value in record:
                decoded = _decode(value)
                if isinstance(decoded, dict):
                    total += sum(
                        len(partners) if isinstance(partners, dict) else 0
                        for partners in decoded.values()
                    )
            values.append(float(total))
        else:
            values.append(sum(_container_size(value) for value in record))
    return values


def _coverage(frame):
    values = []
    denominator = max(len(frame), 1)
    for column in frame.columns:
        if column == "protein_key" or any(x in column.lower() for x in PROVENANCE_WORDS):
            continue
        present = sum(_is_present(value) for value in frame[column])
        values.append((column, present, 100.0 * present / denominator))
    return sorted(values, key=lambda item: (-item[2], item[0]))


def _plot_reports(family, coverage, metric, report_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = os.path.join(report_dir, "figures")
    os.makedirs(figure_dir, exist_ok=True)

    coverage_path = os.path.join(figure_dir, f"{family}_coverage.png")
    shown = coverage[:12]
    fig, axis = plt.subplots(figsize=(8, max(3.2, 0.36 * len(shown))))
    axis.barh([x[0] for x in shown][::-1], [x[2] for x in shown][::-1], color="#315f72")
    axis.set_xlim(0, 100)
    axis.set_xlabel("Rows with a populated value (%)")
    axis.set_title(f"{family}: column coverage")
    fig.tight_layout()
    fig.savefig(coverage_path, dpi=160)
    plt.close(fig)

    distribution_path = os.path.join(figure_dir, f"{family}_distribution.png")
    finite = [x for x in metric if math.isfinite(x)]
    fig, axis = plt.subplots(figsize=(7, 4.2))
    if finite:
        upper = statistics.quantiles(finite, n=100)[98] if len(finite) >= 100 else max(finite)
        clipped = [min(value, upper) for value in finite]
        axis.hist(clipped, bins=min(40, max(5, int(math.sqrt(len(clipped))))), color="#bd5d38")
        if max(finite) > upper:
            axis.set_xlabel(f"Per-row signal count/value (values above p99 clipped to {upper:g})")
        else:
            axis.set_xlabel("Per-row signal count/value")
    axis.set_ylabel("Rows")
    axis.set_title(f"{family}: primary feature distribution")
    fig.tight_layout()
    fig.savefig(distribution_path, dpi=160)
    plt.close(fig)
    return coverage_path, distribution_path


def write_family_report(family, catalog, sidecar, report_dir):
    frame = _read_table(sidecar)
    if "protein_key" not in frame.columns:
        raise ValueError(f"{sidecar} has no protein_key column")
    catalog_keys = set(catalog["protein_key"].astype(str))
    sidecar_keys = set(frame["protein_key"].astype(str))
    coverage = _coverage(frame)
    metric_columns = _metric_columns(family, frame)
    metric = _row_metric(family, frame, metric_columns)
    sanity_checks = _family_sanity_checks(family, catalog, frame)
    coverage_path, distribution_path = _plot_reports(
        family, coverage, metric, report_dir
    )
    positive = sum(value > 0 for value in metric)
    numeric = [value for value in metric if math.isfinite(value)]
    summary = {
        "catalog_rows": len(catalog),
        "sidecar_rows": len(frame),
        "sidecar_columns_excluding_key": len(frame.columns) - 1,
        "duplicate_protein_keys": int(frame["protein_key"].duplicated().sum()),
        "catalog_rows_missing_from_sidecar": len(catalog_keys.difference(sidecar_keys)),
        "sidecar_rows_not_in_catalog": len(sidecar_keys.difference(catalog_keys)),
        "primary_columns": metric_columns,
        "rows_with_nonzero_primary_signal": positive,
        "primary_signal_median": statistics.median(numeric) if numeric else None,
        "primary_signal_mean": statistics.fmean(numeric) if numeric else None,
        "primary_signal_max": max(numeric) if numeric else None,
        "sanity_checks_run": len(sanity_checks),
        "sanity_checks_failed": sum(
            check["status"] == "FAIL" for check in sanity_checks
        ),
    }
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(report_dir, f"{family}.md")
    relative_coverage = os.path.relpath(coverage_path, report_dir).replace("\\", "/")
    relative_distribution = os.path.relpath(distribution_path, report_dir).replace("\\", "/")
    lines = [
        f"# {family} feature report",
        "",
        "This file is generated by `generate_feature_reports.py`; rerun it after changing inputs.",
        "",
        "## Sanity-check summary",
        "",
        "| Check | Value |",
        "|---|---:|",
    ]
    for key, value in summary.items():
        shown_value = ", ".join(value) if isinstance(value, list) else value
        if isinstance(shown_value, float):
            shown_value = f"{shown_value:.3f}"
        lines.append(f"| {key.replace('_', ' ')} | {shown_value} |")
    lines.extend([
        "",
        "A valid sidecar should have zero duplicate keys, zero missing catalog keys, and zero keys outside the catalog.",
        "",
        "## Validation checks",
        "",
        "| Status | Check | Failed rows/items | Checked rows/items | Expectation |",
        "|---|---|---:|---:|---|",
    ])
    for check in sanity_checks:
        lines.append(
            f"| {check['status']} | {check['check']} | {check['failed_rows']:,} | "
            f"{check['checked_rows']:,} | {check['expectation']} |"
        )
    lines.extend([
        "",
        "## Visual checks",
        "",
        f"![Column coverage]({relative_coverage})",
        "",
        f"![Primary distribution]({relative_distribution})",
        "",
        "## Column coverage",
        "",
        "| Column | Populated rows | Coverage |",
        "|---|---:|---:|",
    ])
    for column, present, percentage in coverage:
        lines.append(f"| `{column}` | {present:,} | {percentage:.1f}% |")
    lines.extend([
        "",
        f"Input sidecar: `{os.path.basename(sidecar)}`",
        "",
    ])
    with open(report_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
    return report_path, summary


def write_identifier_report(catalog, report_dir):
    """Write catalog row/identifier coverage as the base-family QC report."""
    identifier_columns = [
        "uniprot_parent_ids", "refseq_protein_ids", "refseq_transcript_ids",
        "ncbi_gene_ids", "ensembl_gene_ids", "ensembl_transcript_ids",
        "ensembl_protein_ids", "hgnc_ids",
    ]
    coverage = []
    for column in identifier_columns:
        if column not in catalog:
            continue
        populated = int(sum(_is_present(value) for value in catalog[column]))
        coverage.append((column, populated, 100.0 * populated / max(len(catalog), 1)))
    row_kinds = catalog["row_kind"].value_counts(dropna=False).to_dict()
    canonical = int(row_kinds.get("swissprot_canonical", 0))
    isoforms = int(row_kinds.get("ncbi_isoform", 0))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = os.path.join(report_dir, "figures")
    os.makedirs(figure_dir, exist_ok=True)
    coverage_path = os.path.join(figure_dir, "identifiers_coverage.png")
    fig, axis = plt.subplots(figsize=(8, 4.5))
    axis.barh(
        [item[0] for item in coverage][::-1],
        [item[2] for item in coverage][::-1],
        color="#315f72",
    )
    axis.set_xlim(0, 100)
    axis.set_xlabel("Rows with at least one identifier (%)")
    axis.set_title("Catalog identifier coverage")
    fig.tight_layout()
    fig.savefig(coverage_path, dpi=160)
    plt.close(fig)

    row_kind_path = os.path.join(figure_dir, "identifiers_row_kinds.png")
    fig, axis = plt.subplots(figsize=(6, 4))
    bars = axis.bar(
        ["Swiss-Prot canonical", "sequence-unique NP isoform"],
        [canonical, isoforms],
        color=["#315f72", "#bd5d38"],
    )
    axis.bar_label(bars, fmt="%d")
    axis.set_ylabel("Rows")
    axis.set_title("Catalog row composition")
    axis.tick_params(axis="x", rotation=12)
    fig.tight_layout()
    fig.savefig(row_kind_path, dpi=160)
    plt.close(fig)

    summary = {
        "catalog_rows": len(catalog),
        "unique_protein_keys": int(catalog["protein_key"].nunique()),
        "duplicate_protein_keys": int(catalog["protein_key"].duplicated().sum()),
        "swissprot_canonical_rows": canonical,
        "sequence_unique_np_isoform_rows": isoforms,
    }
    sanity_checks = _identifier_sanity_checks(catalog)
    summary["sanity_checks_run"] = len(sanity_checks)
    summary["sanity_checks_failed"] = sum(
        check["status"] == "FAIL" for check in sanity_checks
    )
    lines = [
        "# Identifier and catalog report", "",
        "## Sanity-check summary", "",
        "| Check | Value |", "|---|---:|",
    ]
    for key, value in summary.items():
        lines.append(f"| {key.replace('_', ' ')} | {value:,} |")
    lines.extend([
        "", "## Validation checks", "",
        "| Status | Check | Failed rows/items | Checked rows/items | Expectation |",
        "|---|---|---:|---:|---|",
    ])
    for check in sanity_checks:
        lines.append(
            f"| {check['status']} | {check['check']} | {check['failed_rows']:,} | "
            f"{check['checked_rows']:,} | {check['expectation']} |"
        )
    lines.extend([
        "", "## Identifier coverage", "",
        "| Identifier column | Populated rows | Coverage |", "|---|---:|---:|",
    ])
    for column, populated, percentage in coverage:
        lines.append(f"| `{column}` | {populated:,} | {percentage:.1f}% |")
    lines.extend([
        "", "## Visual checks", "",
        "![Identifier coverage](figures/identifiers_coverage.png)", "",
        "![Catalog row composition](figures/identifiers_row_kinds.png)", "",
    ])
    report_path = os.path.join(report_dir, "identifiers.md")
    with open(report_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
    return report_path, summary


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work-dir", required=True)
    p.add_argument(
        "--features", default="all",
        help="all, default, or comma-separated feature families",
    )
    p.add_argument("--output-dir", help="default: WORK_DIR/reports")
    p.add_argument("--skip-unavailable", action="store_true")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    config = rebuild_from_scratch.layout(args.work_dir)
    catalog = _read_table(config["catalog"])
    if "protein_key" not in catalog.columns:
        raise ValueError("catalog has no protein_key column")
    families = rebuild_from_scratch.parse_features(args.features)
    report_dir = os.path.abspath(args.output_dir or os.path.join(args.work_dir, "reports"))
    generated = []
    summaries = {}
    identifier_report, identifier_summary = write_identifier_report(
        catalog, report_dir
    )
    generated.append(identifier_report)
    summaries["identifiers"] = identifier_summary
    print(f"wrote {identifier_report}")
    for family in families:
        sidecar = os.path.join(config["features"], family + ".parquet")
        if not os.path.exists(sidecar):
            if args.skip_unavailable:
                print(f"SKIP {family}: {sidecar} does not exist")
                continue
            raise FileNotFoundError(sidecar)
        report, summary = write_family_report(family, catalog, sidecar, report_dir)
        generated.append(report)
        summaries[family] = summary
        print(f"wrote {report}")
    with open(os.path.join(report_dir, "report_manifest.json"), "w", encoding="utf-8") as handle:
        json.dump(summaries, handle, indent=2)
    return generated


if __name__ == "__main__":
    main()
