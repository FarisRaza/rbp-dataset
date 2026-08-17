"""Post-translational modification (PTM) feature family.

The project snapshot ``df_ptm.csv`` contains canonical UniProt sites in eleven
PTM classes. Coordinates are zero-based residue indexes. In the clean master
build, canonical rows are joined directly and sequence-distinct NCBI isoform
rows remain null because these observations belong to the canonical UniProt
sequence. Alignment projection remains available as an explicit opt-in utility
for analyses that want inferred rather than directly sourced isoform sites.
Direct canonical sites are residue-validated against the current Swiss-Prot
sequence; stale snapshot coordinates are dropped and counted in provenance.

The companion ``ptm.txt`` is a raw tab-separated site export with one-based
site labels such as ``K117``.  :func:`load_raw_site_table` can parse that form,
but the local copy is truncated part-way through N-glycosylation, so the full
wide snapshot is the default source for this project.
"""

from __future__ import annotations

import ast
import csv
import re
from collections import defaultdict


PTM_TYPES = [
    "acetylation",
    "n_glycosylation",
    "o_glycosylation",
    "c_glycosylation",
    "s_glycosylation",
    "methylation",
    "myristoylation",
    "phosphorylation",
    "sumoylation",
    "ubiquitination",
    "s_nitrosylation",
]

RAW_NAME_MAP = {
    "ACETYLATION": "acetylation",
    "N-GLYCOSYLATION": "n_glycosylation",
    "O-GLYCOSYLATION": "o_glycosylation",
    "C-GLYCOSYLATION": "c_glycosylation",
    "S-GLYCOSYLATION": "s_glycosylation",
    "METHYLATION": "methylation",
    "MYRISTOYLATION": "myristoylation",
    "PHOSPHORYLATION": "phosphorylation",
    "SUMOYLATION": "sumoylation",
    "UBIQUITINATION": "ubiquitination",
    "S-NITROSYLATION": "s_nitrosylation",
}

PTM_COLUMNS = [
    column
    for ptm_type in PTM_TYPES
    for column in (
        f"ptm_{ptm_type}",
        f"ptm_{ptm_type}_positions",
        f"ptm_{ptm_type}_residues",
    )
]

PROVENANCE_COLUMNS = [
    "ptm_projection_source_uniprot_ids",
    "ptm_projection_methods",
    "ptm_projection_dropped_sites",
    "ptm_coordinate_system",
]


def _parse_list(value):
    if value in (None, ""):
        return []
    if isinstance(value, list):
        return value
    try:
        parsed = ast.literal_eval(str(value))
    except (SyntaxError, ValueError):
        return []
    return list(parsed) if isinstance(parsed, (list, tuple)) else []


def empty_annotation():
    row = {}
    for ptm_type in PTM_TYPES:
        row[f"ptm_{ptm_type}"] = 0
        row[f"ptm_{ptm_type}_positions"] = []
        row[f"ptm_{ptm_type}_residues"] = []
    row.update(
        ptm_projection_source_uniprot_ids=[],
        ptm_projection_methods=[],
        ptm_projection_dropped_sites=0,
        ptm_coordinate_system="0-based amino-acid index",
    )
    return row


def load_wide_csv(path):
    """Load the project's canonical-UniProt PTM snapshot."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        for source in csv.DictReader(handle):
            accession = source.get("uniprot_id")
            if not accession:
                continue
            row = empty_annotation()
            for ptm_type in PTM_TYPES:
                positions = [
                    int(value)
                    for value in _parse_list(source.get(f"ptm_{ptm_type}_positions"))
                ]
                residues = [
                    str(value)
                    for value in _parse_list(source.get(f"ptm_{ptm_type}_residues"))
                ]
                sites = sorted(set(zip(positions, residues)))
                row[f"ptm_{ptm_type}"] = int(bool(sites))
                row[f"ptm_{ptm_type}_positions"] = [position for position, _ in sites]
                row[f"ptm_{ptm_type}_residues"] = [residue for _, residue in sites]
            out[accession] = row
    return out


def load_raw_site_table(path):
    """Parse raw one-site-per-line PTM data into the same canonical index.

    Expected fields are ``PTM_TYPE, source, UniProt, gene, organism, site``.
    Site labels are one-based (``K117``); output positions are zero-based.
    """
    sites = defaultdict(lambda: defaultdict(set))
    site_pattern = re.compile(r"^([A-Z])([1-9][0-9]*)$")
    with open(path, encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            ptm_type = RAW_NAME_MAP.get(fields[0].strip().upper())
            accession = fields[2].strip()
            match = site_pattern.match(fields[5].strip())
            if not ptm_type or not accession or not match:
                continue
            sites[accession][ptm_type].add((int(match.group(2)) - 1, match.group(1)))
    out = {}
    for accession, grouped in sites.items():
        row = empty_annotation()
        for ptm_type, values in grouped.items():
            ordered = sorted(values)
            row[f"ptm_{ptm_type}"] = int(bool(ordered))
            row[f"ptm_{ptm_type}_positions"] = [x[0] for x in ordered]
            row[f"ptm_{ptm_type}_residues"] = [x[1] for x in ordered]
        out[accession] = row
    return out


def _coordinate_map(source_sequence, target_sequence):
    """Return ``(source-index -> target-index, method)`` for one isoform pair."""
    if source_sequence == target_sequence:
        return {i: i for i in range(len(source_sequence))}, "exact_sequence"
    start = target_sequence.find(source_sequence)
    if start >= 0:
        return {
            i: start + i for i in range(len(source_sequence))
        }, "canonical_is_contiguous_subsequence_of_isoform"
    start = source_sequence.find(target_sequence)
    if start >= 0:
        return {
            start + i: i for i in range(len(target_sequence))
        }, "isoform_is_contiguous_subsequence_of_canonical"

    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(source_sequence, target_sequence)[0]
    mapping = {}
    source_blocks, target_blocks = alignment.aligned
    for (source_start, source_end), (target_start, target_end) in zip(
        source_blocks, target_blocks
    ):
        # Aligned blocks have equal lengths; mismatches are intentionally mapped
        # here and rejected later unless the modified residue is conserved.
        length = min(source_end - source_start, target_end - target_start)
        for offset in range(length):
            mapping[int(source_start + offset)] = int(target_start + offset)
    return mapping, "global_protein_alignment"


def project_annotation(annotation, source_sequence, target_sequence):
    """Project one canonical annotation to an isoform sequence."""
    out = empty_annotation()
    mapping, method = _coordinate_map(source_sequence, target_sequence)
    dropped = 0
    for ptm_type in PTM_TYPES:
        positions = annotation.get(f"ptm_{ptm_type}_positions") or []
        residues = annotation.get(f"ptm_{ptm_type}_residues") or []
        kept = set()
        for source_position, recorded_residue in zip(positions, residues):
            target_position = mapping.get(int(source_position))
            if (
                target_position is None
                or source_position >= len(source_sequence)
                or target_position >= len(target_sequence)
                or source_sequence[source_position] != recorded_residue
                or target_sequence[target_position] != recorded_residue
            ):
                dropped += 1
                continue
            kept.add((target_position, recorded_residue))
        ordered = sorted(kept)
        out[f"ptm_{ptm_type}"] = int(bool(ordered))
        out[f"ptm_{ptm_type}_positions"] = [x[0] for x in ordered]
        out[f"ptm_{ptm_type}_residues"] = [x[1] for x in ordered]
    out["ptm_projection_methods"] = [method]
    out["ptm_projection_dropped_sites"] = dropped
    return out


def _merge_annotations(annotations):
    out = empty_annotation()
    for annotation in annotations:
        for ptm_type in PTM_TYPES:
            sites = set(
                zip(
                    out[f"ptm_{ptm_type}_positions"],
                    out[f"ptm_{ptm_type}_residues"],
                )
            )
            sites.update(
                zip(
                    annotation[f"ptm_{ptm_type}_positions"],
                    annotation[f"ptm_{ptm_type}_residues"],
                )
            )
            ordered = sorted(sites)
            out[f"ptm_{ptm_type}"] = int(bool(ordered))
            out[f"ptm_{ptm_type}_positions"] = [x[0] for x in ordered]
            out[f"ptm_{ptm_type}_residues"] = [x[1] for x in ordered]
        out["ptm_projection_methods"] = sorted(
            set(out["ptm_projection_methods"] + annotation["ptm_projection_methods"])
        )
        out["ptm_projection_dropped_sites"] += annotation[
            "ptm_projection_dropped_sites"
        ]
    return out


def annotate_rows(rows, canonical_annotations, project_isoforms=True):
    """Yield sidecar rows keyed by ``protein_key`` for a clean catalog.

    Set ``project_isoforms=False`` when canonical UniProt annotations must not
    be transferred to sequence-distinct RefSeq products. In that mode every
    PTM cell on a noncanonical row is null rather than a misleading measured 0.
    """
    canonical_sequences = {
        row["uniprot_id"]: row["sequence"]
        for row in rows
        if row.get("row_kind") == "swissprot_canonical" and row.get("uniprot_id")
    }
    for row in rows:
        canonical = row.get("row_kind") == "swissprot_canonical"
        if not canonical and not project_isoforms:
            yield {
                "protein_key": row["protein_key"],
                **{column: None for column in PTM_COLUMNS + PROVENANCE_COLUMNS},
                "ptm_annotation_scope": (
                    "not applicable: source PTM coordinates are canonical-specific"
                ),
            }
            continue
        sources = []
        projected = []
        if canonical:
            accession = row.get("uniprot_id")
            annotation = canonical_annotations.get(accession)
            if not annotation:
                yield {
                    "protein_key": row["protein_key"],
                    **{column: None for column in PTM_COLUMNS + PROVENANCE_COLUMNS},
                    "ptm_annotation_scope": (
                        "canonical UniProt accession absent from supplied PTM snapshot"
                    ),
                }
                continue
            # The archived PTM snapshot can predate the current Swiss-Prot
            # sequence release. Validate every stored residue before presenting
            # its coordinate as a current canonical annotation; stale sites are
            # dropped and counted in the existing audit field.
            direct = project_annotation(
                annotation, row["sequence"], row["sequence"]
            )
            direct["ptm_projection_methods"] = [
                "direct_canonical_UniProt_join_residue_validated"
            ]
            projected.append(direct)
            sources.append(accession)
        else:
            for accession in row.get("uniprot_parent_ids") or []:
                annotation = canonical_annotations.get(accession)
                source_sequence = canonical_sequences.get(accession)
                if not annotation or not source_sequence:
                    continue
                projected.append(
                    project_annotation(annotation, source_sequence, row["sequence"])
                )
                sources.append(accession)
        result = _merge_annotations(projected) if projected else empty_annotation()
        result["ptm_projection_source_uniprot_ids"] = sorted(set(sources))
        yield {
            "protein_key": row["protein_key"],
            **result,
            "ptm_annotation_scope": (
                "direct canonical UniProt PTM coordinates"
                if canonical else
                "residue-conserving projection from canonical UniProt"
            ),
        }


def main(argv=None):
    import argparse
    from catalog_io import read_rows, write_feature_rows

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="clean catalog CSV/JSONL/Parquet")
    parser.add_argument("--ptm-csv", required=True, help="canonical df_ptm.csv snapshot")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    rows = read_rows(args.input)
    source = load_wide_csv(args.ptm_csv)
    columns = ["protein_key"] + PTM_COLUMNS + PROVENANCE_COLUMNS + [
        "ptm_annotation_scope"
    ]
    write_feature_rows(annotate_rows(rows, source), args.output, columns)


if __name__ == "__main__":
    main()
