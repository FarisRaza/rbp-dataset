"""RCSB PDB entry and experimental secondary-structure annotations.

This feature family answers two related questions for every reviewed human
UniProt protein:

* which experimental PDB entries map to the protein; and
* which regular secondary-structure elements (alpha helices and beta strands)
  occur in every mapped PDB chain.

The output is deliberately split in two.  ``SUMMARY_COLUMNS`` are compact
protein-level values suitable for broadcasting onto the isoform table.  The
normalized element table has one row per UniProt/PDB-chain/element observation
and preserves all coordinates; putting millions of elements into a single CSV
cell would be difficult to query and would multiply the data across isoforms.

Sources
-------
PDB-to-UniProt residue mapping comes from the weekly SIFTS
``pdb_chain_uniprot.tsv.gz`` release maintained by PDBe and UniProt.  Regular
secondary structures come from the wwPDB/PDBe entry endpoint, which exposes the
PDB archive's helix and strand annotations in compact, batchable JSON.  PDB IDs
link directly to RCSB entry pages and identify the same wwPDB archive entries.

Coordinates
-----------
``pdb_seq_*`` is the 1-based inclusive polymer/entity sequence numbering
returned as ``residue_number`` by the secondary-structure endpoint.  SIFTS
``RES_BEG/RES_END`` uses that same sequence coordinate system.  ``auth_seq_*``
preserves depositor (author) numbering, including insertion codes.  UniProt
ranges are 1-based inclusive and may contain more than one interval when an
element crosses a mapping boundary.  No approximate coordinate is fabricated
when a SIFTS segment has unequal PDB and UniProt lengths.
"""

from __future__ import annotations

import csv
import datetime as _dt
import gzip
import hashlib
import json
import os
import shutil
import threading
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass


SIFTS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/"
    "pdb_chain_uniprot.tsv.gz"
)
SECONDARY_STRUCTURE_URL = (
    "https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure"
)
RCSB_ENTRY_URL = "https://www.rcsb.org/structure/{}"
USER_AGENT = "Kappel-human-proteome-pipeline/1.0"


# Compact columns appended to the protein/isoform table.  Counts are chain-level
# observations: two conformationally distinct chains in one entry remain two
# observations rather than being silently collapsed.
SUMMARY_COLUMNS = [
    "RCSB_PDB_IDs",
    "RCSB_PDB_count",
    "RCSB_PDB_chains",
    "RCSB_PDB_entries_with_secondary_structure_count",
    "RCSB_PDB_entries_without_secondary_structure_count",
    "RCSB_secondary_structure_observation_count",
    "RCSB_helix_observation_count",
    "RCSB_beta_strand_observation_count",
    "RCSB_secondary_structure_complete_mapping_count",
    "RCSB_secondary_structure_partial_mapping_count",
]


ELEMENT_COLUMNS = [
    "uniprot_id",
    "pdb_id",
    "entity_id",
    "label_asym_id",
    "auth_asym_id",
    "element_uid",
    "element_type",
    "element_index",
    "sheet_id",
    "pdb_seq_start",
    "pdb_seq_end",
    "auth_seq_start",
    "auth_seq_start_insertion_code",
    "auth_seq_end",
    "auth_seq_end_insertion_code",
    "element_length",
    "uniprot_ranges",
    "uniprot_mapped_residue_count",
    "uniprot_mapping_status",
]


@dataclass(frozen=True, order=True)
class MappingSegment:
    """One linear SIFTS mapping segment, all coordinates 1-based inclusive."""

    pdb_start: int
    pdb_end: int
    uniprot_start: int
    uniprot_end: int


@dataclass
class SiftsIndex:
    """Only the SIFTS rows relevant to the requested UniProt accessions."""

    release: str
    accessions: tuple[str, ...]
    pdb_ids: tuple[str, ...]
    by_accession: dict
    by_chain: dict
    mapping_segment_count: int
    invalid_segment_count: int


def utc_now():
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def human_swissprot_accessions(fasta_path):
    """Return reviewed human accessions from a complete Swiss-Prot FASTA.

    The project FASTA contains all organisms, so both the ``sp`` namespace and
    ``OX=9606`` are required.  Only headers are inspected; sequences are not
    loaded into memory.
    """
    out = set()
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">sp|") or " OX=9606" not in line:
                continue
            parts = line.split("|", 2)
            if len(parts) >= 2 and parts[1]:
                out.add(parts[1])
    if not out:
        raise ValueError(f"no reviewed human accessions found in {fasta_path}")
    return tuple(sorted(out))


def download_sifts(out_path, url=SIFTS_URL, force=False, timeout=300):
    """Download the current SIFTS mapping atomically and return ``out_path``."""
    if os.path.exists(out_path) and not force:
        return out_path

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    part = out_path + ".part"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response, \
                open(part, "wb") as fh:
            shutil.copyfileobj(response, fh, length=1024 * 1024)
        # A corrupt/truncated gzip must not replace a usable prior release.
        with gzip.open(part, "rb") as fh:
            if not fh.readline().startswith(b"#"):
                raise ValueError("SIFTS download lacks the expected release header")
        os.replace(part, out_path)
    finally:
        if os.path.exists(part):
            os.remove(part)
    return out_path


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def load_sifts(sifts_path, accessions):
    """Index SIFTS by UniProt accession and by ``(PDB, author chain)``.

    ``by_accession[accession]`` contains PDB ids and author chains for the
    compact summary. ``by_chain[(pdb_id, auth_chain)][accession]`` contains the
    residue mapping segments used to project each secondary-structure element
    onto UniProt coordinates.  Chimeric chains can therefore map to more than
    one accession without assigning one protein's element to another.
    """
    wanted = set(accessions)
    by_accession = defaultdict(lambda: {"pdb_ids": set(), "chains": defaultdict(set)})
    staged_chain = defaultdict(lambda: defaultdict(set))
    pdb_ids = set()
    invalid = 0

    with gzip.open(sifts_path, "rt", encoding="utf-8", errors="replace") as fh:
        release_line = fh.readline().strip()
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "PDB", "CHAIN", "SP_PRIMARY", "RES_BEG", "RES_END", "SP_BEG", "SP_END"
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SIFTS file is missing columns: {sorted(missing)}")

        for row in reader:
            accession = row["SP_PRIMARY"]
            if accession not in wanted:
                continue
            pdb_id = row["PDB"].upper()
            chain = row["CHAIN"]
            pdb_ids.add(pdb_id)
            by_accession[accession]["pdb_ids"].add(pdb_id)
            by_accession[accession]["chains"][pdb_id].add(chain)

            values = [_as_int(row[c]) for c in ("RES_BEG", "RES_END", "SP_BEG", "SP_END")]
            if any(v is None for v in values):
                invalid += 1
                continue
            segment = MappingSegment(*values)
            if segment.pdb_end < segment.pdb_start or \
                    segment.uniprot_end < segment.uniprot_start:
                invalid += 1
                continue
            staged_chain[(pdb_id, chain)][accession].add(segment)

    frozen_chain = {
        key: {acc: tuple(sorted(segments)) for acc, segments in mappings.items()}
        for key, mappings in staged_chain.items()
    }
    frozen_accession = {}
    for accession in accessions:
        record = by_accession.get(accession)
        if record is None:
            frozen_accession[accession] = {"pdb_ids": (), "chains": {}}
        else:
            frozen_accession[accession] = {
                "pdb_ids": tuple(sorted(record["pdb_ids"])),
                "chains": {
                    pdb: tuple(sorted(chains))
                    for pdb, chains in sorted(record["chains"].items())
                },
            }

    return SiftsIndex(
        release=release_line.lstrip("# "),
        accessions=tuple(accessions),
        pdb_ids=tuple(sorted(pdb_ids)),
        by_accession=frozen_accession,
        by_chain=frozen_chain,
        mapping_segment_count=sum(
            len(segments)
            for mappings in frozen_chain.values()
            for segments in mappings.values()
        ),
        invalid_segment_count=invalid,
    )


def _chunks(values, size):
    if size < 1:
        raise ValueError("batch size must be positive")
    for start in range(0, len(values), size):
        yield start // size, tuple(values[start:start + size])


def _cache_path(cache_dir, batch_number, pdb_ids):
    digest = hashlib.sha256(",".join(pdb_ids).encode("ascii")).hexdigest()[:16]
    return os.path.join(cache_dir, f"batch_{batch_number:04d}_{digest}.json.gz")


def _request_secondary_structures(pdb_ids, timeout=300, max_attempts=5):
    """POST one batch, retrying transient failures and splitting as fallback."""
    body = json.dumps(",".join(pdb_id.lower() for pdb_id in pdb_ids)).encode("utf-8")
    request = urllib.request.Request(
        SECONDARY_STRUCTURE_URL,
        data=body,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    last_error = None
    for attempt in range(max_attempts):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                raw = response.read()
            parsed = json.loads(raw)
            if not isinstance(parsed, dict):
                raise ValueError("secondary-structure endpoint returned non-object JSON")
            return parsed
        except (OSError, ValueError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and 400 <= exc.code < 500 \
                    and exc.code not in (408, 413, 429):
                break
            if attempt + 1 < max_attempts:
                time.sleep(min(2 ** attempt, 16))

    # A single problematic or over-large request should not lose the other IDs.
    if len(pdb_ids) > 1:
        middle = len(pdb_ids) // 2
        left = _request_secondary_structures(
            pdb_ids[:middle], timeout=timeout, max_attempts=max_attempts
        )
        right = _request_secondary_structures(
            pdb_ids[middle:], timeout=timeout, max_attempts=max_attempts
        )
        left.update(right)
        return left
    raise RuntimeError(f"failed to retrieve {pdb_ids[0]}: {last_error}")


def _write_batch_cache(path, pdb_ids, data):
    payload = {
        "endpoint": SECONDARY_STRUCTURE_URL,
        "retrieved_at_utc": utc_now(),
        "requested_ids": list(pdb_ids),
        "returned_ids": sorted(str(k).upper() for k in data),
        "data": data,
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # Include the process/thread identity so a resumed launch can never contend
    # with an older orphaned worker over the same temporary filename.  The
    # final cache path remains deterministic and atomic.
    part = f"{path}.{os.getpid()}.{threading.get_ident()}.part"
    try:
        with gzip.open(part, "wt", encoding="utf-8", compresslevel=5) as fh:
            json.dump(payload, fh, separators=(",", ":"))
        for attempt in range(6):
            try:
                os.replace(part, path)
                break
            except PermissionError:
                if attempt == 5:
                    raise
                time.sleep(0.25 * (attempt + 1))
    finally:
        if os.path.exists(part):
            try:
                os.remove(part)
            except PermissionError:
                pass
    return {
        "path": path,
        "requested": len(pdb_ids),
        "returned": len(data),
        "cached": False,
    }


def _validate_cache(path, expected_ids):
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            payload = json.load(fh)
        if tuple(payload.get("requested_ids", ())) != tuple(expected_ids):
            return False
        if not isinstance(payload.get("data"), dict):
            return False
        return True
    except (OSError, ValueError, json.JSONDecodeError):
        return False


def fetch_secondary_structure_batches(
    pdb_ids,
    cache_dir,
    batch_size=500,
    workers=4,
    force=False,
    timeout=300,
    progress=True,
):
    """Fetch every PDB entry into deterministic, gzip-compressed batch caches.

    Completed batches survive interruption.  Cache names include a hash of the
    exact PDB ids, so changing the SIFTS release or batch size cannot silently
    reuse a mismatched response.
    """
    os.makedirs(cache_dir, exist_ok=True)
    batches = [
        (number, ids, _cache_path(cache_dir, number, ids))
        for number, ids in _chunks(tuple(pdb_ids), batch_size)
    ]
    results = []
    to_fetch = []
    for number, ids, path in batches:
        if not force and os.path.exists(path) and _validate_cache(path, ids):
            results.append({
                "path": path,
                "requested": len(ids),
                "returned": None,
                "cached": True,
            })
        else:
            to_fetch.append((number, ids, path))

    if progress:
        print(
            f"secondary structure: {len(pdb_ids):,} PDB entries in "
            f"{len(batches):,} batches; {len(results):,} cached, "
            f"{len(to_fetch):,} to fetch"
        )

    def one(item):
        _, ids, path = item
        data = _request_secondary_structures(ids, timeout=timeout)
        return _write_batch_cache(path, ids, data)

    completed = 0
    if to_fetch:
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = {pool.submit(one, item): item for item in to_fetch}
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                completed += 1
                if progress:
                    print(
                        f"  fetched {completed:>4}/{len(to_fetch)} batches: "
                        f"{result['requested']} requested, {result['returned']} returned",
                        flush=True,
                    )

    order = {path: i for i, (_, _, path) in enumerate(batches)}
    results.sort(key=lambda item: order[item["path"]])
    return results


def _merge_ranges(ranges):
    """Merge overlapping/adjacent 1-based inclusive coordinate intervals."""
    merged = []
    for start, end in sorted(set(ranges)):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def _union_length(ranges):
    return sum(end - start + 1 for start, end in _merge_ranges(ranges))


def map_element_to_uniprot(start, end, segments):
    """Map one PDB sequence interval to UniProt without approximating gaps.

    Returns ``(uniprot_ranges, mapped_residues, status)``.  ``status`` is one of
    ``complete``, ``partial``, ``outside_uniprot_segment`` or
    ``unmapped_length_mismatch``.
    """
    mapped_uniprot = []
    mapped_pdb = []
    mismatch_overlap = False
    any_overlap = False
    for segment in segments:
        left = max(start, segment.pdb_start)
        right = min(end, segment.pdb_end)
        if left > right:
            continue
        any_overlap = True
        pdb_length = segment.pdb_end - segment.pdb_start + 1
        uniprot_length = segment.uniprot_end - segment.uniprot_start + 1
        if pdb_length != uniprot_length:
            mismatch_overlap = True
            continue
        uni_left = segment.uniprot_start + (left - segment.pdb_start)
        uni_right = segment.uniprot_start + (right - segment.pdb_start)
        mapped_uniprot.append((uni_left, uni_right))
        mapped_pdb.append((left, right))

    ranges = _merge_ranges(mapped_uniprot)
    mapped = _union_length(mapped_pdb)
    length = end - start + 1
    if mapped == length and not mismatch_overlap:
        status = "complete"
    elif mapped:
        status = "partial"
    elif mismatch_overlap:
        status = "unmapped_length_mismatch"
    elif any_overlap:
        status = "partial"
    else:
        status = "outside_uniprot_segment"
    return ranges, mapped, status


def _load_cache(path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def _element_record(
    accession,
    pdb_id,
    entity_id,
    chain,
    element,
    element_type,
    index,
    segments,
):
    start_info = element.get("start") or {}
    end_info = element.get("end") or {}
    start = _as_int(start_info.get("residue_number"))
    end = _as_int(end_info.get("residue_number"))
    if start is None or end is None or end < start:
        return None

    ranges, mapped, status = map_element_to_uniprot(start, end, segments)
    # A helix in an affinity tag or another part of a chimeric construct does
    # not belong to this UniProt protein and must not be attributed to it.
    if status == "outside_uniprot_segment":
        return None

    label_chain = chain.get("struct_asym_id") or ""
    auth_chain = chain.get("chain_id") or ""
    uid = f"{pdb_id}.{label_chain}:{element_type}:{start}-{end}:{index}"
    return {
        "uniprot_id": accession,
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "label_asym_id": label_chain,
        "auth_asym_id": auth_chain,
        "element_uid": uid,
        "element_type": element_type,
        "element_index": index,
        "sheet_id": element.get("sheet_id", "") if element_type == "beta_strand" else "",
        "pdb_seq_start": start,
        "pdb_seq_end": end,
        "auth_seq_start": start_info.get("author_residue_number", ""),
        "auth_seq_start_insertion_code": start_info.get("author_insertion_code") or "",
        "auth_seq_end": end_info.get("author_residue_number", ""),
        "auth_seq_end_insertion_code": end_info.get("author_insertion_code") or "",
        "element_length": end - start + 1,
        "uniprot_ranges": json.dumps(ranges, separators=(",", ":")),
        "uniprot_mapped_residue_count": mapped,
        "uniprot_mapping_status": status,
    }


def empty_summary():
    """Serialized table-compatible values for a protein with no PDB mapping."""
    return {
        "RCSB_PDB_IDs": "[]",
        "RCSB_PDB_count": "0",
        "RCSB_PDB_chains": "{}",
        "RCSB_PDB_entries_with_secondary_structure_count": "0",
        "RCSB_PDB_entries_without_secondary_structure_count": "0",
        "RCSB_secondary_structure_observation_count": "0",
        "RCSB_helix_observation_count": "0",
        "RCSB_beta_strand_observation_count": "0",
        "RCSB_secondary_structure_complete_mapping_count": "0",
        "RCSB_secondary_structure_partial_mapping_count": "0",
    }


def write_sifts_mapping_summary(sifts_path, accessions, summary_path):
    """Write every SIFTS-mapped PDB entry for each supplied UniProt accession.

    This is the bounded, fully automatic structure path used by the master
    rebuild. It captures all PDB entry IDs and mapped chains from the current
    weekly SIFTS release without requiring thousands of secondary-structure API
    requests. Secondary-structure fields remain blank because they were not
    measured by this fast mapping-only route; :mod:`stage_rcsb` can enrich the
    same schema when those observations are wanted.
    """
    sifts = load_sifts(sifts_path, sorted(set(accessions)))
    os.makedirs(os.path.dirname(os.path.abspath(summary_path)), exist_ok=True)
    temporary = summary_path + ".part"
    with open(temporary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["uniprot_id"] + SUMMARY_COLUMNS)
        writer.writeheader()
        for accession in sifts.accessions:
            mapping = sifts.by_accession[accession]
            pdb_ids = list(mapping["pdb_ids"])
            chains = {pdb: list(values) for pdb, values in mapping["chains"].items()}
            writer.writerow({
                "uniprot_id": accession,
                "RCSB_PDB_IDs": repr(pdb_ids),
                "RCSB_PDB_count": len(pdb_ids),
                "RCSB_PDB_chains": repr(chains),
                "RCSB_PDB_entries_with_secondary_structure_count": "",
                "RCSB_PDB_entries_without_secondary_structure_count": "",
                "RCSB_secondary_structure_observation_count": "",
                "RCSB_helix_observation_count": "",
                "RCSB_beta_strand_observation_count": "",
                "RCSB_secondary_structure_complete_mapping_count": "",
                "RCSB_secondary_structure_partial_mapping_count": "",
            })
    os.replace(temporary, summary_path)
    return {
        "sifts_release": sifts.release,
        "reviewed_uniprot_accessions": len(sifts.accessions),
        "proteins_with_pdb": sum(
            bool(sifts.by_accession[accession]["pdb_ids"])
            for accession in sifts.accessions
        ),
        "unique_pdb_entries": len(sifts.pdb_ids),
        "mapping_segments": sifts.mapping_segment_count,
        "summary": os.path.abspath(summary_path),
    }


def consolidate(
    sifts,
    batch_results,
    summary_path,
    elements_path,
    metadata_path=None,
    progress_every=10,
):
    """Write the per-protein summary and normalized element table from caches."""
    os.makedirs(os.path.dirname(os.path.abspath(summary_path)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(elements_path)), exist_ok=True)

    counts = defaultdict(Counter)
    entries_with_sse = defaultdict(set)
    returned_pdb_ids = set()
    element_rows = 0
    skipped_outside = 0
    skipped_invalid = 0

    elements_part = elements_path + ".part"
    opener = gzip.open if elements_path.endswith(".gz") else open
    with opener(elements_part, "wt", encoding="utf-8", newline="", compresslevel=1) \
            if opener is gzip.open else opener(elements_part, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ELEMENT_COLUMNS)
        writer.writeheader()

        for batch_index, result in enumerate(batch_results, 1):
            payload = _load_cache(result["path"])
            data = payload["data"]
            returned_pdb_ids.update(str(pdb_id).upper() for pdb_id in data)

            for raw_pdb_id, entry in data.items():
                pdb_id = str(raw_pdb_id).upper()
                for molecule in entry.get("molecules", []):
                    entity_id = molecule.get("entity_id", "")
                    for chain in molecule.get("chains", []):
                        auth_chain = chain.get("chain_id") or ""
                        accession_segments = sifts.by_chain.get((pdb_id, auth_chain))
                        if not accession_segments:
                            continue
                        secondary = chain.get("secondary_structure") or {}
                        families = (
                            ("helix", secondary.get("helices") or []),
                            ("beta_strand", secondary.get("strands") or []),
                        )
                        for element_type, elements in families:
                            for index, element in enumerate(elements, 1):
                                for accession, segments in accession_segments.items():
                                    record = _element_record(
                                        accession, pdb_id, entity_id, chain, element,
                                        element_type, index, segments,
                                    )
                                    if record is None:
                                        start = _as_int((element.get("start") or {}).get("residue_number"))
                                        end = _as_int((element.get("end") or {}).get("residue_number"))
                                        if start is None or end is None or end < start:
                                            skipped_invalid += 1
                                        else:
                                            skipped_outside += 1
                                        continue
                                    writer.writerow(record)
                                    element_rows += 1
                                    entries_with_sse[accession].add(pdb_id)
                                    counts[accession]["all"] += 1
                                    counts[accession][element_type] += 1
                                    status = record["uniprot_mapping_status"]
                                    if status == "complete":
                                        counts[accession]["complete"] += 1
                                    else:
                                        counts[accession]["partial"] += 1

            if progress_every and batch_index % progress_every == 0:
                print(
                    f"  consolidated {batch_index:>4}/{len(batch_results)} batches; "
                    f"{element_rows:,} element rows",
                    flush=True,
                )

    os.replace(elements_part, elements_path)

    summary_part = summary_path + ".part"
    with open(summary_part, "w", encoding="utf-8", newline="") as fh:
        fieldnames = ["uniprot_id"] + SUMMARY_COLUMNS
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for accession in sifts.accessions:
            mapping = sifts.by_accession[accession]
            pdb_ids = list(mapping["pdb_ids"])
            chains = {pdb: list(values) for pdb, values in mapping["chains"].items()}
            with_sse = entries_with_sse.get(accession, set())
            row_counts = counts[accession]
            writer.writerow({
                "uniprot_id": accession,
                "RCSB_PDB_IDs": repr(pdb_ids),
                "RCSB_PDB_count": len(pdb_ids),
                "RCSB_PDB_chains": repr(chains),
                "RCSB_PDB_entries_with_secondary_structure_count": len(with_sse),
                "RCSB_PDB_entries_without_secondary_structure_count": max(
                    len(pdb_ids) - len(with_sse), 0
                ),
                "RCSB_secondary_structure_observation_count": row_counts["all"],
                "RCSB_helix_observation_count": row_counts["helix"],
                "RCSB_beta_strand_observation_count": row_counts["beta_strand"],
                "RCSB_secondary_structure_complete_mapping_count": row_counts["complete"],
                "RCSB_secondary_structure_partial_mapping_count": row_counts["partial"],
            })
    os.replace(summary_part, summary_path)

    metadata = {
        "created_at_utc": utc_now(),
        "scope": "Reviewed human Swiss-Prot accessions (FASTA headers with OX=9606)",
        "coordinate_conventions": {
            "pdb_seq": "1-based inclusive polymer/entity sequence coordinates",
            "auth_seq": "depositor numbering; insertion codes stored separately",
            "uniprot_ranges": "JSON list of 1-based inclusive intervals",
        },
        "sources": {
            "sifts_mapping": SIFTS_URL,
            "secondary_structure": SECONDARY_STRUCTURE_URL,
            "rcsb_entry_url_template": RCSB_ENTRY_URL,
        },
        "sifts_release": sifts.release,
        "counts": {
            "human_swissprot_accessions": len(sifts.accessions),
            "proteins_with_pdb": sum(
                bool(sifts.by_accession[acc]["pdb_ids"]) for acc in sifts.accessions
            ),
            "unique_pdb_entries_requested": len(sifts.pdb_ids),
            "unique_pdb_entries_returned": len(returned_pdb_ids),
            "pdb_entries_without_endpoint_record": len(set(sifts.pdb_ids) - returned_pdb_ids),
            "sifts_mapping_segments": sifts.mapping_segment_count,
            "invalid_sifts_segments": sifts.invalid_segment_count,
            "secondary_structure_element_rows": element_rows,
            "elements_outside_uniprot_segments_not_attributed": skipped_outside,
            "elements_with_invalid_coordinates_skipped": skipped_invalid,
        },
        "files": {
            "summary": os.path.abspath(summary_path),
            "elements": os.path.abspath(elements_path),
        },
        "definitions": {
            "observation_count": (
                "Chain-specific helix/beta-strand observations; equivalent chains "
                "are retained because their assigned structure can differ."
            ),
            "entries_without_secondary_structure": (
                "Mapped PDB entries for which no helix or beta-strand overlapped "
                "the UniProt-mapped part of any mapped chain."
            ),
            "partial_mapping": (
                "An element overlaps the UniProt mapping but is not fully covered, "
                "or an overlapping SIFTS segment has unequal coordinate lengths."
            ),
        },
    }
    if metadata_path:
        part = metadata_path + ".part"
        with open(part, "w", encoding="utf-8") as fh:
            json.dump(metadata, fh, indent=2, sort_keys=True)
            fh.write("\n")
        os.replace(part, metadata_path)
    return metadata


def load_summary(path):
    """Read the compact summary into ``{uniprot_id: {column: cell}}``."""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        missing = set(["uniprot_id"] + SUMMARY_COLUMNS).difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"RCSB summary is missing columns: {sorted(missing)}")
        return {
            row["uniprot_id"]: {column: row.get(column, "") for column in SUMMARY_COLUMNS}
            for row in reader
            if row.get("uniprot_id")
        }
