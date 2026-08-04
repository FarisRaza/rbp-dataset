"""Domains feature family -- named folded domains and zinc fingers.

Source
------
The **UniProt/Swiss-Prot flat file** ``uniprot_sprot.dat``, parsed with
``Bio.SwissProt``. Despite the documentation naming PROSITE and Pfam, and
despite an ``df_np_unique.interpro.tsv`` sitting in the project folder, no
InterPro/InterProScan output reaches these columns. The annotations are the
curated ``DOMAIN`` and ``ZN_FING`` features of the Swiss-Prot entry itself.

That matters here: it means the 2,901 proteins we are adding -- all reviewed
human Swiss-Prot accessions -- can be annotated fully offline from a file
already on disk, with no web service and no RefSeq accession involved.

Produces schema.DOMAIN_GEOMETRY_COLUMNS (7 columns). The companion
schema.DOMAIN_CIDER_COLUMNS (20 columns) are produced by ``cider.per_domain``
from this module's ``Domains_discrete_seq`` output.

Coordinates
-----------
Biopython exposes ``feature.location.start`` as 0-based and ``location.end`` as
the exclusive end, so region sequences are plain ``sequence[start:end]`` slices
and the stored ranges are 0-based half-open. Verified against Q9NQ94 (A1CF):
stored ``{'RRM': [(55,134),(135,218),(230,303)]}`` corresponds to UniProt's
1-based RRM 1 = 56..134, RRM 2 = 136..218, RRM 3 = 231..303.

Granularity
-----------
Domain annotation exists per *protein*, not per isoform. In the master table it
is attached only to rows with ``Dominant_Isoform == 1``; every other isoform row
carries NaN in all 27 domain columns. Rows added by this pipeline are canonical
Swiss-Prot sequences and so are all dominant, and all get domain values.
"""

import re
from collections import defaultdict

FEATURE_TYPES = ("DOMAIN", "ZN_FING")


def normalize_domain(name):
    """Collapse repeat-numbered domain names onto a single vocabulary entry.

    ``"RRM 1"``, ``"RRM 2"``, ``"RRM 3"`` all become ``"RRM"``.

    This reproduces the original regex exactly, quirks included, because the
    already-annotated proteins were built with it and a different normalisation
    would make the new rows incomparable. (That set is the 17,519 rows of
    ``Compressed_Domains.csv`` -- fewer than the master's 17,582 proteins,
    because domain annotation needs a canonical Swiss-Prot sequence to attach
    coordinates to.) The quirks:

      * only a *space-separated* trailing number is stripped, so ``SH2``,
        ``C2``, ``W2`` and ``E1`` survive intact;
      * qualifier suffixes are not normalised, so ``"C2H2-type"`` and
        ``"C2H2-type; degenerate"`` remain distinct domain names.
    """
    return re.sub(r"\s\d+.*$", "", name).strip()


def _feature_bounds(feature):
    """Return (start, end) as ints, or None for a fuzzy/unparseable location.

    Swiss-Prot marks uncertain endpoints with ``<``, ``>`` or ``?``. The
    original pipeline called int() unguarded and would crash on those; skipping
    the feature is the safer behaviour and affects only a handful of entries.
    """
    try:
        return int(str(feature.location.start)), int(str(feature.location.end))
    except (ValueError, TypeError):
        return None


def domains_for_record(record):
    """Extract the domain geometry of one ``Bio.SwissProt`` record.

    Returns a dict over the first five schema.DOMAIN_GEOMETRY_COLUMNS
    (``Domains`` .. ``Domains_range``). Sequence-bearing columns are added by
    `attach_sequences`, which needs the protein sequence.

    A protein with no DOMAIN/ZN_FING features yields ``[]`` and ``{}`` -- which
    is what the master table stores for annotated-but-domainless proteins, as
    distinct from the NaN it stores for non-dominant isoform rows.
    """
    stats = defaultdict(lambda: {"count": 0, "total_size": 0, "ranges": []})

    for feature in record.features:
        if getattr(feature, "type", None) not in FEATURE_TYPES:
            continue
        bounds = _feature_bounds(feature)
        if bounds is None:
            continue
        start, end = bounds
        # ZN_FING features occasionally carry no note; fall back to the feature
        # type, which is what produced the literal 'ZN_FING' domain name in the
        # existing vocabulary.
        note = feature.qualifiers.get("note", feature.type)
        name = normalize_domain(note)
        stats[name]["count"] += 1
        stats[name]["total_size"] += end - start
        stats[name]["ranges"].append((start, end))

    names = sorted(stats)
    return {
        "Domains": names,
        "Domains_count": {n: stats[n]["count"] for n in names},
        "Domains_avg_size": {
            n: stats[n]["total_size"] / stats[n]["count"] for n in names
        },
        "Domains_total_size": {n: stats[n]["total_size"] for n in names},
        "Domains_range": {n: stats[n]["ranges"] for n in names},
    }


def attach_sequences(geometry, sequence):
    """Add ``Domains_discrete_seq`` and ``Domains_concat_seq`` to `geometry`.

    `sequence` must be the same sequence the coordinates were annotated
    against -- i.e. the canonical Swiss-Prot sequence, not an isoform.

    Selenocysteine is rewritten to cysteine here rather than at scoring time.
    localcider raises on ``U``, and the existing table stores the already
    substituted residue, so doing it once at extraction keeps the stored
    sequence and the scored sequence identical.
    """
    discrete, concat = {}, {}
    for name, ranges in geometry["Domains_range"].items():
        pieces = [sequence[a:b].replace("U", "C") for a, b in ranges]
        discrete[name] = pieces
        concat[name] = "".join(pieces)
    out = dict(geometry)
    out["Domains_discrete_seq"] = discrete
    out["Domains_concat_seq"] = concat
    return out


def parse(dat_path, wanted_accessions, sequences):
    """Stream ``uniprot_sprot.dat`` and extract domains for `wanted_accessions`.

    `wanted_accessions` is a set of UniProt accessions; `sequences` maps
    accession -> canonical sequence (used for the two sequence columns).

    Returns ``{accession: {all 7 domain columns}}``. Accessions not found in the
    flat file are simply absent from the result; the caller decides whether that
    is a null or an empty annotation.

    Matching considers every accession on the entry's ``AC`` line, not just the
    primary one, so an entry reachable only by a secondary accession is still
    found -- but the result is keyed by the accession *we asked for*, so a
    merged/demerged accession cannot silently land under a different id.
    """
    from Bio import SwissProt

    wanted = set(wanted_accessions)
    found = {}
    with open(dat_path, "r", encoding="utf-8", errors="replace") as handle:
        for record in SwissProt.parse(handle):
            hits = wanted.intersection(record.accessions)
            if not hits:
                continue
            geometry = domains_for_record(record)
            for accession in hits:
                seq = sequences.get(accession, record.sequence)
                found[accession] = attach_sequences(geometry, seq)
            wanted -= hits
            if not wanted:
                break
    return found


EMPTY = {
    "Domains": [],
    "Domains_count": {},
    "Domains_avg_size": {},
    "Domains_total_size": {},
    "Domains_range": {},
    "Domains_discrete_seq": {},
    "Domains_concat_seq": {},
}
