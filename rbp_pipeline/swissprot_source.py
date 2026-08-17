"""Download and parse the reviewed human Swiss-Prot canonical proteome.

The source is the UniProt REST ``stream`` endpoint queried for taxon 9606 and
``reviewed:true``.  The downloaded file is standard UniProt text format, so it
contains canonical sequences, GeneID/RefSeq/Ensembl/HGNC cross-references, GO,
and curated feature coordinates in one versioned artifact.

Public entry points
-------------------
``download_human_reviewed``
    Download the current reviewed human set and a SHA256/release manifest.
``iter_records``
    Yield ``(SwissProtMeta, Bio.SwissProt.Record)`` pairs.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass, field


STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    "format=txt&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29"
)


@dataclass
class SwissProtMeta:
    accession: str
    secondary_accessions: list[str]
    entry_name: str
    gene_symbol: str | None
    gene_synonyms: list[str]
    description: str
    sequence: str
    gene_ids: list[str] = field(default_factory=list)
    hgnc_ids: list[str] = field(default_factory=list)
    refseq: list[dict] = field(default_factory=list)
    ensembl: list[dict] = field(default_factory=list)


def _sha256(path, chunk=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def download_human_reviewed(out_path, manifest_path=None, force=False, timeout=900):
    """Download reviewed human UniProtKB entries in Swiss-Prot text format.

    The write is atomic.  ``manifest_path`` defaults to ``out_path +
    '.manifest.json'`` and records the UniProt release header, request URL,
    timestamp, byte count and SHA256.
    """
    if os.path.exists(out_path) and not force:
        return out_path
    import requests

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=os.path.basename(out_path) + ".", suffix=".part",
        dir=os.path.dirname(os.path.abspath(out_path)),
    )
    os.close(fd)
    try:
        with requests.get(STREAM_URL, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            with open(temporary, "wb") as handle:
                for block in response.iter_content(1024 * 1024):
                    if block:
                        handle.write(block)
            headers = {key.lower(): value for key, value in response.headers.items()}
        os.replace(temporary, out_path)
    except Exception:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise

    manifest = {
        "source": "UniProtKB/Swiss-Prot REST stream",
        "query": "organism_id:9606 AND reviewed:true",
        "url": STREAM_URL,
        "uniprot_release": headers.get("x-uniprot-release"),
        "uniprot_release_date": headers.get("x-uniprot-release-date"),
        "downloaded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "bytes": os.path.getsize(out_path),
        "sha256": _sha256(out_path),
    }
    manifest_path = manifest_path or out_path + ".manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return out_path


def _parse_gene_name(text):
    # Biopython 1.85+ exposes GN lines as a list of dictionaries, while older
    # releases exposed the original flat-file string. Support both forms so the
    # rebuild is not tied to one parser version.
    if isinstance(text, list):
        primary = None
        synonyms = []
        for record in text:
            if not isinstance(record, dict):
                continue
            primary = primary or record.get("Name")
            values = record.get("Synonyms") or []
            synonyms.extend(values if isinstance(values, list) else [values])
        return primary, sorted({str(value) for value in synonyms if value})
    name = re.search(r"(?:^|;)\s*Name=([^;{]+)", text or "")
    synonyms = re.search(r"(?:^|;)\s*Synonyms=([^;{]+)", text or "")
    syns = [] if not synonyms else [x.strip() for x in synonyms.group(1).split(",")]
    return (name.group(1).strip() if name else None), syns


def _description(text):
    match = re.search(r"RecName:\s*Full=([^;{]+)", text or "")
    if not match:
        match = re.search(r"SubName:\s*Full=([^;{]+)", text or "")
    return match.group(1).strip() if match else (text or "").strip()


def _isoform_tag(value):
    if not value:
        return None
    match = re.search(r"\[([^]]+)\]", value)
    return match.group(1) if match else None


def _xref_accession(value):
    """Read the leading accession from a cross-reference value.

    Biopython keeps the UniProt isoform tag in the same tuple field, e.g.
    ``'NM_001407583.1. [P38398-7]'``.  The trailing punctuation is formatting,
    not part of the accession.
    """
    if not value:
        return None
    return str(value).split()[0].rstrip(".;")


def metadata(record):
    """Convert a Biopython SwissProt record into stable mapping metadata."""
    gene_symbol, synonyms = _parse_gene_name(record.gene_name)
    gene_ids, hgnc_ids, refseq, ensembl = [], [], [], []
    for xref in record.cross_references:
        database = xref[0]
        values = list(xref[1:])
        if database == "GeneID" and values:
            gene_ids.append(_xref_accession(values[0]))
        elif database == "HGNC" and values:
            hgnc_ids.append(values[0])
        elif database == "RefSeq" and values:
            refseq.append({
                "protein": _xref_accession(values[0] if len(values) > 0 else None),
                "transcript": _xref_accession(values[1] if len(values) > 1 else None),
                "uniprot_isoform": _isoform_tag(" ".join(values)),
            })
        elif database == "Ensembl" and values:
            ensembl.append({
                "transcript": _xref_accession(values[0] if len(values) > 0 else None),
                "protein": _xref_accession(values[1] if len(values) > 1 else None),
                "gene": _xref_accession(values[2] if len(values) > 2 else None),
                "uniprot_isoform": _isoform_tag(" ".join(values)),
            })
    return SwissProtMeta(
        accession=record.accessions[0],
        secondary_accessions=list(record.accessions[1:]),
        entry_name=record.entry_name,
        gene_symbol=gene_symbol,
        gene_synonyms=synonyms,
        description=_description(record.description),
        sequence=str(record.sequence),
        gene_ids=sorted(set(gene_ids)),
        hgnc_ids=sorted(set(hgnc_ids)),
        refseq=refseq,
        ensembl=ensembl,
    )


def iter_records(path):
    """Yield reviewed human ``(metadata, record)`` pairs from ``path``."""
    from Bio import SwissProt
    with open(path, encoding="utf-8", errors="replace") as handle:
        for record in SwissProt.parse(handle):
            # The downloaded source is already human-only, but this check makes
            # the parser safe when given a full Swiss-Prot release.
            taxids = {str(x) for x in getattr(record, "taxonomy_id", [])}
            if taxids and "9606" not in taxids:
                continue
            yield metadata(record), record


def load_manifest(path):
    manifest = path + ".manifest.json"
    if not os.path.exists(manifest):
        return {}
    with open(manifest, encoding="utf-8") as handle:
        return json.load(handle)
