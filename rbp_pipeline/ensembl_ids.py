"""Optional exact-sequence mapping to current Ensembl peptide identifiers.

UniProt and the NCBI product report supply most ENSG/ENST/ENSP links.  This
module is a conservative fallback: it attaches Ensembl identifiers only when
the complete amino-acid sequence matches exactly, and restricts candidates to
the row's known ENSG when one is available.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import os
import re
import tempfile


PEPTIDE_FASTA_URL = (
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/"
    "Homo_sapiens.GRCh38.pep.all.fa.gz"
)


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def download_current_peptides(out_path, force=False, timeout=900):
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
        with requests.get(PEPTIDE_FASTA_URL, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            with open(temporary, "wb") as handle:
                for block in response.iter_content(1024 * 1024):
                    if block:
                        handle.write(block)
            final_url = response.url
        os.replace(temporary, out_path)
    except Exception:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise

    manifest = {
        "source": "Ensembl current human peptide FASTA",
        "url": PEPTIDE_FASTA_URL,
        "resolved_url": final_url,
        "downloaded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "bytes": os.path.getsize(out_path),
    }
    with open(out_path + ".manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return out_path


def _field(header, name):
    match = re.search(rf"(?:^|\s){re.escape(name)}:([^\s]+)", header)
    return match.group(1) if match else None


def build_sequence_index(fasta_gz):
    """Map sequence SHA256 to exact Ensembl gene/transcript/protein triples."""
    from Bio import SeqIO

    index = {}
    with gzip.open(fasta_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq).upper().replace("*", "")
            if not sequence:
                continue
            mapping = {
                "gene": _field(record.description, "gene"),
                "transcript": _field(record.description, "transcript"),
                "protein": record.id,
            }
            key = sequence_sha256(sequence)
            if mapping not in index.setdefault(key, []):
                index[key].append(mapping)
    return index


def exact_mappings(sequence, index, ensembl_gene_ids=None):
    mappings = list(index.get(sequence_sha256(sequence), []))
    genes = {str(x).split(".")[0] for x in ensembl_gene_ids or []}
    if genes:
        restricted = [
            mapping for mapping in mappings
            if str(mapping.get("gene") or "").split(".")[0] in genes
        ]
        if restricted:
            mappings = restricted
    return mappings

