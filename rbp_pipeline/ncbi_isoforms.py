"""Download and parse the NCBI human RefSeq protein/isoform package.

NCBI Datasets is used in bulk, once for taxon 9606.  The product report is
essential: FASTA headers identify protein accessions, while the report links
those proteins to NCBI Gene, transcript, isoform-name and Ensembl identifiers.

This module does *not* decide how many table rows to create.  It exposes the
source records; :mod:`isoform_catalog` subsequently groups them by
``(NCBI Gene ID, protein sequence)`` so many transcripts encoding the same
protein do not inflate the table.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile
from dataclasses import dataclass, field


DATASETS_URLS = {
    "Windows": "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe",
    "Linux": "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets",
    "Darwin": "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets",
}


@dataclass
class GeneInfo:
    gene_id: str
    symbol: str | None = None
    description: str | None = None
    synonyms: list[str] = field(default_factory=list)
    tax_id: str = "9606"
    hgnc_ids: list[str] = field(default_factory=list)
    ensembl_gene_ids: list[str] = field(default_factory=list)
    swiss_prot_accessions: list[str] = field(default_factory=list)
    annotation_names: list[str] = field(default_factory=list)


@dataclass
class ProductInfo:
    gene_id: str
    gene_symbol: str | None
    gene_description: str | None
    refseq_protein: str
    refseq_transcript: str | None
    isoform_name: str | None
    ensembl_transcript: str | None
    ensembl_protein: str | None


def _sha256(path, chunk=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_datasets_cli(destination=None):
    """Return a usable NCBI ``datasets`` executable, downloading if needed."""
    existing = shutil.which("datasets")
    if existing:
        return existing

    system = platform.system()
    if system not in DATASETS_URLS:
        raise RuntimeError(
            f"automatic NCBI Datasets installation is unsupported on {system}; "
            "install `datasets` and place it on PATH"
        )
    suffix = ".exe" if system == "Windows" else ""
    destination = destination or os.path.join(
        tempfile.gettempdir(), "ncbi-datasets-tools", "datasets" + suffix
    )
    if os.path.exists(destination):
        return destination

    os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)
    temporary = destination + ".part"
    urllib.request.urlretrieve(DATASETS_URLS[system], temporary)
    os.replace(temporary, destination)
    if system != "Windows":
        os.chmod(destination, 0o755)
    return destination


def download_human_gene_package(
    out_zip,
    datasets_exe=None,
    include_predicted=False,
    force=False,
):
    """Download human protein FASTA plus gene and product JSONL reports.

    RefSeq packages may contain both curated ``NP_`` and predicted ``XP_``
    products.  The package itself is the same either way; ``include_predicted``
    is recorded in the manifest and enforced when parsing.
    """
    if os.path.exists(out_zip) and not force:
        return out_zip
    datasets_exe = ensure_datasets_cli(datasets_exe)
    os.makedirs(os.path.dirname(os.path.abspath(out_zip)), exist_ok=True)
    temporary = out_zip + ".part"
    if os.path.exists(temporary):
        os.remove(temporary)
    command = [
        datasets_exe,
        "download",
        "gene",
        "taxon",
        "human",
        "--include",
        "protein,product-report",
        "--filename",
        temporary,
        "--no-progressbar",
    ]
    subprocess.run(command, check=True)
    os.replace(temporary, out_zip)
    manifest = {
        "source": "NCBI Datasets gene data package",
        "taxon": "Homo sapiens (9606)",
        "command": command,
        "include_predicted_requested": bool(include_predicted),
        "downloaded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "bytes": os.path.getsize(out_zip),
        "sha256": _sha256(out_zip),
    }
    with open(out_zip + ".manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return out_zip


def extract_gene_package(package_zip, out_dir, force=False):
    """Extract a Datasets gene package and return its three required paths."""
    required = {
        "protein_fasta": os.path.join(out_dir, "ncbi_dataset", "data", "protein.faa"),
        "gene_report": os.path.join(out_dir, "ncbi_dataset", "data", "data_report.jsonl"),
        "product_report": os.path.join(
            out_dir, "ncbi_dataset", "data", "product_report.jsonl"
        ),
    }
    if not force and all(os.path.exists(path) for path in required.values()):
        return required
    os.makedirs(out_dir, exist_ok=True)
    with zipfile.ZipFile(package_zip) as archive:
        archive.extractall(out_dir)
    missing = [name for name, path in required.items() if not os.path.exists(path)]
    if missing:
        raise FileNotFoundError(
            f"NCBI package is missing required members: {', '.join(missing)}"
        )
    return required


def iter_jsonl(path):
    with open(path, encoding="utf-8-sig", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSON at {path}:{line_number}") from exc


def load_gene_reports(path):
    """Return ``NCBI Gene ID -> GeneInfo`` from ``data_report.jsonl``."""
    out = {}
    for report in iter_jsonl(path):
        gene_id = str(report.get("gene_id") or "")
        if not gene_id:
            continue
        nomenclature = report.get("nomenclature_authority") or {}
        hgnc = nomenclature.get("identifier")
        out[gene_id] = GeneInfo(
            gene_id=gene_id,
            symbol=report.get("symbol"),
            description=report.get("description"),
            synonyms=sorted(set(report.get("synonyms") or [])),
            tax_id=str(report.get("tax_id") or "9606"),
            hgnc_ids=sorted({hgnc} if hgnc else set()),
            ensembl_gene_ids=sorted(set(report.get("ensembl_gene_ids") or [])),
            swiss_prot_accessions=sorted(
                set(report.get("swiss_prot_accessions") or [])
            ),
            annotation_names=sorted(
                {
                    annotation.get("annotation_name")
                    for annotation in report.get("annotations") or []
                    if annotation.get("annotation_name")
                }
            ),
        )
    return out


def load_product_reports(path, include_predicted=False):
    """Return one :class:`ProductInfo` per transcript/product association."""
    products = []
    for report in iter_jsonl(path):
        gene_id = str(report.get("gene_id") or "")
        if not gene_id:
            continue
        for transcript in report.get("transcripts") or []:
            protein = transcript.get("protein") or {}
            accession = protein.get("accession_version")
            if not accession:
                continue
            if not include_predicted and not accession.startswith("NP_"):
                continue
            products.append(
                ProductInfo(
                    gene_id=gene_id,
                    gene_symbol=report.get("symbol"),
                    gene_description=report.get("description"),
                    refseq_protein=accession,
                    refseq_transcript=transcript.get("accession_version"),
                    isoform_name=protein.get("isoform_name"),
                    ensembl_transcript=transcript.get("ensembl_transcript"),
                    ensembl_protein=protein.get("ensembl_protein"),
                )
            )
    return products


def load_protein_sequences(path, include_predicted=False):
    """Return ``RefSeq accession.version -> amino-acid sequence``."""
    from Bio import SeqIO

    out = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            accession = record.id
            if not include_predicted and not accession.startswith("NP_"):
                continue
            sequence = str(record.seq).upper().replace("*", "")
            if sequence:
                out[accession] = sequence
    return out
