"""Download and parse the current human RefSeq protein release.

This source is a resilient alternative to the NCBI Datasets taxon-wide gene
package.  It downloads the official split human protein FASTA release and the
GRCh38/T2T GFF annotations.  The FASTA files define the complete protein
accession/sequence universe; GFF ``CDS`` records link each protein accession to
NCBI Gene and RefSeq transcript identifiers.

Only curated ``NP_`` products are exposed to the catalog builder.  Predicted
``XP_`` records are deliberately ignored.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import os
import re
import subprocess
import tempfile
import time
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass, field


PROTEIN_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/"
ASSEMBLIES = (
    "GCF_000001405.40_GRCh38.p14",
    "GCF_009914755.1_T2T-CHM13v2.0",
)
ASSEMBLY_BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/"
    "Homo_sapiens/latest_assembly_versions"
)


@dataclass
class GffProduct:
    gene_ids: set[str] = field(default_factory=set)
    gene_symbols: set[str] = field(default_factory=set)
    transcript_ids: set[str] = field(default_factory=set)
    isoform_names: set[str] = field(default_factory=set)
    assemblies: set[str] = field(default_factory=set)
    ensembl_gene_ids: set[str] = field(default_factory=set)
    ensembl_transcript_ids: set[str] = field(default_factory=set)
    ensembl_protein_ids: set[str] = field(default_factory=set)


@dataclass
class RefSeqProtein:
    accession: str
    sequence: str
    description: str


def versionless(accession):
    return re.sub(r"\.\d+$", "", str(accession or ""))


def _sha256(path, chunk=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(url, destination, force=False, timeout=600):
    if os.path.exists(destination) and not force:
        return destination, {}
    import requests

    os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=os.path.basename(destination) + ".", suffix=".part",
        dir=os.path.dirname(os.path.abspath(destination)),
    )
    os.close(fd)
    headers = {}
    try:
        with requests.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            headers = {key.lower(): value for key, value in response.headers.items()}
            with open(temporary, "wb") as handle:
                for block in response.iter_content(1024 * 1024):
                    if block:
                        handle.write(block)
        os.replace(temporary, destination)
    except Exception:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise
    return destination, headers


def _protein_filenames(timeout=120):
    import requests

    response = requests.get(PROTEIN_BASE_URL, timeout=timeout)
    response.raise_for_status()
    names = set(re.findall(
        r'href="(human\.\d+\.protein\.faa\.gz)"', response.text
    ))
    if not names:
        raise RuntimeError("NCBI RefSeq directory listed no human protein FASTA files")
    return sorted(names, key=lambda value: int(value.split(".")[1]))


def download_accession_product_reports(
    protein_dir,
    gff_paths,
    out_path,
    datasets_exe=None,
    force=False,
    batch_size=100,
):
    """Resolve NP_ products absent from assembly GFFs through NCBI Datasets.

    The RefSeq protein release can be newer than the latest genome annotation
    GFF.  For those accessions only, the lightweight accession summary endpoint
    retrieves current product reports.  Reports are deduplicated by GeneID and
    saved as JSON Lines so the catalog stage remains offline and auditable.
    """
    if os.path.exists(out_path) and not force:
        return out_path
    import ncbi_isoforms

    gff_products, _release = load_gff_products(gff_paths)
    unresolved = [
        protein.accession for protein in iter_np_proteins(protein_dir)
        if versionless(protein.accession) not in gff_products
    ]
    executable = ncbi_isoforms.ensure_datasets_cli(datasets_exe)
    reports_by_gene = {}
    for start in range(0, len(unresolved), batch_size):
        batch = unresolved[start:start + batch_size]
        command = [
            executable, "summary", "gene", "accession", *batch,
            "--report", "product", "--as-json-lines",
        ]
        completed = None
        for attempt in range(3):
            try:
                completed = subprocess.run(
                    command, check=True, capture_output=True, text=True
                )
                break
            except subprocess.CalledProcessError:
                if attempt == 2:
                    raise
                time.sleep(attempt + 1)
        for line in completed.stdout.splitlines():
            if not line.strip():
                continue
            report = json.loads(line)
            gene_id = str(report.get("gene_id") or "")
            if gene_id:
                reports_by_gene[gene_id] = report

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=os.path.basename(out_path) + ".", suffix=".part",
        dir=os.path.dirname(os.path.abspath(out_path)),
    )
    os.close(fd)
    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            for gene_id in sorted(reports_by_gene, key=lambda value: int(value)):
                handle.write(json.dumps(reports_by_gene[gene_id], separators=(",", ":")))
                handle.write("\n")
        os.replace(temporary, out_path)
    except Exception:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise
    return out_path


def download_human_np_release(root, force=False, datasets_exe=None):
    """Download human protein FASTAs and GRCh38/T2T GFFs under ``root``."""
    root = os.path.abspath(root)
    protein_dir = os.path.join(root, "protein")
    gff_dir = os.path.join(root, "gff")
    os.makedirs(protein_dir, exist_ok=True)
    os.makedirs(gff_dir, exist_ok=True)

    artifacts = []
    for name in _protein_filenames():
        path = os.path.join(protein_dir, name)
        _download(PROTEIN_BASE_URL + name, path, force=force)
        artifacts.append({
            "kind": "RefSeq human protein FASTA",
            "url": PROTEIN_BASE_URL + name,
            "path": path,
            "bytes": os.path.getsize(path),
            "sha256": _sha256(path),
        })

    gff_paths = []
    for assembly in ASSEMBLIES:
        name = assembly + "_genomic.gff.gz"
        url = f"{ASSEMBLY_BASE_URL}/{assembly}/{name}"
        path = os.path.join(gff_dir, name)
        _download(url, path, force=force)
        gff_paths.append(path)
        artifacts.append({
            "kind": "NCBI RefSeq genome annotation GFF",
            "assembly": assembly,
            "url": url,
            "path": path,
            "bytes": os.path.getsize(path),
            "sha256": _sha256(path),
        })

    accession_report = os.path.join(root, "accession_product_report.jsonl")
    download_accession_product_reports(
        protein_dir,
        gff_paths,
        accession_report,
        datasets_exe=datasets_exe,
        force=force,
    )
    artifacts.append({
        "kind": "NCBI Datasets accession-level product report for NP_ not in GFF",
        "path": accession_report,
        "bytes": os.path.getsize(accession_report),
        "sha256": _sha256(accession_report),
    })

    manifest = {
        "source": "NCBI RefSeq Homo sapiens FTP release",
        "scope": "curated NP_ proteins; GRCh38 and T2T GeneID mappings",
        "downloaded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "artifacts": artifacts,
    }
    manifest_path = os.path.join(root, "manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return {
        "root": root,
        "protein_dir": protein_dir,
        "gff_paths": gff_paths,
        "accession_product_report": accession_report,
        "manifest": manifest_path,
    }


def existing_release_paths(root):
    root = os.path.abspath(root)
    protein_dir = os.path.join(root, "protein")
    gff_dir = os.path.join(root, "gff")
    proteins = sorted(
        os.path.join(protein_dir, name)
        for name in os.listdir(protein_dir)
        if re.fullmatch(r"human\.\d+\.protein\.faa\.gz", name)
    ) if os.path.isdir(protein_dir) else []
    gffs = [
        os.path.join(gff_dir, assembly + "_genomic.gff.gz")
        for assembly in ASSEMBLIES
    ]
    missing = [path for path in gffs if not os.path.exists(path)]
    if not proteins or missing:
        raise FileNotFoundError(
            "incomplete NCBI RefSeq FTP source; run the download stage first"
        )
    return {
        "root": root,
        "protein_dir": protein_dir,
        "protein_paths": proteins,
        "gff_paths": gffs,
        "accession_product_report": (
            os.path.join(root, "accession_product_report.jsonl")
            if os.path.exists(os.path.join(root, "accession_product_report.jsonl"))
            else None
        ),
        "manifest": os.path.join(root, "manifest.json"),
    }


def iter_np_proteins(protein_dir):
    """Yield one current curated human ``NP_`` accession and sequence."""
    from Bio import SeqIO

    seen = set()
    paths = sorted(
        os.path.join(protein_dir, name)
        for name in os.listdir(protein_dir)
        if re.fullmatch(r"human\.\d+\.protein\.faa\.gz", name)
    )
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                accession = record.id
                if not accession.startswith("NP_"):
                    continue
                if accession in seen:
                    continue
                seen.add(accession)
                sequence = str(record.seq).upper().replace("*", "")
                if sequence:
                    yield RefSeqProtein(
                        accession=accession,
                        sequence=sequence,
                        description=record.description,
                    )


def _attributes(text):
    out = {}
    for item in text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        out[key] = urllib.parse.unquote(value)
    return out


def load_gff_products(paths):
    """Return versionless ``NP_`` accession -> Gene/transcript annotations."""
    products = defaultdict(GffProduct)
    release_lines = []
    for path in paths:
        assembly = os.path.basename(path).removesuffix("_genomic.gff.gz")
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#!annotation-"):
                    release_lines.append(line[2:].strip())
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                attrs = _attributes(fields[8])
                accession = attrs.get("protein_id")
                if not accession or not accession.startswith("NP_"):
                    continue
                product = products[versionless(accession)]
                dbxrefs = attrs.get("Dbxref", "").split(",")
                product.gene_ids.update(
                    value.split(":", 1)[1]
                    for value in dbxrefs if value.startswith("GeneID:")
                )
                for value in dbxrefs:
                    if not value.startswith("Ensembl:"):
                        continue
                    identifier = value.split(":", 1)[1]
                    if identifier.startswith("ENSP"):
                        product.ensembl_protein_ids.add(identifier)
                    elif identifier.startswith("ENST"):
                        product.ensembl_transcript_ids.add(identifier)
                    elif identifier.startswith("ENSG"):
                        product.ensembl_gene_ids.add(identifier)
                if attrs.get("gene"):
                    product.gene_symbols.add(attrs["gene"])
                parent = attrs.get("Parent") or attrs.get("transcript_id")
                if parent:
                    for value in parent.split(","):
                        product.transcript_ids.add(value.removeprefix("rna-"))
                if attrs.get("product"):
                    product.isoform_names.add(attrs["product"])
                product.assemblies.add(assembly)
    return dict(products), sorted(set(release_lines))


def load_accession_product_reports(path):
    """Return NP_ -> Gene/transcript mappings from accession lookup JSONL."""
    products = defaultdict(GffProduct)
    if not path or not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            report = json.loads(line)
            gene_id = str(report.get("gene_id") or "")
            symbol = report.get("symbol")
            for transcript in report.get("transcripts") or []:
                protein = transcript.get("protein") or {}
                accession = protein.get("accession_version")
                if not accession or not accession.startswith("NP_"):
                    continue
                product = products[versionless(accession)]
                if gene_id:
                    product.gene_ids.add(gene_id)
                if symbol:
                    product.gene_symbols.add(symbol)
                if transcript.get("accession_version"):
                    product.transcript_ids.add(transcript["accession_version"])
                if protein.get("isoform_name"):
                    product.isoform_names.add(protein["isoform_name"])
                if transcript.get("ensembl_transcript"):
                    product.ensembl_transcript_ids.add(
                        transcript["ensembl_transcript"]
                    )
                if protein.get("ensembl_protein"):
                    product.ensembl_protein_ids.add(protein["ensembl_protein"])
                product.assemblies.add("NCBI Datasets accession product report")
    return dict(products)
