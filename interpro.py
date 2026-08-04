"""InterPro domain annotation -- provenance and parsing.

This family is **not currently in the 157-column table**. The table's `Domains_*`
columns come from UniProt's own DOMAIN/ZN_FING features (see `domains.py`). The
InterProScan results are a parallel, richer annotation that was produced but
never merged.

How `df_np_unique.interpro.tsv` was actually produced
-----------------------------------------------------
Recovered from the original working session; no script for it existed in the
project folder, which is why this module exists.

Input: ``df_np_unique.fasta`` -- **55,670** sequences, one per unique
``isoform_key`` beginning ``NP_``, written by *Generate FASTA for InterPro.ipynb*
with a bare RefSeq accession as the header and the sequence unwrapped on one
line.

Tool: **InterProScan 5.77-108.0**, run on UCLA's Hoffman2 cluster (SGE) out of
``/u/project/kappel/fraza/InterPro/``. Two earlier attempts on a laptop were
killed by the OOM killer -- InterProScan needs roughly 1-2 GB of RAM per core,
and the 6.6 GB download unpacks to over 20 GB.

    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz
    tar -pxvzf interproscan-5.77-108.0-64-bit.tar.gz

The FASTA was split into ~400 chunks of 140 sequences and submitted as an array
of SGE jobs (`h_rt=24:00:00`, `h_data=8G`, `-pe shared 4`), each running::

    module load java/jdk-11.0.14        # NOT java/11; the cluster default is 8
    ./interproscan.sh \\
        -i chunks/fastas/chunk_${N}.fasta \\
        -f tsv \\
        -o chunks/results/chunk_${N}.interpro.tsv \\
        --cpu 4 -goterms -iprlookup

then concatenated::

    cat chunks/results/*.interpro.tsv > df_np_unique.interpro.tsv

``-goterms`` and ``-iprlookup`` add GO terms and InterPro entry accessions. Both
are pure lookups against already-computed hits, so they cost essentially nothing
beyond the domain scan itself.

`write_chunked_jobs` below regenerates that submission setup, so the run is
reproducible rather than only documented.

Member databases in the output
------------------------------
A default full run: MobiDBLite, Pfam, SMART, Gene3D, FunFam, PRINTS,
ProSiteProfiles, ProSitePatterns, SUPERFAMILY, CDD, PANTHER, Coils, PIRSF,
NCBIfam, Hamap, SFLD.

A parsing bug to avoid
----------------------
The historical notebook reads this file with 14 column names for a 15-column
file, so pandas silently consumes the first column as the index and every label
shifts left by one. That is why it filters with ``tsv["length"] == "Pfam"`` --
the column labelled ``length`` actually holds the member-database name.
`read_tsv` below names all 15 columns correctly.
"""

import csv
import os

csv.field_size_limit(1 << 30)

DEFAULT_TSV = "df_np_unique.interpro.tsv"
VERSION = "5.77-108.0"

#: InterProScan 5 TSV columns, in file order. The last two are present only
#: because the run used -iprlookup and -goterms; a bare run emits 11.
TSV_COLUMNS = [
    "protein_id",        # FASTA header, here a RefSeq NP_ accession
    "md5",
    "length",            # sequence length in residues
    "analysis",          # member database, e.g. Pfam, SMART, PANTHER
    "signature_acc",     # e.g. PF00076
    "signature_desc",    # e.g. RNA recognition motif
    "start",             # 1-based, inclusive
    "stop",              # 1-based, inclusive
    "score",             # e-value or score, analysis-dependent
    "status",            # T when the match is significant
    "date",
    "interpro_acc",      # -iprlookup
    "interpro_desc",     # -iprlookup
    "go_terms",          # -goterms, "|"-separated
    "pathways",
]

#: Member databases that report actual domain models, as opposed to disorder,
#: coils or family-level classification. Use these when comparing against the
#: table's UniProt-derived Domains_* columns.
DOMAIN_DATABASES = {
    "Pfam", "SMART", "ProSiteProfiles", "ProSitePatterns",
    "Gene3D", "SUPERFAMILY", "CDD", "PIRSF", "NCBIfam", "Hamap", "SFLD",
}


def read_tsv(path=DEFAULT_TSV, databases=None, significant_only=True):
    """Stream the InterProScan TSV, yielding one dict per hit.

    `databases` filters to particular member databases (see DOMAIN_DATABASES);
    None keeps everything. `significant_only` keeps rows whose status is ``T``.

    Streaming rather than loading: the file is ~197 MB and pandas' default
    dtype inference on it is slower than the parse itself.
    """
    wanted = set(databases) if databases else None
    minimum = TSV_COLUMNS.index("status") + 1

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            # Split manually rather than with csv.reader: signature descriptions
            # contain bare double quotes (e.g. 5'-3'``exonuclease``), which the
            # csv module reads as quoting and uses to swallow the following
            # fields -- producing short rows that lack even an 'analysis' column.
            fields = line.rstrip("\n").split("\t")
            if len(fields) < minimum:
                continue
            if len(fields) < len(TSV_COLUMNS):
                fields = fields + [""] * (len(TSV_COLUMNS) - len(fields))
            row = dict(zip(TSV_COLUMNS, fields))
            if wanted is not None and row["analysis"] not in wanted:
                continue
            if significant_only and row["status"] not in ("T", ""):
                continue
            yield row


def group_by_protein(path=DEFAULT_TSV, databases=None):
    """Collapse hits into ``{protein_id: {domain_name: [(start, stop), ...]}}``.

    Coordinates are converted to the 0-based half-open convention the table's
    `Domains_range` uses, so the result can be handed to
    `domains.attach_sequences` and `cider.per_domain` unchanged.
    """
    out = {}
    for row in read_tsv(path, databases=databases):
        try:
            start = int(row["start"]) - 1
            stop = int(row["stop"])
        except (TypeError, ValueError):
            continue
        name = row["signature_desc"] or row["interpro_desc"] or row["signature_acc"]
        out.setdefault(row["protein_id"], {}).setdefault(name, []).append((start, stop))
    return out


#: Columns this family contributes, mirroring the shape of the UniProt-derived
#: ``Domains_*`` block so the two are directly comparable.
COLUMNS = [
    "InterPro_domains",
    "InterPro_count",
    "InterPro_range",
    "InterPro_accessions",
    "InterPro_databases",
    "InterPro_go_terms",
    "InterPro_n_hits",
]

EMPTY = {
    "InterPro_domains": [],
    "InterPro_count": {},
    "InterPro_range": {},
    "InterPro_accessions": {},
    "InterPro_databases": {},
    "InterPro_go_terms": {},
    "InterPro_n_hits": 0,
}


def index_by_accession(path=DEFAULT_TSV, databases=None):
    """Build ``{NP_accession: {all InterPro columns}}``.

    `databases` defaults to DOMAIN_DATABASES, which excludes MobiDBLite (a
    disorder predictor, already covered by metapredict), Coils and PANTHER
    (family-level rather than domain-level). Pass None to keep everything.
    """
    if databases is None:
        databases = DOMAIN_DATABASES

    staging = {}
    for row in read_tsv(path, databases=databases):
        try:
            start = int(row["start"]) - 1
            stop = int(row["stop"])
        except (TypeError, ValueError):
            continue
        name = row["signature_desc"] or row["interpro_desc"] or row["signature_acc"]
        if not name:
            continue
        bucket = staging.setdefault(row["protein_id"], {})
        entry = bucket.setdefault(
            name, {"ranges": [], "accessions": set(), "databases": set(), "go": set()}
        )
        entry["ranges"].append((start, stop))
        if row["interpro_acc"] and row["interpro_acc"] != "-":
            entry["accessions"].add(row["interpro_acc"])
        entry["databases"].add(row["analysis"])
        for term in (row["go_terms"] or "").split("|"):
            term = term.split("(")[0].strip()
            if term.startswith("GO:"):
                entry["go"].add(term)

    out = {}
    for accession, domains in staging.items():
        names = sorted(domains)
        out[accession] = {
            "InterPro_domains": names,
            "InterPro_count": {n: len(domains[n]["ranges"]) for n in names},
            "InterPro_range": {n: sorted(domains[n]["ranges"]) for n in names},
            "InterPro_accessions": {n: sorted(domains[n]["accessions"]) for n in names},
            "InterPro_databases": {n: sorted(domains[n]["databases"]) for n in names},
            "InterPro_go_terms": {n: sorted(domains[n]["go"]) for n in names},
            "InterPro_n_hits": sum(len(domains[n]["ranges"]) for n in names),
        }
    return out


def columns_for(protein_hgvs, by_accession):
    """InterPro columns for one table row, keyed via its ``ProteinHGVS`` cell.

    ``ProteinHGVS`` holds one or more comma-joined RefSeq protein accessions;
    the InterProScan run was over exactly those accessions. The first one with
    a hit is used -- rows listing several accessions list them because they are
    the same protein sequence under different RefSeq records.
    """
    if not protein_hgvs or protein_hgvs.startswith("No NP"):
        return dict(EMPTY)
    for accession in str(protein_hgvs).split(","):
        record = by_accession.get(accession.strip())
        if record is not None:
            return record
    return dict(EMPTY)


def write_chunked_jobs(fasta, out_dir, chunk_size=140, iprscan_dir=None,
                       email=None, java_module="java/jdk-11.0.14"):
    """Write the SGE split-and-submit scripts used for the original run.

    Produces ``split_and_submit.sh`` and ``interproscan_job.sh`` in `out_dir`.
    Reproduces the original parameters, including the Java module -- the cluster
    default is Java 8 and InterProScan refuses to start on it.
    """
    os.makedirs(out_dir, exist_ok=True)
    iprscan_dir = iprscan_dir or f"$(dirname {out_dir})/interproscan-{VERSION}"
    notify = f"#$ -M {email}\n#$ -m bea\n" if email else ""

    split_path = os.path.join(out_dir, "split_and_submit.sh")
    with open(split_path, "w", newline="\n") as fh:
        fh.write(f"""#!/bin/bash
# Split {os.path.basename(fasta)} into chunks of {chunk_size} sequences and
# submit one InterProScan job per chunk.
set -euo pipefail

FASTA={fasta}
OUTDIR={out_dir}

mkdir -p "$OUTDIR/fastas" "$OUTDIR/results"

awk -v outdir="$OUTDIR/fastas" '
  /^>/{{n++}}
  {{print > outdir"/chunk_" int((n-1)/{chunk_size}) ".fasta"}}
' "$FASTA"

N=$(ls "$OUTDIR/fastas" | wc -l)
echo "split into $N chunks"

for i in $(seq 0 $((N-1))); do
    qsub "$OUTDIR/interproscan_job.sh" "$i"
done
""")

    job_path = os.path.join(out_dir, "interproscan_job.sh")
    with open(job_path, "w", newline="\n") as fh:
        fh.write(f"""#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=24:00:00,h_data=8G
#$ -pe shared 4
{notify}
CHUNK=$1
OUTDIR={out_dir}
IPRSCAN={iprscan_dir}

echo "Job $JOB_ID chunk $CHUNK started on $(hostname -s) at $(date)"

. /u/local/Modules/default/init/modules.sh
module load {java_module}

"$IPRSCAN/interproscan.sh" \\
    -i "$OUTDIR/fastas/chunk_${{CHUNK}}.fasta" \\
    -f tsv \\
    -o "$OUTDIR/results/chunk_${{CHUNK}}.interpro.tsv" \\
    --cpu 4 \\
    -goterms \\
    -iprlookup

echo "Job $JOB_ID chunk $CHUNK ended at $(date)"
""")

    return split_path, job_path
