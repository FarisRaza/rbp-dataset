"""Stage 1 -- choose the proteins to annotate, and resolve their identifiers.

Selects a *target set* and writes it, plus its identifier cross-references, for
the later stages to consume. Everything downstream reads only
``missing_proteins.json``, so changing the selection here is the single lever
that decides what gets annotated:

    python stage_find_missing.py                 # default: absent from the master
    python stage_find_missing.py --all           # every reviewed human protein
    python stage_find_missing.py --accessions ids.txt   # an explicit list
    python stage_find_missing.py --fasta my.fasta       # sequences from a FASTA

``--all`` is the option to use for "regenerate every feature for every protein"
rather than only filling gaps. The other stages are indifferent to which was
chosen; they simply annotate whatever set lands in the intermediate.

By default the target set is the reviewed human Swiss-Prot accessions in
``uniprot_sprot.fasta`` that do not appear in the ``uniprot_id`` column of the
master table.

Writes three intermediates:
  missing_proteins.json  accessions, FASTA headers and sequences
  missing_idmap.json     each accession's UniProt cross-references
  missing_ensg.json      resolved gene symbol, ENSG and ENSP per accession

ENSG resolution goes through three sources in order of authority, because
UniProt's own mapping alone covers only two thirds of these proteins:

  1. ``HUMAN_9606_idmapping.dat``  -- UniProt's Ensembl cross-reference
  2. ``hgnc_complete_set``         -- by gene symbol
  3. ``target_full.csv``           -- by Open Targets' approved symbol

Together they reach roughly 90%. The remainder are mostly immunoglobulin and
T-cell-receptor segments and uncharacterised ORFs that Ensembl does not carry as
distinct genes, so they legitimately have no ENSG.
"""

import collections
import csv
import json
import re
import sys

import paths

csv.field_size_limit(1 << 30)

CROSS_REFERENCES = {
    "Gene_Name", "Ensembl", "Ensembl_PRO", "GeneID",
    "STRING", "RefSeq", "OpenTargets", "HGNC", "UniProtKB-ID",
}


def human_swissprot(fasta_path):
    """Parse reviewed human entries -> ({accession: sequence}, {accession: header})."""
    sequences, headers = {}, {}
    current, buffer = None, []
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(buffer)
                if "OS=Homo sapiens" in line and "OX=9606" in line:
                    current = line.split("|")[1]
                    headers[current] = line
                else:
                    current = None
                buffer = []
            elif current:
                buffer.append(line)
    if current:
        sequences[current] = "".join(buffer)
    return sequences, headers


def master_accessions(master_path):
    with open(master_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        index = header.index("uniprot_id")
        return {row[index] for row in reader}


def cross_references(idmapping_path, wanted):
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(idmapping_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            accession = parts[0].split("-")[0]
            if accession in wanted and parts[1] in CROSS_REFERENCES:
                out[accession][parts[1]].append(parts[2])
    return {k: dict(v) for k, v in out.items()}


def resolve_gene_ids(accessions, headers, idmap):
    """Resolve gene symbol, ENSG and ENSP for each accession."""
    import pandas as pd

    symbols = {}
    for accession in accessions:
        symbol = (idmap.get(accession, {}).get("Gene_Name") or [None])[0]
        if not symbol:
            match = re.search(r"GN=(\S+)", headers[accession])
            symbol = match.group(1) if match else None
        symbols[accession] = symbol

    ensg = {}
    for accession in accessions:
        raw = (idmap.get(accession, {}).get("Ensembl") or [None])[0]
        ensg[accession] = raw.split(".")[0] if raw else None

    hgnc = pd.read_csv(paths.HGNC, sep="\t", low_memory=False)
    by_symbol = dict(zip(hgnc["symbol"], hgnc["ensembl_gene_id"]))
    for accession in accessions:
        if not ensg[accession]:
            candidate = by_symbol.get(symbols[accession])
            if isinstance(candidate, str):
                ensg[accession] = candidate

    approved = {}
    with open(paths.KAPPEL + r"\target_full.csv", newline="",
              encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        i_id, i_symbol = header.index("id"), header.index("approvedSymbol")
        for row in reader:
            approved[row[i_symbol]] = row[i_id]
    for accession in accessions:
        if not ensg[accession] and symbols[accession] in approved:
            ensg[accession] = approved[symbols[accession]]

    ensp = {
        a: ((idmap.get(a, {}).get("Ensembl_PRO") or [None])[0]) for a in accessions
    }
    return symbols, ensg, ensp


def select(argv):
    """Resolve the target set. Returns (accessions, sequences, headers)."""
    import argparse

    parser = argparse.ArgumentParser(description="choose proteins to annotate")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--all", action="store_true",
                       help="every reviewed human Swiss-Prot protein")
    group.add_argument("--accessions", metavar="FILE",
                       help="file of UniProt accessions, one per line")
    group.add_argument("--fasta", metavar="FILE",
                       help="annotate the sequences in this FASTA instead")
    args = parser.parse_args(argv)

    if args.fasta:
        sequences, headers = human_swissprot(args.fasta)
        if not sequences:  # not UniProt-formatted; fall back to raw parsing
            sequences, headers = any_fasta(args.fasta)
        accessions = sorted(sequences)
        print(f"target set: {len(accessions)} sequences from {args.fasta}")
        return accessions, sequences, headers

    sequences, headers = human_swissprot(paths.UNIPROT_FASTA)
    print(f"reviewed human Swiss-Prot entries: {len(sequences)}")

    if args.accessions:
        with open(args.accessions) as fh:
            wanted = {line.strip() for line in fh if line.strip()}
        accessions = sorted(wanted & set(sequences))
        absent = sorted(wanted - set(sequences))
        print(f"target set: {len(accessions)} of {len(wanted)} requested accessions")
        if absent:
            print(f"  {len(absent)} not in Swiss-Prot, skipped: {absent[:5]}")
        return accessions, sequences, headers

    if args.all:
        accessions = sorted(sequences)
        print(f"target set: all {len(accessions)} reviewed human proteins")
        return accessions, sequences, headers

    present = master_accessions(paths.MASTER_TABLE)
    print(f"proteins already in the master table: {len(present)}")
    accessions = sorted(set(sequences) - present)
    surplus = sorted(present - set(sequences))
    print(f"target set: {len(accessions)} absent from the master")
    print(f"in master but not in human Swiss-Prot (left alone): {len(surplus)}")
    return accessions, sequences, headers


def any_fasta(path):
    """Parse a FASTA whose headers are not UniProt-formatted.

    The accession is taken to be the first whitespace-delimited token of the
    header, with a ``sp|X|Y``-style middle field preferred when present.
    """
    sequences, headers = {}, {}
    current, buffer = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(buffer)
                token = line[1:].split()[0] if len(line) > 1 else ""
                parts = token.split("|")
                current = parts[1] if len(parts) >= 2 else token
                headers[current] = line
                buffer = []
            elif current:
                buffer.append(line)
    if current:
        sequences[current] = "".join(buffer)
    return sequences, headers


def main(argv=None):
    accessions, sequences, headers = select(
        sys.argv[1:] if argv is None else argv)
    missing = accessions

    json.dump(
        {"missing": missing,
         "headers": {a: headers[a] for a in missing},
         "seqs": {a: sequences[a] for a in missing}},
        open(paths.scratch("missing_proteins.json"), "w"),
    )

    idmap = cross_references(paths.IDMAPPING, set(missing))
    json.dump(idmap, open(paths.scratch("missing_idmap.json"), "w"))

    symbols, ensg, ensp = resolve_gene_ids(missing, headers, idmap)
    json.dump({"symbol": symbols, "ensg": ensg, "ensp": ensp},
              open(paths.scratch("missing_ensg.json"), "w"))

    print(f"  gene symbol resolved: {sum(1 for v in symbols.values() if v)}")
    print(f"  ENSG resolved:        {sum(1 for v in ensg.values() if v)}")
    print(f"  ENSP resolved:        {sum(1 for v in ensp.values() if v)}")


if __name__ == "__main__":
    main()
