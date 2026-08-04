"""Identity columns -- accession, sequence, and cross-references.

Produces schema.IDENTITY_COLUMNS (10 columns) plus the two bookkeeping columns
``isoform_number`` and ``ID_list``.

Where each value comes from, and how faithfully it can be reproduced:

    uniprot_id       the Swiss-Prot accession                    exact
    sequence         the canonical Swiss-Prot sequence           exact
    Dominant_Isoform 1, by construction                          exact
    Name             gene symbol                                 exact
    ID               Ensembl gene id(s), ";"-joined              reconstructed
    ENSP             Ensembl protein id, versioned               reconstructed
    ProteinHGVS      RefSeq protein accession(s), ","-joined     reconstructed
    UNIQUE           gene symbol, RBP-census rows only           exact (usually null)
    Description      gene description, RBP-census rows only      exact (usually null)
    HGVSDescription  NCBI protein titles                         left empty

``Dominant_Isoform`` is defined as "this row's sequence equals the canonical
Swiss-Prot sequence". Rows added here *are* the canonical sequence, so it is
always 1. That matters downstream: the GO, Domains and CD-CODE columns are
blanked on rows where it is 0, so these rows keep their annotations.

``UNIQUE`` and ``Description`` come from ``RBPs.xlsx``, the 1,393-gene published
RBP census that seeded the project. They are populated for only 7.9% of proteins
already in the table; a newly-added protein gets them only if its gene symbol is
in that census, and is otherwise null -- which is how the other 92% of existing
rows already look, not a new kind of gap.

``HGVSDescription`` is free-text NCBI protein titles. It is empty on 64,332 of
the master's 65,744 rows (97.9%), and there is no offline source for it, so it
is left empty here. That is consistent with the overwhelming majority of
existing rows rather than a degradation.

``ENSP``, ``ID`` and ``ProteinHGVS`` were originally produced by a step whose
code is not in the project folder. They are rebuilt here from
``HUMAN_9606_idmapping.dat`` (UniProt's own cross-reference export), falling
back to HGNC and Open Targets for the gene-id lookup. Values will be correct but
are not guaranteed byte-identical to what the vanished step would have emitted --
for instance where UniProt lists several Ensembl protein ids, the first is taken.
"""

import re

# Sentinel the historical pipeline used when no RefSeq protein accession existed.
NO_REFSEQ = "No NP found"


def parse_fasta_header(header):
    """Pull the protein name and gene symbol out of a UniProt FASTA header.

    ``>sp|Q9NQ94|A1CF_HUMAN APOBEC1 complementation factor OS=Homo sapiens ...``
    -> ``("APOBEC1 complementation factor", "A1CF")``
    """
    name_match = re.search(r"^>\w+\|[^|]+\|\S+\s+(.+?)\s+OS=", header)
    gene_match = re.search(r"GN=(\S+)", header)
    return (
        name_match.group(1) if name_match else None,
        gene_match.group(1) if gene_match else None,
    )


def load_rbp_census(path):
    """Read ``RBPs.xlsx`` -> ``{gene_symbol: (UNIQUE, Description)}``.

    The second spreadsheet row holds dataset provenance rather than data, hence
    ``skiprows=[1]``.
    """
    import pandas as pd

    table = pd.read_excel(path, header=0, skiprows=[1])
    census = {}
    for _, row in table.iterrows():
        symbol = row.get("Name")
        if isinstance(symbol, str):
            census[symbol] = (row.get("UNIQUE"), row.get("Description"))
    return census


def build_row(accession, sequence, header, idmap, ensg, census):
    """Assemble the identity columns for one newly-added protein.

    `idmap` is that accession's slice of ``HUMAN_9606_idmapping.dat`` --
    a dict of database name -> list of identifiers. `ensg` is the resolved
    Ensembl gene id (or None). `census` is `load_rbp_census`'s result.
    """
    _protein_name, gene_symbol = parse_fasta_header(header)
    if not gene_symbol:
        gene_symbol = _first(idmap.get("Gene_Name"))

    refseq = idmap.get("RefSeq") or []
    ensp = _first(idmap.get("Ensembl_PRO"))

    unique_value, description = census.get(gene_symbol, (None, None))

    ensg_list = [ensg] if ensg else []

    return {
        "uniprot_id": accession,
        "Dominant_Isoform": 1,
        "sequence": sequence,
        "UNIQUE": unique_value,
        "ProteinHGVS": ",".join(refseq) if refseq else NO_REFSEQ,
        "HGVSDescription": None,
        "ENSP": ensp,
        "ID": ";".join(ensg_list) if ensg_list else None,
        "Name": gene_symbol,
        # Strictly the census value, even though the FASTA protein name is right
        # there in `_protein_name`. Description is null on the 92% of existing
        # proteins outside the census, so filling it here would make one column
        # mean two different things depending on which batch a row came from.
        # Use `parse_fasta_header` directly if you want the UniProt name.
        "Description": description,
        # bookkeeping
        "isoform_number": 1,
        "ID_list": ensg_list,
    }


def _first(values):
    """First element of a list, or None. Mirrors the original 'take the first'."""
    if not values:
        return None
    return values[0]
