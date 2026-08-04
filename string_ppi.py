"""STRING feature family -- protein-protein interaction partners.

Source
------
``9606.protein.links.v12.0.txt`` -- STRING v12.0, human, full links file
(~13 M undirected edges, space-delimited, columns
``protein1 protein2 combined_score``). Identifiers carry a ``9606.`` taxon
prefix that is stripped.

Produces schema.STRING_COLUMNS (5 columns):

  ENSP_clean                        this protein's Ensembl protein id, unversioned
  PPI_ENSP_Partners                 {partner ENSP: combined_score}, all partners
  PPI_UniProt_Partners              the same, with partners translated to UniProt
  PPI_ENSP_Partners_in_Dataframe    restricted to partners present in this table
  PPI_UniProt_Partners_in_Dataframe the same, in UniProt space

No confidence cutoff is applied here. STRING's own file floor is
``combined_score >= 150``, so that is the effective threshold, inherited rather
than chosen.

Granularity
-----------
Unlike GO, Domains and CD-CODE, these columns are NOT blanked on non-dominant
isoform rows -- every row carrying an ENSP gets partners.

Which columns depend on the table's contents
--------------------------------------------
**Three** of the five, not two. ``PPI_UniProt_Partners`` looks absolute but is
not: it is built by translating partner ENSPs through a map that only covers
proteins in the table, so an untranslatable partner is dropped. The master's own
cardinalities show it -- for P49588, ``PPI_ENSP_Partners`` holds 1,403 entries
while ``PPI_UniProt_Partners`` and both ``_in_Dataframe`` columns hold 881.

Only ``ENSP_clean`` and ``PPI_ENSP_Partners`` are independent of the table.

Adding proteins therefore changes those three columns for **existing** rows as
well. The protein set only ever grows, so an existing row's dicts can only gain
entries, never lose them. `scan` exploits that: instead of rebuilding the whole
13 M-edge graph (which does not fit in memory here), it collects only the edges
that touch a newly-added protein, which is all the information an existing row
needs in order to be updated.
"""

LINKS_FILE = "9606.protein.links.v12.0.txt"
TAXON_PREFIX = "9606."


def strip_taxon(identifier):
    """``'9606.ENSP00000363105'`` -> ``'ENSP00000363105'``."""
    if identifier.startswith(TAXON_PREFIX):
        return identifier[len(TAXON_PREFIX):]
    return identifier


def clean_ensp(ensp):
    """Drop the version suffix: ``'ENSP00000363105.1'`` -> ``'ENSP00000363105'``."""
    if not ensp:
        return None
    return str(ensp).split(".")[0]


def scan(links_path, focus_ensps, also_index_against=None):
    """One pass over the links file.

    `focus_ensps` is the set of ENSPs whose *complete* partner lists are needed
    -- in practice, the newly-added proteins.

    Returns ``(adjacency, reverse)`` where

      adjacency[e] = {partner: score}   for every e in `focus_ensps`
      reverse[o]   = {f: score}         for every ENSP `o` in
                                        `also_index_against` that partners with
                                        some f in `focus_ensps`

    `reverse` is what existing table rows need in order to extend their
    ``_in_Dataframe`` dicts, and it is small: only edges touching a new protein
    are retained, rather than the full graph.
    """
    focus = set(focus_ensps)
    against = set(also_index_against) if also_index_against else None

    adjacency = {e: {} for e in focus}
    reverse = {}

    with open(links_path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        if "protein1" not in header:
            fh.seek(0)
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            a = strip_taxon(parts[0])
            b = strip_taxon(parts[1])
            a_hit = a in focus
            b_hit = b in focus
            if not (a_hit or b_hit):
                continue
            score = int(parts[2])
            if a_hit:
                adjacency[a][b] = score
                if against is not None and b in against:
                    reverse.setdefault(b, {})[a] = score
            if b_hit:
                adjacency[b][a] = score
                if against is not None and a in against:
                    reverse.setdefault(a, {})[b] = score

    return adjacency, reverse


def columns_for(ensp_versioned, adjacency, ensp_to_uniprot, table_ensps, table_uniprots):
    """Build the 5 STRING columns for one protein.

    `ensp_to_uniprot` should be built over the *expanded* table so that partner
    translation sees the newly-added proteins too.
    """
    clean = clean_ensp(ensp_versioned)
    if not clean:
        # Matches the historical encoding: no ENSP means PPI_ENSP_Partners is
        # empty (NaN in the source table) while the other three are empty dicts.
        return {
            "ENSP_clean": None,
            "PPI_ENSP_Partners": None,
            "PPI_UniProt_Partners": {},
            "PPI_ENSP_Partners_in_Dataframe": {},
            "PPI_UniProt_Partners_in_Dataframe": {},
        }

    partners = adjacency.get(clean)
    if partners is None:
        return {
            "ENSP_clean": clean,
            "PPI_ENSP_Partners": None,
            "PPI_UniProt_Partners": {},
            "PPI_ENSP_Partners_in_Dataframe": {},
            "PPI_UniProt_Partners_in_Dataframe": {},
        }

    uniprot_partners = {
        ensp_to_uniprot[p]: s for p, s in partners.items() if p in ensp_to_uniprot
    }
    return {
        "ENSP_clean": clean,
        "PPI_ENSP_Partners": partners,
        "PPI_UniProt_Partners": uniprot_partners,
        "PPI_ENSP_Partners_in_Dataframe": {
            p: s for p, s in partners.items() if p in table_ensps
        },
        "PPI_UniProt_Partners_in_Dataframe": {
            u: s for u, s in uniprot_partners.items() if u in table_uniprots
        },
    }


def extend_table_relative(uniprot_partners, ensp_in_dataframe, uniprot_in_dataframe,
                          gained, ensp_to_uniprot):
    """Fold newly-eligible partners into an existing row's three table-relative dicts.

    `gained` is ``reverse[this_row_ensp]`` from `scan` -- the partners that only
    became eligible because they were just added to the table.

    Returns the three updated dicts in the order
    ``(PPI_UniProt_Partners, PPI_ENSP_Partners_in_Dataframe,
    PPI_UniProt_Partners_in_Dataframe)``. ``PPI_ENSP_Partners`` is not passed in
    because it lists every partner regardless of table membership and so never
    changes.
    """
    uniprot_out = dict(uniprot_partners or {})
    ensp_in_out = dict(ensp_in_dataframe or {})
    uniprot_in_out = dict(uniprot_in_dataframe or {})

    for partner, score in (gained or {}).items():
        ensp_in_out[partner] = score
        mapped = ensp_to_uniprot.get(partner)
        if mapped is not None:
            uniprot_out[mapped] = score
            uniprot_in_out[mapped] = score

    return uniprot_out, ensp_in_out, uniprot_in_out
