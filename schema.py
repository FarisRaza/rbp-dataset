"""Canonical column schema of the master isoform table.

The reference file is ``Isoform_Post_Merge_PSLab_OpenTargets.csv`` (157 columns).
Every row appended by this pipeline must present exactly these columns, in this
order, so that the expanded table can be read by the existing downstream
notebooks without modification.

The column order below is not arbitrary -- it is the order in which the columns
were accreted by the historical pipeline:

    base identity                      cols   0-9
    CIDER on the whole sequence        cols  10-28
    Metapredict IDR / FOLD geometry    cols  29-40
    CIDER applied per IDR              cols  41-59
    UniProt domains + CIDER per domain cols  60-86
    STRING protein-protein interactions cols 87-91
    CD-CODE condensates                cols  92-101
    Gene Ontology                      cols 102-110
    bookkeeping + Open Targets         cols 111-145
    PSLab phase-separation predictions cols 146-156

``ID_list`` (col 112) is a derived helper column: the ``ID`` field split on
";". It is retained because the Open Targets columns are dicts keyed by the
ENSG values it contains.
"""

# --------------------------------------------------------------------------
# Column groups. Keep these in sync with COLUMNS below.
# --------------------------------------------------------------------------

IDENTITY_COLUMNS = [
    "uniprot_id",
    "Dominant_Isoform",
    "sequence",
    "UNIQUE",
    "ProteinHGVS",
    "HGVSDescription",
    "ENSP",
    "ID",
    "Name",
    "Description",
]

# localcider getters applied to the full sequence.
CIDER_COLUMNS = [
    "FCR",
    "NCPR",
    "isoelectric_point",
    "molecular_weight",
    "countNeg",
    "countPos",
    "countNeut",
    "fraction_negative",
    "fraction_positive",
    "fraction_expanding",
    "amino_acid_fractions",
    "fraction_disorder_promoting",
    "mean_net_charge",
    "mean_hydropathy",
    "uversky_hydropathy",
    "PPII_propensity",
    "kappa",
    "delta",
    "deltaMax",
]

IDR_GEOMETRY_COLUMNS = [
    "IDR_count",
    "IDR_avg_size",
    "IDR_total_size",
    "IDR_range",
    "IDR_discrete_seq",
    "IDR_concat_seq",
    "FOLD_count",
    "FOLD_avg_size",
    "FOLD_total_size",
    "FOLD_range",
    "FOLD_discrete_seq",
    "FOLD_concat_seq",
]

# The same localcider getters as CIDER_COLUMNS, evaluated once per IDR and
# stored as a list ordered to match IDR_discrete_seq.
#
# This is deliberately spelled out rather than derived as "IDR_" + CIDER_COLUMNS:
# the two blocks disagree on where kappa sits. Whole-sequence order ends
# ... PPII_propensity, kappa, delta, deltaMax, whereas the IDR (and Domains)
# blocks put kappa immediately after fraction_disorder_promoting.
IDR_CIDER_COLUMNS = [
    "IDR_FCR",
    "IDR_NCPR",
    "IDR_isoelectric_point",
    "IDR_molecular_weight",
    "IDR_countNeg",
    "IDR_countPos",
    "IDR_countNeut",
    "IDR_fraction_negative",
    "IDR_fraction_positive",
    "IDR_fraction_expanding",
    "IDR_amino_acid_fractions",
    "IDR_fraction_disorder_promoting",
    "IDR_kappa",
    "IDR_mean_net_charge",
    "IDR_mean_hydropathy",
    "IDR_uversky_hydropathy",
    "IDR_PPII_propensity",
    "IDR_delta",
    "IDR_deltaMax",
]

DOMAIN_GEOMETRY_COLUMNS = [
    "Domains",
    "Domains_count",
    "Domains_avg_size",
    "Domains_total_size",
    "Domains_range",
    "Domains_discrete_seq",
    "Domains_concat_seq",
]

# Note the extra "Omega" entry: the per-domain CIDER block carries one more
# metric than the whole-sequence and per-IDR blocks do.
DOMAIN_CIDER_COLUMNS = [
    "Domains_FCR",
    "Domains_NCPR",
    "Domains_isoelectric_point",
    "Domains_molecular_weight",
    "Domains_countNeg",
    "Domains_countPos",
    "Domains_countNeut",
    "Domains_fraction_negative",
    "Domains_fraction_positive",
    "Domains_fraction_expanding",
    "Domains_amino_acid_fractions",
    "Domains_fraction_disorder_promoting",
    "Domains_kappa",
    "Domains_Omega",
    "Domains_mean_net_charge",
    "Domains_mean_hydropathy",
    "Domains_uversky_hydropathy",
    "Domains_PPII_propensity",
    "Domains_delta",
    "Domains_deltaMax",
]

STRING_COLUMNS = [
    "ENSP_clean",
    "PPI_ENSP_Partners",
    "PPI_UniProt_Partners",
    "PPI_ENSP_Partners_in_Dataframe",
    "PPI_UniProt_Partners_in_Dataframe",
]

CDCODE_COLUMNS = [
    "Condensate Name",
    "UID",
    "Condensate Type",
    "Species Tax Id",
    "Proteins",
    "DNA",
    "RNA",
    "C-mods",
    "Condensatopathy",
    "Confidence Score",
]

GO_COLUMNS = [
    "C_ids",
    "C_descriptions",
    "C_evidence",
    "P_ids",
    "P_descriptions",
    "P_evidence",
    "F_ids",
    "F_descriptions",
    "F_evidence",
]

BOOKKEEPING_COLUMNS = [
    "isoform_number",
    "ID_list",
]

# Open Targets. The first five come from the association and expression
# datasets; the remaining 28 are every column of the target dataset except
# its "id" key.
OPENTARGETS_COLUMNS = [
    "diseaseId",
    "datatypeId",
    "score",
    "evidenceCount",
    "tissues",
    "approvedSymbol",
    "biotype",
    "transcriptIds",
    "canonicalTranscript",
    "canonicalExons",
    "genomicLocation",
    "alternativeGenes",
    "approvedName",
    "go",
    "hallmarks",
    "synonyms",
    "symbolSynonyms",
    "nameSynonyms",
    "functionDescriptions",
    "subcellularLocations",
    "targetClass",
    "obsoleteSymbols",
    "obsoleteNames",
    "constraint",
    "tep",
    "proteinIds",
    "dbXrefs",
    "chemicalProbes",
    "homologues",
    "tractability",
    "safetyLiabilities",
    "pathways",
    "tss",
]

PSLAB_COLUMNS = [
    "mean_lambda",
    "faro",
    "shd",
    "ncpr",
    "fcr",
    "scd",
    "ah_ij",
    "nu_svr",
    "Delta G [kT]",
    "Saturation concentration [mg/mL]",
    "Saturation concentration [uM]",
]

COLUMNS = (
    IDENTITY_COLUMNS
    + CIDER_COLUMNS
    + IDR_GEOMETRY_COLUMNS
    + IDR_CIDER_COLUMNS
    + DOMAIN_GEOMETRY_COLUMNS
    + DOMAIN_CIDER_COLUMNS
    + STRING_COLUMNS
    + CDCODE_COLUMNS
    + GO_COLUMNS
    + BOOKKEEPING_COLUMNS
    + OPENTARGETS_COLUMNS
    + PSLAB_COLUMNS
)

assert len(COLUMNS) == 157, len(COLUMNS)
assert len(set(COLUMNS)) == 157, "duplicate column name"


def validate_against(csv_path):
    """Raise if `csv_path`'s header differs from COLUMNS.

    Used as a guard before appending generated rows to the master table: a
    silent column-order drift would corrupt every appended row.
    """
    import csv

    csv.field_size_limit(1 << 30)
    with open(csv_path, newline="", encoding="utf-8", errors="replace") as fh:
        header = next(csv.reader(fh))
    if header != COLUMNS:
        only_file = [c for c in header if c not in set(COLUMNS)]
        only_here = [c for c in COLUMNS if c not in set(header)]
        raise ValueError(
            f"header mismatch for {csv_path}\n"
            f"  in file but not in schema: {only_file}\n"
            f"  in schema but not in file: {only_here}\n"
            f"  file has {len(header)} cols, schema has {len(COLUMNS)}"
        )
    return True
