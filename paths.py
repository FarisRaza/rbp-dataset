"""Locations, interpreters and source files.

Nothing here is machine-specific by default: the data directory is taken to be
the parent of this package, and intermediates land in a temp directory (they are
bulky and regenerable, and the project folder is OneDrive-synced). Point the pipeline somewhere else with environment variables rather than by
editing this file:

    KAPPEL_DIR          directory holding the master table and source files
    EXPANSION_SCRATCH   where intermediates are written
    PSPRED_REPO         clone of KULL-Centre/_2024_buelow_PSpred
    PYTHON_METAPREDICT  interpreter with metapredict installed
    PYTHON_PSLAB        interpreter with scikit-learn 1.6 + MDAnalysis

`setup_environment.py` creates the two interpreters and the PSpred clone in the
default locations, so in the normal case none of these need to be set.

Call `check()` to find out what is missing before running anything expensive.
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

#: Directory holding the master table and every source file.
KAPPEL = os.environ.get("KAPPEL_DIR", os.path.dirname(HERE))

#: Where per-family intermediates are written. These run to a few hundred MB and
#: are fully regenerable, so they default to the system temp directory rather
#: than next to the code -- the project folder is OneDrive-synced and there is no
#: reason to upload them.
SCRATCH = os.environ.get(
    "EXPANSION_SCRATCH",
    os.path.join(tempfile.gettempdir(), "isoform_expansion_intermediates"),
)

#: Clone of https://github.com/KULL-Centre/_2024_buelow_PSpred
PSPRED_REPO = os.environ.get("PSPRED_REPO", os.path.join(HERE, "vendor", "PSpred"))

#: Per-environment interpreters. metapredict (torch) and PSLab (scikit-learn 1.6)
#: cannot share one environment.
#:
#: These deliberately live under a SHORT path rather than inside the package.
#: torch ships files nested ~160 characters deep, and Windows' 260-character
#: MAX_PATH is measured from the drive letter -- installing into a project that
#: already sits under a long path (a synced OneDrive folder, say) fails with
#: "WinError 206: The filename or extension is too long". Hence "isoenvs/mp"
#: rather than "<package>/envs/metapredict".
ENV_ROOT = os.environ.get(
    "EXPANSION_ENVS", os.path.join(tempfile.gettempdir(), "isoenvs")
)
_BIN = "Scripts" if os.name == "nt" else "bin"
_EXE = "python.exe" if os.name == "nt" else "python"

PYTHON_METAPREDICT = os.environ.get(
    "PYTHON_METAPREDICT", os.path.join(ENV_ROOT, "mp", _BIN, _EXE)
)
PYTHON_PSLAB = os.environ.get(
    "PYTHON_PSLAB", os.path.join(ENV_ROOT, "ps", _BIN, _EXE)
)

# ---------------------------------------------------------------- outputs ---
MASTER_TABLE = os.path.join(KAPPEL, "Isoform_Post_Merge_PSLab_OpenTargets.csv")
NEW_ROWS = os.path.join(KAPPEL, "Missing_Proteins_Rows.csv")
EXPANDED_TABLE = os.path.join(
    KAPPEL, "Isoform_Post_Merge_PSLab_OpenTargets_Expanded.csv"
)
EXTENDED_TABLE = os.path.join(
    KAPPEL, "Isoform_Post_Merge_PSLab_OpenTargets_Extended.csv"
)
RCSB_ENRICHED_TABLE = os.path.join(
    KAPPEL, "Isoform_Post_Merge_PSLab_OpenTargets_Extended_RCSB.csv"
)
RCSB_SUMMARY = os.path.join(KAPPEL, "Human_Proteome_RCSB_PDB_Summary.csv")
RCSB_ELEMENTS = os.path.join(
    KAPPEL, "Human_Proteome_RCSB_Secondary_Structure.csv.gz"
)
RCSB_METADATA = os.path.join(
    KAPPEL, "Human_Proteome_RCSB_Secondary_Structure.metadata.json"
)

# ---------------------------------------------------------------- sources ---
UNIPROT_FASTA = os.path.join(KAPPEL, "uniprot_sprot.fasta")
UNIPROT_DAT = os.path.join(KAPPEL, "uniprot_sprot.dat")
IDMAPPING = os.path.join(KAPPEL, "HUMAN_9606_idmapping.dat")
HGNC = os.path.join(KAPPEL, "hgnc_complete_set_2024-10-01.txt")
RBP_CENSUS = os.path.join(KAPPEL, "RBPs.xlsx")
STRING_LINKS = os.path.join(KAPPEL, "9606.protein.links.v12.0.txt")
CDCODE_DIR = os.path.join(KAPPEL, "CD-CODE_Files")
CONDENSATE_TABLE = os.path.join(KAPPEL, "condensates.csv")
CONDENSATE_EXTRA = os.path.join(KAPPEL, "MISSING.csv")
OT_ASSOCIATIONS = os.path.join(KAPPEL, "association_by_datatype_direct_full.csv")
OT_EXPRESSION = os.path.join(KAPPEL, "expression_all.csv")
OT_TARGETS = os.path.join(KAPPEL, "target_full.csv")
OT_DISEASES = os.path.join(KAPPEL, "opentargets_disease_25.12.parquet")
RCSB_SIFTS = os.path.join(SCRATCH, "pdb_chain_uniprot.tsv.gz")

#: Which source files each feature family needs. Used by `check()` so a run
#: fails immediately with a useful message instead of hours in.
REQUIREMENTS = {
    "identity": [UNIPROT_FASTA, IDMAPPING, HGNC, RBP_CENSUS, OT_TARGETS],
    "cider": [],
    "idr": [],
    "domains": [UNIPROT_DAT],
    "go": [UNIPROT_DAT],
    "string": [STRING_LINKS, MASTER_TABLE],
    "cdcode": [CDCODE_DIR, CONDENSATE_TABLE, CONDENSATE_EXTRA],
    "opentargets": [OT_ASSOCIATIONS, OT_EXPRESSION, OT_TARGETS],
    "pslab": [PSPRED_REPO],
    "rcsb": [UNIPROT_FASTA],
}

#: Where to obtain each source file, for the error message when one is absent.
PROVENANCE = {
    UNIPROT_FASTA: "UniProt: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    UNIPROT_DAT: "UniProt: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
    IDMAPPING: "UniProt: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
    HGNC: "HGNC: https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2024-10-01.txt",
    STRING_LINKS: "STRING v12.0: https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
    CDCODE_DIR: "CD-CODE (https://cd-code.org) -- per-condensate member tables, downloaded by hand",
    CONDENSATE_TABLE: "CD-CODE (https://cd-code.org) -- the condensate summary table",
    RBP_CENSUS: "project file, NOT downloadable: the 1,393-gene published RBP "
                "census (Baltz 2012 / Castello 2012 / Kramer 2014 / Beckmann "
                "2015). Only supplies UNIQUE and Description; without it those "
                "two columns are null, which is already true of 92% of rows.",
    CONDENSATE_EXTRA: "project file, NOT downloadable: a single hand-made row "
                      "(CCNT1 / O60563) inserted at condensate index 228. "
                      "Without it every CD-CODE assignment after 228 shifts.",
    OT_ASSOCIATIONS: "Open Targets Platform export (associationByDatatypeDirect)",
    OT_EXPRESSION: "Open Targets Platform export (expression)",
    OT_TARGETS: "Open Targets Platform export (target)",
    OT_DISEASES: "Open Targets Platform 25.12 disease metadata: "
                 "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/disease/disease.parquet",
    PSPRED_REPO: "git clone https://github.com/KULL-Centre/_2024_buelow_PSpred  (or run setup_environment.py)",
    RCSB_SIFTS: "SIFTS weekly mapping: https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz",
}


def scratch(name):
    """Absolute path to an intermediate, creating the directory on first use."""
    os.makedirs(SCRATCH, exist_ok=True)
    return os.path.join(SCRATCH, name)


def check(families=None, verbose=True):
    """Report which inputs and interpreters are present.

    Returns a list of human-readable problems; empty means ready to run.
    Pass `families` to check only some of them, e.g. ``check(["domains"])``.
    """
    wanted = families or list(REQUIREMENTS)
    problems = []

    if verbose:
        print(f"data directory : {KAPPEL}")
        print(f"intermediates  : {SCRATCH}")
        print()

    seen = set()
    for family in wanted:
        for item in REQUIREMENTS.get(family, []):
            if item in seen:
                continue
            seen.add(item)
            ok = os.path.exists(item)
            if verbose:
                print(f"  [{'x' if ok else ' '}] {os.path.basename(item)}")
            if not ok:
                problems.append(
                    f"missing {item}\n      obtain from: "
                    f"{PROVENANCE.get(item, 'unknown source')}"
                )

    if verbose:
        print()
    for label, interpreter, needed_by in (
        ("metapredict interpreter", PYTHON_METAPREDICT, "idr"),
        ("PSLab interpreter", PYTHON_PSLAB, "pslab"),
    ):
        if needed_by not in wanted:
            continue
        ok = os.path.exists(interpreter)
        if verbose:
            print(f"  [{'x' if ok else ' '}] {label}: {interpreter}")
        if not ok:
            problems.append(
                f"missing {label} at {interpreter}\n"
                f"      create it with: python setup_environment.py"
            )

    return problems


if __name__ == "__main__":
    problems = check()
    if problems:
        print(f"{len(problems)} problem(s):\n")
        for problem in problems:
            print(f"  - {problem}")
        sys.exit(1)
    print("all inputs present")
