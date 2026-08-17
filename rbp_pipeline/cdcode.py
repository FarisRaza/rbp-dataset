"""CD-CODE feature family -- biomolecular condensate membership.

Source
------
Two hand-assembled sets of files downloaded from https://cd-code.org:

  * ``condensates.csv``     -- 327 rows, one per condensate, holding the ten
                               attributes that become the ten table columns.
  * ``CD-CODE_Files/download-data*.csv`` -- one member table per condensate,
                               listing the UniProt accessions in it.

There is no API call or download script anywhere in the project; these files
were fetched one condensate at a time through the web UI.

Produces schema.CDCODE_COLUMNS (10 columns). Each is a list with one element per
condensate the protein belongs to, and the ten lists are positionally parallel --
element *i* of ``Condensate Name``, ``UID``, ``Condensate Type`` etc. all
describe the same condensate.

The fragile part
----------------
A condensate's identity is carried *only* by the position of its member file in
a sorted listing -- the files do not record which condensate they belong to.
Reproducing the original ordering therefore requires reproducing two hand-made
index surgeries verbatim; `condensate_frames` does that and
`verify_alignment` checks the result against each condensate's declared member
count. Run the verification before trusting any output: renaming, adding or
re-downloading a single file silently shifts every condensate assignment after
it.

Note this fragility is about *building the condensate index*. Once built, the
protein lookup itself is a proper dictionary keyed on UniProt accession, so
adding new proteins is safe.

Numeric serialisation
---------------------
The historical columns were written from pandas scalars and so contain values
like ``[np.int64(9606)]``, which ``ast.literal_eval`` cannot parse. This module
casts to plain Python ints, so newly written rows are literal_eval-safe. That is
a deliberate deviation: it makes new rows *more* readable than old ones, and
does not change the values.
"""

import os
import re

MEMBER_DIR = "CD-CODE_Files"
CONDENSATE_TABLE = "condensates.csv"
EXTRA_ROW = "MISSING.csv"

# The ten columns of condensates.csv, in file order. "Name" is renamed to
# "Condensate Name" on the way into the table, because the table already has a
# "Name" column holding the gene symbol.
SOURCE_COLUMNS = [
    "Name", "UID", "Condensate Type", "Species Tax Id", "Proteins",
    "DNA", "RNA", "C-mods", "Condensatopathy", "Confidence Score",
]
RENAME = {"Name": "Condensate Name"}


def condensate_frames(root):
    """Reproduce the historical ordering of the 327 member tables.

    Files named ``download-data (N).csv`` sort by integer N; the remaining
    timestamp-named files sort lexicographically, which for this naming scheme
    is also chronological. Two manual moves then follow, exactly as performed
    originally:

      1. the two frames at indices 325..326 are lifted out and reinserted at 225
      2. a hand-made single-row frame (``MISSING.csv``, the CCNT1/O60563 entry)
         is inserted at index 228

    Returns a list of DataFrames whose index matches the row index of
    ``condensates.csv``.
    """
    import pandas as pd

    path = os.path.join(root, MEMBER_DIR)
    files = [f for f in os.listdir(path)
             if f.startswith("download") and f.endswith(".csv")]

    numbered = re.compile(r"download-data \((\d+)\)")
    with_number = [f for f in files if numbered.match(f)]
    with_number.sort(key=lambda f: int(numbered.match(f).group(1)))
    without_number = sorted(f for f in files if f not in set(with_number))
    ordered = with_number + without_number

    frames = [pd.read_csv(os.path.join(path, f)) for f in ordered]

    moved = frames[325:327]
    frames = frames[:325] + frames[327:]
    frames = frames[:225] + moved + frames[225:]

    frames.insert(228, pd.read_csv(os.path.join(root, EXTRA_ROW)))

    # The insert leaves one more frame than there are condensates: there are 327
    # member files and 327 rows in condensates.csv, so the hand-added row makes
    # 328. The surplus sits at the end and is a duplicate download -- the last
    # file is byte-identical to an earlier one, covering the same condensate.
    #
    # Left in place it is not harmless: build_index would map that condensate's
    # proteins (P17600/SYN1 and Q92777/SYN2) to index 327, and every lookup of
    # them would raise IndexError against the 327-row attribute table.
    #
    # Dropping a duplicate loses nothing, but dropping a *distinct* frame would
    # silently discard a real condensate, so that is checked rather than assumed.
    expected = len(pd.read_csv(os.path.join(root, CONDENSATE_TABLE)))
    while len(frames) > expected:
        surplus = frames[-1]
        if not any(_same_members(surplus, other) for other in frames[:-1]):
            raise ValueError(
                f"{len(frames)} member tables for {expected} condensates, and the "
                f"surplus one is not a duplicate of any other -- the file ordering "
                f"in {MEMBER_DIR} no longer matches {CONDENSATE_TABLE}. Every "
                f"condensate assignment after the first divergence would be wrong."
            )
        frames = frames[:-1]

    return frames


def _same_members(left, right):
    """Do two member tables list exactly the same proteins?"""
    if "UniProt ID" not in left.columns or "UniProt ID" not in right.columns:
        return False
    return set(left["UniProt ID"].dropna()) == set(right["UniProt ID"].dropna())


def verify_alignment(root, frames=None):
    """Check each member table's row count against its declared ``Proteins``.

    Returns a list of (index, declared, actual) triples for every disagreement;
    an empty list means the ordering reproduces the original. This is the only
    available check that the positional identity is intact.
    """
    import pandas as pd

    if frames is None:
        frames = condensate_frames(root)
    names = pd.read_csv(os.path.join(root, CONDENSATE_TABLE))

    problems = []

    # Check the counts agree first. Iterating only over `names` -- as an earlier
    # version did -- cannot see a surplus trailing frame, which is precisely the
    # failure that breaks lookups.
    if len(frames) != len(names):
        problems.append(
            (None, f"{len(names)} condensates", f"{len(frames)} member tables")
        )

    for i in range(len(names)):
        if i >= len(frames):
            problems.append((i, int(names.iloc[i]["Proteins"]), None))
            continue
        declared = int(names.iloc[i]["Proteins"])
        actual = len(frames[i])
        if declared != actual:
            problems.append((i, declared, actual))
    return problems


def build_index(root):
    """Return (protein_to_condensates, condensate_attributes).

    ``protein_to_condensates`` maps a UniProt accession to the list of
    condensate row-indices it appears in. ``condensate_attributes`` is the
    ``condensates.csv`` table as a list of plain-Python dicts.
    """
    import pandas as pd

    frames = condensate_frames(root)
    names = pd.read_csv(os.path.join(root, CONDENSATE_TABLE))

    protein_to_condensates = {}
    for index, frame in enumerate(frames):
        if "UniProt ID" not in frame.columns:
            continue
        for accession in frame["UniProt ID"].dropna().unique():
            protein_to_condensates.setdefault(accession, []).append(index)

    attributes = []
    for _, row in names.iterrows():
        attributes.append({c: _plain(row[c]) for c in SOURCE_COLUMNS})

    return protein_to_condensates, attributes


def _plain(value):
    """Cast a numpy scalar to its Python equivalent so repr() stays parseable."""
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, AttributeError):
            return value
    return value


def lookup(accession, protein_to_condensates, attributes):
    """Return the 10 CD-CODE columns for one accession.

    A protein in no condensate gets ten empty lists, which is what the table
    stores for annotated-but-condensate-free proteins.
    """
    indices = protein_to_condensates.get(accession, [])
    out = {}
    for source in SOURCE_COLUMNS:
        column = RENAME.get(source, source)
        out[column] = [attributes[i][source] for i in indices]
    return out


EMPTY = {RENAME.get(c, c): [] for c in SOURCE_COLUMNS}
