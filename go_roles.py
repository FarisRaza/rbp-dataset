"""GO functional-role flags -- transcription, translation, mRNA stability.

Derived columns, not a new data source: these read the ``P_descriptions`` and
``F_descriptions`` already in the table (see `go.py`) and mark whether each
protein carries a term of the given kind. No download, no join key -- every row
that has GO terms can be flagged, existing rows included.

Origin
------
The keyword patterns come from the original functional-role session, which asked
two questions of the RBP table: which RBPs play no role in transcription, and
which control mRNA translation/stability. Both were answered by substring search
over the GO description columns, checking **both** biological process (``P_``)
and molecular function (``F_``) -- molecular-function terms like "translation
repressor activity" would be missed by looking at process alone.

That session produced one combined ``role_in_translation_stability`` column.
Translation regulation and mRNA stability are distinct biology, so they are
split here into separate columns, with the original combined flag retained so
earlier analyses stay reproducible.

Columns produced
----------------
    role_in_transcription           1 if any transcription-related P or F term
    role_in_translation             1 if any translation-regulation term
    role_in_mrna_stability          1 if any mRNA stabilisation/destabilisation term
    role_in_translation_stability   1 if either of the previous two -- the
                                    original session's combined definition

A caveat inherited from the method: this is a **substring match over term
names**, not a traversal of the GO DAG. A protein annotated only to a child term
whose name does not contain the keyword will be missed, and the flags describe
"has a term mentioning X", not "is established to regulate X". `matched_terms`
exists so the vocabulary driving each flag can be audited rather than trusted.
"""

import ast
import re

#: Any term naming transcription. Deliberately broad, matching the original.
#: Note this also catches "reverse transcription" and "transcription factor
#: binding" -- decide with `matched_terms` whether that is wanted.
TRANSCRIPTION = r"transcription"

#: Regulation of translation, in either direction.
TRANSLATION = (
    r"regulation of translation|"
    r"translational repress|"
    r"translational activat|"
    r"translation repressor activity|"
    r"translation activator activity"
)

#: mRNA stabilisation and destabilisation. "mRNA stab" also matches
#: "mRNA stability", "3'-UTR-mediated mRNA stabilization" and similar.
MRNA_STABILITY = r"mRNA stab|mRNA destabiliz"

ROLE_COLUMNS = [
    "role_in_transcription",
    "role_in_translation",
    "role_in_mrna_stability",
    "role_in_translation_stability",
]

#: The GO columns searched. Both aspects, for every flag -- the original session
#: switched to this after starting with process only.
SEARCH_COLUMNS = ["P_descriptions", "F_descriptions"]

_COMPILED = {
    "role_in_transcription": re.compile(TRANSCRIPTION, re.IGNORECASE),
    "role_in_translation": re.compile(TRANSLATION, re.IGNORECASE),
    "role_in_mrna_stability": re.compile(MRNA_STABILITY, re.IGNORECASE),
}


def parse_terms(cell):
    """Read one GO description cell into a list of term names.

    These columns hold a Python repr of a list (``"['mRNA binding', ...]"``).
    A cell that will not parse is returned as a single-element list, so a
    plain string still works.
    """
    if cell is None:
        return []
    if isinstance(cell, list):
        return cell
    text = str(cell).strip()
    if not text or text in ("[]", "nan"):
        return []
    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return [text]
    if isinstance(parsed, list):
        return [t for t in parsed if isinstance(t, str)]
    return [text] if isinstance(parsed, str) else []


def flags_for(row):
    """Return the four role columns for one row.

    `row` maps column name -> cell, and needs only the SEARCH_COLUMNS.
    Matching is done per parsed term rather than against the raw cell so that
    list punctuation cannot create a false match.
    """
    terms = []
    for column in SEARCH_COLUMNS:
        terms.extend(parse_terms(row.get(column)))

    out = {}
    for name, pattern in _COMPILED.items():
        out[name] = int(any(pattern.search(term) for term in terms))

    out["role_in_translation_stability"] = int(
        out["role_in_translation"] or out["role_in_mrna_stability"]
    )
    return out


def matched_terms(rows, which="role_in_transcription"):
    """Count the distinct GO terms driving one flag, most frequent first.

    Use this before trusting a flag: it shows exactly which term names the
    pattern is catching, which is how an over-broad keyword gets spotted.

    `rows` is an iterable of dicts (or of pandas rows).
    """
    import collections

    pattern = _COMPILED[which]
    counter = collections.Counter()
    for row in rows:
        for column in SEARCH_COLUMNS:
            for term in parse_terms(row.get(column)):
                if pattern.search(term):
                    counter[term] += 1
    return counter.most_common()
