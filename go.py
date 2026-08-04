"""Gene Ontology feature family.

Source
------
The **GO cross-references inside ``uniprot_sprot.dat``** -- *not* ``goa_human.gaf``,
despite that file sitting in the project folder. The GAF is loaded by one of the
historical notebooks and then never used; a separate, newer pipeline consumes it
to build differently-named columns (``go_cc_ids`` / ``go_cc_names``) that do not
appear in this table.

Each Swiss-Prot ``DR   GO;`` line looks like::

    DR   GO; GO:0005634; C:nucleus; IDA:UniProtKB.

which ``Bio.SwissProt`` exposes as a cross_reference tuple::

    ("GO", "GO:0005634", "C:nucleus", "IDA:UniProtKB")

Produces schema.GO_COLUMNS (9 columns): for each of the three GO aspects
(C = cellular component, P = biological process, F = molecular function) a list
of term ids, a list of term names, and a list of evidence codes. The three lists
within an aspect are positionally parallel -- index *i* of ``C_ids``,
``C_descriptions`` and ``C_evidence`` all describe the same annotation.

Evidence strings keep their full ``CODE:SOURCE`` form (``'IDA:UniProtKB'``,
``'IBA:GO_Central'``). No evidence-code filtering is applied, matching the
historical pipeline -- IEA (electronic, uncurated) annotations are included.

Granularity
-----------
Per protein. In the master table these columns are blanked to NaN on rows with
``Dominant_Isoform == 0``. Rows added by this pipeline are all dominant.
"""

ASPECTS = {"C": "C", "P": "P", "F": "F"}

EMPTY = {
    "C_ids": [], "C_descriptions": [], "C_evidence": [],
    "P_ids": [], "P_descriptions": [], "P_evidence": [],
    "F_ids": [], "F_descriptions": [], "F_evidence": [],
}


def go_for_record(record):
    """Extract the 9 GO columns from a ``Bio.SwissProt`` record."""
    out = {k: list(v) for k, v in EMPTY.items()}

    for ref in record.cross_references:
        if not ref or ref[0] != "GO":
            continue
        go_id = ref[1] if len(ref) > 1 else None
        aspect_and_name = ref[2] if len(ref) > 2 else ""
        evidence = ref[3] if len(ref) > 3 else None

        aspect, _, description = aspect_and_name.partition(":")
        aspect = aspect.strip()
        if aspect not in ASPECTS:
            continue

        # The historical code used aspect_and_name.split(":")[1], which silently
        # truncates any GO term name containing a colon. partition() keeps the
        # whole name, which is the correct reading of the flat-file format.
        out[f"{aspect}_ids"].append(go_id)
        out[f"{aspect}_descriptions"].append(description.strip())
        out[f"{aspect}_evidence"].append(evidence)

    return out
