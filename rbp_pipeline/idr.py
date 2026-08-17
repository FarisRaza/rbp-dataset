"""IDR feature family -- intrinsically disordered / folded region geometry.

Source
------
`metapredict` v3 (Emenecker, Griffith & Holehouse), network version ``V3``.
Nothing is downloaded; disorder is predicted from sequence by a bundled model.

Produces schema.IDR_GEOMETRY_COLUMNS (12 columns): the count, average size,
total size, residue ranges, per-region sequences and concatenated sequence for
both the disordered (``IDR_*``) and the folded (``FOLD_*``) partition of each
protein.

Call parameters
---------------
The historical table was built with::

    metapredict.predict_disorder(seqs, version="V3", return_domains=True,
                                 disorder_threshold=0.5, force_disable_batch=True)

``disorder_threshold=0.5`` is the load-bearing choice -- it is *not* the
metapredict default of 0.42, and using the default shifts almost every boundary.
The remaining domain-carving parameters (minimum_IDR_size=12,
minimum_folded_domain=50, gap_closure=10) are left at their defaults.

Verified
--------
Recomputing IDR_range and FOLD_range for the first 8 rows of the master table
reproduces the stored boundaries exactly.

Windows note
------------
metapredict ships a Cython accelerator (``backend.cython.domain_definition``)
that needs MSVC to build. Where it is absent, `install_python_fallback` reroutes
the call to metapredict's own pure-Python implementation of the same routine.
metapredict's test suite (``tests/test_cython_v_python.py``) asserts the two
produce identical output, so this changes speed, not results.
"""

import numpy as np

# metapredict's alphabet is the 20 canonical amino acids. Unlike the CIDER path
# -- which deletes unknown residues -- every non-canonical character is
# SUBSTITUTED here, so that the sequence length is preserved and the returned
# residue indices still address the original sequence.
#
# This is the substitution table the historical notebooks used.
METAPREDICT_SUBSTITUTIONS = {
    "B": "N",  # Asx -> Asn
    "U": "C",  # selenocysteine -> cysteine
    "X": "G",  # unknown -> glycine
    "Z": "Q",  # Glx -> Gln
}

DISORDER_THRESHOLD = 0.5
NETWORK_VERSION = "V3"


def install_python_fallback():
    """Route metapredict's Cython entry point at its pure-Python twin.

    No-op when the compiled extension imported cleanly. Call once, before the
    first prediction.
    """
    from metapredict.backend import domain_definition as dd

    compiled = getattr(dd, "CYTHON_build_domains_from_values", None)
    # When the extension is missing, metapredict installs a stub that raises
    # ModuleNotFoundError; detect it by module of origin rather than by name.
    if compiled is not None and getattr(compiled, "__module__", "") not in (
        "metapredict.backend.domain_definition",
    ):
        return False  # real compiled extension present, leave it alone
    dd.CYTHON_build_domains_from_values = getattr(dd, "__build_domains_from_values")
    return True


def sanitize_for_metapredict(seq):
    """Substitute non-canonical residues, preserving length and indexing."""
    for src, dst in METAPREDICT_SUBSTITUTIONS.items():
        seq = seq.replace(src, dst)
    return seq


def _summarise(boundaries, sequence, prefix):
    """Turn a list of [start, end] boundaries into the six stored columns."""
    ranges = [[int(a), int(b)] for a, b in boundaries]
    seqs = [sequence[a:b] for a, b in ranges]
    sizes = [b - a for a, b in ranges]
    return {
        f"{prefix}_count": len(ranges),
        f"{prefix}_avg_size": float(np.mean(sizes)) if sizes else 0.0,
        f"{prefix}_total_size": int(sum(sizes)),
        f"{prefix}_range": ranges,
        f"{prefix}_discrete_seq": seqs,
        f"{prefix}_concat_seq": "".join(seqs),
    }


def predict(sequences, batch_size=200, show_progress=False):
    """Predict IDR/FOLD geometry for many sequences.

    `sequences` is an iterable of raw protein sequences. Yields one dict per
    input, in input order, holding all 12 schema.IDR_GEOMETRY_COLUMNS.

    Region sequences are sliced from the *sanitised* sequence, which is what
    the historical table stores -- so an X-containing protein will show a G at
    that position in IDR_discrete_seq. Boundaries index the original sequence
    identically, because sanitisation is length-preserving.
    """
    import metapredict as meta

    sequences = list(sequences)
    for start in range(0, len(sequences), batch_size):
        chunk = sequences[start : start + batch_size]
        clean = [sanitize_for_metapredict(s) for s in chunk]
        results = meta.predict_disorder(
            clean,
            version=NETWORK_VERSION,
            return_domains=True,
            disorder_threshold=DISORDER_THRESHOLD,
            force_disable_batch=True,
            show_progress_bar=show_progress,
        )
        for seq, res in zip(clean, results):
            row = {}
            row.update(_summarise(res.disordered_domain_boundaries, seq, "IDR"))
            row.update(_summarise(res.folded_domain_boundaries, seq, "FOLD"))
            yield row
