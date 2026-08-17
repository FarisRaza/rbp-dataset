"""CIDER feature family -- sequence-composition and charge-patterning metrics.

Source
------
`localcider` (Holehouse et al.), class ``localcider.sequenceParameters.SequenceParameters``.
Nothing is downloaded; every value is computed from the amino-acid sequence alone.

This module reproduces three blocks of the master table:

  * whole-sequence CIDER              -> schema.CIDER_COLUMNS         (19 cols)
  * CIDER evaluated per IDR           -> schema.IDR_CIDER_COLUMNS     (19 cols)
  * CIDER evaluated per domain region -> schema.DOMAIN_CIDER_COLUMNS  (20 cols)

The three blocks are not the same metric list. The per-domain block additionally
carries ``Omega``; the whole-sequence block does not. The per-IDR and per-domain
blocks also order ``kappa`` differently from the whole-sequence block. Those
differences are historical, but they are load-bearing for column alignment, so
each block gets its own explicit getter list below.

Verified
--------
Recomputing the whole-sequence block for Q9NQ94 (row 0 of the master table)
reproduces all 18 scalar columns to floating-point equality.

Sequence sanitisation
---------------------
localcider only accepts the 20 canonical amino acids. The historical notebooks
mapped ``U`` (selenocysteine) to ``C`` and deleted ``X``. That exact convention
is reproduced here by `sanitize_for_cider`; see `idr.py` for the *different*
convention metapredict was run under, which is deliberately not shared.
"""

from localcider.sequenceParameters import SequenceParameters

# Whole-sequence block. Note kappa/delta/deltaMax trail PPII_propensity here.
WHOLE_SEQUENCE_GETTERS = [
    ("FCR", "get_FCR"),
    ("NCPR", "get_NCPR"),
    ("isoelectric_point", "get_isoelectric_point"),
    ("molecular_weight", "get_molecular_weight"),
    ("countNeg", "get_countNeg"),
    ("countPos", "get_countPos"),
    ("countNeut", "get_countNeut"),
    ("fraction_negative", "get_fraction_negative"),
    ("fraction_positive", "get_fraction_positive"),
    ("fraction_expanding", "get_fraction_expanding"),
    ("amino_acid_fractions", "get_amino_acid_fractions"),
    ("fraction_disorder_promoting", "get_fraction_disorder_promoting"),
    ("mean_net_charge", "get_mean_net_charge"),
    ("mean_hydropathy", "get_mean_hydropathy"),
    ("uversky_hydropathy", "get_uversky_hydropathy"),
    ("PPII_propensity", "get_PPII_propensity"),
    ("kappa", "get_kappa"),
    ("delta", "get_delta"),
    ("deltaMax", "get_deltaMax"),
]

# Per-IDR block: same 19 metrics, but kappa sits directly after
# fraction_disorder_promoting.
PER_IDR_GETTERS = [
    ("FCR", "get_FCR"),
    ("NCPR", "get_NCPR"),
    ("isoelectric_point", "get_isoelectric_point"),
    ("molecular_weight", "get_molecular_weight"),
    ("countNeg", "get_countNeg"),
    ("countPos", "get_countPos"),
    ("countNeut", "get_countNeut"),
    ("fraction_negative", "get_fraction_negative"),
    ("fraction_positive", "get_fraction_positive"),
    ("fraction_expanding", "get_fraction_expanding"),
    ("amino_acid_fractions", "get_amino_acid_fractions"),
    ("fraction_disorder_promoting", "get_fraction_disorder_promoting"),
    ("kappa", "get_kappa"),
    ("mean_net_charge", "get_mean_net_charge"),
    ("mean_hydropathy", "get_mean_hydropathy"),
    ("uversky_hydropathy", "get_uversky_hydropathy"),
    ("PPII_propensity", "get_PPII_propensity"),
    ("delta", "get_delta"),
    ("deltaMax", "get_deltaMax"),
]

# Per-domain block: the per-IDR list plus Omega, inserted after kappa.
PER_DOMAIN_GETTERS = (
    PER_IDR_GETTERS[:13]
    + [("Omega", "get_Omega")]
    + PER_IDR_GETTERS[13:]
)

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sanitize_for_cider(seq):
    """Coerce `seq` into localcider's 20-letter alphabet.

    ``U`` (selenocysteine) is substituted with ``C``, its closest chemical
    analogue; every other non-canonical character (``X``, ``B``, ``Z``, ``O``)
    is deleted. Deleting rather than substituting keeps composition fractions
    honest -- an unknown residue should not be counted as a specific one -- at
    the cost of shifting downstream residue indices, which is why this function
    is only ever applied to sequences that are about to be scored, never to
    sequences whose coordinates are still needed.
    """
    if not seq:
        return ""
    return "".join(c for c in seq.replace("U", "C") if c in CANONICAL_AA)


def compute_block(seq, getters):
    """Evaluate `getters` against one sequence.

    Returns a dict keyed by the short metric name (``"FCR"``, not
    ``"get_FCR"``). Returns None for every metric when the sanitised sequence
    is empty or is a single residue, since localcider's patterning metrics
    (kappa, delta, Omega) are undefined there and raise.
    """
    clean = sanitize_for_cider(seq)
    if len(clean) < 2:
        return {name: None for name, _ in getters}

    sp = SequenceParameters(clean)
    out = {}
    for name, getter in getters:
        try:
            out[name] = getattr(sp, getter)()
        except Exception:
            # kappa/Omega are undefined for sequences with no charged residues
            # or no proline; localcider signals this by raising. A null is the
            # honest encoding, and matches what the historical table stores.
            out[name] = None
    return out


def whole_sequence(seq):
    """CIDER for a full protein sequence -> dict over schema.CIDER_COLUMNS."""
    return compute_block(seq, WHOLE_SEQUENCE_GETTERS)


def per_region(seqs, getters):
    """CIDER for an ordered collection of subsequences.

    Returns a dict of metric-name -> list, with one list entry per input
    sequence, preserving input order. This is the shape the ``IDR_*`` columns
    use: a list parallel to ``IDR_discrete_seq``.
    """
    blocks = [compute_block(s, getters) for s in seqs]
    return {name: [b[name] for b in blocks] for name, _ in getters}


def per_idr(idr_seqs):
    """CIDER per IDR -> dict of lists over schema.IDR_CIDER_COLUMNS."""
    return per_region(idr_seqs, PER_IDR_GETTERS)


def per_domain(domain_seq_map):
    """CIDER per domain region.

    `domain_seq_map` maps a domain name to the ordered list of that domain's
    discrete region sequences, e.g. ``{"RRM": [seq1, seq2, seq3]}``.

    Returns a dict of metric-name -> ``{domain_name: [value_per_region]}``,
    matching how the ``Domains_*`` columns are stored: a dict keyed by domain
    name whose values are lists parallel to ``Domains_discrete_seq[name]``.
    """
    out = {name: {} for name, _ in PER_DOMAIN_GETTERS}
    for domain_name, seqs in domain_seq_map.items():
        block = per_region(seqs, PER_DOMAIN_GETTERS)
        for metric, values in block.items():
            out[metric][domain_name] = values
    return out
