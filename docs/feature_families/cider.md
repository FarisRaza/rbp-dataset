# CIDER

`annotate_cider.py` calculates whole-sequence localCIDER charge, composition,
hydropathy, and patterning metrics. It also applies the same calculations to
each metapredict IDR and each curated UniProt domain.

Inputs are the catalog sequence plus the IDR and domain sidecars; there is no
external data download. Noncanonical residues are sanitized only during
scoring, and the policy is recorded in the sidecar.

Sanity checks: one unique sidecar key per catalog row; scalar whole-sequence
metrics on every scorable sequence; per-IDR and per-domain list lengths agree
with their geometry sidecars; FCR lies in `[0, 1]` and NCPR in `[-1, 1]`.

The generated report plots column coverage and the distribution of whole-chain
NCPR. Inspect extreme values and any missing whole-sequence metrics.
