# PSLab phase-separation predictions

`annotate_pslab.py` applies the published PSpred models to every metapredict IDR
and stores sequence features, predicted free energy, and saturation
concentrations in the same order as `IDR_discrete_seq`.

The family requires the isolated scikit-learn 1.6 environment and model files
installed by `setup_environment.py`; predictions are sequence- and IDR-specific.

Sanity checks: every output list length equals `IDR_count`; concentrations are
nonnegative when finite; rows with no IDR contain empty lists; model-loading
version constraints are satisfied.

The report plots output coverage and the number/distribution of available
per-IDR predictions. Compare free energies only under the same model release.
