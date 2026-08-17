"""PSLab feature family -- predicted phase-separation propensity per IDR.

Source
------
The predictor from S. von Bülow, G. Tesei, F. K. Zaidi, T. Mittag and
K. Lindorff-Larsen, "Prediction of phase-separation propensities of disordered
proteins from sequence", PNAS (2025). Code and trained models:
https://github.com/KULL-Centre/_2024_buelow_PSpred

The project's ``isoform_idr_pslab_predictions.csv`` was produced by that repo's
``PSLab.ipynb`` Colab notebook, whose batch-prediction cell emits exactly the
13 columns of that file. This module reproduces that cell locally so the new
rows do not require a Colab round-trip.

Produces schema.PSLAB_COLUMNS (11 columns). Eight are sequence features
computed analytically; the last three come from two MLP regressors trained on
CALVADOS 2 slab simulations, at the model's fixed conditions of T = 293 K and
ionic strength 150 mM.

    mean_lambda   mean per-residue "stickiness"
    faro          fraction of aromatic residues
    shd           sequence hydropathy decoration
    ncpr          net charge per residue
    fcr           fraction of charged residues
    scd           sequence charge decoration
    ah_ij         summary of pairwise Ashbaugh-Hatch interaction potentials
    nu_svr        SVR-predicted Flory scaling exponent
    Delta G [kT]                       transfer free energy, dilute -> dense
    Saturation concentration [mg/mL]   exp(predicted log c_sat)
    Saturation concentration [uM]      the same, divided by molecular weight

Granularity
-----------
One prediction per IDR, not per protein. Each stored column is a list ordered to
match ``IDR_discrete_seq``. A protein with no IDRs stores an empty list in all
11 columns.

Verified
--------
Re-running the predictor on IDR ``Q9NQ94_1_1`` ("MESNHKSGDGLSGT") reproduces all
11 stored values, including the two ML predictions, to floating-point equality.

Environment
-----------
The shipped ``.joblib`` models were pickled under scikit-learn 1.6 and reference
classes that lived in ``__main__`` at pickling time; `load_models` re-establishes
those names before unpickling. This module therefore needs its own interpreter
(scikit-learn==1.6, MDAnalysis, numba, biopython, joblib) separate from the one
running metapredict.
"""

import os
import sys

FEATURES = ["mean_lambda", "faro", "shd", "ncpr", "fcr", "scd", "ah_ij", "nu_svr"]

# Model conditions are fixed by training; exposed here so callers can record them.
TEMPERATURE_K = 293
IONIC_STRENGTH_MM = 150
CHARGE_TERMINI = True


def _bootstrap(repo_root):
    """Put the predictor's own modules on the path and cd into their directory.

    ``sequence.py`` and ``predictor.py`` import each other by bare module name
    and read ``residues.csv`` / ``svr_model_nu.joblib`` by relative path, so the
    process must run from ``scripts_colab`` with those two data files alongside.
    """
    scripts = os.path.join(repo_root, "scripts_colab")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    for name in ("residues.csv",):
        target = os.path.join(scripts, name)
        if not os.path.exists(target):
            import shutil

            shutil.copy(os.path.join(repo_root, "data", name), target)
    nu = os.path.join(scripts, "svr_model_nu.joblib")
    if not os.path.exists(nu):
        import shutil

        shutil.copy(os.path.join(repo_root, "models", "svr_model_nu.joblib"), nu)
    os.chdir(scripts)
    return scripts


def _apply_patches():
    """Two fixes to the upstream ``sequence`` module, neither changing results.

    1. ``SeqFeatures.__init__`` eagerly computes ``kappa``, which divides by the
       count of charged residues and so raises ZeroDivisionError on an IDR that
       has none. ``kappa`` is not one of the eight model features, so the whole
       prediction was being lost to an unused attribute. Guarding it lets those
       IDRs produce their real values.

    2. ``SeqFeatures`` calls ``joblib.load(nu_file)`` on *every* sequence,
       re-reading the nu SVR from disk each time. Caching the load leaves the
       model identical and removes the dominant cost.
    """
    import functools
    import joblib
    import sequence as _sequence

    original_kappa = _sequence.calc_kappa_manual

    def safe_kappa(seq):
        try:
            return original_kappa(seq)
        except ZeroDivisionError:
            return float("nan")

    _sequence.calc_kappa_manual = safe_kappa
    _sequence.load = functools.lru_cache(maxsize=8)(joblib.load)


def load_models(repo_root):
    """Load the two MLP regressors, the nu SVR and the residue table.

    Returns (models, residues, nu_file) ready to pass to `predict`.
    """
    _bootstrap(repo_root)

    import __main__
    import joblib
    import pandas as pd
    import predictor as _predictor

    _apply_patches()

    # The pickles name ``Model`` and ``AttrSetter`` as ``__main__`` attributes.
    __main__.Model = _predictor.Model
    __main__.AttrSetter = _predictor.AttrSetter

    models = {
        "dG": joblib.load(os.path.join(repo_root, "models/idrome90/mlp/dG/model.joblib")),
        "logcdil_mgml": joblib.load(
            os.path.join(repo_root, "models/idrome90/mlp/logcdil_mgml/model.joblib")
        ),
    }
    residues = pd.read_csv("residues.csv").set_index("one")
    return models, residues, "svr_model_nu.joblib"


def predict(seqs, models, residues, nu_file):
    """Predict the 11 PSLab values for each sequence in `seqs`.

    Yields one dict per input sequence, in input order. Sequences shorter than
    a few residues make the decoration metrics degenerate; the predictor is
    still evaluated, matching the historical behaviour, which ran every IDR
    metapredict emitted (minimum IDR size 12).
    """
    import numpy as np
    import sequence as _sequence
    from predictor import X_from_seq

    for seq in seqs:
        try:
            feats = _sequence.SeqFeatures(
                seq, residues=residues, charge_termini=CHARGE_TERMINI, nu_file=nu_file
            )
            X = X_from_seq(
                seq,
                FEATURES,
                residues=residues,
                charge_termini=CHARGE_TERMINI,
                nu_file=nu_file,
            )
            row = {f: float(getattr(feats, f)) for f in FEATURES}

            dG = float(np.mean(models["dG"].predict(X)))
            log_csat = float(np.mean(models["logcdil_mgml"].predict(X)))
            csat_mgml = float(np.exp(log_csat))

            row["Delta G [kT]"] = dG
            row["Saturation concentration [mg/mL]"] = csat_mgml
            row["Saturation concentration [uM]"] = csat_mgml / feats.mw * 1e6
        except Exception:
            # One pathological IDR must not abort a batch of thousands; a null
            # row is recorded and counted by the caller.
            row = {c: None for c in FEATURES}
            row.update({
                "Delta G [kT]": None,
                "Saturation concentration [mg/mL]": None,
                "Saturation concentration [uM]": None,
            })
        yield row


EMPTY = {
    "mean_lambda": [],
    "faro": [],
    "shd": [],
    "ncpr": [],
    "fcr": [],
    "scd": [],
    "ah_ij": [],
    "nu_svr": [],
    "Delta G [kT]": [],
    "Saturation concentration [mg/mL]": [],
    "Saturation concentration [uM]": [],
}
