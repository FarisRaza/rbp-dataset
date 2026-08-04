"""Fidelity check for the PSLab family. Runs under paths.PYTHON_PSLAB.

Re-scores IDRs taken from ``isoform_idr_pslab_predictions.csv`` -- the output of
the original Colab run -- and compares against the values stored there. This
exercises the whole path including both MLP ensembles, so a match means the
local predictor reproduces the Colab predictor exactly.

Invoked by ``validate.py --all``.
"""

import csv
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import paths
import pslab

csv.field_size_limit(1 << 30)

PREDICTIONS = os.path.join(paths.KAPPEL, "isoform_idr_pslab_predictions.csv")
HOW_MANY = 25
TOLERANCE = 1e-6


def main():
    with open(PREDICTIONS, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        sample = [next(reader) for _ in range(HOW_MANY)]

    models, residues, nu_file = pslab.load_models(paths.PSPRED_REPO)
    columns = pslab.FEATURES + [
        "Delta G [kT]",
        "Saturation concentration [mg/mL]",
        "Saturation concentration [uM]",
    ]

    passed = failed = 0
    sequences = [row["Sequence"] for row in sample]
    for row, actual in zip(sample, pslab.predict(sequences, models, residues, nu_file)):
        for column in columns:
            stored = float(row[column])
            got = actual[column]
            if got is not None and abs(got - stored) <= TOLERANCE * max(1.0, abs(stored)):
                passed += 1
            else:
                failed += 1
                print(f"    FAIL {row['Name']}.{column}")
                print(f"      stored     {stored}")
                print(f"      recomputed {got}")

    print(f"  PSLab: {passed} matched, {failed} differed  "
          f"[{'ok' if not failed else 'FAILED'}]")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
