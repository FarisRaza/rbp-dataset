"""Fidelity check for the IDR family. Runs under paths.PYTHON_METAPREDICT.

Recomputes IDR_range and FOLD_range for proteins already in the master table and
compares against the stored boundaries. Invoked by ``validate.py --all``.
"""

import ast
import csv
import sys
import warnings

warnings.filterwarnings("ignore")

import idr
import paths
from validate import SAMPLE, load_sample

csv.field_size_limit(1 << 30)


def main():
    idr.install_python_fallback()
    rows, index = load_sample(SAMPLE)

    sequences = [rows[a][index["sequence"]] for a in rows]
    accessions = list(rows)

    passed = failed = 0
    for accession, result in zip(accessions, idr.predict(sequences)):
        for column in ("IDR_range", "FOLD_range"):
            stored = ast.literal_eval(rows[accession][index[column]])
            actual = result[column]
            if stored == actual:
                passed += 1
            else:
                failed += 1
                print(f"    FAIL {accession}.{column}")
                print(f"      stored     {stored}")
                print(f"      recomputed {actual}")

    print(f"  IDRs (metapredict): {passed} matched, {failed} differed  "
          f"[{'ok' if not failed else 'FAILED'}]")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
