"""Stage 6 -- PSLab phase-separation prediction, one row per IDR.

Runs under paths.PYTHON_PSLAB. IDRs are named with the historical scheme
``<uniprot>_<isoform_number>_<idr_number>``; every newly-added protein is a
single canonical isoform, so the isoform number is always 1.

Throughput is bounded by the model itself -- each prediction averages a
50-member cross-validation ensemble, twice (once for dG, once for log c_sat).
"""

import json
import time
import warnings

warnings.filterwarnings("ignore")

import paths
import pslab


def main():
    idr_rows = json.load(open(paths.scratch("missing_idr.json")))

    names, sequences = [], []
    for accession, row in idr_rows.items():
        for i, sequence in enumerate(row["IDR_discrete_seq"], start=1):
            names.append(f"{accession}_1_{i}")
            sequences.append(sequence)
    print(f"{len(sequences)} IDRs from {len(idr_rows)} proteins")

    models, residues, nu_file = pslab.load_models(paths.PSPRED_REPO)
    print("models loaded")

    out = {}
    start = time.time()
    for i, (name, row) in enumerate(
        zip(names, pslab.predict(sequences, models, residues, nu_file))
    ):
        out[name] = row
        if (i + 1) % 500 == 0:
            rate = (i + 1) / (time.time() - start)
            remaining = (len(sequences) - i - 1) / rate / 60
            print(f"  {i+1}/{len(sequences)}  {rate:.1f} idr/s  "
                  f"eta {remaining:.1f} min", flush=True)

    json.dump(out, open(paths.scratch("missing_pslab.json"), "w"))

    failed = sum(1 for v in out.values() if v["Delta G [kT]"] is None)
    print(f"done: {len(out)} IDRs in {time.time()-start:.0f}s")
    if failed:
        print(f"  {failed} IDRs could not be scored and are stored as nulls")

    # `predict` nulls a row rather than aborting a batch of thousands for one
    # pathological sequence. That is right for a handful, but a *wholesale*
    # failure -- wrong scikit-learn version, missing residues.csv, unreadable
    # models -- would otherwise exit 0 with every PSLab column silently empty.
    if out and failed / len(out) > 0.5:
        raise SystemExit(
            f"{failed}/{len(out)} IDRs failed to score. That is not a few odd "
            f"sequences; check the PSLab environment (scikit-learn must be 1.6) "
            f"and that {paths.PSPRED_REPO} is complete."
        )


if __name__ == "__main__":
    main()
