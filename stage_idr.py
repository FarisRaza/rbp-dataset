"""Stage 2 -- metapredict IDR/FOLD geometry. Runs under paths.PYTHON_METAPREDICT."""

import json
import time
import warnings

warnings.filterwarnings("ignore")

import idr
import paths


def main():
    if idr.install_python_fallback():
        print("using metapredict's pure-Python domain builder "
              "(Cython extension unavailable)")

    data = json.load(open(paths.scratch("missing_proteins.json")))
    accessions = data["missing"]
    sequences = [data["seqs"][a] for a in accessions]

    out = {}
    start = time.time()
    for i, (accession, row) in enumerate(
        zip(accessions, idr.predict(sequences, batch_size=100))
    ):
        out[accession] = row
        if (i + 1) % 500 == 0:
            rate = (i + 1) / (time.time() - start)
            print(f"  {i+1}/{len(accessions)}  {rate:.1f} seq/s", flush=True)

    json.dump(out, open(paths.scratch("missing_idr.json"), "w"))

    with_idr = sum(1 for v in out.values() if v["IDR_count"])
    total = sum(v["IDR_count"] for v in out.values())
    print(f"done: {len(out)} proteins in {time.time()-start:.0f}s")
    print(f"  {with_idr} have at least one IDR; {total} IDRs in total")


if __name__ == "__main__":
    main()
