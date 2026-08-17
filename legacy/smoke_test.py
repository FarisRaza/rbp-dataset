"""End-to-end acceptance test on a handful of proteins.

Runs every stage for a small target set in an isolated scratch directory, then
checks the assembled rows. Takes a couple of minutes rather than the ~40 that a
full run needs, so it is the fastest way to confirm a fresh install works.

    python smoke_test.py

It writes nothing into the project folder and does not disturb a full run's
intermediates: both the scratch directory and the output file are temporary.
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile

csv.field_size_limit(1 << 30)

HERE = os.path.dirname(os.path.abspath(__file__))

# A deliberately mixed set: a well-studied RBP with domains, condensates and
# many interaction partners; two proteins that are in the master table; and one
# selenoprotein, which exercises the non-canonical-residue path.
ACCESSIONS = ["Q9NQ94", "P49588", "P31946", "P18283", "O75526"]

EXPECT_NON_EMPTY = [
    "uniprot_id", "sequence", "FCR", "kappa",
    "IDR_count", "IDR_range", "C_ids", "Name",
]


def run(command, env, label):
    print(f"\n--- {label}")
    result = subprocess.run(command, cwd=HERE, env=env)
    if result.returncode != 0:
        raise SystemExit(f"FAILED: {label}")


def main():
    scratch = tempfile.mkdtemp(prefix="isoform_smoke_")
    out_path = os.path.join(scratch, "smoke_rows.csv")
    ids_path = os.path.join(scratch, "accessions.txt")
    with open(ids_path, "w") as fh:
        fh.write("\n".join(ACCESSIONS) + "\n")

    env = dict(os.environ, EXPANSION_SCRATCH=scratch)
    import paths

    print(f"scratch: {scratch}")
    print(f"target : {len(ACCESSIONS)} accessions")

    try:
        run([sys.executable, "-u", "stage_find_missing.py", "--accessions", ids_path],
            env, "stage 1: target set and identifiers")
        run([paths.PYTHON_METAPREDICT, "-u", "stage_idr.py"], env, "stage 2: IDRs")
        run([sys.executable, "-u", "stage_swissprot.py"], env, "stage 3: domains + GO")
        run([sys.executable, "-u", "stage_string.py"], env, "stage 4: STRING")
        run([sys.executable, "-u", "stage_opentargets.py"], env, "stage 5: Open Targets")
        run([paths.PYTHON_PSLAB, "-u", "stage_pslab.py"], env, "stage 6: PSLab")
        run([sys.executable, "-u", "build_rows.py", scratch, out_path],
            env, "assembly")

        print("\n--- checking output")
        import schema

        with open(out_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            rows = list(reader)

        problems = []
        if header != schema.COLUMNS:
            problems.append(f"header has {len(header)} columns, expected 157")
        if len(rows) != len(ACCESSIONS):
            problems.append(f"{len(rows)} rows, expected {len(ACCESSIONS)}")
        for row in rows:
            if len(row) != len(schema.COLUMNS):
                problems.append(f"row width {len(row)}, expected {len(schema.COLUMNS)}")

        by_accession = {r[0]: dict(zip(header, r)) for r in rows}
        for accession in ACCESSIONS:
            record = by_accession.get(accession)
            if record is None:
                problems.append(f"{accession} missing from output")
                continue
            for column in EXPECT_NON_EMPTY:
                if not record[column].strip():
                    problems.append(f"{accession}.{column} is empty")

        print(f"  {len(rows)} rows x {len(header)} columns")
        for accession in ACCESSIONS:
            record = by_accession.get(accession, {})
            print(f"  {accession}  {record.get('Name', '?'):10s} "
                  f"IDRs={record.get('IDR_count', '?'):>3s}  "
                  f"domains={record.get('Domains', '?')[:34]}")

        if problems:
            print(f"\n{len(problems)} PROBLEM(S):")
            for problem in problems:
                print(f"  - {problem}")
            return 1

        print("\nSMOKE TEST PASSED")
        return 0
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
