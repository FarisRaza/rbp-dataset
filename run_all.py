"""Driver: run every stage of the expansion in dependency order.

Each stage writes an intermediate into `paths.SCRATCH` and is skipped if that
intermediate already exists, so a failed run can be resumed without repeating
the expensive steps. Pass ``--force`` to rebuild everything.

Stages run under different interpreters because metapredict (torch) and PSLab
(scikit-learn 1.6) cannot share an environment; see `paths.py`.

    python run_all.py            # resume, skipping completed stages
    python run_all.py --force    # rebuild from scratch
    python run_all.py --list     # show stages and their status

Any other argument is forwarded to stage 1, which chooses the target set:

    python run_all.py --all                  # every reviewed human protein
    python run_all.py --accessions ids.txt   # a specific list
    python run_all.py --fasta my.fasta       # sequences from a FASTA

Passing a selector implies --force, since every later stage derives from stage
1's output and would otherwise be reused against the wrong protein set.
"""

import os
import subprocess
import sys

import paths

HERE = os.path.dirname(os.path.abspath(__file__))

# (name, interpreter, script, produces)
STAGES = [
    ("find missing proteins", sys.executable, "stage_find_missing.py",
     ["missing_proteins.json", "missing_idmap.json", "missing_ensg.json"]),
    ("IDRs (metapredict)", paths.PYTHON_METAPREDICT, "stage_idr.py",
     ["missing_idr.json"]),
    ("domains + GO (Swiss-Prot)", sys.executable, "stage_swissprot.py",
     ["missing_domains.json", "missing_go.json"]),
    ("STRING", sys.executable, "stage_string.py",
     ["string_new.pkl", "string_reverse.pkl"]),
    ("Open Targets", sys.executable, "stage_opentargets.py",
     ["missing_opentargets.pkl"]),
    ("PSLab", paths.PYTHON_PSLAB, "stage_pslab.py",
     ["missing_pslab.json"]),
]


def done(products):
    return all(os.path.exists(paths.scratch(p)) for p in products)


def main(argv):
    if "-h" in argv or "--help" in argv:
        # Without this, --help would be treated as a selector: it implies
        # --force, gets forwarded to stage 1, whose argparse prints help and
        # exits 0 -- after which the driver cheerfully rebuilds everything.
        print(__doc__)
        return 0

    force = "--force" in argv
    listing = "--list" in argv

    # Anything else is forwarded to stage 1, which is what decides the target
    # set (--all, --accessions FILE, --fasta FILE). Selecting a different set
    # implies rebuilding, since every later stage keys off stage 1's output.
    selector = [a for a in argv if a not in ("--force", "--list")]
    if selector:
        force = True

    if listing:
        for name, _, script, products in STAGES:
            print(f"  [{'x' if done(products) else ' '}] {name:28s} {script}")
        return 0

    for name, interpreter, script, products in STAGES:
        if done(products) and not force:
            print(f"skip  {name} (already built)")
            continue
        print(f"run   {name}  [{os.path.basename(interpreter)}]")
        arguments = selector if script == "stage_find_missing.py" else []
        if not os.path.exists(interpreter):
            # subprocess would raise FileNotFoundError before a return code
            # exists, so the failure message below would never be reached.
            print(f"FAILED at stage: {name}\n"
                  f"  interpreter not found: {interpreter}\n"
                  f"  create it with: python setup_environment.py")
            return 1
        result = subprocess.run(
            [interpreter, "-u", os.path.join(HERE, script)] + arguments, cwd=HERE
        )
        if result.returncode != 0:
            print(f"FAILED at stage: {name}")
            return result.returncode

    print("\nassembling new rows")
    if subprocess.run(
        [sys.executable, "-u", os.path.join(HERE, "build_rows.py"), paths.SCRATCH],
        cwd=HERE,
    ).returncode != 0:
        return 1

    print("\nwriting expanded master table")
    return subprocess.run(
        [sys.executable, "-u", os.path.join(HERE, "append_to_master.py"),
         paths.SCRATCH],
        cwd=HERE,
    ).returncode


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
