"""Stage 5 -- Open Targets rows for the ENSGs of the new proteins.

The three source files total roughly 2.9 GB and one field (``tissues``) can
exceed 50 KB per gene, so each is streamed and filtered to the wanted ENSGs
rather than loaded whole.
"""

import json
import pickle
import time

import opentargets as ot
import paths


def main():
    ids = json.load(open(paths.scratch("missing_ensg.json")))
    wanted = {e for e in ids["ensg"].values() if e}
    print(f"{len(wanted)} distinct ENSGs wanted")

    start = time.time()
    associations = ot.load_associations(
        paths.OT_ASSOCIATIONS, wanted)
    print(f"associations:      {len(associations):5d} genes  ({time.time()-start:.0f}s)")

    start = time.time()
    targets, target_columns = ot.load_targets(
        paths.OT_TARGETS, wanted)
    print(f"target annotation: {len(targets):5d} genes, {len(target_columns)} columns"
          f"  ({time.time()-start:.0f}s)")

    start = time.time()
    expression = ot.load_expression(
        paths.OT_EXPRESSION, wanted)
    print(f"expression:        {len(expression):5d} genes  ({time.time()-start:.0f}s)")

    with open(paths.scratch("missing_opentargets.pkl"), "wb") as fh:
        pickle.dump({
            "associations": associations,
            "expression": expression,
            "targets": targets,
            "target_columns": target_columns,
        }, fh, protocol=4)
    print("wrote missing_opentargets.pkl")


if __name__ == "__main__":
    main()
