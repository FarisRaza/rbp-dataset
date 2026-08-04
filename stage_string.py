"""Stage 4 -- STRING partner data for the new proteins, plus existing-row deltas.

Produces two intermediates:

  string_new.pkl      complete partner dicts for the new proteins, and the
                      ENSP/UniProt sets of the *expanded* table
  string_reverse.pkl  for each existing ENSP, the newly-eligible partners it
                      gains -- what `append_to_master.py` needs to bring the
                      two ``_in_Dataframe`` columns up to date

The full graph is never held in memory. Only edges touching a new protein are
retained, which is sufficient for both jobs and keeps the footprint to a few
hundred megabytes rather than several gigabytes.
"""

import csv
import json
import pickle
import time

import paths
import string_ppi

csv.field_size_limit(1 << 30)


def main():
    start = time.time()
    # ENSP_clean varies *between isoform rows of the same protein*, so a
    # first-one-wins map would silently drop most of them -- about 28% of the
    # table's distinct ENSPs. Every one matters: they define which partners
    # count as "in the dataframe", and which existing rows can be updated.
    existing = {}          # uniprot_id -> set of ENSP_clean
    existing_ensps = set()
    with open(paths.MASTER_TABLE, newline="", encoding="utf-8",
              errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        i_uniprot = header.index("uniprot_id")
        i_ensp = header.index("ENSP_clean")
        for row in reader:
            accession = row[i_uniprot]
            ensp = row[i_ensp] or None
            bucket = existing.setdefault(accession, set())
            if ensp:
                bucket.add(ensp)
                existing_ensps.add(ensp)
    multi = sum(1 for v in existing.values() if len(v) > 1)
    print(f"existing proteins: {len(existing)}  ({time.time()-start:.0f}s)")
    print(f"  distinct ENSP_clean: {len(existing_ensps)}  "
          f"({multi} proteins carry more than one)")

    ids = json.load(open(paths.scratch("missing_ensg.json")))
    new_ensp = {a: string_ppi.clean_ensp(v) for a, v in ids["ensp"].items()}
    new_ensps = {e for e in new_ensp.values() if e}
    print(f"new proteins: {len(new_ensp)}, of which {len(new_ensps)} have an ENSP")

    all_ensps = existing_ensps | new_ensps
    all_uniprots = set(existing) | set(new_ensp)
    print(f"expanded table: {len(all_uniprots)} proteins, {len(all_ensps)} ENSPs")

    start = time.time()
    adjacency, reverse = string_ppi.scan(
        paths.STRING_LINKS, new_ensps, also_index_against=existing_ensps
    )
    print(f"STRING scan: {time.time()-start:.0f}s")
    print(f"  {sum(len(v) for v in adjacency.values())} edges for the new proteins")
    print(f"  {len(reverse)} existing ENSPs gain at least one partner")

    # Every ENSP a protein carries maps back to it, not just the first. Where
    # two proteins somehow share an ENSP the later one wins, matching the
    # original `set_index(...).to_dict()` behaviour.
    ensp_to_uniprot = {}
    for accession, ensps in existing.items():
        for ensp in ensps:
            ensp_to_uniprot[ensp] = accession
    for accession, ensp in new_ensp.items():
        if ensp:
            ensp_to_uniprot[ensp] = accession

    with open(paths.scratch("string_new.pkl"), "wb") as fh:
        pickle.dump({
            "adjacency": adjacency,
            "new_ensp": new_ensp,
            "ensp_to_uniprot": ensp_to_uniprot,
            "all_ensps": all_ensps,
            "all_uniprots": all_uniprots,
        }, fh, protocol=4)
    with open(paths.scratch("string_reverse.pkl"), "wb") as fh:
        pickle.dump(reverse, fh, protocol=4)
    print("wrote string_new.pkl and string_reverse.pkl")


if __name__ == "__main__":
    main()
