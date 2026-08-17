"""Stage 3 -- domains and GO, from one streaming pass over uniprot_sprot.dat.

Both families read the same 4 GB flat file, so they share a single pass. The
scan stops as soon as every wanted accession has been seen, which in practice is
about two thirds of the way through.

Also asserts that the sequence in the flat file matches the one in the FASTA for
every accession. If those disagreed, the domain coordinates -- which index the
flat-file sequence -- would not address the sequence stored in the table.
"""

import json
import time
import warnings

warnings.filterwarnings("ignore")

import domains as domains_module
import go as go_module
import paths


def main():
    from Bio import SwissProt

    data = json.load(open(paths.scratch("missing_proteins.json")))
    sequences = data["seqs"]
    remaining = set(data["missing"])
    total_wanted = len(remaining)

    domain_rows, go_rows, disagreements = {}, {}, []
    scanned = 0
    start = time.time()

    with open(paths.UNIPROT_DAT, "r", encoding="utf-8", errors="replace") as fh:
        for record in SwissProt.parse(fh):
            scanned += 1
            if scanned % 200000 == 0:
                print(f"  scanned {scanned} records, {len(remaining)} still wanted",
                      flush=True)
            hits = remaining.intersection(record.accessions)
            if not hits:
                continue
            geometry = domains_module.domains_for_record(record)
            terms = go_module.go_for_record(record)
            for accession in hits:
                sequence = sequences.get(accession, record.sequence)
                if record.sequence != sequence:
                    disagreements.append(accession)
                domain_rows[accession] = domains_module.attach_sequences(
                    geometry, sequence)
                go_rows[accession] = terms
            remaining -= hits
            if not remaining:
                break

    json.dump(domain_rows, open(paths.scratch("missing_domains.json"), "w"))
    json.dump(go_rows, open(paths.scratch("missing_go.json"), "w"))

    print(f"done: scanned {scanned} records in {time.time()-start:.0f}s")
    print(f"  found {len(domain_rows)}/{total_wanted}; "
          f"{len(remaining)} absent from the flat file")
    print(f"  {sum(1 for v in domain_rows.values() if v['Domains'])} have "
          f"at least one annotated domain")
    print(f"  {sum(1 for v in go_rows.values() if v['C_ids'] or v['P_ids'] or v['F_ids'])}"
          f" have at least one GO term")
    print(f"  sequence disagreement between .dat and .fasta: {len(disagreements)}")
    if disagreements:
        print("   !! domain coordinates may not address the stored sequence:",
              disagreements[:10])


if __name__ == "__main__":
    main()
