"""Build current RCSB PDB mappings and secondary structures for human Swiss-Prot.

The expensive network phase is resumable: every deterministic batch is cached
as gzip JSON under ``EXPANSION_SCRATCH/rcsb_secondary_batches``.  Re-running the
command fetches only absent or corrupt batches, then rewrites the final summary
and normalized element table atomically.

    python rbp_pipeline/stage_rcsb.py
    python rbp_pipeline/stage_rcsb.py --workers 6 --batch-size 500
    python rbp_pipeline/stage_rcsb.py --refresh-sifts
"""

import argparse
import json
import os

import paths
import rcsb


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", default=paths.UNIPROT_FASTA)
    p.add_argument("--sifts", default=paths.RCSB_SIFTS)
    p.add_argument("--summary", default=paths.RCSB_SUMMARY)
    p.add_argument("--elements", default=paths.RCSB_ELEMENTS)
    p.add_argument("--metadata", default=paths.RCSB_METADATA)
    p.add_argument(
        "--cache-dir", default=paths.scratch("rcsb_secondary_batches")
    )
    p.add_argument("--batch-size", type=int, default=500)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--refresh-sifts", action="store_true")
    p.add_argument("--force-fetch", action="store_true")
    p.add_argument(
        "--limit-accessions",
        type=int,
        help="testing only: keep the first N accessions (use explicit output paths)",
    )
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    rcsb.download_sifts(args.sifts, force=args.refresh_sifts)

    accessions = rcsb.human_swissprot_accessions(args.fasta)
    if args.limit_accessions:
        accessions = accessions[:args.limit_accessions]
    print(f"reviewed human proteins: {len(accessions):,}")

    sifts = rcsb.load_sifts(args.sifts, accessions)
    proteins_with_pdb = sum(
        bool(sifts.by_accession[acc]["pdb_ids"]) for acc in accessions
    )
    print(f"SIFTS release: {sifts.release}")
    print(f"proteins with PDB: {proteins_with_pdb:,}")
    print(f"unique PDB entries: {len(sifts.pdb_ids):,}")
    print(f"mapping segments: {sifts.mapping_segment_count:,}")
    if sifts.invalid_segment_count:
        print(f"invalid/incomplete mapping segments retained only for ID lists: "
              f"{sifts.invalid_segment_count:,}")

    batches = rcsb.fetch_secondary_structure_batches(
        sifts.pdb_ids,
        args.cache_dir,
        batch_size=args.batch_size,
        workers=args.workers,
        force=args.force_fetch,
        timeout=args.timeout,
    )
    metadata = rcsb.consolidate(
        sifts,
        batches,
        args.summary,
        args.elements,
        metadata_path=args.metadata,
    )
    print(json.dumps(metadata["counts"], indent=2, sort_keys=True))
    print(f"summary  -> {os.path.abspath(args.summary)}")
    print(f"elements -> {os.path.abspath(args.elements)}")
    print(f"metadata -> {os.path.abspath(args.metadata)}")


if __name__ == "__main__":
    main()
