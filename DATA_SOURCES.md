# Data sources

Source data and generated tables are intentionally excluded from Git. Put
manually prepared files in one directory and pass it with `--source-dir`.

```text
sources-manual/
|-- 9606.protein.links.v12.0.txt
|-- df_np_unique.interpro.tsv
|-- df_ptm.csv
|-- rbp_master_with_eclip.csv
|-- condensates.csv
|-- MISSING.csv
`-- CD-CODE_Files/
    `-- download-data*.csv
```

```bash
python build_dataset.py all \
  --work-dir /data/human-proteome-build \
  --source-dir /data/sources-manual \
  --features all
```

`SOURCE_MANIFEST.json` is the machine-readable companion to this page.

## Downloaded automatically

| Source | Used for |
|---|---|
| Reviewed human UniProtKB/Swiss-Prot | Canonical sequences, identifiers, GO, curated domains |
| Human RefSeq FTP release | Current curated `NP_` sequences and GRCh38/T2T GeneID/transcript mappings |
| Current Ensembl GRCh38 peptide FASTA | Conservative exact-sequence ENSG/ENST/ENSP fallback |
| Open Targets Platform 25.12 | Tissue expression and disease/condition associations |
| Weekly SIFTS UniProt/PDB mapping | Every PDB entry and chain mapped to each canonical UniProt protein |

Downloads are cached under `WORK_DIR/sources` and receive provenance manifests
where the upstream format permits it. Existing files are kept unless `--force`
is used.

The default NCBI route is NP-only. The alternative NCBI Datasets route supports
predicted products with `--ncbi-source datasets --include-predicted`.

## eCLIP: ENCODE and ENCORI

Required master input: `rbp_master_with_eclip.csv`

Override: `--eclip-table /path/rbp_master_with_eclip.csv`

The source-building workflow is in
[`eclip_source_pipeline/README.md`](eclip_source_pipeline/README.md). It keeps
ENCODE and ENCORI as separate evidence families to avoid treating overlapping
data as independent replication. POSTAR3 and optional Skipper summaries are
also retained. The resulting annotation is gene-level and is broadcast to
mapped protein isoforms with an explicit scope column.

## InterProScan

Required master input: `df_np_unique.interpro.tsv`

Override: `--interpro-tsv /path/results.interpro.tsv`

Generate an accession-stable RefSeq FASTA from the catalog and run InterProScan
with InterPro and GO lookups:

```bash
python rbp_pipeline/export_catalog_fasta.py \
  --catalog /data/build/catalog/human_protein_isoforms.parquet \
  --row-kind isoform \
  --identifier refseq \
  --output /data/interpro/refseq_proteins.fasta

interproscan.sh \
  -i /data/interpro/refseq_proteins.fasta \
  -f tsv \
  -o /data/sources-manual/df_np_unique.interpro.tsv \
  --cpu 4 -goterms -iprlookup
```

The archived project snapshot used InterProScan 5.77-108.0. A current run may
use a newer release, but its version should be recorded because member-database
content changes over time. `rbp_pipeline/interpro.py` also contains the chunked
SGE job generator used for the original large run.

## STRING

Required master input: `9606.protein.links.v12.0.txt`

Override: `--string-links /path/9606.protein.links.v12.0.txt`

Download and decompress the human full-links file from the STRING v12 download
site. The feature code streams it and does not load the complete network into
memory.

## PTMs

Required master input: `df_ptm.csv`

Override: `--ptm-csv /path/df_ptm.csv`

This project snapshot is a canonical-UniProt wide table covering eleven PTM
classes. Canonical rows are joined directly. Sequence-distinct NCBI isoforms
remain null because the observations use canonical UniProt coordinates.

The available raw `ptm.txt` is truncated and cannot recreate the complete wide
snapshot. Preserve `df_ptm.csv` in archival storage, or create a new complete
iPTMNet export and record its retrieval date. Alignment projection code remains
available in `rbp_pipeline/ptm.py` for explicit inferred analyses, but the
master table does not present projected sites as directly sourced data.

## CD-CODE

Required inputs: `condensates.csv`, `MISSING.csv`, and `CD-CODE_Files/`

Override: `--cdcode-root /path/source-directory`

The historical project used one member download per condensate and one manual
positional correction. The code validates source/member ordering before
annotation and fails rather than silently shifting assignments. Preserve the
validated bundle externally because no stable keyed bulk export was used.

## RCSB/PDB

No manual source is needed for PDB identifiers. The master downloads current
SIFTS data and records all PDB IDs and mapped chains for canonical UniProt rows.

An optional enriched summary containing PDBe secondary-structure observations
can be generated and supplied:

```bash
python rbp_pipeline/export_catalog_fasta.py \
  --catalog /data/build/catalog/human_protein_isoforms.parquet \
  --row-kind canonical \
  --identifier uniprot \
  --output /data/rcsb/reviewed_human_swissprot.fasta

python rbp_pipeline/stage_rcsb.py \
  --fasta /data/rcsb/reviewed_human_swissprot.fasta \
  --summary /data/rcsb/Human_Proteome_RCSB_PDB_Summary.csv \
  --elements /data/rcsb/Human_Proteome_RCSB_Secondary_Structure.csv.gz
```

Pass that file with `--rcsb-summary`. The automatic route is faster and still
contains every SIFTS-mapped PDB entry; only the secondary-structure fields are
omitted.

## Model dependency acquired by setup

PSLab uses the published `_2024_buelow_PSpred` repository and trained models.
`python setup_environment.py` clones it and creates a scikit-learn 1.6
environment. Metapredict 3.0.2 receives a separate environment. Neither the
environments nor `vendor/` should be committed.

## Archive with a data release

For exact reproducibility, archive checksums and retrieval versions for:

- `df_ptm.csv`;
- the complete CD-CODE bundle;
- `rbp_master_with_eclip.csv` and its upstream export dates;
- the InterProScan TSV and software/database version;
- the STRING v12 links file;
- run manifests, catalog audit, feature reports, and the final Parquet/CSV.

These artifacts belong in a versioned data repository or release attachment,
not in the source-code repository.
