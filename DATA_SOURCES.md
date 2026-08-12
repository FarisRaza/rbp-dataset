# Data sources

Large data files are intentionally not stored in this repository. The master
pipeline records source paths, releases, and checksums in its manifests. Put
manually acquired files in one directory and pass that directory with
`--source-dir`, or use the individual override shown below.

```text
sources-manual/
├── 9606.protein.links.v12.0.txt
├── df_np_unique.interpro.tsv
├── df_ptm.csv
├── rbp_master_with_eclip.csv
├── Human_Proteome_RCSB_PDB_Summary.csv       # optional
├── condensates.csv
├── MISSING.csv
└── CD-CODE_Files/
    └── download-data*.csv
```

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/human-proteome-build \
  --source-dir /data/sources-manual
```

`SOURCE_MANIFEST.json` is the machine-readable companion to this page.

## Downloaded automatically by the master pipeline

| Source | Used for | Acquisition |
|---|---|---|
| Reviewed human UniProtKB/Swiss-Prot | Canonical rows, identifiers, GO, curated domains | UniProt REST query `organism_id:9606 AND reviewed:true` |
| NCBI Datasets human gene package | RefSeq protein sequences, transcripts, isoform names, Gene IDs | NCBI Datasets CLI; downloaded and extracted automatically |
| Current Ensembl GRCh38 peptide FASTA | Conservative exact-sequence ENSG/ENST/ENSP fallback | Ensembl FTP `current_fasta` |
| Open Targets Platform release 25.12 | Tissue expression and disease/condition associations | Typed Parquet release datasets; release can be changed with `--opentargets-release` |

Each download is cached under `WORK_DIR/sources`. Existing downloads are kept
unless `--force` is used. UniProt, NCBI, and Ensembl downloads receive local
manifest files containing provenance information.

Official documentation:

- [UniProt REST API](https://www.uniprot.org/help/api_queries)
- [NCBI Datasets gene packages](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-gene-data-package/)
- [Ensembl FTP downloads](https://www.ensembl.org/info/data/ftp/index.html)
- [Open Targets datasets](https://platform-docs.opentargets.org/data-access/datasets)

## Public sources requiring a separate preparation step

### eCLIP: ENCODE, ENCORI, POSTAR3, and optional Skipper

Required master input: `rbp_master_with_eclip.csv`

Override: `--eclip-table /path/rbp_master_with_eclip.csv`

The complete source-building workflow is in
[`eclip_source_pipeline/README.md`](eclip_source_pipeline/README.md). It
automatically obtains ENCODE GRCh38 peaks, but the following files must be
downloaded/exported by the user because their upstream interfaces do not expose
one stable bulk artifact used by the historical analysis:

- a current GENCODE GRCh38 GTF, named `gencode.gtf`;
- ENCORI sites, named `sites_annotated.bed`;
- POSTAR3 sites, named `postar3_sites.bed`;
- the RBP census/base table, named `rbp_dataframe.csv`;
- optional Skipper published references and window archives.

The join deliberately retains ENCODE and ENCORI as separate evidence families.
The output is gene-level and is broadcast to protein isoforms with an explicit
scope column.

### InterProScan

Required master input: `df_np_unique.interpro.tsv`

Override: `--interpro-tsv /path/results.interpro.tsv`

Run InterProScan on the RefSeq protein sequences represented in the clean
catalog. Preserve the RefSeq protein accession as each FASTA record's ID.

```bash
python export_catalog_fasta.py \
  --catalog /data/build/catalog/human_protein_isoforms.parquet \
  --row-kind isoform --identifier refseq \
  --output /data/interpro/refseq_proteins.fasta

interproscan.sh \
  -i /data/interpro/refseq_proteins.fasta \
  -f tsv \
  -o df_np_unique.interpro.tsv \
  --cpu 4 -goterms -iprlookup
```

The historical snapshot used InterProScan 5.77-108.0. A current rebuild may use
a newer release, but should preserve its version in the run notes and should
not compare hit counts across releases without checking member-database changes.
For large input, `interpro.write_chunked_jobs(...)` writes the same chunked SGE
job structure used for the historical run. See the [InterProScan download and
installation documentation](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html).

### STRING v12

Required master input: `9606.protein.links.v12.0.txt`

Override: `--string-links /path/9606.protein.links.v12.0.txt`

Download the human organism-specific **full links** file from the
[STRING v12 download page](https://version-12-0.string-db.org/cgi/download),
decompress it, and retain the filename above. The code streams this large file;
it is not loaded entirely into memory.

### RCSB/PDB (optional)

Required input when `rcsb` is selected: `Human_Proteome_RCSB_PDB_Summary.csv`

Override: `--rcsb-summary /path/Human_Proteome_RCSB_PDB_Summary.csv`

This source can be regenerated from SIFTS and PDBe APIs:

```bash
python export_catalog_fasta.py \
  --catalog /data/build/catalog/human_protein_isoforms.parquet \
  --row-kind canonical --identifier uniprot \
  --output /data/rcsb/reviewed_human_swissprot.fasta

python stage_rcsb.py \
  --fasta /data/rcsb/reviewed_human_swissprot.fasta \
  --summary /path/Human_Proteome_RCSB_PDB_Summary.csv \
  --elements /path/Human_Proteome_RCSB_Secondary_Structure.csv.gz
```

RCSB is optional because it is slower and experimental structure coverage is
canonical-sequence-specific. Select it with `--features all` or include `rcsb`
in an explicit list.

## Project snapshots requiring archival or a new upstream rebuild

### PTMs

Required input: `df_ptm.csv`

Override: `--ptm-csv /path/df_ptm.csv`

This is a canonical-UniProt, wide site table covering eleven PTM classes. The
feature code joins canonical sites directly and conservatively projects only
residue-preserving sites to NCBI isoforms by sequence alignment.

The local `ptm.txt` appears to be an iPTMNet-style site export based on its
schema and source labels, but it is truncated part-way through the export and
cannot reproduce `df_ptm.csv`. Therefore either:

1. archive the complete `df_ptm.csv` outside Git (institutional storage or a
   versioned data release), or
2. create a new complete export from the [iPTMNet API](https://research.bioinformatics.udel.edu/iptmnet/api/doc/), convert it to the schema documented in `ptm.py`, and record the retrieval date/version.

This is a known source-reproducibility gap, not a code gap: `annotate_ptm.py`
and the isoform projection logic are fully included.

### CD-CODE

Required inputs: `condensates.csv`, `MISSING.csv`, and `CD-CODE_Files/`

Override: `--cdcode-root /path/directory-containing-these-files`

The historical source consisted of one member download per condensate from
[CD-CODE](https://cd-code.org). No stable bulk API/download was used, and the
files require a documented positional correction. The code verifies membership
counts before annotation and refuses to run if the ordering no longer aligns.
For exact reproduction, preserve this source bundle in external archival
storage. A future replacement should use a stable bulk export keyed by
condensate UID.

## Model/code dependency acquired by setup

PSLab uses the published
[`_2024_buelow_PSpred`](https://github.com/KULL-Centre/_2024_buelow_PSpred)
repository and its trained models. `python setup_environment.py` clones it and
creates a separate scikit-learn 1.6 environment. IDR prediction similarly gets
its own metapredict 3.0.2 environment. Neither environment nor vendored model
repository should be committed to Git.

## What should be archived outside Git

At minimum, archive the following with checksums if the exact meeting-era table
must remain reproducible:

- `df_ptm.csv`;
- `condensates.csv`, `MISSING.csv`, and all of `CD-CODE_Files/`;
- the final `rbp_master_with_eclip.csv` and its upstream export dates;
- the InterProScan TSV and InterProScan version;
- the STRING v12 links file;
- any final Parquet/CSV table release.

The Git repository should contain code and documentation only. Data release
artifacts can be attached separately through Zenodo, Figshare, an institutional
repository, or GitHub Releases when size permits.
