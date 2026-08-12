# Human proteome and isoform feature pipeline

This repository reconstructs a protein-level dataset containing:

- one guaranteed canonical row for every reviewed human UniProtKB/Swiss-Prot
  entry;
- sequence-deduplicated human NCBI RefSeq protein isoforms;
- explicit UniProt, NCBI Gene, HGNC, RefSeq, ENSG, ENST, and ENSP mappings;
- independently generated feature sidecars that are joined by `protein_key`;
- clean tissue-expression and disease/condition annotations from Open Targets.

The pipeline is designed so each scientific feature family can be rerun alone,
and the master command can rebuild a selected subset or the complete table.
Large source data, generated Parquet files, virtual environments, and model
checkouts are deliberately excluded from Git.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python setup_environment.py
```

`setup_environment.py` creates isolated environments for metapredict and PSLab
because their model dependencies conflict. It also clones the published PSLab
predictor and model files. Do not commit `.venv`, those environments, or
`vendor/`.

## Data sources

The master downloads UniProt, NCBI, Ensembl, and Open Targets automatically.
Several feature families need a manually downloaded, computed, or archived
source. See [`DATA_SOURCES.md`](DATA_SOURCES.md) before the first full run and
use [`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json) for automation.

The most important known gaps are:

- PTM annotation currently requires the complete archived `df_ptm.csv` or a
  newly built iPTMNet export; the local historical raw export was truncated.
- Exact historical CD-CODE reproduction requires the externally archived
  `condensates.csv`, `MISSING.csv`, and `CD-CODE_Files/` bundle.

These data files should live outside the repository.

## Quick start

Build the default feature set:

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/human-proteome-build \
  --source-dir /data/human-proteome-sources
```

The default features are `idr`, `domains`, `go`, `eclip`, `interpro`, `ptm`,
`opentargets`, `cdcode`, `string`, `go_roles`, `cider`, and `pslab`. RCSB is
optional because its network stage is slower.

Build every family, including RCSB:

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/human-proteome-build \
  --source-dir /data/human-proteome-sources \
  --features all
```

Stages are resumable. Existing downloads, catalogs, and feature sidecars are
kept unless `--force` is supplied.

## Select only particular feature families

`--features` accepts `default`, `all`, or a comma-separated list. Only requested
sidecars are assembled, reducing computation and final-table storage.

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/rna-binding-subset \
  --source-dir /data/human-proteome-sources \
  --features idr,cider,eclip,interpro,ptm,go
```

Dependencies are added and ordered automatically. For example, selecting only
`cider` computes the IDR and domain sidecars it needs, but the final table still
contains only the requested CIDER family. Similarly, `pslab` computes IDRs and
`go_roles` computes GO annotations as hidden dependencies.

To tolerate unavailable manual sources during an exploratory run:

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/minimal-build \
  --features idr,domains,go,cider \
  --skip-unavailable
```

## Select only particular proteins

Selectors are applied after the complete canonical-plus-isoform catalog is
constructed and before any feature family runs. Different selectors should use
different work directories, or `--force` when intentionally replacing a
catalog.

By identifiers:

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/rbp-build \
  --proteins "P09651,Q01844,ENSG00000116044" \
  --features idr,cider,eclip,interpro,ptm
```

By a text file containing one identifier per line:

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/rbp-build \
  --protein-list rbp_accessions.txt
```

By exact amino-acid sequences in FASTA:

```bash
python rebuild_from_scratch.py all \
  --work-dir /data/fasta-build \
  --protein-fasta proteins.fasta \
  --strict-selection
```

Identifier matching accepts protein keys, gene symbols, NCBI Gene/HGNC IDs,
UniProt accessions, RefSeq protein/transcript accessions, and ENSG/ENST/ENSP
IDs. Versioned and versionless RefSeq/Ensembl IDs are recognized. A canonical
UniProt accession selects its canonical row and mapped NCBI isoform rows. FASTA
selection is exact-sequence-specific. Multiple selector types are combined by
union.

The catalog audit records matched and unmatched requests. Use
`--strict-selection` to fail rather than continue when anything is unmatched.

## Outputs

```text
WORK_DIR/
├── sources/                              # downloads; not Git content
├── catalog/
│   ├── human_protein_isoforms.parquet
│   └── catalog.audit.json
├── features/
│   ├── idr.parquet
│   ├── eclip.parquet
│   └── ... one keyed sidecar per selected family
├── human_proteome_isoforms_features.parquet
├── human_proteome_isoforms_features.parquet.manifest.json
└── run_manifest.json
```

The assembled table can also be CSV:

```bash
python rebuild_from_scratch.py assemble \
  --work-dir /data/human-proteome-build \
  --features idr,cider,eclip \
  --output /data/human-proteome-build/subset.csv
```

In Python, the requested dataframe is simply:

```python
import pandas as pd
df = pd.read_parquet("/data/human-proteome-build/human_proteome_isoforms_features.parquet")
```

## Run one family directly

Each `annotate_<family>.py` program reads the clean catalog and writes one keyed
sidecar. For example:

```bash
python annotate_interpro.py \
  --input /data/build/catalog/human_protein_isoforms.parquet \
  --interpro-tsv /data/sources/df_np_unique.interpro.tsv \
  --output /data/build/features/interpro.parquet
```

Available families:

| Family | Main source/method |
|---|---|
| `cider` | localCIDER sequence metrics; also per IDR/domain |
| `idr` | metapredict V3 disorder/fold geometry |
| `domains` | reviewed UniProt DOMAIN and ZN_FING features |
| `go` | UniProt GO cross-references |
| `eclip` | separate ENCODE, ENCORI, POSTAR3, and Skipper summaries |
| `interpro` | InterProScan domain models on RefSeq proteins |
| `ptm` | canonical PTM sites with residue-conserving isoform projection |
| `opentargets` | normalized tissue expression and disease associations |
| `cdcode` | CD-CODE condensate membership |
| `string` | STRING v12 protein interactions |
| `go_roles` | derived transcription/translation/mRNA-stability flags |
| `pslab` | phase-separation predictions per IDR |
| `rcsb` | optional PDB/SIFTS and secondary-structure summary |

Implementation details and provenance live in each family module's docstring.
See the individual [`docs/feature_families`](docs/feature_families/README.md)
guides for interpretation and family-specific sanity checks.
The historical compatibility workflow remains documented in the source files,
while [`REBUILD_FROM_SCRATCH.md`](REBUILD_FROM_SCRATCH.md) describes the clean
catalog design in more depth.

`export_catalog_fasta.py` exports canonical, NCBI-isoform, or all catalog
sequences with stable UniProt, RefSeq, or `protein_key` headers. It is the
documented preparation step for InterProScan and the optional RCSB workflow.

## Exploratory and sanity-check reports

After a run, generate one Markdown report and two PNG figures per selected
feature family:

```bash
python generate_feature_reports.py \
  --work-dir /data/human-proteome-build \
  --features all \
  --skip-unavailable
```

Reports are written to `WORK_DIR/reports/<family>.md`. Each report checks key
uniqueness and catalog coverage, summarizes populated columns and the primary
feature distribution, and embeds coverage/distribution figures. This provides
a fast sanity check without reading the very wide assembled table.

To keep a snapshot of reports in the GitHub repository, direct the output to a
tracked documentation directory after a finalized run:

```bash
python generate_feature_reports.py \
  --work-dir /data/human-proteome-build \
  --features all \
  --output-dir docs/generated_reports
git add docs/generated_reports
```

These reports and PNGs are small derived documentation; the underlying data
sidecars and final Parquet table should remain outside Git.

## Tests

```bash
python -m unittest -v test_rebuild.py
python -m py_compile *.py eclip_source_pipeline/*.py
```

Before publishing a data release, preserve `run_manifest.json`,
`catalog.audit.json`, every source manifest/checksum, and the final table's
assembly manifest alongside the data artifact.
