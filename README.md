# Human proteome and RNA-binding-protein dataset

This repository rebuilds the Kappel Lab protein/isoform table from source data.
It produces one canonical row for every reviewed human UniProtKB/Swiss-Prot
entry plus sequence-unique curated NCBI RefSeq `NP_` isoforms mapped to those
proteins.

The final table follows the feature-block layout of
`Isoform_Post_Merge_PSLab_OpenTargets.csv`, with three intentional changes:

- explicit UniProt, NCBI Gene ID, RefSeq, ENSG, ENST, and ENSP identifiers;
- revised Open Targets tissue-expression and disease/condition annotations;
- `dominant_isoform = 1` only when the row sequence exactly matches its mapped
  canonical UniProt sequence.

UniProt-coordinate annotations are null on sequence-distinct isoforms. Features
computed directly from amino-acid sequence, including CIDER, metapredict IDRs,
IDR-CIDER, and PSLab, are computed for every sequence.

## Repository layout

```text
build_dataset.py          one master command
setup_environment.py      installs metapredict and PSLab environments
rbp_pipeline/             current catalog, feature-family, and assembly code
docs/feature_families/    interpretation and sanity-check documentation
eclip_source_pipeline/    raw ENCODE/ENCORI/POSTAR/Skipper preparation
legacy/                   obsolete append-to-an-existing-CSV workflow
tests/                    fast offline tests
```

Large biological sources, generated tables, virtual environments, and model
checkouts are excluded from Git.

## Install

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python setup_environment.py
python -m unittest discover -s tests -v
```

`setup_environment.py` creates separate metapredict and PSLab environments
because their model dependencies conflict. It also downloads the published
PSpred models.

## Prepare source data

The master command automatically downloads current reviewed human Swiss-Prot,
the human RefSeq `NP_` release, Ensembl peptides, Open Targets, and the current
SIFTS PDB mapping. Some archived or separately prepared inputs cannot be stored
in Git:

- `rbp_master_with_eclip.csv`
- `df_np_unique.interpro.tsv`
- `df_ptm.csv`
- `9606.protein.links.v12.0.txt`
- `condensates.csv`, `MISSING.csv`, and `CD-CODE_Files/`

Place them together in a source directory. Acquisition and regeneration details
are in [DATA_SOURCES.md](DATA_SOURCES.md).

## Build the complete table

```bash
python build_dataset.py all \
  --work-dir /data/human-proteome-build \
  --source-dir /data/human-proteome-sources \
  --features all
```

The workflow is resumable. Existing downloads, catalogs, and feature sidecars
are reused unless `--force` is supplied. The main output is:

```text
WORK_DIR/human_proteome_isoforms_features.parquet
```

Use `--output /path/table.csv` when CSV is required. Parquet is much smaller
and preserves typed identifier lists.

## Choose column families

`--features` selects complete feature-family blocks. Available values are:

```text
idr, domains, go, eclip, interpro, ptm, opentargets, cdcode,
string, go_roles, cider, pslab, rcsb
```

For example, build only sequence/RNA-binding-oriented annotations:

```bash
python build_dataset.py all \
  --work-dir /data/rbp-feature-build \
  --source-dir /data/human-proteome-sources \
  --features idr,cider,pslab,eclip,interpro,go,go_roles
```

Dependencies are automatic. Requesting `cider` also computes IDR and domain
sidecars needed for per-region CIDER, but only requested family blocks are
included in the final table.

To keep only exact final columns, add `--columns`. `protein_key` is always kept:

```bash
python build_dataset.py assemble \
  --work-dir /data/rbp-feature-build \
  --features idr,cider,eclip \
  --columns uniprot_id,dominant_isoform,Name,ENSG,sequence,IDR_count,FCR,encori_published_n_regions
```

Run `python build_dataset.py --help` for all options.

## Choose rows

Selectors are applied before feature computation, so they reduce both runtime
and sidecar size. Multiple selector types are combined by union.

Select UniProt, Ensembl, RefSeq, NCBI Gene/HGNC IDs, protein keys, or symbols:

```bash
python build_dataset.py all \
  --work-dir /data/selected-proteins \
  --source-dir /data/human-proteome-sources \
  --proteins "P09651,Q01844,ENSG00000116044" \
  --features idr,cider,eclip,go
```

Select identifiers from a text file:

```bash
python build_dataset.py all \
  --work-dir /data/selected-proteins \
  --protein-list rbp_accessions.txt
```

Select exact amino-acid sequences from FASTA:

```bash
python build_dataset.py all \
  --work-dir /data/fasta-selection \
  --protein-fasta proteins.fasta \
  --strict-selection
```

A UniProt accession retains its canonical row and mapped NCBI isoform rows.
FASTA selection is exact-sequence-specific. `--strict-selection` fails if any
request is unmatched.

## Feature scope

| Family | Scope in the final table |
|---|---|
| CIDER, metapredict IDR, IDR-CIDER, PSLab | every amino-acid sequence |
| UniProt domains and domain-CIDER | canonical UniProt sequence only |
| UniProt GO and derived GO RNA-role flags | canonical UniProt sequence only |
| PTM snapshot | canonical UniProt sequence only |
| CD-CODE | canonical UniProt protein only |
| RCSB/PDB | canonical UniProt sequence only; all current SIFTS-mapped PDB IDs |
| InterPro | matching RefSeq protein sequence |
| Open Targets | ENSG/gene-level, broadcast to mapped isoforms |
| eCLIP | gene-level, broadcast to mapped isoforms; ENCODE and ENCORI kept separate |
| STRING | represented ENSP identifiers |

Each sidecar includes a scope/provenance field where the distinction matters.
Missing eCLIP data means “not measured in the supplied compilation,” not a
negative binding result.

## Outputs

```text
WORK_DIR/
|-- sources/                                  downloaded source snapshots
|-- catalog/
|   |-- human_protein_isoforms.parquet
|   |-- catalog.audit.json
|   `-- reports/                              mapping-gap reports
|-- features/                                 one keyed Parquet per family
|-- human_proteome_isoforms_features_clean.parquet
|-- human_proteome_isoforms_features.parquet  sorted compatibility output
`-- run_manifest.json
```

The clean assembly retains native identifiers and provenance. The final output
adds the familiar legacy aliases and feature-family order, replaces outdated
Open Targets columns, and sorts by UniProt with the canonical sequence first.

## Run one family directly

Every `rbp_pipeline/annotate_<family>.py` file is independently runnable. For
example:

```bash
python rbp_pipeline/annotate_interpro.py \
  --input /data/build/catalog/human_protein_isoforms.parquet \
  --interpro-tsv /data/sources/df_np_unique.interpro.tsv \
  --output /data/build/features/interpro.parquet
```

Scientific implementation and source notes live in the corresponding module,
such as `rbp_pipeline/interpro.py` or `rbp_pipeline/eclip.py`.

## Sanity-check reports

```bash
python rbp_pipeline/generate_feature_reports.py \
  --work-dir /data/human-proteome-build \
  --features all \
  --skip-unavailable
```

This writes one Markdown summary plus coverage/distribution figures per family.
Generated data remain outside Git; manifests and report summaries should be
archived with any released table.
