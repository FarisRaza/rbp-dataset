# Rebuild the human protein/isoform dataset from scratch

This is the clean, current-source pipeline. It is separate from `run_all.py`,
whose purpose is to append missing canonical proteins to the historical table.
The clean rebuild prioritizes biological identity and provenance over
bit-for-bit reproduction of old CSV cells.

## One-command run

```bash
python -m pip install -r requirements.txt
python setup_environment.py
python rebuild_from_scratch.py all --work-dir D:/human_isoform_rebuild
```

On Windows, replace the work directory with any short path that has enough free
space. The run is resumable; existing downloads, catalog and feature sidecars
are kept unless `--force` is supplied.

The final output is:

```text
<work-dir>/human_proteome_isoforms_features.parquet
```

Intermediate outputs are intentional deliverables:

```text
sources/                                  downloaded, versioned source artifacts
catalog/human_protein_isoforms.parquet    identity + sequence table
catalog/catalog.audit.json                row and deduplication counts
features/<family>.parquet                 one independently rerunnable family
features/opentargets_expression_long.parquet  one row per ENSG/tissue
features/opentargets_disease_long.parquet     one row per ENSG/condition/datatype
features/opentargets_source_index.sqlite      bounded-memory source index
run_manifest.json                         requested families and paths
```

Use `--output final.csv` if a CSV is required. Parquet is strongly recommended:
sequences and nested feature values make the CSV much larger and slower.

## What is a row?

1. Every reviewed human Swiss-Prot accession creates exactly one canonical row.
   Swiss-Prot proteins are therefore not lost merely because NCBI lacks a
   corresponding RefSeq product.
2. NCBI RefSeq products create one row per unique `(NCBI Gene ID, amino-acid
   sequence)`.
3. If several transcripts or RefSeq protein accessions encode the same protein
   sequence, all identifiers are retained in lists on one row.
4. If an NCBI sequence is identical to the Swiss-Prot canonical sequence at the
   same gene, its identifiers enrich the canonical row instead of creating a
   duplicate row.
5. Curated `NP_` products are included by default. Add `--include-predicted` to
   include `XP_` products as well.
6. NCBI proteins with no reviewed Swiss-Prot parent are excluded by default;
   add `--include-orphan-refseq` only if that broader NCBI-defined universe is
   wanted.

This is why transcript-rich genes such as BRCA1 do not create hundreds of
duplicate protein rows.

The number of “missing proteins” depends on which historical file is compared.
The older table was about 5,000 reviewed proteins short; the more recent master
had already gained roughly 2,600, leaving about 2,900. The clean rebuild avoids
hard-coding either number and starts from the current reviewed human Swiss-Prot
release each time.

## Identifier mapping

Identifier columns are lists, not overloaded comma-separated strings:

- `uniprot_parent_ids`, `uniprot_isoform_ids`
- `refseq_protein_ids`, `refseq_transcript_ids`
- `ensembl_gene_ids`, `ensembl_transcript_ids`, `ensembl_protein_ids`
- `hgnc_ids`, `ncbi_gene_id`
- `gene_synonyms`, `ncbi_gene_ids`, `uniprot_secondary_accessions`

Mappings are collected in this order:

1. UniProt GeneID, RefSeq, Ensembl and HGNC cross-references;
2. NCBI gene and product reports;
3. optional exact full-sequence matching to the current Ensembl human peptide
   FASTA, restricted to a known ENSG when possible.

`identifier_mapping_methods` states which evidence populated a row.
`identifier_ambiguity` records multiple parents or multiple exact Ensembl
matches. The code never chooses one ENSP and broadcasts it to every isoform.

The stable key is `sp:<accession>` for canonical proteins and
`refseq:<GeneID>:<sequence-hash-prefix>` for NCBI sequence-unique rows.

## Run one feature family

Every annotation program consumes the clean catalog and writes a sidecar keyed
by `protein_key`. Typical commands are below; every program also supports
`--help`.

| Family | Command | Primary input or method |
|---|---|---|
| CIDER | `python annotate_cider.py --input catalog.parquet --output cider.parquet --idr-features idr.parquet --domain-features domains.parquet` | amino-acid sequence; optional IDR/domain sidecars |
| IDR/fold geometry | `%PYTHON_METAPREDICT% annotate_idr.py --input catalog.parquet --output idr.parquet` | metapredict V3, threshold 0.5 |
| Curated domains | `python annotate_domains.py --input catalog.parquet --output domains.parquet --swissprot human.dat` | Swiss-Prot `DOMAIN` and `ZN_FING` features |
| Gene Ontology | `python annotate_go.py --input catalog.parquet --output go.parquet --swissprot human.dat` | UniProt GO cross-references |
| GO role flags | `python annotate_go_roles.py --input catalog.parquet --output go_roles.parquet --go-features go.parquet` | derived from GO term names |
| eCLIP/CLIP | `python annotate_eclip.py --input catalog.parquet --output eclip.parquet --eclip-table rbp_master_with_eclip.csv` | ENCODE, ENCORI, POSTAR3 and Skipper compilation |
| InterPro | `python annotate_interpro.py --input catalog.parquet --output interpro.parquet --interpro-tsv df_np_unique.interpro.tsv` | InterProScan TSV joined by RefSeq protein |
| PTM | `python annotate_ptm.py --input catalog.parquet --output ptm.parquet --ptm-csv df_ptm.csv` | canonical UniProt sites; alignment projection to isoforms |
| Open Targets | `python annotate_opentargets.py ...` | tissue expression plus named disease/condition associations, keyed by ENSG |
| CD-CODE | `python annotate_cdcode.py ... --cdcode-root <source-dir>` | condensate membership by UniProt parent |
| STRING | `python annotate_string.py ... --string-links 9606.protein.links.v12.0.txt` | all represented ENSP queries |
| PSLab | `%PYTHON_PSLAB% annotate_pslab.py ... --idr-features idr.parquet --pspred-repo <repo>` | one prediction per metapredict IDR |
| RCSB/PDBe | `python annotate_rcsb.py ... --rcsb-summary Human_Proteome_RCSB_PDB_Summary.csv` | canonical SIFTS/PDB summary |

Scientific logic and source notes live in the correspondingly named modules:
`cider.py`, `idr.py`, `domains.py`, `go.py`, `eclip.py`, `interpro.py`, `ptm.py`,
`opentargets.py`, `cdcode.py`, `string_ppi.py`, `pslab.py`, `go_roles.py` and
`rcsb.py`. `feature_runner.py` only standardizes catalog input and keyed output.
The raw-source ENCODE/ENCORI/POSTAR/Skipper workflow is preserved under
`eclip_source_pipeline/`; see its README for execution order and criteria.
Use `summarize_eclip.py` to regenerate the eCLIP summary statistics and figures.

## Scope matters

- Sequence-derived CIDER and IDR features are computed separately for every
  sequence row.
- UniProt domain coordinates and RCSB residue mappings apply only to canonical
  Swiss-Prot rows. They are not copied to sequence-distinct isoforms.
- GO, CD-CODE and eCLIP are protein- or gene-level measurements and may be
  inherited/broadcast to isoforms. Their sidecars include an explicit scope
  column.
- eCLIP is gene-level. A row with no CLIP values means “not measured in the
  supplied compilation,” not “measured and no RNA binding was found.”
- InterPro is sequence-specific and joins through each row's RefSeq protein
  accessions.
- PTMs begin as canonical UniProt coordinates. For an NCBI isoform, a site is
  retained only if sequence alignment maps the position and the modified amino
  acid is conserved. Projection method, parent accession and dropped-site count
  are stored.
- Open Targets annotations are gene-level and are broadcast to canonical and
  isoform rows that share an ENSG. Missing data means "not present in the
  release", not a measured negative result.

## Open Targets: clean expression and condition columns

The clean rebuild does **not** import the broad historical `target_full.csv`
dump. It downloads three typed, pinned Open Targets datasets and writes only
the requested information:

- `opentargets_tissue_expression`: flat records containing ENSG, tissue EFO or
  UBERON ID/name, organ/system labels, RNA value/unit/z-score/level, protein
  level/reliability and cell-type protein levels;
- `opentargets_expression_tissue_count`;
- `opentargets_disease_associations`: named conditions, descriptions,
  therapeutic areas and evidence broken down by datatype;
- `opentargets_disease_count`, `opentargets_disease_names` and
  `opentargets_therapeutic_areas`;
- `opentargets_release` and `opentargets_annotation_scope`.

`max_datatype_score` is explicitly the maximum of the per-datatype scores. It
is not mislabeled as the Open Targets overall association score. The complete
per-datatype scores and evidence counts remain in `evidence_by_datatype`.

The master intentionally defaults to the reproducible `25.12` release, which
matches the compact tissue-expression schema used by the local historical
export. Choose another release that still publishes the compatible
`expression` dataset with `--opentargets-release YY.MM`. Open Targets changed
that product to the substantially different `baseline_expression` dataset in
26.06, so 26.06 is not silently treated as the same measurement. The
`download` and `all` commands retrieve `expression`,
`association_by_datatype_direct` and `disease` Parquet files into
`<work-dir>/sources/opentargets/<release>/`.

Run only this family with already-downloaded sources:

```bash
python annotate_opentargets.py \
  --input catalog/human_protein_isoforms.parquet \
  --output features/opentargets.parquet \
  --opentargets-expression sources/opentargets/25.12/expression \
  --opentargets-associations sources/opentargets/25.12/association_by_datatype_direct \
  --opentargets-diseases sources/opentargets/25.12/disease.parquet \
  --opentargets-release 25.12 \
  --opentargets-cache features/opentargets_source_index.sqlite \
  --opentargets-expression-long-output features/opentargets_expression_long.parquet \
  --opentargets-disease-long-output features/opentargets_disease_long.parquet
```

The two `*-long-output` arguments are optional. They provide analysis-ready
tables without requiring users to parse nested cells from the assembled
protein table. The SQLite cache keeps the full-proteome run bounded in memory;
if `--opentargets-cache` is omitted, an equivalent temporary index is deleted
after the sidecar is written. The loaders also accept the older local CSV
association and multiline expression exports for backward compatibility.

## Source availability

The master script downloads Swiss-Prot, NCBI Datasets, Ensembl and the pinned
Open Targets expression/association/disease datasets automatically.
Several historical feature sources are local snapshots and must be placed in
`KAPPEL_DIR` (the parent data directory by default):

- `rbp_master_with_eclip.csv`
- `df_np_unique.interpro.tsv`
- `df_ptm.csv`
- STRING links and the CD-CODE source directory

Run with `--skip-unavailable` to build every family whose source is present, or
select a subset, for example:

```bash
python rebuild_from_scratch.py features \
  --work-dir D:/human_isoform_rebuild \
  --features domains,go,idr,cider,ptm
```

RCSB is optional in the master default because `stage_rcsb.py` first needs to
generate its current normalized structure files. Select it with
`--features all` after that summary exists.

## Assemble selected sidecars

The master does this automatically. It can also be called directly:

```bash
python assemble_features.py \
  --catalog catalog/human_protein_isoforms.parquet \
  --feature features/cider.parquet \
  --feature features/ptm.parquet \
  --output final.parquet
```

Assembly rejects duplicate `protein_key` values and duplicate column names, then
verifies that the output row count exactly equals the catalog row count.

## Verification

Fast offline tests cover the row invariants, sequence deduplication, PTM
projection and keyed assembly:

```bash
python -m unittest -v test_rebuild.py
```
