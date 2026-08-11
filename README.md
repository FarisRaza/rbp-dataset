# Isoform table — feature generation

Code to compute every column of `Isoform_Post_Merge_PSLab_OpenTargets.csv` for
any set of human proteins, consistent with how the existing rows were made.

For a clean rebuild from current sources—including one guaranteed row per
reviewed human Swiss-Prot protein, sequence-deduplicated NCBI RefSeq isoforms,
revamped ENSG/ENST/ENSP mapping, independent feature sidecars, and PTMs—start
with [`REBUILD_FROM_SCRATCH.md`](REBUILD_FROM_SCRATCH.md) and run:

```bash
python rebuild_from_scratch.py all --work-dir /path/to/rebuild
```

The older `run_all.py` workflow below remains available when exact compatibility
with the historical 157-column table is the goal.

The table's ten column families are covered by nine modules — `cider.py` serves
three of them (whole sequence, per IDR, per domain). Each module is standalone
with one documented entry point, so a family can be regenerated on its own, or
all of them run together to produce complete 157-column rows.

---

## Quick start

```bash
python setup_environment.py
```

```bash
python paths.py
```

```bash
python run_all.py --all
```

`setup_environment.py` builds the two extra interpreters and vendors the PSLab
predictor. `paths.py` reports which source files are present and where to get
any that are missing. `run_all.py` runs every stage and writes complete rows.

Choose what to annotate:

```bash
python run_all.py --all
```

```bash
python run_all.py --accessions my_ids.txt
```

```bash
python run_all.py --fasta my_proteins.fasta
```

With no selector it annotates exactly the proteins absent from the master table.

---

## Regenerating one family at a time

Every family reads the target set written by stage 1 and writes its own
intermediate, so families are independent and can be run, re-run or replaced
individually. Run stage 1 first, then whichever families you need:

| Family | Columns | Module | Command | Writes |
|---|---|---|---|---|
| *(target set)* | — | `identity.py` | `python stage_find_missing.py [--all]` | `missing_proteins.json`, `missing_idmap.json`, `missing_ensg.json` |
| IDRs | 12 | `idr.py` | `$PYTHON_METAPREDICT stage_idr.py` | `missing_idr.json` |
| Domains | 7 | `domains.py` | `python stage_swissprot.py` | `missing_domains.json` |
| GO | 9 | `go.py` | *(same run as domains)* | `missing_go.json` |
| STRING | 5 | `string_ppi.py` | `python stage_string.py` | `string_new.pkl`, `string_reverse.pkl` |
| Open Targets | 33 | `opentargets.py` | `python stage_opentargets.py` | `missing_opentargets.pkl` |
| PSLab | 11 | `pslab.py` | `$PYTHON_PSLAB stage_pslab.py` | `missing_pslab.json` |
| CIDER ×3 | 58 | `cider.py` | *(computed during assembly)* | — |
| CD-CODE | 10 | `cdcode.py` | *(computed during assembly)* | — |

Then assemble:

```bash
python build_rows.py
```

CIDER and CD-CODE have no stage of their own because both are cheap and need no
large scan — CIDER is pure computation from sequence, CD-CODE a dictionary
lookup — so `build_rows.py` calls them inline.

Domains and GO share one stage because both are read from the same 4 GB
Swiss-Prot flat file and a single pass serves both.

## The ten column families

Every module takes protein sequences (and, where needed, identifiers) and
returns a dict of column-name → value, in the exact on-disk encoding the master
table uses. Column names and order live in `schema.py`, which asserts itself
against the real header.

### 1. CIDER — `cider.py` (19 columns)

Sequence composition and charge patterning from **localcider**
(`SequenceParameters`). No downloads; computed from sequence alone.

```python
import cider
cider.whole_sequence(sequence)          # -> FCR, NCPR, kappa, delta, ...
```

Also supplies the per-IDR and per-domain variants (below). The three blocks are
not interchangeable: the domain block carries an extra `Omega`, and the IDR and
domain blocks order `kappa` differently from the whole-sequence block. Each has
its own explicit getter list.

Sanitisation: `U`→`C`, other non-canonical residues deleted.

### 2. IDRs — `idr.py` (12 columns)

Disordered/folded region geometry from **metapredict** V3.

```python
import idr
idr.install_python_fallback()            # only needed where Cython is unavailable
for row in idr.predict(sequences):
    row["IDR_range"], row["FOLD_range"], ...
```

The load-bearing parameter is `disorder_threshold=0.5` — *not* metapredict's
0.42 default. Using the default shifts nearly every boundary.

Sanitisation differs from CIDER's on purpose: `B`→`N`, `U`→`C`, `X`→`G`, `Z`→`Q`,
substituting rather than deleting so that length is preserved and the returned
indices still address the original sequence.

### 3. CIDER on IDRs — `cider.py` (19 columns)

```python
cider.per_idr(idr_discrete_seqs)         # -> {"FCR": [v1, v2, ...], ...}
```

One value per IDR, ordered to match `IDR_discrete_seq`.

### 4. Domains — `domains.py` (7 columns)

**Not InterPro, not Pfam, not PROSITE** — despite what the documentation
spreadsheet says, and despite `df_np_unique.interpro.tsv` sitting in the project
folder. These are the curated `DOMAIN` and `ZN_FING` features of the
**Swiss-Prot flat file** `uniprot_sprot.dat`, read with `Bio.SwissProt`.

```python
import domains
geometry = domains.domains_for_record(record)        # a Bio.SwissProt record
row = domains.attach_sequences(geometry, sequence)
```

Coordinates are 0-based half-open, so region sequences are plain
`sequence[start:end]` slices. `normalize_domain()` collapses `RRM 1`/`RRM 2`
onto `RRM` and reproduces the original regex quirks exactly — a different
normalisation would make new rows incomparable to the existing 17,582.

### 5. CIDER on Domains — `cider.py` (20 columns)

```python
cider.per_domain({"RRM": [seq1, seq2, seq3]})
# -> {"FCR": {"RRM": [v1, v2, v3]}, ..., "Omega": {...}}
```

Computed per domain *occurrence*, not on the concatenation.

### 6. STRING PPI — `string_ppi.py` (5 columns)

Partners and combined scores from `9606.protein.links.v12.0.txt`, keyed on the
unversioned Ensembl protein id.

```python
import string_ppi
adjacency, reverse = string_ppi.scan(links_path, focus_ensps, also_index_against=existing)
string_ppi.columns_for(ensp, adjacency, ensp_to_uniprot, table_ensps, table_uniprots)
```

No confidence cutoff is applied; STRING's own file floor of 150 is the effective
threshold. The full 13M-edge graph is never held in memory — only edges touching
a protein of interest are retained.

**Three of the five columns are table-relative**, not two. Besides the obvious
`_in_Dataframe` pair, `PPI_UniProt_Partners` is too — translating a partner ENSP
to UniProt requires that partner to be in the table, so untranslatable partners
are silently dropped. The master's own counts show it: for P49588,
`PPI_ENSP_Partners` holds 1,403 entries while the other three hold 881.

So adding proteins changes all three for *existing* rows. The protein set only
grows, so the update is purely additive; `append_to_master.py` applies it in the
same streaming pass that copies rows. Only `ENSP_clean` and
`PPI_ENSP_Partners` are independent of what else is in the table.

### 7. CD-CODE — `cdcode.py` (10 columns)

Condensate membership from 327 hand-downloaded member tables plus
`condensates.csv`.

```python
import cdcode
assert not cdcode.verify_alignment(data_dir)         # run this first
protein_to_condensates, attributes = cdcode.build_index(data_dir)
cdcode.lookup(accession, protein_to_condensates, attributes)
```

A condensate's identity is carried *only* by its member file's position in a
sorted listing, reconstructed via two hardcoded index surgeries. **Always run
`verify_alignment()` first** — it checks every member table's row count against
its declared protein count, and currently reports 0 mismatches. Re-downloading
or renaming any file in `CD-CODE_Files/` would silently corrupt every assignment
after it.

Numeric values are cast to plain Python ints, so new rows are
`ast.literal_eval`-safe; existing rows contain `np.int64(9606)` reprs that are
not.

### 8. Gene Ontology — `go.py` (9 columns)

**Not `goa_human.gaf`** — that file is loaded by one historical notebook and then
never used. These come from the `DR   GO;` cross-references inside
`uniprot_sprot.dat`.

```python
import go
go.go_for_record(record)   # -> C_ids/C_descriptions/C_evidence, P_*, F_*
```

The three lists within an aspect are positionally parallel. No evidence-code
filtering: IEA annotations are included, matching the existing rows.

### 9. Historical Open Targets compatibility (33 columns)

The older `run_all.py` compatibility workflow preserves the original disease,
expression and broad target annotation dump keyed on ENSG:

```python
import opentargets as ot
associations = ot.load_associations(path, wanted_ensgs)
targets, target_columns = ot.load_targets(path, wanted_ensgs)
expression = ot.load_expression(path, wanted_ensgs)
ot.columns_for(ensg_list, associations, expression, targets, target_columns)
```

Every column is a **dict keyed by ENSG**, because one protein can map to several
(`ID` is semicolon-separated). Sources total ~2.9 GB and one field can exceed
50 KB per gene, so all three are streamed and filtered rather than loaded.

The clean `rebuild_from_scratch.py` workflow replaces this block with eight
explicit columns containing normalized tissue expression and named
disease/condition associations. It does not read `target_full.csv`, and it also
writes optional long-form expression and disease Parquet tables. See
[`REBUILD_FROM_SCRATCH.md`](REBUILD_FROM_SCRATCH.md#open-targets-clean-expression-and-condition-columns).

### 10. PSLab — `pslab.py` (11 columns)

Predicted phase-separation propensity per IDR, from
[KULL-Centre/_2024_buelow_PSpred](https://github.com/KULL-Centre/_2024_buelow_PSpred)
(von Bülow, Tesei, Zaidi, Mittag & Lindorff-Larsen, *PNAS* 2025). The original
values came from that repo's Colab notebook; this runs the same models locally.

```python
import pslab
models, residues, nu_file = pslab.load_models(repo_root)
for row in pslab.predict(idr_sequences, models, residues, nu_file):
    row["Delta G [kT]"], row["Saturation concentration [uM]"], ...
```

Eight features are analytic; `Delta G` and the two saturation concentrations
come from MLP ensembles trained on CALVADOS 2 slab simulations, at fixed
T = 293 K and I = 150 mM. One prediction per IDR, stored as lists ordered to
match `IDR_discrete_seq`.

### Not yet in the table — `interpro.py`

InterProScan results exist (`df_np_unique.interpro.tsv`, 55,418 proteins) but
**no InterPro column reaches the 157-column table** — `Domains_*` comes from
UniProt. This module records how that TSV was produced, since no script for it
survived in the project folder:

InterProScan **5.77-108.0** on UCLA's Hoffman2 cluster, over `df_np_unique.fasta`
(55,670 RefSeq `NP_` sequences), split into ~400 chunks of 140 and submitted as
SGE jobs running `./interproscan.sh -i chunk_N.fasta -f tsv -o chunk_N.interpro.tsv
--cpu 4 -goterms -iprlookup` under `module load java/jdk-11.0.14` (the cluster
default is Java 8, which InterProScan refuses). `write_chunked_jobs()`
regenerates that submission setup.

`read_tsv()` / `group_by_protein()` parse the results with all 15 columns named
correctly. The historical notebook passes 14 names for a 15-column file, so
every label shifts by one — which is why it filters on `tsv["length"] == "Pfam"`.
`group_by_protein()` emits 0-based half-open ranges, so its output can be handed
straight to `domains.attach_sequences` and `cider.per_domain`.

### Identity columns — `identity.py` (10 + 2 columns)

`uniprot_id`, `sequence`, `Dominant_Isoform` and `Name` are exact. `ID`, `ENSP`
and `ProteinHGVS` are **reconstructed** from `HUMAN_9606_idmapping.dat` with
HGNC and Open Targets fallbacks — the step that originally produced them is not
in the project folder, so values are correct but not guaranteed byte-identical.
`HGVSDescription` is left empty (no offline source; already empty on 97.9% of
existing rows). `UNIQUE` and `Description` come from the 1,393-gene RBP census
in `RBPs.xlsx` and are null otherwise, as they are for 92% of existing proteins.

---

## Verifying fidelity

Every claim above was checked by recomputing values for proteins **already in
the master table** and diffing against what is stored:

```bash
python validate.py --all
```

| Family | Check | Result |
|---|---|---|
| CIDER | 19 columns × 6 proteins | 114/114 exact |
| Domains + GO | 16 columns × 6 proteins | 96/96 exact |
| STRING | full partner dicts (1,403 and 1,296 partners) | exact |
| CD-CODE | 10 columns + all 327 member counts | 48/48; 0 alignment mismatches |
| IDRs | `IDR_range` + `FOLD_range` × 8 proteins | 8/8 exact |
| PSLab | 11 values × 25 IDRs from the original Colab output | exact |

Re-run this after changing any module or refreshing any source file. A
regression means new rows would be subtly incomparable to existing ones.

---

## Corrections to the documentation spreadsheet

- **Domains are from UniProt, not PROSITE/Pfam** (family 4 above).
- **GO is from UniProt, not `goa_human.gaf`** (family 8 above).
- **`Domains_*` are dicts, not lists**: `Domains_count` is `{'RRM': 3}`,
  `Domains_range` is `{'RRM': [(55, 134), (135, 218), (230, 303)]}`.
- **`tss` is not PSLab's total sequence stickiness.** It is the last column of
  Open Targets' `target_full.csv`, a transcription start coordinate. PSLab does
  emit a `tss` feature but it was dropped before the merge.
- **Counts are stale.** The spreadsheet says 14,956 unique proteins / 55,670
  isoforms; the current table has **17,582 / 65,744**.

## How many proteins are missing

**2,901**, not the ~5,000 in the older notebooks.

| Table | Unique `uniprot_id` | Gap vs. Swiss-Prot |
|---|---|---|
| `Isoform_IDRs_Localization_Condensates_Updated.csv` (older base) | 14,945 | 5,476 |
| `Isoform_Post_Merge_PSLab_OpenTargets.csv` (current) | 17,582 | **2,901** |

Human Swiss-Prot holds 20,420 reviewed accessions. The 5,476 figure was correct
for the table it was computed against; the master has since gained ~2,600
proteins. 63 accessions are in the master but not in human Swiss-Prot and are
left alone.

Of the 2,901, **92 are RNA-binding** — 68 by GO molecular function, 43 by a
classical RBD (18 with RRM) — including RBMY1B/C/D/E, HNRNPA1L3, HNRNPCL3/CL4,
RBMXL2 and ZC3H11B/C.

---

## RCSB PDB secondary structures

`rcsb.py` and `stage_rcsb.py` add a current experimental-structure feature
family for the full reviewed human proteome. The mapping is not inferred from
gene symbols: the current weekly SIFTS `pdb_chain_uniprot.tsv.gz` release maps
each PDB author chain and residue segment to a canonical UniProt accession.
Regular secondary structures are the PDB archive's alpha-helix and beta-strand
ranges exposed by the batch PDBe API; the PDB IDs link to the corresponding
RCSB entry pages.

```bash
python stage_rcsb.py
```

The network phase is batched, gzip-cached and resumable under
`EXPANSION_SCRATCH/rcsb_secondary_batches`. A current full run covers 20,420
reviewed human proteins, of which 8,982 have at least one mapped experimental
PDB entry (77,110 distinct entries in the 2026-08-03 SIFTS release).

Outputs:

- `Human_Proteome_RCSB_PDB_Summary.csv` — one row per reviewed protein, with
  PDB IDs, author chains and ten compact coverage/element counts;
- `Human_Proteome_RCSB_Secondary_Structure.csv.gz` — one row per
  UniProt/PDB-chain/helix-or-strand observation, retaining PDB polymer
  coordinates, depositor coordinates and mapped UniProt ranges;
- `Human_Proteome_RCSB_Secondary_Structure.metadata.json` — releases, source
  URLs, coordinate conventions and QC totals.

Elements in affinity tags or other construct segments outside the SIFTS
UniProt mapping are not attributed to the protein. Elements that cross a
mapping boundary are retained and flagged `partial`; unequal-length mapping
segments are never forced into approximate coordinates.

Append the compact summary to the existing 242-column extended table with one
streaming pass:

```bash
python append_rcsb_to_extended.py
```

The normalized element table remains separate and joins on `uniprot_id`; this
avoids placing millions of structure observations into single spreadsheet
cells and duplicating them across every isoform row.

---

## Layout

| File | Role |
|---|---|
| `schema.py` | the 157 columns, in order; validates itself against the master header |
| `paths.py` | locations, interpreters, source-file provenance; `python paths.py` to check |
| `setup_environment.py` | builds the two venvs and vendors the PSLab predictor |
| `cider.py` `idr.py` `domains.py` `go.py` `cdcode.py` `string_ppi.py` `opentargets.py` `pslab.py` `identity.py` `rcsb.py` | one per feature family |
| `stage_*.py` | run one family over the whole target set, writing an intermediate |
| `build_rows.py` | assemble intermediates into complete 157-column rows |
| `append_to_master.py` | stream master + new rows into the expanded table |
| `run_all.py` | drive every stage in order, under the right interpreter |
| `validate.py` `validate_idr.py` `validate_pslab.py` `validate_rcsb.py` | fidelity and coordinate-mapping checks |

Stages are resumable: each is skipped when its intermediate already exists.

## Environments

Three, because metapredict and PSLab pin incompatible dependencies:

- **base** — localcider, pandas, numpy, biopython, openpyxl. Everything except IDRs and PSLab.
- **`envs/metapredict`** — metapredict 3.0.2 + CPU torch.
- **`envs/pslab`** — scikit-learn **1.6 exactly** (the shipped `.joblib` models will not unpickle cleanly on newer), MDAnalysis, numba.

`setup_environment.py` creates both, under a deliberately short path
(`%TEMP%\isoenvs` by default, override with `EXPANSION_ENVS`) rather than inside
this package. That is not tidiness: torch nests files about 160 characters deep,
Windows caps full paths at 260 measured from the drive letter, and a package
living in a synced OneDrive folder is already ~55 characters in. Installing into
`<package>/envs/` fails with *WinError 206: The filename or extension is too
long*.

Two further Windows-only metapredict problems, both handled automatically:

- its source distribution contains test fixtures with `|` in the filename, which
  Windows forbids — the script unpacks the archive itself, skipping them;
- it wants MSVC to compile a Cython accelerator — the script drops the extension
  from the build, and `idr.install_python_fallback()` routes to metapredict's own
  pure-Python implementation of the same routine, which metapredict's test suite
  (`tests/test_cython_v_python.py`) asserts is identical. Only speed changes.

On Linux and macOS the wheel installs normally and none of this applies.

## Known limitations

- `ENSP`, `ID`, `ProteinHGVS` are reconstructed, not recovered (see above).
- Coverage is bounded by what the proteins have: of the 2,901, 90.1% resolve to
  an ENSG and 67.0% to an ENSP. The rest are largely immunoglobulin and
  T-cell-receptor segments, olfactory receptors and uncharacterised ORFs that
  Ensembl does not carry as distinct genes, so their STRING and Open Targets
  columns are legitimately empty.
- CD-CODE membership depends on file ordering (see family 7).
- Non-canonical residues are substituted, not flagged — but only 4 of the 2,901
  sequences contain any, all selenoproteins (13 `U` positions, no `X`/`B`/`Z`).
  Listed in `Missing_Proteins_NonCanonical_Residues.csv`. The "346 sequences,
  3,176 X positions" warning in the old notebook described the older 5,476-set;
  those proteins have since been added to the table.
