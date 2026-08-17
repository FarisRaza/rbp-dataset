# Column reference

This reference integrates the definitions from `Isoform Table Documentation
(1).xlsx` with the current clean rebuild schema. Wording has been normalized,
typos have been corrected, and current names are used where the rebuild made
an ambiguous historical field explicit.

The workbook's metadata counts and its broad historical Open Targets target
dump are not release specifications. Catalog size is determined by the pinned
source releases, and Open Targets now contains normalized tissue-expression and
disease/condition annotations only.

Nested lists and dictionaries are native values in the clean Parquet catalog.
Feature sidecars store nested cells as strict JSON strings so they have the same
representation in CSV and Parquet. Coordinates are zero-based and half-open
unless a field says otherwise.

## Identifiers and row provenance

The first table lists the compact, analysis-facing aliases added to the final
table. The clean catalog retains the more explicit fields in the second table.

| Final column | Description |
|---|---|
| `uniprot_id` | Reviewed UniProtKB/Swiss-Prot parent accession. |
| `dominant_isoform` | `1` only when this row's amino-acid sequence exactly matches the canonical Swiss-Prot sequence of a mapped UniProt parent; otherwise `0`. |
| `sequence` | Amino-acid sequence represented by the row. |
| `UNIQUE` | RBP-census gene label retained from the eCLIP compilation when available. This is not the stable row key. |
| `ProteinHGVS` | Comma-separated RefSeq `NP_` protein accessions represented by the sequence-deduplicated row. |
| `HGVSDescription` | Available NCBI isoform/product names. Despite the historical name, this is not a fully normalized HGVS expression. |
| `ENSG` | Semicolon-separated Ensembl gene identifiers. |
| `ENST` | Semicolon-separated Ensembl transcript identifiers. |
| `ENSP` | Semicolon-separated Ensembl protein identifiers. |
| `ncbi_gene_id` | Primary NCBI Gene identifier. |
| `protein_key` | Stable row key used to join every feature sidecar. |
| `ID` | Compatibility alias for `ENSG`. |
| `Name` | Gene symbol. |
| `Description` | Gene/protein description. |
| `isoform_number` | Row number within a UniProt group after canonical-first sorting. |
| `ID_list` | JSON representation of all ENSG values for compatibility with historical analyses. |

| Clean catalog column(s) | Description |
|---|---|
| `row_kind` | `swissprot_canonical` or sequence-distinct `ncbi_isoform`. |
| `length_aa`, `sequence_sha256`, `sequence_source` | Sequence length, exact-sequence hash, and source provenance. |
| `tax_id` | NCBI taxonomy identifier; human rows use 9606. |
| `gene_symbol`, `gene_synonyms`, `gene_description` | Gene-level names and description. |
| `ncbi_gene_ids`, `hgnc_ids` | All mapped NCBI Gene and HGNC identifiers. |
| `uniprot_secondary_accessions`, `uniprot_entry_name` | UniProt secondary accessions and entry name. |
| `uniprot_parent_ids`, `uniprot_isoform_ids` | Mapped canonical parent accessions and UniProt isoform identifiers. |
| `swissprot_canonical_accessions`, `is_swissprot_canonical` | Canonical-accession evidence and row flag. |
| `refseq_protein_ids`, `refseq_transcript_ids`, `ncbi_isoform_names` | RefSeq proteins, coding transcripts, and NCBI product names represented by the row. |
| `ensembl_gene_ids`, `ensembl_transcript_ids`, `ensembl_protein_ids` | Explicit Ensembl identifier lists. |
| `identifier_mapping_methods`, `identifier_ambiguity`, `canonical_match_method` | How mappings were made, any unresolved ambiguity, and how canonical equivalence was established. |
| `swissprot_release`, `ncbi_annotation_release`, `ensembl_release`, `build_timestamp_utc` | Source-release and build provenance. |

## CIDER: whole sequence

These localCIDER values are computed for every amino-acid sequence.

| Column | Description |
|---|---|
| `FCR` | Fraction of charged residues. |
| `NCPR` | Net charge per residue. |
| `isoelectric_point` | Predicted isoelectric point. |
| `molecular_weight` | Predicted molecular weight in daltons. |
| `countNeg` | Number of negatively charged D/E residues. |
| `countPos` | Number of positively charged R/K residues. |
| `countNeut` | Number of residues treated as neutral. |
| `fraction_negative` | Fraction of negatively charged residues. |
| `fraction_positive` | Fraction of positively charged residues. |
| `fraction_expanding` | Fraction of E/D/R/K/P residues predicted to promote chain expansion. |
| `amino_acid_fractions` | Mapping from amino-acid code to sequence fraction. |
| `fraction_disorder_promoting` | Fraction of residues classified as disorder-promoting. |
| `mean_net_charge` | Absolute mean net charge. |
| `mean_hydropathy` | Mean hydropathy on the skewed Kyte-Doolittle scale used by localCIDER. |
| `uversky_hydropathy` | Mean hydropathy on the normalized Uversky/Kyte-Doolittle scale. |
| `PPII_propensity` | Overall polyproline-II propensity. |
| `kappa` | Charge-patterning/segregation score. |
| `delta` | Sequence charge-patterning delta score. |
| `deltaMax` | Maximum possible delta for the sequence composition. |
| `cider_sequence_sanitization` | Exact policy used for noncanonical residues while scoring. |

## Metapredict IDRs and folded regions

| Column | Description |
|---|---|
| `IDR_count` | Number of predicted intrinsically disordered regions. |
| `IDR_avg_size` | Mean IDR length in amino acids. |
| `IDR_total_size` | Total number of residues assigned to IDRs. |
| `IDR_range` | Ordered `[start, end]` ranges for IDRs. |
| `IDR_discrete_seq` | Ordered list containing each IDR sequence separately. |
| `IDR_concat_seq` | Concatenation of all IDR sequences. |
| `FOLD_count` | Number of complementary non-IDR/folded regions. |
| `FOLD_avg_size` | Mean folded-region length in amino acids. |
| `FOLD_total_size` | Total number of residues assigned to folded regions. |
| `FOLD_range` | Ordered `[start, end]` ranges for folded regions. |
| `FOLD_discrete_seq` | Ordered list containing each folded-region sequence separately. |
| `FOLD_concat_seq` | Concatenation of all folded-region sequences. |
| `idr_method` | Predictor/version and threshold used. |

## CIDER on IDRs

Each `IDR_*` value below is a list ordered to match `IDR_discrete_seq`.

| Column | Per-IDR value |
|---|---|
| `IDR_FCR` | Fraction of charged residues. |
| `IDR_NCPR` | Net charge per residue. |
| `IDR_isoelectric_point` | Predicted isoelectric point. |
| `IDR_molecular_weight` | Molecular weight in daltons. |
| `IDR_countNeg`, `IDR_countPos`, `IDR_countNeut` | Negative, positive, and neutral residue counts. |
| `IDR_fraction_negative`, `IDR_fraction_positive` | Negative and positive residue fractions. |
| `IDR_fraction_expanding` | Chain-expanding residue fraction. |
| `IDR_amino_acid_fractions` | Amino-acid composition mapping. |
| `IDR_fraction_disorder_promoting` | Disorder-promoting residue fraction. |
| `IDR_kappa` | Charge-segregation score. |
| `IDR_mean_net_charge` | Absolute mean net charge. |
| `IDR_mean_hydropathy`, `IDR_uversky_hydropathy` | Mean hydropathy on the localCIDER scales. |
| `IDR_PPII_propensity` | Polyproline-II propensity. |
| `IDR_delta`, `IDR_deltaMax` | Delta and maximum-delta charge-patterning scores. |

## Curated UniProt domains

These fields come from Swiss-Prot `DOMAIN` and `ZN_FING` features and are null
for sequence-distinct NCBI isoforms.

| Column | Description |
|---|---|
| `Domains` | Ordered names of curated folded domains. |
| `Domains_count` | Number of occurrences grouped by domain name. |
| `Domains_avg_size` | Mean region length grouped by domain name. |
| `Domains_total_size` | Total annotated length grouped by domain name. |
| `Domains_range` | Ranges grouped by domain name. |
| `Domains_discrete_seq` | Separate sequence for every region, grouped by domain name. |
| `Domains_concat_seq` | Concatenated region sequence for each domain name. |
| `uniprot_domain_annotation_scope` | Whether coordinates are direct canonical annotations or not applicable to the row. |

## CIDER on curated domains

The following are dictionaries of lists keyed by domain name, with inner values
ordered like `Domains_discrete_seq`: `Domains_FCR`, `Domains_NCPR`,
`Domains_isoelectric_point`, `Domains_molecular_weight`, `Domains_countNeg`,
`Domains_countPos`, `Domains_countNeut`, `Domains_fraction_negative`,
`Domains_fraction_positive`, `Domains_fraction_expanding`,
`Domains_amino_acid_fractions`, `Domains_fraction_disorder_promoting`,
`Domains_kappa`, `Domains_mean_net_charge`, `Domains_mean_hydropathy`,
`Domains_uversky_hydropathy`, `Domains_PPII_propensity`, `Domains_delta`, and
`Domains_deltaMax`. Their definitions are the same as the whole-sequence CIDER
metrics. `Domains_Omega` is the localCIDER Omega patterning value for each
domain region.

## STRING protein interactions

| Current column | Description |
|---|---|
| `string_query_ensp_ids` | Version-normalized ENSP identifiers queried for the row. |
| `string_partners_ensp_by_query` | Mapping from query ENSP to partner ENSP and STRING confidence score. |
| `string_partners_uniprot_by_query` | Conservative translation of partner ENSPs to unambiguous UniProt parents. |
| `string_partners_ensp_in_catalog_by_query` | Partner ENSPs represented in the selected catalog. |
| `string_version` | STRING release used. |

The final table also exposes compatibility aliases `ENSP_clean`,
`PPI_ENSP_Partners`, `PPI_UniProt_Partners`,
`PPI_ENSP_Partners_in_Dataframe`, and
`PPI_UniProt_Partners_in_Dataframe` for the corresponding values above.

## CD-CODE condensates

All ten values are positionally parallel lists. Element *i* in every column
describes the same condensate.

| Column | Description |
|---|---|
| `Condensate Name` | Condensate names. |
| `UID` | CD-CODE condensate identifiers. |
| `Condensate Type` | Biomolecular or synthetic classification. |
| `Species Tax Id` | Taxonomy identifier for the condensate record. |
| `Proteins` | Number/indicator describing protein involvement in the source record. |
| `DNA` | DNA involvement in condensate formation. |
| `RNA` | RNA involvement in condensate formation. |
| `C-mods` | Condensate modification metadata supplied by CD-CODE; the historical workbook left this definition unresolved. |
| `Condensatopathy` | Disease/condensatopathy implication. |
| `Confidence Score` | Source evidence-confidence score. |
| `cdcode_annotation_scope` | Direct canonical membership or non-applicable isoform scope. |

## Gene Ontology

| Column | Description |
|---|---|
| `C_ids`, `C_descriptions`, `C_evidence` | Parallel lists of cellular-component GO IDs, names, and evidence codes. |
| `P_ids`, `P_descriptions`, `P_evidence` | Parallel lists of biological-process GO IDs, names, and evidence codes. |
| `F_ids`, `F_descriptions`, `F_evidence` | Parallel lists of molecular-function GO IDs, names, and evidence codes. |
| `go_annotation_scope` | Direct canonical UniProt annotation or non-applicable isoform scope. |
| `go_source_uniprot_ids` | UniProt accession(s) from which the GO annotations came. |

## Open Targets expression and disease associations

These fields replace the historical mixture of `diseaseId`, `datatypeId`,
`score`, `evidenceCount`, `tissues`, and 28 largely undocumented target-profile
columns. They are gene-level and may repeat across rows sharing an ENSG.
Version suffixes are removed from catalog ENSGs before joining because Open
Targets keys targets by versionless Ensembl gene identifier.

| Column | Description |
|---|---|
| `opentargets_tissue_expression` | Normalized tissue records including ENSG, tissue identifiers/names, organ and anatomical-system labels, RNA value/unit/z-score/level, and protein level/reliability/cell types. |
| `opentargets_expression_tissue_count` | Number of unique ENSG/tissue records. |
| `opentargets_disease_associations` | Normalized target-condition records retaining disease ID/name/description, therapeutic areas, datatype-specific scores, and evidence counts. |
| `opentargets_disease_count` | Number of unique associated disease/condition identifiers. |
| `opentargets_disease_names` | Deduplicated disease ID/name pairs. |
| `opentargets_therapeutic_areas` | Deduplicated therapeutic-area ID/name pairs. |
| `opentargets_release` | Pinned Open Targets Platform release. |
| `opentargets_annotation_scope` | ENSG join behavior and interpretation of absent records. |

The feature build also writes normalized long tables for tissue records and
disease/datatype records so individual observations can be analyzed without
unpacking nested cells.

## PSLab phase-separation predictions

Every column is a list ordered to match `IDR_discrete_seq`.

| Column | Description |
|---|---|
| `mean_lambda` | Mean per-residue stickiness. |
| `faro` | Aromatic-residue fraction used by the model. |
| `shd` | Sequence hydropathy decoration. |
| `ncpr` | Net charge per residue used by PSLab. |
| `fcr` | Fraction of charged residues used by PSLab. |
| `scd` | Sequence charge decoration. |
| `ah_ij` | Pairwise interaction-potential summary. |
| `nu_svr` | Predicted polymer scaling exponent. |
| `Delta G [kT]` | Predicted free energy for dilute-to-dense phase transition. |
| `Saturation concentration [mg/mL]` | Predicted saturation concentration in mg/mL. |
| `Saturation concentration [uM]` | Predicted saturation concentration in micromolar. |
| `pslab_annotation_scope` | States that predictions are made once per metapredict IDR in the same order. |

The workbook contained a provisional `tss` definition (“total sequence
stickiness”), but `tss` is not emitted by the current PSLab feature block.

## eCLIP and related CLIP evidence

The exact eCLIP columns follow the supplied compilation. Names are retained so
ENCODE, ENCORI, POSTAR3, and Skipper remain analytically separate.

| Column pattern | Description |
|---|---|
| `has_*` | Whether the corresponding source/evidence subset measured the gene. |
| `*_n_*` | Source-specific site, peak, dataset, or region counts; these are provenance/depth indicators and should not be compared directly across sources. |
| `*_frac_<region>` | Fraction of source observations in 3' UTR, 5' UTR, CDS, noncoding exon, intron, or intergenic regions. |
| `*_dominant` | Region with the largest source-specific fraction. |
| `skipper_enrich_*` | Skipper region enrichment relative to its background model. |
| `n_eclip_sources` | Number of represented CLIP evidence sources. |
| `rbp_census_unique` | RBP-census label carried by the source compilation. |
| `eclip_annotation_scope` | Gene-level broadcast and measured/absent interpretation. |

## InterProScan

| Column | Description |
|---|---|
| `InterPro_domains` | Retained signature/domain names. |
| `InterPro_count` | Hit counts grouped by signature/domain name. |
| `InterPro_range` | Sequence ranges grouped by signature/domain name. |
| `InterPro_accessions` | InterPro accessions grouped by signature/domain name. |
| `InterPro_databases` | Member databases supporting each retained hit. |
| `InterPro_go_terms` | GO terms supplied by InterProScan for retained hits. |
| `InterPro_n_hits` | Total number of retained hit ranges. |
| `interpro_source_refseq_protein` | RefSeq protein accession used for the direct join. |
| `interpro_version` | InterProScan/parser release. |
| `interpro_annotation_scope` | Direct match or explanation of missing/indistinguishable no-hit state. |

## Post-translational modifications

For each PTM type—acetylation, N/O/C/S-glycosylation, methylation,
myristoylation, phosphorylation, sumoylation, ubiquitination, and
S-nitrosylation—the table contains three columns:

| Pattern | Description |
|---|---|
| `ptm_<type>` | Binary indicator that at least one canonical site of this type exists. |
| `ptm_<type>_positions` | Zero-based canonical-sequence residue positions. |
| `ptm_<type>_residues` | Residue codes parallel to the position list. |
| `ptm_projection_source_uniprot_ids` | Canonical UniProt source accessions. |
| `ptm_projection_methods` | Direct/projection method provenance. The master build uses direct canonical annotation only. |
| `ptm_projection_dropped_sites` | Number of sites dropped because their archived residue/coordinate no longer matches the current canonical sequence, or during an explicit opt-in isoform projection. |
| `ptm_coordinate_system` | Coordinate convention. |
| `ptm_annotation_scope` | Direct canonical annotation or non-applicable isoform scope. |

## GO-derived regulatory roles

| Column | Description |
|---|---|
| `role_in_transcription` | `1` when a GO biological-process or molecular-function name matches the documented transcription rule. |
| `role_in_translation` | `1` when a term matches translation-regulation rules. |
| `role_in_mrna_stability` | `1` when a term matches mRNA stabilization/destabilization rules. |
| `role_in_translation_stability` | Logical OR of translation and mRNA-stability flags, retained for compatibility. |
| `go_role_annotation_scope` | Notes that labels are derived from canonical GO term names. |

## RCSB/PDB structures

| Column | Description |
|---|---|
| `RCSB_PDB_IDs` | All SIFTS-mapped PDB entry identifiers for the canonical UniProt protein. |
| `RCSB_PDB_count` | Number of unique mapped PDB entries. |
| `RCSB_PDB_chains` | PDB-chain mapping records. |
| `RCSB_PDB_entries_with_secondary_structure_count` | Mapped entries with retrieved secondary-structure annotations. |
| `RCSB_PDB_entries_without_secondary_structure_count` | Mapped entries lacking retrieved secondary-structure annotations. |
| `RCSB_secondary_structure_observation_count` | Total retained regular secondary-structure observations. |
| `RCSB_helix_observation_count` | Retained helix observations. |
| `RCSB_beta_strand_observation_count` | Retained beta-strand observations. |
| `RCSB_secondary_structure_complete_mapping_count` | Observations with complete SIFTS-to-UniProt coordinate mapping. |
| `RCSB_secondary_structure_partial_mapping_count` | Observations with partial mapping. |
| `rcsb_annotation_scope` | Direct canonical mapping, absent canonical mapping, or non-inheritance to an isoform. |
