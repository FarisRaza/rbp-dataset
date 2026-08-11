# Build the ENCODE / ENCORI / POSTAR3 / Skipper CLIP source table

These scripts were recovered from the completed eCLIP analysis session. They
build the gene-level `rbp_master_with_eclip.csv` source consumed by
`../annotate_eclip.py`. Run them from one working directory because their
documented defaults intentionally exchange named intermediate files there.

## Recommended current workflow

1. Download a current GENCODE GRCh38 GTF as `gencode.gtf` and run
   `build_gencode_regions.py`. It writes merged, strand-aware 3'UTR, 5'UTR,
   CDS, noncoding-exon and intron tracks under `gencode_regions/`.
2. Run `download_encode_peaks_grch38.py`. It keeps all GRCh38 peak files,
   distinguishes individual replicates from IDR sets, and writes a manifest.
   `download_encode_peaks.py` is retained only for provenance and is superseded.
3. Build the two ENCODE families:

   ```bash
   python annotate_peaks.py --mode precedence
   python encode_reproducible.py \
     --min-neglog10p 3 --min-log2fc 3 --min-experiments 1 \
     --out encode_region_counts_published.csv
   python encode_reproducible.py \
     --min-neglog10p 3 --min-experiments 2 \
     --out encode_region_counts_matched.csv
   ```

4. Place the ENCORI site export at `sites_annotated.bed`. Build the published
   and matched-criterion families:

   ```bash
   python encori_published_criterion.py \
     --out encori_region_counts_published.csv
   python encori_reproducible.py \
     --max-log10p -3 --min-datasets 2 \
     --out encori_region_counts_matched.csv
   ```

5. Place the POSTAR3 site export at `postar3_sites.bed` and run:

   ```bash
   python annotate_db_sites.py --source postar3
   ```

6. Optional: place the published Skipper enrichment references and window
   archives under `skipper_published/`, then run `build_skipper_table.py`.
7. Place the RBP census/base table at `rbp_dataframe.csv` and run
   `join_eclip_to_master.py`. It writes `rbp_master_with_eclip.csv` plus an
   unmatched-name report.

The exact expected filenames, criteria and caveats are in each script's module
docstring. `compare_criteria.py`, `encori_calibrate.py`, and `compute_fdr.py` are
diagnostic/audit scripts; they are not required for the final join.

## Important interpretation

- ENCODE and ENCORI are kept as separate column families. Their experiments
  overlap, so pooling them as independent evidence would be pseudo-replication.
- Region fractions/enrichments are analysis variables. Raw site/peak counts are
  sequencing depth and dataset-availability provenance and are not comparable
  across sources.
- The final table is gene-level. `annotate_eclip.py` explicitly broadcasts a
  measured gene's values to all of its protein isoform rows and labels that
  scope.
- POSTAR3's score is not a usable significance statistic; its family is
  therefore marked unfiltered rather than given a fabricated cutoff.

Dependencies for this source-building layer are `requests`, `numpy`, `pandas`,
and optionally `statsmodels` for the FDR implementation check. The historical
`encode_region_pipeline.py` also uses `pyranges`; the recommended scripts above
do not require that legacy pipeline.

After the join, `../summarize_eclip.py` produces source coverage, region
composition, ENCODE–ENCORI agreement tables and meeting-ready PNG figures.
