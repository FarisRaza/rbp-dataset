# Intrinsically disordered regions

`annotate_idr.py` runs metapredict V3 with the documented 0.5 disorder
threshold and stores disordered/folded region counts, coordinates, sequences,
and total/average lengths.

Column definitions: [metapredict IDRs and folded regions](../COLUMN_REFERENCE.md#metapredict-idrs-and-folded-regions).

The annotation is sequence-specific, so canonical and NCBI isoform rows are
predicted independently. No identifier inheritance is used.

Sanity checks: `IDR_total_size + FOLD_total_size` equals protein length;
coordinates are nonoverlapping and within sequence bounds; sequence-list lengths
equal their counts; and every catalog key occurs exactly once.

The generated report plots IDR column coverage and the per-row IDR-count
distribution. In the meeting-era expanded table, 18,053 of 20,483 unique
UniProt proteins had at least one IDR in any represented row; that historical
number is context, not a fixed expectation for a current rebuild.

Run `python qc/check_idr.py --work-dir WORK_DIR`.
