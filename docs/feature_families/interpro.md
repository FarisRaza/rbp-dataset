# InterProScan domains

`annotate_interpro.py` parses correctly aligned 15-column InterProScan TSV
output and joins hits through RefSeq protein accessions. It retains signature
names, ranges, InterPro accessions, member databases, GO terms, and hit counts.

The default parser retains domain-model databases and excludes disorder/coils
and family-only classifications already represented elsewhere.

Sanity checks: all ranges are zero-based, half-open, and within sequence
bounds; `InterPro_n_hits` equals the number of retained ranges; source RefSeq
accessions are recorded; InterProScan version is preserved.

The generated report plots annotation coverage and retained hit counts. Release
or member-database changes can shift counts and should be noted when comparing
runs.
