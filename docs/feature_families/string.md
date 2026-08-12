# STRING protein interactions

`annotate_string.py` streams the human STRING v12 full-links file and retains
partners for ENSP identifiers represented in the selected catalog. Partner
scores remain grouped by query ENSP, with a conservative UniProt translation
when one parent mapping is unambiguous.

Because the set of represented proteins affects which partners can be
translated into in-dataset UniProt IDs, subset and full-proteome runs need not
produce identical translated-partner coverage.

Sanity checks: query ENSPs are version-normalized; scores are numeric and within
STRING bounds; duplicate query-partner pairs are absent; sidecar keys are
unique.

The report plots STRING field coverage and the number of retained ENSP partners
per row. Record STRING version with every data release.
