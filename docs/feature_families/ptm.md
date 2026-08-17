# Post-translational modifications

`annotate_ptm.py` covers acetylation; N-, O-, C-, and S-glycosylation;
methylation; myristoylation; phosphorylation; sumoylation; ubiquitination; and
S-nitrosylation.

Canonical UniProt sites are joined directly. Sequence-distinct NCBI isoforms
remain null because the source observations are defined on canonical UniProt
coordinates. This avoids presenting inferred alignment projections as sourced
experimental annotations.

Sanity checks: indicator columns agree with nonempty position/residue lists;
positions and residues have equal lengths; coordinates are zero-based and in
bounds; canonical source and annotation scope are recorded. Alignment
projection remains available in `rbp_pipeline/ptm.py` as an explicit analysis
utility, but it is disabled in the master-table build.

The generated report compares PTM-class coverage and the number of PTM classes
per row. See `DATA_SOURCES.md` for the known gap in recreating the complete
historical wide PTM source.
