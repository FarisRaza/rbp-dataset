# Post-translational modifications

`annotate_ptm.py` covers acetylation; N-, O-, C-, and S-glycosylation;
methylation; myristoylation; phosphorylation; sumoylation; ubiquitination; and
S-nitrosylation.

Canonical UniProt sites are joined directly. For NCBI isoforms, sequence
alignment projects a site only when the modified residue is conserved; deleted
or changed residues are counted as dropped rather than silently copied.

Sanity checks: indicator columns agree with nonempty position/residue lists;
positions and residues have equal lengths; coordinates are zero-based and in
bounds; projection source and method are recorded.

The generated report compares PTM-class coverage and the number of PTM classes
per row. See `DATA_SOURCES.md` for the known gap in recreating the complete
historical wide PTM source.
