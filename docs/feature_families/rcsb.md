# RCSB/PDB structural coverage

`annotate_rcsb.py` joins compact SIFTS/PDBe summaries including PDB and chain
counts, mapped residue coverage, and secondary-structure composition.

This optional annotation is applied only to canonical Swiss-Prot rows because
experimental structures are mapped to those sequences; it is not inherited to
sequence-distinct NCBI isoforms.

Sanity checks: PDB and chain counts are nonnegative; coverage fractions are
bounded; secondary-structure residue counts do not exceed mapped residues;
canonical-only scope is respected.

The report plots structural-field coverage and PDB mappings per row. Sparse
coverage is biologically expected and should not be treated as pipeline failure.
