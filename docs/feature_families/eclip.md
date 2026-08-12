# eCLIP and CLIP-derived RNA-binding evidence

`annotate_eclip.py` joins the gene-level compilation built by
`eclip_source_pipeline/`: ENCODE, ENCORI, POSTAR3, and optional Skipper are kept
as separate evidence families to avoid treating overlapping experiments as
independent replication.

The source table is joined by gene symbol and broadcast to protein isoforms.
The scope column distinguishes measured genes from proteins lacking a supplied
CLIP record.

Sanity checks: source-specific availability flags agree with site/peak counts;
region fractions are finite and bounded; enrichment values are interpreted only
when their source has data; ENCODE and ENCORI columns are not pooled.

The generated report shows per-column/source coverage and a per-row CLIP signal
distribution. `summarize_eclip.py` additionally creates region-composition and
ENCODE-versus-ENCORI comparison figures for the historical source table.
