"""eCLIP feature family -- where each RBP binds on RNA.

Source
------
``rbp_master_with_eclip.csv``, built by the project's ``join_eclip_to_master.py``
from five independent CLIP resources. This module does not recompute any of
that; it exposes the finished 72 columns for joining onto the isoform table.

Granularity mismatch -- read this before using it
--------------------------------------------------
The eCLIP table is **per gene**, keyed on HGNC-style symbol, and covers only the
1,392 genes of the RBP census. The isoform table is per isoform and covers
20,483 proteins. So:

  * a value attaches to every isoform row of a gene, not to a specific isoform --
    eCLIP peaks are genomic and were never resolved to individual protein
    isoforms;
  * roughly 92% of rows get nothing, because their gene has no CLIP data. That
    is absence of measurement, not absence of binding, and the ``has_*`` flags
    exist to keep the two distinguishable.

Why the columns are split by source rather than merged
------------------------------------------------------
ENCODE eCLIP is a subset of both ENCORI and POSTAR3 -- POSTAR3 site ids
literally contain ``human_RBP_eCLIP_ENCODE_*`` -- so pooling them would be
pseudo-replication. Keeping them apart also preserves real disagreements: AGO2
is 0.52 3'UTR in ENCORI but 0.25 in POSTAR3, and that gap survives significance
filtering.

Which columns to analyse
------------------------
Use ``*_frac_*`` and ``skipper_enrich_*``. **Do not compare ``*_n_*`` counts
across sources.** ``encori_n_regions`` mostly tracks how many datasets exist for
an RBP (TARDBP has 100 ENCORI datasets and 1.24M sites; a once-profiled RBP has
a few thousand), and ``encode_n_peaks`` tracks sequencing depth. The same
nominal p <= 1e-3 threshold retains 9.1% of ENCODE but 97.3% of ENCORI, because
each source arrives pre-filtered to a different depth. The counts are provenance.

The `enrich_*` columns matter for a specific reason: raw overlap alone shows
nearly every RBP touching a 3'UTR somewhere, which is a statement about genomic
real estate rather than biology. Enrichment against background is what makes the
fraction interpretable.
"""

import csv
import os

csv.field_size_limit(1 << 30)

DEFAULT_TABLE = "rbp_master_with_eclip.csv"

#: Join key: the isoform table's ``Name`` against the eCLIP table's ``Name``,
#: both HGNC-style gene symbols. 1,391 of the 1,392 eCLIP genes match.
JOIN_COLUMN = "Name"

SOURCE_PREFIXES = ("encode_", "encori_", "postar_", "skipper_")

REGIONS = ["utr3", "utr5", "cds", "ncrna_exon", "intron", "intergenic"]


def _is_eclip_column(name):
    return name.startswith(SOURCE_PREFIXES) or name.startswith("has_") \
        or name == "n_eclip_sources"


def load(path=None, kappel_dir=None):
    """Read the eCLIP table -> ``({gene_symbol: {column: value}}, column_names)``.

    Only the 72 eCLIP-derived columns are kept; the other 877 columns of that
    file are the RBP census and are already represented in the isoform table.
    """
    if path is None:
        path = os.path.join(kappel_dir or ".", DEFAULT_TABLE)

    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        columns = [c for c in reader.fieldnames if _is_eclip_column(c)]
        by_gene = {}
        for row in reader:
            symbol = row.get(JOIN_COLUMN)
            if symbol:
                by_gene[symbol] = {c: row.get(c, "") for c in columns}
    return by_gene, columns


def columns_for(gene_symbol, by_gene, columns):
    """The eCLIP columns for one protein, keyed by its gene symbol.

    A gene with no CLIP data gets an empty string in every column, and its
    ``has_*`` flags stay empty too -- deliberately distinct from a measured 0,
    which would claim the RBP was profiled and found not to bind.
    """
    record = by_gene.get(gene_symbol)
    if record is None:
        return {c: "" for c in columns}
    return {c: record.get(c, "") for c in columns}


def dominant_region(record, source="encori_published"):
    """Convenience: the region holding the largest fraction of an RBP's sites.

    Returns None when that source has no data for the gene. The stored
    ``*_dominant`` column already carries this; this recomputes it from the
    fractions so a different source or region subset can be asked for.
    """
    best, best_value = None, -1.0
    for region in REGIONS:
        raw = record.get(f"{source}_frac_{region}", "")
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if value > best_value:
            best, best_value = region, value
    return best
