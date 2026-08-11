import pandas as pd
import pyranges as pr
import glob
import gzip
import os

# ============================================================
# STEP 1: Build CDS and intron beds from GENCODE GTF
# (5UTR/3UTR beds already exist from the ENCORI run - reused as-is)
# ============================================================
print("Step 1: Building CDS and intron beds...")

gtf_cols = ['Chromosome', 'source', 'feature', 'Start', 'End', 'score',
            'Strand', 'frame', 'attributes']
gtf = pd.read_csv('gencode.gtf.gz', sep='\t', header=None, names=gtf_cols,
                   comment='#', compression='gzip')
gtf['transcript_id'] = gtf['attributes'].str.extract(r'transcript_id "([^"]+)"')

# CDS bed
cds = gtf[gtf['feature'] == 'CDS'][['Chromosome', 'Start', 'End', 'Strand']].copy()
cds['Start'] = cds['Start'] - 1
cds.to_csv('gencode_cds.bed', sep='\t', header=False, index=False)
print(f"Wrote gencode_cds.bed: {len(cds)} rows")

# Intron bed: transcript span minus exons, per transcript
transcript = gtf[gtf['feature'] == 'transcript'][['Chromosome', 'Start', 'End', 'Strand']].copy()
transcript['Start'] = transcript['Start'] - 1
exon = gtf[gtf['feature'] == 'exon'][['Chromosome', 'Start', 'End', 'Strand']].copy()
exon['Start'] = exon['Start'] - 1

transcript_gr = pr.PyRanges(transcript)
exon_gr = pr.PyRanges(exon)
intron_gr = transcript_gr.subtract(exon_gr)
intron_gr.df.to_csv('gencode_intron.bed', sep='\t', header=False, index=False)
print(f"Wrote gencode_intron.bed: {len(intron_gr.df)} rows")

# ============================================================
# STEP 2: Load all region beds, merge each
# ============================================================
print("\nStep 2: Loading region beds...")

def load_merged_bed(fname, extra_cols=None):
    cols = ['Chromosome', 'Start', 'End', 'Strand'] if extra_cols is None else extra_cols
    df = pd.read_csv(fname, sep='\t', header=None)
    df.columns = cols[:df.shape[1]]
    return pr.PyRanges(df[['Chromosome', 'Start', 'End', 'Strand']]).merge()

regions = {
    '5UTR': load_merged_bed('gencode_5utr.bed', ['Chromosome','Start','End','utr_type','score','Strand']),
    '3UTR': load_merged_bed('gencode_3utr.bed', ['Chromosome','Start','End','utr_type','score','Strand']),
    'CDS': load_merged_bed('gencode_cds.bed'),
    'intron': load_merged_bed('gencode_intron.bed'),
}

# ============================================================
# STEP 3: Load ENCODE peaks, tag each with RBP name from filename
# ============================================================
print("\nStep 3: Loading ENCODE eCLIP peak files...")

peak_files = glob.glob('encode_peaks/*.bed.gz')
print(f"Found {len(peak_files)} peak files")

all_peaks = []
for fpath in peak_files:
    rbp_name = os.path.basename(fpath).split('_')[0]
    cols = ["Chromosome","Start","End","Name","Score","Strand",
            "signalValue","pValue","qValue","peak"]
    df = pd.read_csv(fpath, sep='\t', header=None, names=cols, compression='gzip')
    df['RBP'] = rbp_name
    all_peaks.append(df[['Chromosome','Start','End','Strand','RBP','signalValue','pValue']])

peaks = pd.concat(all_peaks, ignore_index=True)
print(f"Total peaks loaded: {len(peaks)}")

# ============================================================
# STEP 4: Classify each peak by region type
# ============================================================
print("\nStep 4: Classifying peaks by region...")

# downcast to shrink memory footprint before building PyRanges
peaks['Start'] = peaks['Start'].astype('int32')
peaks['End'] = peaks['End'].astype('int32')
peaks['Chromosome'] = peaks['Chromosome'].astype('category')
peaks['RBP'] = peaks['RBP'].astype('category')
peaks['signalValue'] = peaks['signalValue'].astype('float32')
peaks['pValue'] = peaks['pValue'].astype('float32')

# build once, reuse for all 4 region joins instead of rebuilding each time
peaks_gr = pr.PyRanges(peaks)

annotated_chunks = []
for region_name, region_gr in regions.items():
    print(f"  Joining against {region_name}...")
    hit = peaks_gr.join(region_gr)
    if len(hit) > 0:
        hit_df = hit.df.drop_duplicates(subset=['Chromosome','Start','End','RBP'])
        hit_df['region'] = region_name
        annotated_chunks.append(hit_df[['Chromosome','Start','End','RBP','region','signalValue','pValue']])
    del hit

del peaks_gr

peak_region_table = pd.concat(annotated_chunks, ignore_index=True)
peak_region_table.to_csv('encode_rbp_peak_region_annotation.csv', index=False)
print(f"Wrote encode_rbp_peak_region_annotation.csv: {len(peak_region_table)} rows")

# confidence-filtered version: signalValue > 1 (2-fold+ enrichment over SMInput),
# pValue > 3 (CLIPper reports -log10(p), so this is p < 0.001)
confident_table = peak_region_table[(peak_region_table['signalValue'] > 1) &
                                     (peak_region_table['pValue'] > 3)]
confident_table.to_csv('encode_rbp_peak_region_annotation_confident.csv', index=False)
print(f"Wrote encode_rbp_peak_region_annotation_confident.csv: {len(confident_table)} rows "
      f"({len(confident_table)/len(peak_region_table):.1%} of all classified peaks)")

# ============================================================
# STEP 5: Per-RBP region summary (both unfiltered and confidence-filtered)
# ============================================================
print("\nStep 5: Building per-RBP region summary...")

region_counts = peak_region_table.groupby(['RBP','region']).size().unstack(fill_value=0)
total_per_rbp = peaks.groupby('RBP').size().rename('total_peaks')
region_counts = region_counts.join(total_per_rbp)
region_fracs = region_counts.drop(columns='total_peaks').div(region_counts['total_peaks'], axis=0)
region_fracs['total_peaks'] = region_counts['total_peaks']
region_fracs.to_csv('encode_rbp_region_fractions.csv')

print("Done. Top 20 RBPs by 3'UTR fraction (unfiltered):")
if '3UTR' in region_fracs.columns:
    print(region_fracs.sort_values('3UTR', ascending=False).head(20))

# same summary, but only using confident (signalValue > 1, pValue > 3) peaks
confident_counts = confident_table.groupby(['RBP','region']).size().unstack(fill_value=0)
confident_total = confident_table.groupby('RBP').size().rename('total_peaks')
confident_counts = confident_counts.join(confident_total)
confident_fracs = confident_counts.drop(columns='total_peaks').div(confident_counts['total_peaks'], axis=0)
confident_fracs['total_peaks'] = confident_counts['total_peaks']
confident_fracs.to_csv('encode_rbp_region_fractions_confident.csv')

print("\nDone. Top 20 RBPs by 3'UTR fraction (confidence-filtered):")
if '3UTR' in confident_fracs.columns:
    print(confident_fracs.sort_values('3UTR', ascending=False).head(20))

