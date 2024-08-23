import pandas as pd
import numpy as np
import sys

chr=sys.argv[1]

df=pd.read_csv(f'potato_sift_annots_{chr}.bed',sep='\t',header=None,
               names=['#chr','start','end','ref_allele','alt_allele','ref_amino','alt_amino','variant_type','amino_pos','sift_score','sift_median','num_seqs','sift_prediction'])

df = df.sort_values(by=['#chr', 'start', 'end', 'ref_allele', 'alt_allele', 'sift_score'])

# Find the index of the minimum value in col 10 for each group
idx = df.groupby(['#chr', 'start', 'end', 'ref_allele', 'alt_allele'])['sift_score'].idxmin()
idx = idx.dropna()

# Select the rows with the minimum value in col 10
df_result = df.loc[idx]

df_result.to_csv(f'FINAL_potato_{chr}_filtered.bed.gz',sep='\t',header=True,index=False,compression='gzip')