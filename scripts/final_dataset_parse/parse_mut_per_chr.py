# Imports
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from pard.grantham import grantham
import blosum as bl
import logging
import time
import sys 
import joblib
import subprocess
import random

def is_transition(row):
   return 1 if ((row['REF_ALLELE'] == 'A' and row['ALT_ALLELE'] == 'G') or
                (row['REF_ALLELE'] == 'G' and row['ALT_ALLELE'] == 'A') or
                (row['REF_ALLELE'] == 'C' and row['ALT_ALLELE'] == 'T') or
                (row['REF_ALLELE'] == 'T' and row['ALT_ALLELE'] == 'C')) else 0  

def grantham_score(row):
  # Check for both '.' and '*' in either column
  if '.' in (row['REF_AMINO'], row['ALT_AMINO']) or '*' in (row['REF_AMINO'], row['ALT_AMINO']):
    return 0
  else:
    return grantham(row['REF_AMINO'], row['ALT_AMINO'])
def blosum62_score(row):
  # Check for both '.' and '*' in either column
  if '.' in (row['REF_AMINO'], row['ALT_AMINO']) or '*' in (row['REF_AMINO'], row['ALT_AMINO']):
    return 0
  else:
    blosum62_mat=bl.BLOSUM(62, default=0)
    return (blosum62_mat[row['REF_AMINO']][row['ALT_AMINO']])
# A function to log messages
def log_message(message):
    current_time = time.localtime()
    time_string = time.strftime("%Y-%m-%d %H:%M:%S", current_time)
    logging.info(f'{time_string} : {message}')
    logging.getLogger().handlers[0].flush()


novel_path=str(sys.argv[1]) # Path for input (features with a label)
working_dir=str(sys.argv[2]) # Output dirctory for the parsed dataset
current_chr=str(sys.argv[3]) # Chromosome
scaler_name=str(sys.argv[4]) # Scaler from the training, for normalization
nov_or_neut=str(sys.argv[5]) # Novel or neutral dataset
path_muts='/home/labs/alevy/omerbar/features/arabidopsis/mutation_features'
logging.basicConfig(filename=f"{working_dir}/parse_muts.log", level=logging.INFO, format="%(message)s")
log_message(f'{current_chr}: Starting chr...')

# Neutral dataset
colnames=["#chr", "start", "end", "REF_ALLELE",
          "ALT_ALLELE_SIFT", "REF_AMINO","ALT_AMINO", "VARIANT_TYPE","AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","SIFT_PREDICTION","codon_pos", 
          "genomic_region_CDS", "genomic_region_exon","genomic_region_5UTR", "genomic_region_gene", "genomic_region_mRNA", "genomic_region_mRNA_TE_gene",
          "genomic_region_pseudogene ","genomic_region_3UTR", "genomic_region_tansposable_element_gene", "gene_count", "genomic_region_sRNA","genomic_region_intergenic",
          "kmers_13", "gc_35", "upstream_dist", "downstream_dist", "exon_count", "splice_junc", "repeat", "chromhmm",
          "TF_markers", "methylation_mark", "chrom_access_mark", "DNA_hypersensitivity_mark", "atacseq_mark",
          "H3K4me3", "H3K9me1", "H3K27me3", "H3K4me1", "H3K4me2", "H3K23ac", "H3K9ac", "H3K36ac", "H3K27ac", 
          "gerp_exp", "gerp_sub", "gerp_ta", "lrt_phylop_masked", "phast_est", "exp_median","score_phylop_masked","gerp_phylop_masked",
          "ALT_ALLELE", "label"]

# Novel dataset
df = pd.read_csv(f'{path_muts}/{novel_path}/{nov_or_neut}_features_chr_{current_chr}.bed.gz', sep='\t', names=colnames,compression='gzip')#,skiprows=1) # Removed extra feats that are descriminative.

log_message(f'{current_chr}: Novel set read, {len(df)} lines read.')
# Capitalize ref allele
df['REF_ALLELE']=df['REF_ALLELE'].str.upper()
# Filter out shitty ambiguous ref alleles
df = df[~df['REF_ALLELE'].isin(['K','M','W','S','D','Y','R'])]
# Filtering out '-' or 'N' mutations
df=df[(df.REF_ALLELE!='-') & (df.REF_ALLELE!='N')]

log_message(f'{current_chr}: Combining mutation datasets...')
# Combining datasets

# Removing the ALT_ALLELE_SIFT because we already have the alt allele for non cds as well
df=df.drop(columns=['ALT_ALLELE_SIFT'])
# Changing missing values for sift score to 2 to reflect this as well
#df['SIFT_PREDICTION'].fillna(2, inplace=True)
# Fixing shitty chromhmm
if len(df.loc[df["chromhmm"] == "-1", "chromhmm"])>0:
    log_message(f'mode for chromhmm: {df.loc[df["chromhmm"] != "-1", "chromhmm"].mode()[0]}')
    log_message(f'number of values changed: {len(df.loc[df["chromhmm"] == "-1", "chromhmm"])}')
    df.loc[df["chromhmm"] == "-1", "chromhmm"] = df.loc[df["chromhmm"] != "-1", "chromhmm"].mode()[0]
# Changing chromhmm to categorical to make it a dummy feature (it is discrete)
df["chromhmm"]=df["chromhmm"].astype('category')

# Adding features
# transition and transversion
df['is_transition'] = df.apply(is_transition, axis=1)
# grantham scores
df['Grantham_score'] = df.apply(grantham_score, axis=1)
# blosum 62 scores
df['Blosum62_score'] = df.apply(blosum62_score, axis=1)

# Getting dummies
# Bcause the chromosome names are ChrX, we need to trim them, so that get_dummies won't confuse it.. We will convert it back for the homogenousness of the analysis
df['#chr'] = df['#chr'].str[3:].astype(int)
df=pd.get_dummies(df,dtype=int)

  
# Lump UTRs
df['genomic_region_3UTR'] = (df['genomic_region_3UTR'] | df['genomic_region_5UTR']).astype(int)
# Lump Stop mutations
df['VARIANT_TYPE_STOP-LOSS'] = (df['VARIANT_TYPE_STOP-LOSS'] | df['VARIANT_TYPE_STOP-GAIN']).astype(int)
# Moving the label to the last column
df=df[[col for col in df.columns if col != 'label'] + ['label']]

# Saving not normalized dataset
df.to_csv(f'{working_dir}/df_chr_{current_chr}_not_normed.csv', sep='\t', header=True, index=False)

log_message(f'{current_chr}: Normalizing...')
cols_to_norm = ["AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","codon_pos", "gene_count", "kmers_13", "gc_35", "upstream_dist", "downstream_dist", "exon_count",
                "TF_markers", "methylation_mark", "chrom_access_mark", "DNA_hypersensitivity_mark", "atacseq_mark",
                "H3K4me3", "H3K9me1", "H3K27me3", "H3K4me1", "H3K4me2", "H3K23ac", "H3K9ac", "H3K36ac", "H3K27ac",
                "gerp_exp", "gerp_sub", "gerp_ta", "lrt_phylop_masked", "phast_est", "exp_median","score_phylop_masked","gerp_phylop_masked",
               ]
log_message(f'Columns normalized:{cols_to_norm}')

# Fitting standard scaler used on training data
scaler = joblib.load(scaler_name)

# Applying scalers to standarize and normalize data
df[cols_to_norm]=scaler.transform(df[cols_to_norm]).round(5)

# Adding back the 'Chr' before each chromosome
df['#chr'] = 'Chr' + df['#chr'].astype(str)
# Saving the positions
df_pos=df[['#chr','start','end']]
df.drop(columns=['#chr','start','end'],inplace=True)

log_message(f'{current_chr}: Saving datasets...')
# Saving datasets
df.to_csv(f'{working_dir}/df_chr_{current_chr}_normed.csv',sep='\t',header=True,index=False)
df_pos.to_csv(f'{working_dir}/df_chr_{current_chr}_pos.csv',sep='\t',header=True,index=False)