# Imports
import pandas as pd
import os 
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import scipy.signal as signal
from pard.grantham import grantham
import blosum as bl
import logging
import time
import sys 
import joblib
import subprocess
import random
import string

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
def log_message(message):
    current_time = time.localtime()
    time_string = time.strftime("%Y-%m-%d %H:%M:%S", current_time)
    logging.info(f'{time_string} : {message}')
    logging.getLogger().handlers[0].flush()
def human_readable(input_str):
    number = int(input_str)
    if number >= 1_000_000:
        formatted_number = np.floor(number / 1_000_000 * 100) / 100
        return f"{formatted_number:.2f}M"
    elif number >= 1_000:
        formatted_number = np.floor(number / 1_000) 
        return f"{formatted_number:.0f}K"
    else:
        return str(number)

novel_path=str(sys.argv[1])
neut_path=str(sys.argv[2])
neut_fraction = float(sys.argv[3])
novel_fraction = float(sys.argv[4])
post_fix=str(sys.argv[5])
path_muts='/home/labs/alevy/omerbar/features/arabidopsis/mutation_features/'

# Creating a random thingy to save logs and the dataset
rand_char=random.choice(string.ascii_uppercase)
rand_num=random.randint(10, 99)
rand_phrase=f"{rand_char}{rand_num}"
path=f'/home/labs/alevy/omerbar/NN_general/NEW_AT/dataset_{rand_phrase}_standardized'
if post_fix!='':
    path=f'/home/labs/alevy/omerbar/NN_general/NEW_AT/dataset_{rand_phrase}_standardized_{post_fix}'
subprocess.run(['mkdir',path])
logging.basicConfig(filename=f"{path}/parse_muts.log", level=logging.INFO, format="%(message)s")
log_message('Starting...')

# Neutral dataset
os.chdir(os.path.join(path_muts,neut_path))

colnames=["#chr", "start", "end", "REF_ALLELE",
          "ALT_ALLELE_SIFT", "REF_AMINO","ALT_AMINO", "VARIANT_TYPE","AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","SIFT_PREDICTION","codon_pos", 
          "genomic_region_CDS", "genomic_region_exon","genomic_region_5UTR", "genomic_region_gene", "genomic_region_mRNA", "genomic_region_mRNA_TE_gene",
          "genomic_region_pseudogene ","genomic_region_3UTR", "genomic_region_tansposable_element_gene", "gene_count", "genomic_region_sRNA","genomic_region_intergenic",
          "kmers_13", "gc_35", "upstream_dist", "downstream_dist", "exon_count", "splice_junc", "repeat", "chromhmm",
          "TF_markers", "methylation_mark", "chrom_access_mark", "DNA_hypersensitivity_mark", "atacseq_mark",
          "H3K4me3", "H3K9me1", "H3K27me3", "H3K4me1", "H3K4me2", "H3K23ac", "H3K9ac", "H3K36ac", "H3K27ac", 
          "gerp_exp", "gerp_sub", "gerp_ta", "lrt_phylop_masked", "phast_est", "exp_median","score_phylop_masked","gerp_phylop_masked",
          "ALT_ALLELE", "label"]
df_neutral = pd.read_csv('neutral_features.bed', sep='\t', names=colnames,)#nrows=1000) # Removed extra feats that are descriminative. #skiprows = lambda x: x not in rand_read_ind)
log_message(f'Neutral set read, {len(df_neutral)} lines total, {len(df_neutral)*neut_fraction} used.')

# Subsetting randomly
neut_rand_subset = df_neutral.sample(frac=neut_fraction, replace=False)

# Novel dataset
os.chdir(os.path.join(path_muts,novel_path))
df_novel = pd.read_csv('novel_features.bed', sep='\t', names=colnames,)#nrows=1000)#,skiprows=1) # Removed extra feats that are descriminative.
log_message(f'Novel set read, {len(df_novel)} lines, {len(df_novel)*novel_fraction} used.')

# Subsetting randomly, keeping chr mutation rates
#novel_rand_subset = df_novel.groupby('#chr', group_keys=False).apply(lambda x: x.sample(frac=novel_fraction))
novel_rand_subset = df_novel.sample(frac=novel_fraction, replace=False)

log_message('Combining mutation datasets...')
# Combining datasets
df=pd.concat([novel_rand_subset,neut_rand_subset])
# Capitalize ref allele
df['REF_ALLELE']=df['REF_ALLELE'].str.upper()
# Filter out shitty ambiguous ref alleles
df = df[~df['REF_ALLELE'].isin(['K','M','W','S','D','Y','R'])]
# Filtering out hard masked and indels
df=df[(df.REF_ALLELE!='-') & (df.REF_ALLELE!='N')]
# Removing the ALT_ALLELE_SIFT because we already have the alt allele for non cds as well
df=df.drop(columns=['ALT_ALLELE_SIFT'])
# Changing missing values for sift score to 2 to reflect this as well
#df['SIFT_PREDICTION'].fillna(2, inplace=True)
# Changing chromhmm missing value to most common value - selected positions in the end of chromosomes
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
log_message(df.columns)
df['genomic_region_3UTR'] = (df['genomic_region_3UTR'] | df['genomic_region_5UTR']).astype(int)
# Lump Stop mutations
df['VARIANT_TYPE_STOP-LOSS'] = (df['VARIANT_TYPE_STOP-LOSS'] | df['VARIANT_TYPE_STOP-GAIN']).astype(int)
# Moving the label to the last column
df=df[[col for col in df.columns if col != 'label'] + ['label']]

# Getting the size of this dataset
run_size=human_readable(len(df))

# Saving not normalized dataset
df.to_csv(f'{path}/{run_size}_df_both_not_normalized.csv', sep='\t', header=True, index=False)

log_message('Normalizing...')
cols_to_norm = ["AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","codon_pos", "gene_count", "kmers_13", "gc_35", "upstream_dist", "downstream_dist", "exon_count",
                "TF_markers", "methylation_mark", "chrom_access_mark", "DNA_hypersensitivity_mark", "atacseq_mark",
                "H3K4me3", "H3K9me1", "H3K27me3", "H3K4me1", "H3K4me2", "H3K23ac", "H3K9ac", "H3K36ac", "H3K27ac",
                "gerp_exp", "gerp_sub", "gerp_ta", "lrt_phylop_masked", "phast_est", "exp_median","score_phylop_masked","gerp_phylop_masked",
               ]
log_message(f'Columns normalized:{cols_to_norm}')

# Normalizing
scaler = StandardScaler()
# Applying scalers to standarize and normalize data

df[cols_to_norm]=scaler.fit_transform(df[cols_to_norm]).round(5)

# Save scaler 
joblib.dump(scaler, f'{path}/{run_size}_standard_scaler.joblib')


# Split the datasets
train_df, test_df = train_test_split(df, test_size=0.1,shuffle=True)

# Adding back the 'Chr' before each chromosome
df['#chr'] = 'Chr' + df['#chr'].astype(str)
# Saving the positions
train_df_positions=train_df[['#chr','start','end']]
test_df_positions=test_df[['#chr','start','end']]

# Dropping the positions
train_df=train_df.drop(['#chr','start','end'],axis=1)
test_df=test_df.drop(['#chr','start','end'],axis=1)


log_message('Saving datasets...')
# Saving datasets
train_df.to_csv(f'{path}/{run_size}_train_df_normed.csv',sep='\t',header=True,index=False)
train_df_positions.to_csv(f'{path}/{run_size}_train_df_positions.csv',sep='\t',header=True,index=False)

test_df.to_csv(f'{path}/{run_size}_test_df_normed.csv',sep='\t',header=True,index=False)
test_df_positions.to_csv(f'{path}/{run_size}_test_df_positions.csv',sep='\t',header=True,index=False)