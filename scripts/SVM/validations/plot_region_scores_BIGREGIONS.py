import pandas as pd
import matplotlib.pyplot as plt
import sys
import re
import joblib
from pard.grantham import grantham
import blosum as bl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

# This script needs to involve all of the steps, from getting a region of positions from the feature's dataset:
# 1. Parse it (mut_parse, NN_general) - normalize, pick cols, add specific features added there such as protein substitutions matrices
# 2. Predict upon these positions using the SVM, should be fairly easy
#    This step will include predicting per all mini learners, per row. Each row will also be duplicated 4 times, to include a different mutation identity
# 3. Plot all these predictions

# This is defined to create the hierarchy of the genomic region feature, and to add the intron region
# mRNA, exon and gene are removed regions, as they are substituted by other regions.
def collapse_genomic_regions(df):
    # Remove mrna and exon, they are contained in cds
    df_final=df[[col for col in df.columns if col not in ['genomic_region_exon','genomic_region_mRNA']]]
    # Create the intron region = gene - cds - utr
    df_final['genomic_region_INTRON']=df_final['genomic_region_gene']-df_final['genomic_region_CDS']-df_final['genomic_region_3UTR']-df_final['genomic_region_5UTR']
    df_final['genomic_region_INTRON']= df_final['genomic_region_INTRON'].apply(lambda x: 0 if x < 0 else x)
    df_final=df_final[[col for col in df_final.columns if col not in ['genomic_region_gene']]]
    return df_final
    
# This is defined in order to parse the predictions easily in order to plot
def parse_df(df,dummy_col):
        # Subset genomic regions via regex
        r=re.compile(f"{dummy_col}*")
        feature_cols=list(filter(r.match, df.columns))
        # Subset for these in order to combine them
        subset_df = df[feature_cols]
        # Revert the df back to categorical
        feature_series=subset_df.idxmax(axis=1)
        # Remove prefix
        feature_series=feature_series.str.replace(f'{dummy_col}_', '')
        # Add back to df
        df[f'{dummy_col}'] = feature_series
        # Remove the dummy features
        df=df.drop(columns=df.filter(regex=f'{dummy_col}_[^_]*').columns)
        return df
    
# A feature added in step 1
def is_transition(row):
   return 1 if ((row['REF_ALLELE'] == 'A' and row['ALT_ALLELE'] == 'G') or
                (row['REF_ALLELE'] == 'G' and row['ALT_ALLELE'] == 'A') or
                (row['REF_ALLELE'] == 'C' and row['ALT_ALLELE'] == 'T') or
                (row['REF_ALLELE'] == 'T' and row['ALT_ALLELE'] == 'C')) else 0  
   
# A feature added in step 1
def grantham_score(row):
  # Check for both '.' and '*' in either column
  if '.' in (row['REF_AMINO'], row['ALT_AMINO']) or '*' in (row['REF_AMINO'], row['ALT_AMINO']):
    return 0
  else:
    return grantham(row['REF_AMINO'], row['ALT_AMINO'])

# A feature added in step 1
def blosum62_score(row):
  # Check for both '.' and '*' in either column
  if '.' in (row['REF_AMINO'], row['ALT_AMINO']) or '*' in (row['REF_AMINO'], row['ALT_AMINO']):
    return 0
  else:
    blosum62_mat=bl.BLOSUM(62, default=0)
    return (blosum62_mat[row['REF_AMINO']][row['ALT_AMINO']])

# Scale the scores after predicting, in order to plot them scaled
def scale_score_phred(df,quantiles, colname):
    # Get the quantile the score belongs to
    return np.round((-10*np.log10(1-np.searchsorted(quantiles, df[colname])/len(quantiles))),5)

# Args
df_name=sys.argv[1] # The features of these regions, from the feature file 
scaler_name=sys.argv[2] # The scaler for this batch, to keep true to it
sample_df_name=sys.argv[3]
num_learners=int(sys.argv[4])
is_region_big=bool(sys.argv[5])
novel_or_neut=int(sys.argv[6]) # 0 for neutral, 1 for novel
#model_name=sys.argv[4]  
#cols_to_drop=sys.argv[5:]
quantile_file='../../final_novel_preds/novel_quantile_values.txt'
quantiles=np.loadtxt(quantile_file)


# Taken from ~/NN_general/parse_mut_dataset_standard_norm.py
colnames=["#chr", "start", "end", "REF_ALLELE",
          "SIFT_MUTATION", "REF_AMINO","ALT_AMINO", "VARIANT_TYPE","AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","SIFT_PREDICTION","codon_pos", 
          "genomic_region_CDS", "genomic_region_exon","genomic_region_5UTR", "genomic_region_gene", "genomic_region_mRNA", "genomic_region_mRNA_TE_gene",
          "genomic_region_pseudogene ","genomic_region_3UTR", "genomic_region_tansposable_element_gene", "gene_count", "genomic_region_sRNA","genomic_region_intergenic",
          "kmers_13", "gc_35", "upstream_dist", "downstream_dist", "exon_count", "splice_junc", "repeat", "chromhmm",
          "TF_markers", "methylation_mark", "chrom_access_mark", "DNA_hypersensitivity_mark", "atacseq_mark",
          "H3K4me3", "H3K9me1", "H3K27me3", "H3K4me1", "H3K4me2", "H3K23ac", "H3K9ac", "H3K36ac", "H3K27ac", 
          "gerp_exp", "gerp_sub", "gerp_ta", "lrt_phylop_masked", "phast_est", "exp_median","score_phylop_masked","gerp_phylop_masked",]

# Taken from ~/NN_general/parse_mut_dataset_standard_norm.py
cols_to_norm = ["AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","codon_pos", "gene_count", "kmers_13", "gc_35", "upstream_dist", "downstream_dist", "exon_count",
                "TF_markers", "methylation_mark", "chrom_access_mark", "DNA_hypersensitivity_mark", "atacseq_mark",
                "H3K4me3", "H3K9me1", "H3K27me3", "H3K4me1", "H3K4me2", "H3K23ac", "H3K9ac", "H3K36ac", "H3K27ac",
                "gerp_exp", "gerp_sub", "gerp_ta", "lrt_phylop_masked", "phast_est", "exp_median","score_phylop_masked","gerp_phylop_masked",
               ]
df=pd.read_csv(df_name, sep='\t',header=None,names=colnames)
df['REF_ALLELE']=df['REF_ALLELE'].str.upper()
# Filter out shitty ambiguous ref alleles
df = df[~df['REF_ALLELE'].isin(['K','M','W','S','D','Y','R'])]
# Filtering out '-' or 'N' mutations
df=df[(df.REF_ALLELE!='-') & (df.REF_ALLELE!='N')]
print('Processing Data...')
# Process the data - add features added during parsing stage [stage 1]
# The transition mutation has to be defined after we settle for a mutation, which only happens later on
# grantham scores
df['Grantham_score'] = df.apply(grantham_score, axis=1)
# blosum 62 scores
df['Blosum62_score'] = df.apply(blosum62_score, axis=1)
if len(df.loc[df["chromhmm"] == "-1", "chromhmm"])>0:
    df.loc[df["chromhmm"] == "-1", "chromhmm"] = df.loc[df["chromhmm"] != "-1", "chromhmm"].mode()[0]
# Changing chromhmm to categorical to make it a dummy feature (it is discrete)
df["chromhmm"]=df["chromhmm"].astype('category')
# Normalize - first we normalize, then we remove the features not using for this model's training
scaler = joblib.load(scaler_name)
df[cols_to_norm]=scaler.transform(df[cols_to_norm]).round(5)
# Remove features 
df['#chr'] = df['#chr'].str[3:].astype(int)
df=pd.get_dummies(df,dtype=int)

# Because certain regions will not contain all the possible dummies for a certain feature, we need to add these as 'empty' dummies, all 0
# For that we need a sample df containing all of these, and we use the test df for this model's run
# Adding shitty empty columns for missing features
sample_df = pd.read_csv(f'{sample_df_name}', sep='\t', header=0,nrows=2) # Reading only for the header
new_cols=[col for col in sample_df.columns if col not in df.columns and col!='label' and 'ALT_ALLELE' not in col]
for col in new_cols:
    df[col]=0
# We set the label to 1 because the following mutation we will predict upon are generated randomly
# Ultimately we may hit some neutral mutation by accident, but, well.. I'm not going to make a check for this one in a million chance atm
df['label']=novel_or_neut

# One last thing to do would be to generate the 4 different mutations per base we predict upon
# This will depend on whether the column is a coding one or not, or more precisely, whether we have 'ALT_ALLELE' as '.' or not
# Checkpoint 1
df.to_csv(f'{df_name.split(".")[0]}_checkpoint_1.csv',sep='\t',index=False,header=True)
df_coding=df[df['SIFT_MUTATION_.']==0].copy()
df_noncoding=df[df['SIFT_MUTATION_.']==1].copy()
df_coding['ALT_ALLELE_A']=df_coding['SIFT_MUTATION_A'] # Here we should already have 4 duplicates per row, because we crossed with SIFT features.
df_coding['ALT_ALLELE_T']=df_coding['SIFT_MUTATION_T'] # Here we should already have 4 duplicates per row, because we crossed with SIFT features.
df_coding['ALT_ALLELE_C']=df_coding['SIFT_MUTATION_C'] # Here we should already have 4 duplicates per row, because we crossed with SIFT features.
df_coding['ALT_ALLELE_G']=df_coding['SIFT_MUTATION_G'] # Here we should already have 4 duplicates per row, because we crossed with SIFT features.
df_coding=parse_df(df_coding,'ALT_ALLELE')
df_noncoding_duplicated = pd.concat([df_noncoding]*4, ignore_index=True)
df_noncoding_duplicated['ALT_ALLELE'] = ['A', 'C', 'T', 'G'] * len(df_noncoding)
df=pd.concat([df_noncoding_duplicated,df_coding])
# transition and transversion
df=parse_df(df,'REF_ALLELE')
df=parse_df(df,'SIFT_MUTATION')
df.drop(columns=['SIFT_MUTATION'], inplace=True)
#df=parse_df(df,'ALT_ALLELE')

df['is_transition'] = df.apply(is_transition, axis=1)
df=pd.get_dummies(df,dtype=int)
df.sort_values(by=['start','end']).to_csv(f'{df_name.split(".")[0]}_checkpoint_2.csv',sep='\t',index=False,header=True)

# Reorder and rename to fit the training of the SVM
# Man this is becoming a pain in the ass and it seems to wrong to do this crap manually
#df.drop(columns=['ALT_ALLELE_A','ALT_ALLELE_T','ALT_ALLELE_C','ALT_ALLELE_G'], inplace=True)
df=df[["#chr","start","end","AMINO_POS","SIFT_SCORE","SIFT_MEDIAN","NUM_SEQS","codon_pos","genomic_region_CDS","genomic_region_exon","genomic_region_5UTR",
       "genomic_region_gene","genomic_region_mRNA","genomic_region_mRNA_TE_gene","genomic_region_pseudogene ","genomic_region_3UTR",
       "genomic_region_tansposable_element_gene","gene_count","genomic_region_sRNA","genomic_region_intergenic","kmers_13","gc_35",
       "upstream_dist","downstream_dist","exon_count","splice_junc","repeat","TF_markers","methylation_mark","chrom_access_mark",
       "DNA_hypersensitivity_mark","atacseq_mark","H3K4me3","H3K9me1","H3K27me3","H3K4me1","H3K4me2","H3K23ac","H3K9ac","H3K36ac",
       "H3K27ac","gerp_exp","gerp_sub","gerp_ta","lrt_phylop_masked","phast_est","exp_median","score_phylop_masked","gerp_phylop_masked",
       "is_transition","Grantham_score","Blosum62_score","REF_ALLELE_A","REF_ALLELE_C","REF_ALLELE_G","REF_ALLELE_T","REF_AMINO_*",
       "REF_AMINO_.","REF_AMINO_A","REF_AMINO_C","REF_AMINO_D","REF_AMINO_E","REF_AMINO_F","REF_AMINO_G","REF_AMINO_H","REF_AMINO_I",
       "REF_AMINO_K","REF_AMINO_L","REF_AMINO_M","REF_AMINO_N","REF_AMINO_P","REF_AMINO_Q","REF_AMINO_R","REF_AMINO_S","REF_AMINO_T",
       "REF_AMINO_V","REF_AMINO_W","REF_AMINO_Y","ALT_AMINO_*","ALT_AMINO_.","ALT_AMINO_A","ALT_AMINO_C","ALT_AMINO_D","ALT_AMINO_E",
       "ALT_AMINO_F","ALT_AMINO_G","ALT_AMINO_H","ALT_AMINO_I","ALT_AMINO_K","ALT_AMINO_L","ALT_AMINO_M","ALT_AMINO_N","ALT_AMINO_P",
       "ALT_AMINO_Q","ALT_AMINO_R","ALT_AMINO_S","ALT_AMINO_T","ALT_AMINO_V","ALT_AMINO_W","ALT_AMINO_Y","VARIANT_TYPE_.",
       "VARIANT_TYPE_NONSYNONYMOUS","VARIANT_TYPE_START-LOST","VARIANT_TYPE_STOP-GAIN","VARIANT_TYPE_STOP-LOSS","VARIANT_TYPE_SYNONYMOUS",
       "SIFT_PREDICTION_DELETERIOUS","SIFT_PREDICTION_TOLERATED","chromhmm_E1","chromhmm_E10","chromhmm_E11","chromhmm_E12","chromhmm_E13",
       "chromhmm_E14","chromhmm_E15","chromhmm_E16","chromhmm_E17","chromhmm_E18","chromhmm_E19","chromhmm_E2","chromhmm_E20","chromhmm_E21",
       "chromhmm_E22","chromhmm_E23","chromhmm_E24","chromhmm_E25","chromhmm_E26","chromhmm_E27","chromhmm_E28","chromhmm_E29","chromhmm_E3",
       "chromhmm_E30","chromhmm_E31","chromhmm_E32","chromhmm_E33","chromhmm_E34","chromhmm_E35","chromhmm_E36","chromhmm_E4","chromhmm_E5",
       "chromhmm_E6","chromhmm_E7","chromhmm_E8","chromhmm_E9","ALT_ALLELE_A","ALT_ALLELE_C","ALT_ALLELE_G","ALT_ALLELE_T","label",]]
# Checkpoint 2
df.sort_values(by=['start','end']).to_csv(f'{df_name.split(".")[0]}_checkpoint_3.csv',sep='\t',index=False,header=True)
#df=parse_df(df,'REF_ALLELE') # Reran this to get this nicely for prediction analysis later on

# Predict via the SVM model
# Because we made an ensemble SVM, we need to predict each prediction 50 times, 1 per mini learner, and then accumulate the results
# I wanted to use the other scripts already made, but it seems far easier to do this right here in the loop
# Since we would probably predict upon a few Kbps, it makes sense to just save them in memory instead of saving all of the mini predictions before accumulating them together
# Forgive me father for I have sinned
print('Predicting...')
# init final prediction vector
predictions=[]
# Go over the bases
for row_id in range(len(df)):
    # Get the row from the df
    row_series = df.iloc[row_id].copy()
    # Construct a DataFrame with the same columns as df
    row = pd.DataFrame([row_series], columns=row_series.index)
    # We want to remove the information features, in order to predict. They are only removed from the row to predict upon, not from the actual df
    # Get the info features
    info_features_names = ['#chr', 'start', 'end','label']
    info_features = row[info_features_names]
    # Remove the info features
    row_for_prediction = row.drop(columns=info_features_names,inplace=False)
    # Init the predictions for this row
    row_predictions=[]
    # Predict with all the 50 mini learners
    for mini_learner_id in range(num_learners):
        # Load the learner
        mini_learner = joblib.load(f'../../models/svm_model_{mini_learner_id}.joblib')
        # Reshape it to feed into svm
        row_np = row_for_prediction#.values.reshape(1, -1)
        # Predict and append to the prediction list
        mini_prediction = mini_learner.predict_proba(row_np)[:,1]
        row_predictions.append(mini_prediction)
    # Add this row's prediction into the prediction vector, aggregated
    predictions.append(np.mean(row_predictions))
    print(f'{row_id}/{len(df)}')
# After done with all predictions, append them to the final df
df['Prediction'] = predictions
df['Scaled_score']=scale_score_phred(df,quantiles,'Prediction')
df['#chr'] = 'Chr' + df['#chr'].astype(str)

# Parsing for plotting
df=collapse_genomic_regions(df)
df=parse_df(df,'REF_ALLELE')
df=parse_df(df,'ALT_ALLELE')
df=parse_df(df,'genomic_region')
df=parse_df(df,'VARIANT_TYPE') 
# Getting the genomic region feature to contain all subset mutation variant types
df['parsed_genomic_region']=df['genomic_region'].copy()
# This is to combine both variant type and genomic region, for the CDS positions. We have this issue of the first bp not seen by SIFT, and so it is of VARIANT_TYPE=='.'. 
# This fucks it up as we will always have 'CDS' in the legend as long as we're plotting CDS regions from their start.
df.loc[df['genomic_region'] == 'CDS', 'parsed_genomic_region'] = df.loc[df['genomic_region'] == 'CDS', 'VARIANT_TYPE']
df = df[df['parsed_genomic_region'] != 'CDS']
df = df[df['parsed_genomic_region'] != '.']
#df.loc[(df['genomic_region'] == 'CDS') & (df['VARIANT_TYPE'] != '.'), 'parsed_genomic_region'] = df.loc[(df['genomic_region'] == 'CDS') & (df['VARIANT_TYPE'] != '.'), 'VARIANT_TYPE']
# Definition of regulatory, we're not sure if it's correct
#df.loc[(df['parsed_genomic_region'] == 'intergenic') & ((df['funTFBS'] == 1) | (df['CEmotif'] == 1) | (df['motifTFBS'] == 1)), 'parsed_genomic_region'] = 'regulatory'

# Checkpoint 3
df.loc[np.isinf(df['Scaled_score']), 'Scaled_score'] = 50
df.sort_values(by=['start','end']).to_csv(f'{df_name.split(".")[0]}_checkpoint_3.csv',sep='\t',index=False,header=True)

#df=pd.read_csv('check2_final_regions_predicted.csv',sep='\t',header=0)
print('Plotting...')

# Plot
# Divide the df into different mutation alleles
df_A=df[df['ALT_ALLELE']=='A']
df_C=df[df['ALT_ALLELE']=='C']
df_T=df[df['ALT_ALLELE']=='T']
df_G=df[df['ALT_ALLELE']=='G']
# Set score to 0 for mutations that are the same as the reference allele
df_A.loc[df_A['REF_ALLELE']=='A','Scaled_score']=0
df_C.loc[df_C['REF_ALLELE']=='C','Scaled_score']=0
df_T.loc[df_T['REF_ALLELE']=='T','Scaled_score']=0
df_G.loc[df_G['REF_ALLELE']=='G','Scaled_score']=0

# Dynamically create a color map according to the number of unique genomic regions in this bases range
#unique_regions = df['parsed_genomic_region'].unique()
#num_regions = len(unique_regions)
#color_map = mpl.colormaps['tab10']
#colors_dict = dict(zip(unique_regions, [color_map(i) for i in range(num_regions)]))
colors_dict = {
    'intergenic'    : 'slategrey',
    'CDS'           : 'royalblue',
    '5UTR'          : 'bisque',
    '3UTR'          : 'fuschia',
    'SYNONYMOUS'    : 'darkgreen',
    'NONSYNONYMOUS' : 'tomato',
    'START-LOST'    : 'darkorange',
    'STOP-LOSS'     : 'darkturquoise',
    'STOP-GAIN'     : 'black',
    'INTRON'        : 'peru',
}

# Create the plot
fig = plt.figure(figsize=(12, 10),)# constrained_layout=True)
# Used to get appropriate room for shared y label (leftmost grid) and legend (rightmost grid)
gs = gridspec.GridSpec(5, 3, height_ratios=[20, 20, 20, 20, 1], width_ratios=[1, 20, 3])
ax = [fig.add_subplot(gs[i, 1]) for i in range(4)]
# Init legend and axis lists
legend_handles = []
legend_labels = set() 

# Iterate through the rows and plot each bar with the specified color
for i, df in zip(range(4), [df_A, df_C, df_T, df_G]):
    bar_positions = []  # Store the x-axis positions of the bars
    ref_alleles = []
    for index, row in df.iterrows():
        # Get the genomic region type to match the color
        variant_type = row['parsed_genomic_region']
        curr_color = colors_dict.get(variant_type, 'white')
        print(f'For region: {variant_type}, matched color:{curr_color}')
        # Create the appropriate bar for the row
        bar = ax[i].bar(row['start'], row['Scaled_score'], color=curr_color, label=variant_type)
        bar_positions.append(row['start'])
        ref_alleles.append(row['REF_ALLELE'])
        # Add to the legend
        if variant_type not in legend_labels:
            legend_handles.append(bar)
            legend_labels.add(variant_type)

    # Incase the region is not too big
    if ~is_region_big:
    # Set x-ticks and x-tick labels for each subplot
        ax[i].set_xticks(bar_positions)
        ax[i].set_xticklabels(ref_alleles, rotation=90, ha='center', fontsize=8)
    ax[i].set_ylim(0, 50)
# Set subplot titles
ax[0].set_title('A Mutation')
ax[1].set_title('C Mutation')
ax[2].set_title('T Mutation')
ax[3].set_title('G Mutation')

# Add subplot for shared x label
ax_xlabel = fig.add_subplot(gs[4, 1])
ax_xlabel.set_axis_off()

# Add the legend
plt.figlegend(handles=legend_handles, labels=legend_labels,)
# Add titles
#fig.suptitle(f'Scaled Prediction Scores for mutations\nchr {df_A.iloc[0,0]}, pos {df_A.iloc[0,1]}-{df_A.iloc[-1,1]}', fontsize=16)
fig.text(0.005, 0.5, 'Scaled Score', va='center', rotation='vertical',fontsize=20)  # Adjust position as needed
fig.text(0.47, 0.02, 'Segment', ha='center', fontsize=20)
#plt.subplots_adjust(bottom=1)
plt.tight_layout()
fig.savefig(f'{df_name.split(".")[0]}.svg', dpi=300, bbox_inches='tight')