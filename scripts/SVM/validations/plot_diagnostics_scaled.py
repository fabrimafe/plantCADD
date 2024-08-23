import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import re
import sys
import os
from subprocess import call
from matplotlib.ticker import FixedLocator, FixedFormatter

call(['mkdir','-p','feat_plots'])
# Getting arguements
path=os.getcwd()
pred_file=sys.argv[1]
# A function to parse the dataset - un-dummy the feature for plots
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

def plot_prediction_quantiles_by_feature_barstacked(df,disc_col, sparse=False,prefix=''):
    title=f'Scaled prediction score quantiles colored by {disc_col}'
    figname=f'{path}/feat_plots/2bin_SCALED_bstack_by_{disc_col}.svg'
    if sparse:
       df=df[df['genomic_region']=='CDS']
       figname=f'{path}/feat_plots/bstack_by_{disc_col}_SPARSE.svg'
       title=f'NN prediction score quantiles colored by {disc_col}, CDS regions'
    plt.figure()
    # Create quantiles for the df
    df['Scaled_prediction_quantiles'] = pd.cut(df['Scaled_score'], bins = list(range(0, 49, 2)))
    #df['Scaled_prediction_quantiles'] = df['Scaled_prediction_quantiles'] + 1
    
    # Group the data by quantiles and Genomic_Region
    grouped = df.groupby(['Scaled_prediction_quantiles', disc_col]).size().reset_index(name='count')
    grouped['Proportion']=grouped['count']/grouped['count'].sum()
    grouped['Proportion'] = grouped['Proportion'] / grouped.groupby('Scaled_prediction_quantiles')['Proportion'].transform('sum')
    pivot = grouped.pivot(index='Scaled_prediction_quantiles', columns=disc_col, values='Proportion')
    
    color_dict=['C0','C1','C2','C3','C4','black','C5','C6','C7','C8','C9','lime']
    # Plot each category with custom colors
    ax = pivot.plot(kind='bar', stacked=True, color=color_dict)

        
    custom_tick_positions = [i for i in range(0,24)]
    custom_tick_labels = [i for i in range(2,49,2)]
    ax.xaxis.set_major_locator(FixedLocator(custom_tick_positions))
    ax.xaxis.set_major_formatter(FixedFormatter(custom_tick_labels))
    ax.set_xlabel('Phred Scale Score')
    ax.set_ylabel('Proportion')
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    # Save the plot
    if prefix != '':
        figname=f'{figname[:-4]}_{prefix}.svg'
    plt.savefig(figname, dpi=240, bbox_inches='tight')
    plt.close()
    plt.clf()

def collapse_genomic_regions(df):
    # Remove mrna and exon, they are contained in cds
    df_final=df[[col for col in df.columns if col not in ['genomic_region_exon','genomic_region_mRNA']]]
    # Create the intron region = gene - cds - utr
    df_final['genomic_region_INTRON']=df_final['genomic_region_gene']-df_final['genomic_region_CDS']-df_final['genomic_region_3UTR']-df_final['genomic_region_5UTR']
    df_final['genomic_region_INTRON']= df_final['genomic_region_INTRON'].apply(lambda x: 0 if x < 0 else x)
    df_final=df_final[[col for col in df_final.columns if col not in ['genomic_region_gene']]]
    return df_final

# Save the plot
# Read the dataset of predictions
df=pd.read_csv(os.path.join(path,pred_file), sep='\t', header=0)
df = df.rename(columns={'genomic_region_intergenic': 'genomic_region_INTERGENIC',
                        'genomic_region_pseudogene': 'genomic_region_PSEUDOGENE',
                        'genomic_region_mRNA_TE_gene': 'genomic_region_mRNA-TE-GENE'})

print('Reading done.')
# Unlump UTRs
df['genomic_region_3UTR'] = df['genomic_region_3UTR'] - df['genomic_region_5UTR']
# Lump Stop mutations
df['VARIANT_TYPE_STOP-LOSS'] = df['VARIANT_TYPE_STOP-LOSS'] - df['VARIANT_TYPE_STOP-GAIN']
# Converting dummies to leveled variables
df=collapse_genomic_regions(df)
df=parse_df(df,'genomic_region') 
df=parse_df(df,'VARIANT_TYPE') 
print('Parsing done.')
df['parsed_genomic_region']=df['genomic_region']
df.loc[df['genomic_region'] == 'CDS', 'parsed_genomic_region'] = df.loc[df['genomic_region'] == 'CDS', 'VARIANT_TYPE']
# Dividing the dataset to novel and neutral
df_novel=df[df.label==1]
#df_neutral=df[df.label==0]
#plot_prediction_quantiles_by_feature_barstacked(df,'genomic_region',False)
plot_prediction_quantiles_by_feature_barstacked(df_novel,'parsed_genomic_region',False,'NOVEL')
#plot_prediction_quantiles_by_feature_barstacked(df_neutral,'genomic_region',False,'NEUTRAL')
plot_prediction_quantiles_by_feature_barstacked(df_novel,'label',False)


print('Plotted genomic_regions.')

df_CDS=df[df.genomic_region=='CDS']#.drop(columns=['Scaled_prediction_quantiles'])
df_CDS_novel=df_CDS[df_CDS.label==1]
#df_CDS_neutral=df_CDS[df_CDS.label==0]
plot_prediction_quantiles_by_feature_barstacked(df_CDS_novel,'VARIANT_TYPE',False)
plot_prediction_quantiles_by_feature_barstacked(df_CDS_novel,'VARIANT_TYPE',False,'NOVEL')
#plot_prediction_quantiles_by_feature_barstacked(df_CDS_neutral,'VARIANT_TYPE',False,'NEUTRAL')

print('Plotted variant_types.')


