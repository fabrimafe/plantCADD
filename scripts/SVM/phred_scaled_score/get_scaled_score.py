import pandas as pd
import numpy as np
import sys

def scale_score_phred(df,quantiles):
    # Get the quantile the score belongs to
    return np.round((-10*np.log10(1-np.searchsorted(quantiles, df.Prediction)/len(quantiles))),5)

df_name=sys.argv[1]
quantile_file=sys.argv[2] #Defined by novel predictions

# Load quantiles
quantiles=np.loadtxt(quantile_file)

# Load the dataframe of mutations
df=pd.read_csv(df_name,sep='\t',header=0)

df['Scaled_score']=scale_score_phred(df,quantiles)
# In case we get inf, that is, the score value is higher than the max quantile
if df['Scaled_score'].isna().any():
    # Equate it to max quantile
    df[df['Scaled_score']==np.inf]=df['Scaled_score'][~np.isinf(df['Scaled_score'])].max()

df.to_csv(f'{df_name}',sep='\t',header=True,index=False)