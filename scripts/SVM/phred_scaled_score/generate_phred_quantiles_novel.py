import pandas as pd
import numpy as np
import sys
import os

df_name=sys.argv[1]

# Read novel scores
novel_df=pd.read_csv(f'{df_name}',sep='\t',names=['Prediction'],header=None) # In bash run for i in *; do cat ${i}|tail -n +2 |awk '{print $NF}' >> PREDICTION_VALUES.txt;done
# Create the quantiles, [0,1], resolution: 10e-5
quantiles = np.linspace(0, 1, int(1/1e-5) + 1)
# Get the score values that fit per quantile generated
quantile_values = np.quantile(novel_df.Prediction, quantiles)
# Save this as a dataframe
pd.DataFrame({"Quantile_values":quantile_values}).to_csv('novel_quantile_values.txt',sep='\t',header=False,index=False)