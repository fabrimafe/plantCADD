import numpy as np
import pandas as pd
import sys

# Get data and read it
test_df_name=sys.argv[1]
test_pos_name=sys.argv[2]
df_folder=sys.argv[3]
pred_folder=sys.argv[4]
prefix=sys.argv[5]

test_df=pd.read_csv(f'{df_folder}/{test_df_name}',header=0,sep='\t')
test_pos=pd.read_csv(f'{df_folder}/{test_pos_name}',header=0,sep='\t')

# Add all predictons from mini SVMs
pred_df=pd.DataFrame()
for i in range(49):
    predictions = np.load(f'{pred_folder}/predictions_model_{i}_{test_df_name}.npy')
    pred_df[f'pred_{i}']=predictions

# get average per row from predictions df and add to test df
test_df['Prediction']=pred_df.mean(axis=1)

# Add positions
final_df=pd.concat([test_pos,test_df],axis=1)

# Write it
final_df.to_csv(f'${prefix}_final_predictions_{test_df_name}',sep='\t',index=False,header=True)