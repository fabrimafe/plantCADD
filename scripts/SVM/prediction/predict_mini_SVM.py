import pandas as pd
import numpy as np
import joblib
import sys

model_id=sys.argv[1]
test_df_name=sys.argv[2]
data_folder=sys.argv[3]
save_folder=sys.argv[4]


# Read test data
test_data = pd.read_csv(f'{data_folder}/{test_df_name}',sep='\t',header=0)

# Load the trained SVM model
svm_model = joblib.load(f'models/svm_model_{model_id}.joblib')

# Get the labels for later
test_data.drop(columns=['label'],inplace=True)

# Get predictions and append them
predictions = svm_model.predict_proba(test_data)

# Append positions and save
np.save(f'{save_folder}/predictions_model_{model_id}_{test_df_name}.npy',predictions[:,1])