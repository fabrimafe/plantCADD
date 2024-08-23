import pandas as pd
from sklearn import svm
import joblib
import numpy as np
import logging
import time
import sys

def log_message(message):
    current_time = time.localtime()
    time_string = time.strftime("%Y-%m-%d %H:%M:%S", current_time)
    logging.info(f'{time_string} : {message}')
    logging.getLogger().handlers[0].flush()

# Get df names from user
train_df_name=sys.argv[1]
batch_size=int(sys.argv[2])
model_id=int(sys.argv[3])

# Init logger
curr_path=sys.path[0]
current_time = time.localtime()
time_string = time.strftime("%Y-%m-%d_%H:%M", current_time)
logging.basicConfig(filename=f"{curr_path}/model_{model_id}_train.log", level=logging.INFO, format="%(message)s")

# Read and parse train set
train_df=pd.read_csv(f'bootstrapped_training_sets/{train_df_name}',header=0,sep='\t',compression='gzip')
train_labels=train_df['label']
train_df.drop(columns=['label'],inplace=True)
log_message(f"Training data loaded from {train_df_name}")

# Truncate dataset if necessary
num_samples = len(train_labels)
num_epochs = num_samples // batch_size
if num_epochs * batch_size < num_samples:
    train_df = train_df[:num_epochs * batch_size]
    train_labels = train_labels[:num_epochs * batch_size]

# Initialize SVM model
model = svm.SVC(kernel='linear', probability=True)

# Train the model
for epoch in range(num_epochs):
    log_message(f"Epoch {epoch + 1}/{num_epochs}:")
    start_idx = epoch * batch_size
    end_idx = (epoch + 1) * batch_size
    batch_X = train_df[start_idx:end_idx]
    batch_y = train_labels[start_idx:end_idx]
    model.fit(batch_X, batch_y)

# Save the trained model
joblib.dump(model, f'svm_model_{model_id}.joblib')