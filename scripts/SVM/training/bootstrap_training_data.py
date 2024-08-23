import pandas as pd
import numpy as np
import sys
import subprocess

def bootstrap_sample_pandas(data):
    """
    Generate bootstrap sample with replacement from the given pandas DataFrame.

    Parameters:
        data (DataFrame): The original dataset.

    Returns:
        A bootstrap sample as a pandas DataFrame.
    """
    indices = np.random.choice(data.index, size=len(data), replace=True)
    bootstrap_sample = data.loc[indices]
    return bootstrap_sample

# Example usage:
# Assuming you have a pandas DataFrame called 'df' with your data
train_df_name = sys.argv[1]
bootstrap_index = sys.argv[2]
df = pd.read_csv(train_df_name, header=0, sep='\t')
num_bootstrap_samples = 50
subprocess.run(["mkdir","-p","bootstrapped_training_sets"])
bootstrap_sample = bootstrap_sample_pandas(df)
bootstrap_sample.to_csv(f"bootstrapped_training_sets/bootstrap_sample_{bootstrap_index}.csv.gz", sep='\t', index=False, header=True,compression='gzip')
