import pandas as pd
import sys

df_name = sys.argv[1]
output_name = sys.argv[2]
num_of_muts=int(sys.argv[3])

df = pd.read_csv(df_name, sep='\t', header=0)

# Sort the DataFrame by 'Scaled_score' in descending order and select the top 1000 rows
top_1k_df = df.nlargest(num_of_muts, 'Scaled_score')

top_1k_df[['#chr', 'start', 'end']].sort_values(by=['#chr', 'start']).to_csv(output_name, sep='\t', index=False, header=True)
