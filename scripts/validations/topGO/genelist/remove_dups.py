import pandas as pd
import sys

df_name=sys.argv[1]
output_file_name=sys.argv[2]

# Read the df
df=pd.read_csv(df_name,sep='\t',header=0)

# Sort it according to chr/start and then transcript length
df_sorted = df.sort_values(by=[df.columns[0], df.columns[1], df.columns[-1]]).reset_index(drop=True)

# Remove duplicates (defined by chr/start/end, and choose the one based on the longest transcript length)
df_no_duplicates = df_sorted.groupby(list(df.columns[:3])).last().reset_index()

df_no_duplicates = df_no_duplicates[df.columns[:-6].tolist() + [df.columns[-2]]]

# Save it
df_no_duplicates.to_csv(output_file_name, index=False,header=True,sep='\t')