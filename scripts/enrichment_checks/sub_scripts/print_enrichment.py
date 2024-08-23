import pandas as pd
import sys
df_name=sys.argv[1]
df=pd.read_csv(df_name,sep='\t',header=None)
print(f"Number of positions: {len(df)}")
print(f"Codon enrichment distribution:")
print(df.iloc[:,-1].value_counts()/len(df))
