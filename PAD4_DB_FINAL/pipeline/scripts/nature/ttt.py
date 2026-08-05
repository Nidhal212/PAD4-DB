import pandas as pd
df = pd.read_parquet('data/processed/pad4_compounds.parquet')
# Replace 'mechanism' with whatever column was auto-detected (or just print all columns)
print(df['assay_mechanism'].unique()) # or the specific column name