import pandas as pd
target = 'SOZMHIJABUOUSN-ORMVGQBSSAN-N'  # fill full IK
cliffs = pd.read_parquet('data/processed/activity_cliffs.parquet')
involved = cliffs[
    (cliffs['inchi_key_a'].str.startswith('SOZMHIJABUOUSN')) |
    (cliffs['inchi_key_b'].str.startswith('SOZMHIJABUOUSN'))
]
print(f'Cliff pairs involving SOZMHIJABUOUSN: {len(involved)}')