import json
import pandas as pd
import glob as glob

output_filename = 'standardized_molecules_2.csv'

done_mol_files = glob.glob('sanitized_molecules/*.json')
all_molecules = []

print('reading in all files')
for filename in done_mol_files:
	with open(filename, "r") as fid:
		all_molecules.append(fid.read())

print('done reading in files')
all_molecules = ','.join(all_molecules)
all_molecules = '[%s]'%all_molecules

print('done making json')
df = pd.read_json(all_molecules)

print('done making dataframe')

df.to_csv(output_filename,index=False, encoding = 'utf-8')
print('done')