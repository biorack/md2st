import json
import pandas as pd

output_filename = 'standardized_molecules.csv'

done_mol_files = glob.glob('sanitized_molecules/*.json')
all_molecules = []
for filename in done_mol_files:
	with open(filename, "r") as fid:
	    all_molecules.append(json.load(fid))

df = pd.DataFrame(all_molecules)

df.to_csv(output_filename,index=False)