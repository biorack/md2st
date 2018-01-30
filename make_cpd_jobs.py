
from __future__ import print_function

import sys
import os
import glob
import json

import argparse

parser = argparse.ArgumentParser()

# To make the input integers
parser.add_argument('--num', type=int)

with open("/global/homes/b/bpb/repos/md2st/original_molecules.json", "r") as fid:
    all_molecules = json.load(fid)

python_binary = '/global/common/software/m2650/python-cori/bin/python'
py_file = '/global/homes/b/bpb/repos/md2st/standardize_compounds.py'

args = parser.parse_args()
NUM_AT_A_TIME = args.num

# with open("all_molecules.pkl", "rb") as f:
#     all_molecules = pickle.load(f)

done_mol_files = glob.glob('sanitized_molecules/*.json')
print(len(done_mol_files))
# indices in list of mols in pickle file
all_indices = range(0,len(all_molecules))

#get the indices of mol files that have been completed
done_indices = [int(os.path.basename(f).split('_')[0]) for f in done_mol_files]

#indices to be processed
indices = list(set(all_indices) - set(done_indices))

print(len(indices))

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

lists_of_indices = list(chunks(indices,NUM_AT_A_TIME))

with open('taskfile.sh','w') as fid:
	for task_indices in lists_of_indices:
		fid.write('%s %s --mol_indices '%(python_binary,py_file))
		for idx in task_indices:
			fid.write('%d '%idx)
		fid.write('\n')
