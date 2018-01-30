
from __future__ import print_function

import sys
import os
from rdkit import Chem


import glob

# sys.path.append('/global/homes/b/bpb/repos/metatlas/metatlas/compounds')
from metatlas.compounds import structure_cleaning as sc
from molvs import tautomer
# from molvs import Standardizer

import json
import argparse

parser = argparse.ArgumentParser()

# To make the input integers
parser.add_argument('--mol_indices', nargs='+', type=int)


# with open("/global/homes/b/bpb/repos/magi/workflow/parse_flatfiles/all_molecules.pkl", "rb") as f:
    # all_molecules = pickle.load(f)

with open("/global/homes/b/bpb/repos/magi/workflow/parse_flatfiles/original_molecules.json", "r") as fid:
    all_molecules = json.load(fid)

save_path = '/global/homes/b/bpb/repos/magi/workflow/parse_flatfiles/sanitized_molecules'

def mol_to_neutral_desalted_canonical(mol_entry,filename):
    """
    given an rdkit mol and filename
    convert the mol into a standard form and
    save it as an mdl mol file.
    """
    canon = tautomer.TautomerCanonicalizer()
    mol = Chem.MolFromSmiles(mol_entry['original_smiles'])
    print('%s SMILES: %s'%(mol_entry['name'],mol_entry['original_smiles']))
    if (mol.GetNumAtoms()>0) and (mol is not None) and (not '*' in mol_entry['original_smiles']):
        Chem.SanitizeMol(mol)
        Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        mol, status = sc.NeutraliseCharges(mol)
        mol, status = sc.desalt(mol)
#         try:
#             mol = standardizer.standardize(mol)
#             Chem.SanitizeMol(mol)
#             Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
#         except:
        Chem.SanitizeMol(mol)
        Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        mol = canon.canonicalize(mol)
        Chem.SanitizeMol(mol)
        Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
        inchikey = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
        mol_entry['standardized_smiles'] = smiles
        mol_entry['standardized_inchikey'] = inchikey
    with open(filename, 'w') as fid:
        json.dump(mol_entry, fid)
        # Chem.MolToMolFile(mol,filename,includeStereo=True)



args = parser.parse_args()
for idx in args.mol_indices:
    output_filename = os.path.join(save_path,'%d_%s.json'%(idx,all_molecules[idx]['unique_id']))
    if not os.path.isfile(output_filename):
        mol_to_neutral_desalted_canonical(all_molecules[idx],output_filename)
    
