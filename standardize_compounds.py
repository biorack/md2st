
from __future__ import print_function

import sys
import os
from rdkit import Chem

import glob

from molvs import tautomer
# from molvs import Standardizer

import json
import argparse

parser = argparse.ArgumentParser()

# To make the input integers
parser.add_argument('--mol_indices', nargs='+', type=int)

with open("/global/homes/b/bpb/repos/md2st/original_molecules.json", "r") as fid:
    all_molecules = json.load(fid)

save_path = '/global/homes/b/bpb/repos/md2st/sanitized_molecules'
if not os.path.exists(save_path):
    os.makedirs(save_path)

""" contribution from Hans de Winter """
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y)) for x,y in patts]

_reactions=None
def neutralise_charges(mol, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = Chem.AllChem.ReplaceSubstructs(mol, reactant, product, replaceAll=True)
            mol = rms[0]
    if replaced:
        Chem.SanitizeMol(mol)
        Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        return (mol, True)
    else:
        return (mol, False)

def desalt(mol):
    #input is an rdkit mol
    #returns an rdkit mol keeping the biggest component
    #returns original mol if only one component
    #returns a boolean indicated if cleaning was necessary
    d = Chem.rdmolops.GetMolFrags(mol) #these are atom indices
    if len(d) == 1: #If there are fragments or multiple molecules this will be greater than 1 
        return mol,False
    my_smiles=Chem.MolToSmiles(mol,True)
    parent_atom_count=0;
    disconnected=my_smiles.split('.')
    #With GetMolFrags, we've already established that there is more than one disconnected structure
    status = False
    for s in disconnected:
        little_mol = Chem.MolFromSmiles(s,sanitize=True)
        #Sanitize=True will fail for choline sulfate.  Can't sanitize the radical.
        if little_mol is not None:
            count = little_mol.GetNumAtoms()
            if count > parent_atom_count:
                parent_atom_count = count
                parent_mol = little_mol
                status = True
    return parent_mol,status

def mol_to_neutral_desalted_canonical(mol_entry,filename):
    """
    given an rdkit mol and filename
    convert the mol into a standard form and
    save it as an mdl mol file.
    """
    canon = tautomer.TautomerCanonicalizer()
    if mol_entry['original_smiles'] is not None:
        mol = Chem.MolFromSmiles(mol_entry['original_smiles'])
        print('%s SMILES: %s'%(mol_entry['name'],mol_entry['original_smiles']))
        if (mol.GetNumAtoms()>0) and (mol is not None) and (not '*' in mol_entry['original_smiles']):
            Chem.SanitizeMol(mol)
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            mol, status = neutralise_charges(mol)
            mol, status = desalt(mol)
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
    
