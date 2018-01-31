from __future__ import print_function

import sys
import os
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import PandasTools
import gzip

import pandas as pd
import numpy as np
import re
import glob


from io import open

import uuid




#Specify locations of flatfile databases
biocyc_paths = ['/global/homes/b/bpb/Downloads/metacyc_21.5/data/','../../data/tier1/Eco_21.5/data/','../../data/tier1/Sco_17.5/data/']
path_to_chebi = '/global/homes/b/bpb/Downloads/ChEBI_complete.sdf.gz'
path_to_hmdb = '/global/homes/b/bpb/Downloads/hmdb_structures.sdf'

path_to_chebi_names = '/global/homes/b/bpb/Downloads/chebi_names.tsv'
path_to_lipidmaps = '/global/homes/b/bpb/Downloads/LMSDFDownload12Dec17.tar.gz'
output_file = 'original_molecules.json'

def parse_biocyc_flat_file(my_files,attributes =  ['UNIQUE-ID','COMMON-NAME']):
    """
     We will only use the following

     ['UNIQUE-ID',
     'COMMON-NAME']

    """
    header_lines = []
    all_lines = []
    for my_file in my_files:
        with open(my_file,'r',encoding='utf-8', errors='ignore') as fid:
            for line in fid.readlines():
                if line.startswith('#'):
                    header_lines.append(line)
                else:
                    all_lines.append(line)
    header_lines = '\n'.join(header_lines)
    all_entities = '\n'.join(all_lines).split('//\n')

    # attributes = [a.strip() for a in header_lines.split('Attributes:\n')[-1].replace('#','').split('\n')]
    # attributes = [a for a in attributes if len(a)>0 ]
    
    empty_dict = {}
    for a in attributes:
        empty_dict[a] = ''

    all_dicts = []
    for a in all_entities:
        new_entity = copy.deepcopy(empty_dict)
        for a_attr in a.split('\n'):
            a_attr_key_value = a_attr.split(' - ')
            if a_attr_key_value[0] in new_entity:
                try:
                    new_entity[a_attr_key_value[0]].append(a_attr_key_value[1])
                except:
                    new_entity[a_attr_key_value[0]] = [a_attr_key_value[1]]
        all_dicts.append(new_entity)
    df = pd.DataFrame(all_dicts)
    return df


print('doing biocyc')
molecules = []
for path in biocyc_paths:
    mol_files = glob.glob(os.path.join(path,'MetaCyc-MOLfiles/*.mol'))
    for f in mol_files:
        cpd_id = os.path.basename(f).replace('.mol','')
        with open(f,'r',encoding='utf-8', errors='ignore') as fid:
            t = fid.read()
        name = t.split('\n')[0].strip().strip('"')       
        mol = Chem.MolFromMolFile(f,sanitize=True)
        if mol is not None:
            formula = CalcMolFormula(mol)
            try:
                Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
                smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
                original_smiles = smiles
            except:
                original_smiles = None
            molecules.append({'original_id':str(cpd_id),
                  'name':str(name),
                  'source':str('BioCyc'),
                  'formula':str(formula),
                  'original_smiles':str(original_smiles),
                  'unique_id':str(uuid.uuid4())})
        else:
            molecules.append({'original_id':str(cpd_id),
                  'name':str(name),
                  'source':str('BioCyc'),
                  'formula':None,
                  'original_smiles':None,
                  'unique_id':str(uuid.uuid4())})
            

print('doing chebi')
done_chebi_ids = []
with gzip.open(path_to_chebi) as fid:
    suppl = Chem.rdmolfiles.ForwardSDMolSupplier(fid,sanitize=True)
    for mol in suppl:
        if mol is not None:
            formula = CalcMolFormula(mol)
            try:
                Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
                smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
                original_smiles = smiles
            except:
                original_smiles = None
            d = mol.GetPropsAsDict()
            done_chebi_ids.append(d['ChEBI ID'])
            molecules.append({'original_id':str(d['ChEBI ID']),
                  'name':str(d['ChEBI Name']),
                  'source':str('CHEBI'),
                  'formula':str(formula),
                  'original_smiles':str(original_smiles),
                  'unique_id':str(uuid.uuid4())})

# add chebi compounds that only have a name
chebi_names = pd.read_csv(path_to_chebi_names,sep='\t')
chebi_names = chebi_names[['COMPOUND_ID','NAME']]
chebi_names.columns = ['compound_id','name']
chebi_names.drop_duplicates(subset='compound_id',inplace=True)
chebi_names['compound_id'] = chebi_names['compound_id'].apply(lambda x: 'CHEBI:%d'%x)
missing_chebi_names = chebi_names[~ chebi_names['compound_id'].isin(done_chebi_ids)]
print(len(missing_chebi_names))
for i,row in missing_chebi_names.iterrows():
	print(row['compound_id'],row['name'])
	molecules.append({'original_id':str(row['compound_id']),
	  'name':str(row['name']),
	  'source':str('CHEBI'),
	  'formula':None,
	  'original_smiles':None,
	  'unique_id':str(uuid.uuid4())})

            
"""
LOAD LIPID MAPS

Annotations include:
LM_ID, COMMON_NAME, CATEGORY, MAIN_CLASS, SUB_CLASS, CHEBI_ID, INCHI_KEY,LIPID_MAPS_CMPD_URL,
PUBCHEM_SUBSTANCE_URL, SYSTEMATIC_NAME, SYNONYMS,EXACT_MASS, FORMULA, LIPIDBANK_ID,
PUBCHEM_SID, KEGG_HI KEY, INCHI STRING, STATUS
"""
print('doing lipid maps')
with gzip.open(path_to_lipidmaps) as fid:
    suppl = Chem.rdmolfiles.ForwardSDMolSupplier(fid,sanitize=True)
    for mol in suppl:
        if mol is not None:
            formula = CalcMolFormula(mol)
            try:
                Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
                smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
                original_smiles = smiles
            except:
                original_smiles = None
            d = mol.GetPropsAsDict()
            if 'COMMON_NAME' in d:
                name_str = 'COMMON_NAME'
            else:
                name_str = 'SYSTEMATIC_NAME'
            molecules.append({'original_id':str(d['LM_ID']),
                  'name':str(d[name_str]),
                  'source':str('Lipid Maps'),
                  'formula':str(formula),
                  'original_smiles':str(original_smiles),
                  'unique_id':str(uuid.uuid4())})

# """
# LOAD HMDB
# """
print('doing hmdb')
hmdb_df = Chem.PandasTools.LoadSDF(path_to_hmdb, idName='ID', molColName='ROMol', includeFingerprints=False, isomericSmiles=False, smilesName=None, embedProps=False)
for i,row in hmdb_df.iterrows():
    mol = row['ROMol']
    formula = CalcMolFormula(mol)
    try:
        Chem.SanitizeMol(mol)
        Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
        original_smiles = smiles
    except:
        print('HMDB SANITIZATION ERROR',row['HMDB_ID'],row['GENERIC_NAME'])
        original_smiles = None
    molecules.append({'original_id':row['HMDB_ID'],
                      'name':row['GENERIC_NAME'],
                      'source':'HMDB',
                      'formula':formula,
                      'original_smiles':original_smiles,
                      'unique_id':str(uuid.uuid4())})


import json
with open(output_file, 'w',encoding="latin1") as fid:
    fid.write(unicode(json.dumps(molecules)))#,ensure_ascii=False)))