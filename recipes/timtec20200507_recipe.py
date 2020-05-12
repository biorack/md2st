


""" contribution from Hans de Winter """
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
#         # Salts This was a test to deal with bonded salts in SDF files.  Hopefully unnecessary
#         ('[Na]','[H]'),
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

def count_carbons(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

def desalt(mol):
    #input is an rdkit mol
    #returns an rdkit mol keeping the biggest component
    #returns original mol if only one component
    #returns a boolean indicated if cleaning was necessary
    d = GetMolFrags(mol) #these are atom indices
    if len(d) == 1: #If there are fragments or multiple molecules this will be greater than 1
        return mol,False
    my_smiles=MolToSmiles(mol,True)
    parent_atom_count=0;
    disconnected=my_smiles.split('.')
    #With GetMolFrags, we've already established that there is more than one disconnected structure
    status = False
    for s in disconnected:
        little_mol = MolFromSmiles(s,sanitize=True)
        #Sanitize=True will fail for choline sulfate.  Can't sanitize the radical.
        if little_mol is not None:
            count = little_mol.GetNumAtoms()
            if (count > parent_atom_count):
#             if (count_carbons(little_mol)>0) & (count > parent_atom_count):
                parent_atom_count = count
                parent_mol = little_mol
                status = True
    return parent_mol,status

def force_to_unicode(text):
    "If text is unicode, it is returned as is. If it's str, convert it to Unicode using UTF-8 encoding"
    return text if isinstance(text, unicode) else text.decode('utf8')

def mol_to_neutral_desalted_canonical(mol_entry,filename):
    """
    given an rdkit mol and filename
    convert the mol into a standard form and
    save it as an mdl mol file.
    """
#     canon = tautomer.TautomerCanonicalizer()
    if mol_entry['original_smiles'] is not None:
        mol = Chem.MolFromSmiles(mol_entry['original_smiles'])
        # print('%s SMILES: %s'%(force_to_unicode(mol_entry['name']),mol_entry['original_smiles']))
        print('SMILES: %s'%(mol_entry['original_smiles']))
        if (not '*' in mol_entry['original_smiles']) and (mol is not None) and (mol.GetNumAtoms()>0):
            Chem.SanitizeMol(mol)
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            mol, status = neutralise_charges(mol)
            mol = desalt(mol)
#             try:
#                 mol = standardizer.standardize(mol)
#                 Chem.SanitizeMol(mol)
#                 Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
#             except:
#                 Chem.SanitizeMol(mol)
#                 Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
#             mol = canon.canonicalize(mol)
            Chem.SanitizeMol(mol)
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
            inchikey = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
            mol_entry['standardized_smiles'] = smiles
            mol_entry['standardized_inchikey'] = inchikey
    with open(filename, 'w') as fid:
        json.dump(mol_entry, fid)
        # Chem.MolToMolFile(mol,filename,includeStereo=True)





sdf_files = ['/Users/bpb/Downloads/NDL-3000_POL-TT021818LBL2sdf.sdf',
            '/Users/bpb/Downloads/FL-500_L-TT021818LBL2.sdf',
            '/Users/bpb/Downloads/NPL-800_POL-TT021818LBL2sdf.sdf 2']

maps = ['NDL3000','FL500','NPL800']

sdf = []

for i,my_file in enumerate(sdf_files):
    print(maps[i],my_file)
    if 'NPL-800' in my_file:
        with open(my_file,'rb') as fid:
            my_str = fid.read()
        with open('/Users/bpb/Downloads/fix_file.sdf','wb') as fid:
            fid.write(my_str.decode(encoding='Latin').encode('utf-8'))
        with open('/Users/bpb/Downloads/fix_file.sdf','rb') as fid:
            temp = PandasTools.LoadSDF(fid)
    else:
        temp = PandasTools.LoadSDF(my_file)
    temp.columns = [c.lower() for c in temp.columns]
    drop_cols = ['ACCEPTORS', 'AMOUNT', 'Absorption', 'Acceptors', 'Brutto-formula',
           'CAS NUMBER', 'CHIRAL_CENTERS', 'Chiral', 'DONORS', 'Donors', 'EMAIL', 'LIPINSKY',
           'Lipinsky', 'LogP', 'LogS', 'Molecular weight', 'N+O',
            'PO', 'ROTATION_BONDS','SALTDATA', 'SMILES', 'SUPPLIER', 'WEBSITE','inchi','inchikey','order','plate','well']
    drop_cols = [c.lower() for c in drop_cols]
    drop_cols = [c for c in drop_cols if c in temp.columns]

    temp.drop(columns=drop_cols,inplace=True)
    temp['original_smiles'] = temp['romol'].apply(MolToSmiles)

    if not 'name' in temp.columns:
        print('DOES NOT HAVE NAME')
        print(temp.columns)
        temp['name'] = np.nan
    idx = (pd.notna(temp['iupacname'])) & (pd.isna(temp['name']))
    temp.loc[idx,'name'] = temp.loc[idx,'iupacname']
    temp.loc[pd.isna(temp['name']),'name'] = ''


    temp['code'] = temp.apply(lambda x: '%s-%s'%(x['library'],x['idnumber']),axis=1)

    temp.drop(columns=['romol','iupacname','library','idnumber'],inplace=True)
    temp['library'] = maps[i]

    sdf.append(temp)
    print('\n\n\n')
sdf = pd.concat(sdf,axis=0,sort=True)

sdf.reset_index(drop=True,inplace=True)

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

sdf['smiles'] = ''
sdf['inchi_key'] = ''
sdf['inchi'] = ''

for i,row in sdf.iterrows():
    mol = MolFromSmiles(str(row['original_smiles']))
    #################################################
    # THIS IS THE BIG CHANGE FOR THESE STRUCTURES
    # SEVERAL ENTRIES HAD SALTS BONDED WRONG IN sdf
    # FILES. SEEMS LIKE TO INCHI KNOCKS THE SALTS
    # OFF.
    #################################################
    mol = MolToInchi(mol)
    mol = MolFromInchi(mol)
    SanitizeMol(mol)
    #################################################
    #################################################
    Kekulize(mol, clearAromaticFlags=True)
    mol, status = neutralise_charges(mol)
    mol, status = desalt(mol)
    SanitizeMol(mol)
    mol, status = neutralise_charges(mol)
    SanitizeMol(mol)

    Kekulize(mol, clearAromaticFlags=True)
    new_smiles = MolToSmiles(mol,isomericSmiles=True)
    new_inchikey = MolToInchiKey(mol)
    new_inchi = MolToInchi(mol)
    mw = ExactMolWt(mol)
    sdf.loc[i,'smiles'] = new_smiles
    sdf.loc[i,'inchi_key'] = new_inchikey
    sdf.loc[i,'inchi'] = new_inchi
    sdf.loc[i,'neutral_mass'] = mw
    sdf.loc[i,'formula'] = CalcMolFormula(mol)

setup_cols ={'inchi_key':'metatlas_inchikey',
             'inchi':'metatlas_inchi',
             'formula':'metatlas_formula',
             'neutral_mass':'metatlas_mw',
            'code':'original_id',
            'library':'source'}
sdf.rename(columns=setup_cols).to_csv('/Users/bpb/Downloads/Tim-Tec-Compounds.tab',sep='\t',index=None)


df = pd.merge(sdf[['code','inchi_key','inchi','name','neutral_mass','original_smiles']],c18_df[['code','library','label','mz','adduct']],left_on='code',right_on='code',how='inner')
df.shape
