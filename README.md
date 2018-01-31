# md2st
## Metabolite Databases to Standardized Tautomers

To maintain the original compound found in a metabolite database and a standardized representation for grouping metabolites that you would like to consider the same (e.g. arpartate and aspartic acid) you need a clean code to bring in structures (hopefully from sdf of mol files) and a standardization recipe.

We aren't suggesting that this is the best way to standardize molecules, but its how we do it: neutralize; desalt; and find canonical tautomer.  The canonical tautomer comes from MolVS.  The neutralization and desalting procedures are from various rdkit help forums.

A major limitation of this approach is that stereo-chemistry is lost when creating the standardized forms.  D and L isomers will be consolidated as identical standard forms.

# Procedure

1. original_to_json.py 
    - Import original compounds from a variety of sources and create a json file with database id, name, source, formula and original structure as a kekule form SMILES structural representation.  Each is assigned a unique identifier in case you must go back to an entry.

2. make_cpd_jobs.py
    - For the compounds in the original json file, create a taskfile for the standardization of each. Set NUM_AT_A_TIME to a bigger number initially (less jobs) and gradually make it smaller if you need to debug.

Start out like this:

```bash
python make_cpd_jobs.py --num 100
```

Then, after those finish, the rest can be submitted like this.

```bash
python make_cpd_jobs.py --num 1
```

There is one large conjugated molecule that can not finish in 30 minutes.  That single molecule, I ran on a login node and it finishes in about 45 minutes.

3. standardize_compounds.py
    - for each entry in the json files creates two additional fields. First is a canonical tautomer SMILES string and the second is an inchi key for the inchi of that SMILES.

4. create_compound_reference_table.py
    - a simple script that loads all the json files and makes a nice csv table for future use. It leaves you with original molecules in kekule-SMILES format, names in UTF-8 encoded strings, and (our definition of) standardized structures and inchi keys.

Step 3 (standardize_compounds.py) should be run through a scheduler in a parallel programming environment.  I use TaskFarmer at NERSC.  Each molecule takes a little while to find the canonical tautomer form.  An example SBATCH file is included for those inclined to also do it this way. Using 64 nodes, it takes about an hour to find the standard tautomer for about 200,000 compounds.


# Input Files

## BioCyc

You will need to download the input metabolite databases.  Code here uses  MetaCyc mol files; ChEBI sdf complete; and Lipid Maps sdf complete.

Using Metacyc 21.5 flat files:
Get them here: http://bioinformatics.ai.sri.com/ecocyc/dist/flatfiles-52983746/

Username: biocyc-flatfiles

Password: ##########

Search your email for actual password. Subject = "Pathway Tools and BioCyc distribution available"

## ChEBI and LipidMaps SDF files

findable by web searching.

## HMDB SDF file

This is coming soon.  There is a ton of information for each record and parsing it takes several hours.  Once this process is bug free, these molecules can be added.