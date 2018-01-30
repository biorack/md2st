# md2st
Metabolite Databases to Standardized Tautomers

To maintain the original compound found in a metabolite database and a standardized representation for grouping metabolites that you would like to consider the same (e.g. arpartate and aspartic acid) you need a clean code to bring in structures (hopefully from sdf of mol files) and a standardization recipe.

We aren't suggesting that this is the best way to standardize molecules, but its how we do it: 1) neutralize; 2) desalt; 3) canonical tautomer.  The canonical tautomer comes from MolVS.  The neutralize and desalt codes are from various rdkit help forums.

You will need to download the input metabolite databases.  Code here uses  MetaCyc mol files; ChEBI sdf complete; and Lipid Maps sdf complete.


1. original_to_json.py 
    - Import original compounds from a variety of sources and create a json file with database id, name, source, formula and original structure as a kekule form SMILES structural representation.  Each is assigned a unique identifier in case you must go back to an entry.

2. make_cpd_jobs.py
    - For the compounds in the original json file, create a taskfile for the standardization of each. Set NUM_AT_A_TIME to a bigger number initially (less jobs) and gradually make it smaller if you need to debug.

3. standardize_compounds.py
    - for each entry in the json files creates two additional fields. First is a canonical tautomer SMILES string and the second is an inchi key for the inchi of that SMILES.

4. create_compound_reference_table.py
    - a simple script that loads all the json files and makes a nice csv table for future use.

Step 3 (standardize_compounds.py) should be run through a scheduler in a parallel programming environment.  I use TaskFarmer at NERSC.  Each molecule takes a little while to find the canonical tautomer form.  An example SBATCH file is included for those inclined to also do it this way.


