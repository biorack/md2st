# md2st
metabolite databases to standardized tautomers

1. original_to_json.py 
    - Import original compounds from a variety of sources and create a json file with database id, name, source, formula and original structure as a kekule form SMILES structural representation.  Each is assigned a unique identifier in case you must go back to an entry.

2. make_cpd_jobs.py
    - For the compounds in the original json file, create a taskfile for the standardization of each. Set NUM_AT_A_TIME to a bigger number initially (less jobs) and gradually make it smaller if you need to debug.

3. standardize_compounds.py
    - for each entry in the json files creates two additional fields. First is a canonical tautomer SMILES string and the second is an inchi key for the inchi of that SMILES.

4. create_compound_reference_table.py
    - a simple script that loads all the json files and makes a nice csv table for future use.

Step 3 (standardize_compounds.py) should be run through a scheduler in a parallel programming environment.  I use TaskFarmer at NERSC.  Each molecule takes a little while to find the canonical tautomer form.  An example SBATCH file is included for those inclined to also do it this way.


