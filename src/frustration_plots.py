from Bio.PDB import *
import numpy as np
import warnings
from Bio import BiopythonWarning
import sys
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain
import os
import re
import time



'''
TO DO : 


or isolate the funtions that make the structural alignment
or organise all the plot generation with options of clculate only a part of it of alls of it

'''
'''
This scripts take one PDB File containing multiples monomers (chains) and make a frustration study of the monomers using FrustratometeR 
the results files are added to results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_monomer_for_a_frame 
Usage:

python3 frustration_plots.py chemin/vers/fichier.pdb


'''

# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)
#1.
def extract_info_from_filename(pdb_file):
    """
    Extract the type (MtEnc/TmEnc) and number from the filename.
    Handles both formats:
    - MtEnc_0.pdb
    - MtEnc_monomer_0.pdb
    """
    filename = os.path.basename(pdb_file)

    # Pattern for both formats:
    # 1. MtEnc_0.pdb
    # 2. MtEnc_monomer_0.pdb
    match = re.search(r'(MtEnc|TmEnc)(?:_monomer)?_(\d+)', filename)

    if match:
        enc_type = match.group(1)
        enc_number = match.group(2)
        return enc_type, enc_number
    else:
        raise ValueError(
            f"Filename '{filename}' doesn't match expected patterns. "
            "Expected formats:\n"
            "1. '.../MtEnc_<number>.pdb'\n"
            "2. '.../MtEnc_monomer_<number>.pdb'\n"
            "3. Same with TmEnc\n"
            "where <number> is one or more digits."
        )
#2.
def load_monomers(pdb_file):
    """Load the n monomers, corresponding to the chains in a list"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)

    chains = []
    for model in structure :
        for chain in model:
            chains.append(chain)
    return chains
#2.2
def load_pdb(pdb_file):
    """Load the PDB file and return it"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)
    return structure


#3.
def save_each_monomer_as_pdb(aligned_monomers, output_dir, enc_type1 , enc_number1, enc_type2 = None ,enc_number2 = None ):
    """
    Save each monomer in a file (atom number for 1 file is limited to 99999 by PDBIO)
    Args:
        aligned_monomers: List of chain object (Monomers)
        output_dir
        enc_type: Type of encapsulin (MtEnc/TmEnc) have to be generalised
        enc_number: Number of the encapsulin file

    """
    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"aligned monomere a {len(aligned_monomers)} monomers")
    for i, monomer in enumerate(aligned_monomers, start=1):
        if enc_type2 == None :
            # Create the structure
            structure = Structure.Structure(f"MONOMER_{i}")
            model = Model.Model(0)
            structure.add(model)

            # Create a new chain
            new_chain = Chain.Chain("A")

            # add the residues of the monomer
            for residue in monomer:
                new_chain.add(residue)

            model.add(new_chain)

            # addapt the serial number of atoms
            atom_number = 1
            for atom in structure.get_atoms():
                atom.set_serial_number(atom_number)
                atom_number += 1
            output_file = os.path.join(output_dir, f"{enc_type1}{enc_number1}_monomer{i}.pdb")

            # create de pdb
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_file)


            print(f"Monomer {i} save in {output_file}")
        else :

            structure = Structure.Structure(f"MONOMER_{i}")
            structure.add(monomer)


            enc_type = {1:enc_type1, 2:enc_type2}
            enc_number = {1:enc_number1, 2:enc_number2}
            output_file = os.path.join(output_dir, f"{enc_type[i]}{enc_number[i]}_monomer{i}.pdb")


             # create de pdb
            io = PDBIO()
            io.set_structure(monomer)
            io.save(output_file)
            print(f"Monomer {i} save in {output_file}")

#4 frustratometeR
def calculate_frustration(PDB_directory, output_directory):
    """
    this function call the R package FrustratometeR for the calculation of frustration.
    :param PDB_directory: the directory of the pdb file used for the frustration calculation
    :param output_directory: the directory to save the frustration results
    """



#when only a file given
def main(pdb_file1):
    start_time = time.time()
    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filename
    enc_type, enc_number = extract_info_from_filename(pdb_file1)

    #1.1. create the output directory path
    frustration_dir = f"FRUSTRATION_{enc_type.upper()}"
    capsids_dir = f"{enc_type.upper()}_CAPSIDS"
    monomer_dir = "FRUSTRATION_monomer_for_a_frame"
    results_dir = os.path.join("../results", frustration_dir, capsids_dir,monomer_dir, f"{enc_type}{enc_number}_monomers")


    #1.2 create the repository if it not exist
    os.makedirs(results_dir, exist_ok=True)

    # 2. Charger les monomères
    monomers = load_monomers(pdb_file1)
    print(f"Loaded {len(monomers)} monomers")

     # 3. create PDB files
    save_each_monomer_as_pdb(monomers, results_dir, enc_type, enc_number)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

if __name__ == "__main__":
    if len(sys.argv) != 2 :
        print("Usage: python analyze_encapsulin.py path/to/encapsulin.pdb number_of_wanted_monomer")
        sys.exit(1)
    elif len(sys.argv) ==2:
        try:
            main(sys.argv[1])
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)
