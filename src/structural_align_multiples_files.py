from Bio.PDB import *
import numpy as np
import warnings
from Bio import BiopythonWarning
import sys
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain
import os
import re


'''
This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a wanted monomer, and make a structural alignment of this monomer and mesure the  mean RMSF for this monomers in each file (frames) 
A plot is add to the plot directory, with name of type rmsf_with_std_per_frame_<TmEnc|MmEnc>_monomer_<(i)>.png
The aligned PDB files are added to the result directory. 
Usage:

python3 structural_align_multiples_files.py path/to/directory/ number_of_wanted_monomer

'''

# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)

#1. Extract the type and the names of all the pdb files of a given directory
def extract_info_from_directory(pdb_directory):
    """
    Take a directory path and extract encapsulin type and return files in a dictionary with the frame number as key.

    Args:
        pdb_directory (str): Path to directory containing PDB files

    Returns:
        tuple: (encapsulin_type, dictionary_of_files)
            dictionary_of_files: {enc_number: file_path}

    Raises:
        ValueError: If directory doesn't exist or contains no valid files
    """
    if not os.path.isdir(pdb_directory):
        raise ValueError(f"Directory not found: {pdb_directory}")

    pattern = re.compile(r'(MtEnc|TmEnc)(?:_monomer)?_(\d+)\.pdb$', re.IGNORECASE)
    files_dict = {}
    enc_type = None

    for filename in sorted(os.listdir(pdb_directory)):
        match = pattern.match(filename)
        if match:
            current_type = match.group(1)
            frame_number = match.group(2)

            # Verify consistent encapsulin type
            if enc_type is None:
                enc_type = current_type
            elif current_type != enc_type:
                raise ValueError(f"Inconsistent encapsulin types: found both {enc_type} and {current_type}")

            # Add to dictionary (convert number to int for proper sorting)
            files_dict[int(frame_number)] = os.path.join(pdb_directory, filename)

    if not enc_type or not files_dict:
        raise ValueError(f"No valid encapsulin files found in {pdb_directory}")

    return enc_type, files_dict

#2. Load the wanted monomer of each frame in a list
def load_monomers(files_dict, monomer_number):
    """Take a dico of file path and a wanted monomer and load the monomer of each frame in a list"""
    parser = PDBParser(QUIET=True)
    chains = []
    for frame_number, path in sorted(files_dict.items()):
        structure = parser.get_structure("encapsulin", path)

        for model in structure :
            chain = list(model)[int(monomer_number)]
            chains.append(chain)
    return chains

def main(pdb_directory, monomer_number):

    #1.
    enc_type, files_dict = extract_info_from_directory(pdb_directory)
    print(f"Encapsulin type : {enc_type}")
    print(f"{len(files_dict)} frames :")
    for frame_number, path in sorted(files_dict.items()):
        print(f"frame {frame_number}: {os.path.basename(path)}")
    #2
    print(load_monomers(files_dict, monomer_number))

if __name__ == "__main__":
    if len(sys.argv) != 3 :
        print("Usage: python analyze_encapsulin.py path/to/encapsulin.pdb number_of_wanted_monomer")
        sys.exit(1)
    elif len(sys.argv) ==3:
        try:
            main(sys.argv[1], sys.argv[2])
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)
