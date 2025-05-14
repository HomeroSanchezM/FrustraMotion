from Bio.PDB import *
import numpy as np
import warnings
from Bio import BiopythonWarning
import sys
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain
import os
import re
import gzip
import shutil
import time



'''
This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a wanted monomer, and make a frustration study of the monomers in each file (frames) using FrustratometeR 
the results files are added to results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_frames_for_a_monomer 
A plot is add to the plots/frustration directory, with name of type rmsf_with_std_per_frame_<TmEnc|MmEnc>_monomer_<(i)>.png
Usage:

python3 frustration_plots_multiples_files.py path/to/directory/ number_of_the_wanted_monomer

'''
# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)

#1. Extract the type, the names of all the pdb files, and the id of all the chains of a given directory
def extract_info_from_directory(pdb_directory):
    """
    Take a directory path and extract encapsulin type and return the files path in a dictionary with the frame number as key.

    Args:
        pdb_directory (str): Path to directory containing PDB files

    Returns:
        tuple: (encapsulin_type, dictionary_of_files)
            dictionary_of_files: {enc_number: file_path}
            dico_of_chain_id: {chain_number: chain_id}

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

    parser = PDBParser(QUIET=True)
    chain_dict = {}
    structure = parser.get_structure("encapsulin", files_dict[0])
    for model in structure :
        chains = list(model)
        for i, chain in enumerate(chains, start=1):
            chain_dict[int(i)] = chain



    return enc_type, files_dict, chain_dict

def extract_chain_to_temp(pdb_path, chain_id, temp_dir="../tmp"):
    """Extract a single chain to a temporary file"""
    os.makedirs(temp_dir, exist_ok=True)
    output_path = os.path.join(temp_dir, f"chain_{chain_id}_{os.path.basename(pdb_path)}")

    with open(pdb_path) as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM')):
                # Split using whitespace (multiple spaces/tabs)
                fields = line.split()
                # Le champ chain_id est normalement le 5ème après split()
                # Mais cela dépend du format exact - peut nécessiter ajustement
                if len(fields) >= 5 and fields[4] == chain_id:
                    f_out.write(line)

    return output_path

#2. Load the wanted monomer of each frame in a list
def load_monomers(files_dict,chain_dict, monomer_number):
    """Take a dico of file path and a wanted monomer and load the monomer of each frame in a list"""
    chains = []
    temp_files = []
    chain_id = chain_dict[int(monomer_number)].get_id()
    for frame_number, path in sorted(files_dict.items()):
            print(f"Extracting chain {chain_id} from file {os.path.basename(path)} ")
            temp_path = extract_chain_to_temp(path, chain_id)
            temp_files.append(temp_path)
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("temp", temp_path)
            chains.append(structure[0][chain_id])
    # Clean up temporary files
    for temp_file in temp_files:
        os.remove(temp_file)

    return chains


#4.
def save_each_monomer_as_pdb(aligned_monomers, output_dir, enc_type, file_dict, monomer_number):
    """
    Save each monomer in a file with naming convention: "<enc_type><frame_number>_monomer<monomer_number>.pdb"

    Args:
        aligned_monomers: List of chain objects (Monomers)
        output_dir: Directory to save the PDB files
        enc_type: Type of encapsulin (MtEnc/TmEnc)
        file_dict: Dictionary mapping frame numbers to file paths
        monomer_number: Number of the wanted monomer
    """
    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"Saving {len(aligned_monomers)} aligned monomers")

    # Check we have same number of monomers and frames
    if len(aligned_monomers) != len(file_dict):
        raise ValueError(f"Mismatch between number of monomers ({len(aligned_monomers)}) and frames ({len(file_dict)})")

    # Iterate through aligned monomers and frame numbers together
    for (frame_number, _), monomer in zip(sorted(file_dict.items()), aligned_monomers):
        # Create the structure
        structure = Structure.Structure(f"MONOMER_{frame_number}")
        model = Model.Model(0)
        structure.add(model)

        # Create a new chain
        new_chain = Chain.Chain("A")

        # Add the residues of the monomer
        for residue in monomer:
            new_chain.add(residue)

        model.add(new_chain)

        # Adjust the serial number of atoms
        atom_number = 1
        for atom in structure.get_atoms():
            atom.set_serial_number(atom_number)
            atom_number += 1

        # Create output filename
        output_file = os.path.join(output_dir, f"{enc_type}{frame_number}_monomer{monomer_number}.pdb")

        # Save the PDB file
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)

        print(f"Saved frame {frame_number} as {output_file}")



def main(pdb_directory, monomer_number):
    start_time = time.time()

    #1.
    enc_type, files_dict, chain_dict = extract_info_from_directory(pdb_directory)
    print(f"Encapsulin type : {enc_type}")
    print(f"Number of Monomers by structure : {len(chain_dict)}")
    print(f"Number of frames : {len(files_dict)} ")
    #print(f"{len(files_dict)} frames :")
    #for frame_number, path in sorted(files_dict.items()):
    #    print(f"frame {frame_number}: {os.path.basename(path)}")

    #1.1. create the output directory path
    frustration_dir = f"FRUSTRATION_{enc_type.upper()}"
    capsids_dir = f"{enc_type.upper()}_CAPSIDS"
    frames_dir = "FRUSTRATION_frames_for_a_monomer"
    results_pdb_dir = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir, capsids_dir,frames_dir, f"{enc_type}_monomer_{monomer_number}_monomers")
    results_frustration_dir = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir, capsids_dir,frames_dir, f"{enc_type}_monomer_{monomer_number}_frustration")

    plots_dir = os.path.join("../plots", enc_type, "frustration")

    #1.2 create the repository if it not exist
    os.makedirs(results_pdb_dir, exist_ok=True)
    os.makedirs(results_frustration_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    #2
    monomers = load_monomers(files_dict,chain_dict, monomer_number)
    print(f"Loaded {len(monomers)} monomers")


    # 4. create aligned PDB files
    save_each_monomer_as_pdb(monomers, results_pdb_dir, enc_type, files_dict, monomer_number)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

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
