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

#3.
def structural_alignment(monomers):
    """Align all the monomers to the first one"""
    super_imposer = Superimposer()

     # Handle empty input
    if not monomers:
        return []

    # Determine input type (Chain or Structure)
    reference = monomers[0]
    is_chain = hasattr(reference, 'get_atoms') and not hasattr(reference, 'get_chains')
    is_structure = hasattr(reference, 'get_chains')
    if not (is_chain or is_structure):
        raise ValueError("Input must be either list of Chains or list of Structures")

    aligned = [reference]
    if is_chain :
        # parse the CA atoms of the reference monomer
        ref_atoms = []
        for residue in reference:
            for atom in residue:
                if atom.id == 'CA':
                    ref_atoms.append(atom)

        for monomer in monomers[1:]:
            # parse the CA atoms of the mobile monomers
            mobile_atoms = []
            for residue in monomer:
                for atom in residue:
                    if atom.id == 'CA':
                        mobile_atoms.append(atom)
            # Align to the reference
            super_imposer.set_atoms(ref_atoms, mobile_atoms)
            # apply the regulation to the mobile monomers
            super_imposer.apply(monomer)
            # add the regulated monomers to the list aligned
            aligned.append(monomer)
    else :
        print("on rentre dans le else")
        # parse the CA atoms of the reference monomer
        ref_atoms = []
        for model in reference:
            for monomer in model:
                for residue in monomer:
                    for atom in residue:
                        if atom.id == 'CA':
                            ref_atoms.append(atom)

        for struct in monomers[1:]:
            # parse the CA atoms of the mobile monomers
            mobile_atoms = []
            for model in struct:
                for monomer in model:
                    for residue in monomer:
                        for atom in residue:
                            if atom.id == 'CA':
                                mobile_atoms.append(atom)

            # Align to the reference
            super_imposer.set_atoms(ref_atoms, mobile_atoms)
            # apply the regulation to the mobile monomers
            super_imposer.apply(monomer)
            # add the regulated monomers to the list aligned
            aligned.append(monomer)

    return aligned

#4.
def save_each_monomer_as_pdb(aligned_monomers, output_dir, enc_type1 , file_dict1,monomer_number, enc_type2 = None ,file_dict2 = None ):
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
    for i, monomer in enumerate(aligned_monomers, start=1)and frame_number, path in sorted(file_dict1.items()) :
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
            output_file = os.path.join(output_dir, f"{enc_type1}{frame_number}_monomer{monomer_number}.pdb")

            # create de pdb
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_file)


            print(f"Monomer {i} save in {output_file}")
        else :

            structure = Structure.Structure(f"MONOMER_{i}")
            structure.add(monomer)


            enc_type = {1:enc_type1, 2:enc_type2}
            enc_number = {1:file_dict1, 2:file_dict2}
            #output_file = os.path.join(output_dir, f"{enc_type[i]}{enc_number[i]}_monomer{monomer_number}.pdb")


             # create de pdb
            #io = PDBIO()
            #io.set_structure(monomer)
            #io.save(output_file)
            print(f"Work in process")






def main(pdb_directory, monomer_number):

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
    results_dir = os.path.join("../results", frustration_dir, capsids_dir, f"{enc_type}_monomer_{monomer_number}_aligned_frames")

    plots_dir = os.path.join("../plots", enc_type)

    #1.2 create the repository if it not exist
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    #2
    monomers = load_monomers(files_dict,chain_dict, monomer_number)
    print(f"Loaded {len(monomers)} monomers")

     # 3. Alignement structural
    aligned = structural_alignment(monomers)
    print("Structural alignment completed")

    # 4. create aligned PDB files
    save_each_monomer_as_pdb(aligned, results_dir, enc_type, files_dict, monomer_number)


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
