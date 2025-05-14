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
import subprocess



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



#4 frustratometeR
def calculate_frustration(PDB_directory, output_directory):

    """
    this function call the R package FrustratometeR for the calculation of frustration.
    :param PDB_directory: the directory of the pdb file used for the frustration calculation
    :param output_directory: the directory to save the frustration results
    """

    # check if directories exists
    if not os.path.isdir(PDB_directory):
        raise ValueError(f"PDB directory not found: {PDB_directory}")
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # the temp R sript
    r_script = """
    # FrustratometeR calculation
    library(frustratometeR)

    # Lire le fichier PDB
    pdb_file <- "{pdb_file}"
    output_dir <- "{output_dir}"

    # Calculer la frustration
    results <- calculate_frustration(PdbFile = pdb_file, Mode = "singleresidue", ResultsDir = output_dir, Graphics = FALSE)

    """

    # Pour chaque fichier PDB dans le répertoire
    # TO DO make the frustration calculation in order for the files
    for pdb_file in os.listdir(PDB_directory):
        if pdb_file.endswith('.pdb'):
            full_path = os.path.join(PDB_directory, pdb_file)

            # Personnaliser le script R pour ce fichier
            current_script = r_script.format(
                pdb_file=full_path,
                output_dir=output_directory
            )
            #print(current_script)
            # Écrire le script temporaire
            with open('temp_frustration.R', 'w') as f:
                f.write(current_script)


            # Exécuter Rscript
            try:
                subprocess.run(['Rscript', 'temp_frustration.R'], check=True)
                print(f"Successfully processed {pdb_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {pdb_file}: {e}")
            finally:
                # Nettoyer le fichier temporaire
                if os.path.exists('temp_frustration.R'):
                    os.remove('temp_frustration.R')


#5
def dico_frustIndex(frustration_file):
    """
    Create a dictionary  with frustration indices from FrustratometeR output file.

    Args:
        frustration_file (str): Path to the frustration results file

    Returns:
        dict: {residue_key: FrstIndex} where residue_key is "AAresNum" (e.g. "M4")

    Example output:
        {'M4': 0.759, 'E5': -0.582, 'F6': 1.027, ...}
    """
    frustration_dict = {}

    with open(frustration_file, 'r') as f:
        # Skip header line if present
        header = f.readline()

        for line in f:
            # Skip empty lines
            if not line.strip():
                continue

            parts = line.split()

            # Check we have enough columns (at least 8)
            if len(parts) >= 8:
                res_num = parts[0]  # Residue number (4, 5, 6...)
                aa = parts[3]  # Amino acid (M, E, F...)
                frst_index = float(parts[7])  # FrstIndex value

                # Create key like "M4" and add to dictionary
                key = f"{aa}{res_num}"
                frustration_dict[key] = frst_index

    return frustration_dict


import os
import re


def dico_of_dico_frustIndex(frustration_directory):
    """
    This function takes the full path of the directory where the frustration results files are
    and for each result, builds a dictionary of the residues frustrationIndex, adding it to a main dictionary.
    Extracts the frame numberthat is before "monomer" in filenames (e.g., 0 from "TmEnc0_monomer1.pdb" or "MtEnc0_monomer1.pdb")

    :param frustration_directory: Path to directory containing frustration results
    :return: Dictionary {frame_number: {'M4': 0.759, ...}, ...}
    """
    dico_frames = {}

    for directory in os.listdir(frustration_directory):
        doc_path = os.path.join(frustration_directory, directory, "FrustrationData")

        for file in os.listdir(doc_path):
            if file.endswith("pdb_singleresidue"):
                full_path = os.path.join(doc_path, file)
                name_file = os.path.basename(full_path)

                # regex work with both TmEnc and MtEnc prefixes
                match = re.search(r'(?:Tm|Mt)Enc(\d+)_monomer', name_file)
                if match:
                    frame_number = int(match.group(1))
                    dico_monomer = dico_frustIndex(full_path)
                    dico_frames[frame_number] = dico_monomer
                else:
                    print(f"Warning: Could not parse encounter number from filename: {name_file}")

    return dico_frames


#8. plot mean frustration + std per residue
#8.1. dico of mean frustration for each res
def dico_mean_frustration (dico_monomers) :
    """

    :param dico_monomers: dico with the frustration info
    :return:
    """
     # Initialize dictionary to store sum and count
    dico_mean = {}

    # First pass: collect sum and count for each residue
    for monomer_data in dico_monomers.values():
        for residue, frustration in monomer_data.items():
            if residue not in dico_mean:
                dico_mean[residue] = {'mean': 0.0, 'std': []}
            dico_mean[residue]['mean'] += frustration
            dico_mean[residue]['std'].append(frustration)

    #print(dico_mean)
    #print(residue_stats)
    #print(len(dico_monomers))
    # Second pass: calculate mean for each residue
    dico_mean_and_std = {}
    for residue, stats in dico_mean.items():
        if residue not in dico_mean_and_std:
                dico_mean_and_std[residue] = {'mean': 0.0, 'std': 0.0}
        dico_mean_and_std[residue]['mean'] = stats['mean']/len(dico_monomers)
        dico_mean_and_std[residue]['std'] = np.std(stats['std'])
    return dico_mean_and_std


def plot_frustration_per_res(dico, enc_type, monomer_number, plots_dir):
    """
    Plot mean frustration with standard deviation per residue, connecting dots with lines.
    Saves the plot as a PNG file in the specified directory.

    Parameters:
    - dico (dict): Dictionary with residue IDs as keys and {'mean': float, 'std': float} as values.
                   Example: {'M4': {'mean': -0.135, 'std': 0.668}, ...}
    - enc_type (str): Encounter type (e.g., "H1", "H2").
    - enc_number (int or str): Encounter number or identifier.
    - plots_dir (str): Directory path to save the plot.
    """
    # Extract data from dictionary
    residues = list(dico.keys())
    means = [dico[res]['mean'] for res in residues]
    stds = [dico[res]['std'] for res in residues]

    # Create a compact figure (smaller size)
    plt.figure(figsize=(16, 6))

    # Plot means with error bars for std and connect dots with lines
    plt.errorbar(
        range(len(residues)),  # Use numeric positions for better line connection
        means,
        yerr=stds,
        fmt='-o',  # '-' connects dots, 'o' shows markers
        color='b',
        ecolor='r',
        capsize=3,  # Smaller caps
        capthick=1,  # Thinner caps
        linewidth=0.5,  # Thinner connecting line
        markersize=2,  # Smaller markers
        label='Mean ± Std'
    )

    # Add horizontal line at y=0 for reference
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)

    # Customize plot
    plt.title(f'Mean Frustration per Residue ({enc_type}, Monomer {monomer_number})', fontsize=10)
    plt.xlabel('Residue', fontsize=9)
    plt.ylabel('Mean Frustration', fontsize=9)

    # Customize x-axis: Show every 3rd residue to avoid clutter
    plt.xticks(
        range(len(residues))[::3],
        residues[::3],
        rotation=45,
        fontsize=7,
        ha='right'
    )

    plt.grid(True, linestyle=':', alpha=0.5)  # Lighter grid
    plt.legend(fontsize=8)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save and close
    name_plot = f"frustration_per_frame_{enc_type}_monomer_{monomer_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()



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
    #print(f"Loaded {len(monomers)} monomers")


    # 3. create aligned PDB files
    save_each_monomer_as_pdb(monomers, results_pdb_dir, enc_type, files_dict, monomer_number)

    # 4. calculation of frustration
    calculate_frustration(results_pdb_dir, results_frustration_dir)

    #6
    dico_monomers = dico_of_dico_frustIndex(results_frustration_dir)
    #print(dico_monomers)


    #7 grafical views

    #8.
    dico_mean = dico_mean_frustration(dico_monomers)
    #print(dico_mean)

    #9
    plot_frustration_per_res(dico_mean, enc_type, monomer_number, plots_dir)





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
