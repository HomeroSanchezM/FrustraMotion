from curses.ascii import isdigit

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
import subprocess
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from pandas.core.dtypes.inference import is_integer

'''
TO DO : 
make the frustration calculation in order for the files


'''
'''
This scripts take one PDB File containing multiples monomers (chains) and make a frustration study of the monomers using FrustratometeR 
the results files are added to results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_monomer_for_a_frame 
A plot is add to the plots/frustration directory, with name of type 
frustration_per_res_<TmEnc|MmEnc>_ALL_monomer_<(i)>.png   and 
frustration_per_res_<TmEnc|MmEnc>_monomer_<(i)>.png
Usage:

python3 frustration_plots.py path/to/file.pdb

If 2 file given, frustration will be calculated for the monomers of the 2 files and a scatter-plot of the frustration means for each residue will be made

python3 frustration_plots.py path/to/file1.pdb path/to/file2.pdb 


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
            "where <number> is one or more digits.\n"
            "If you're tring to use an option, make sure is write in the correct format:\n"
            "-option=value"
        )
#2.
def load_monomers(pdb_file):
    """Load the n monomers, corresponding to the chains in a list"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)

    chains = []
    for model in structure :
        for chain in model:
            chain_clean = []
            for residue in chain:
                name = residue.get_resname()
                if name != 'HSD' :
                    chain_clean.append(residue)
                else :
                    residue.resname = 'HIS'
                    chain_clean.append(residue)
            chains.append(chain_clean)

    return chains
#2.2
def load_pdb(pdb_file):
    """Load the PDB file and return it"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)
    return structure
#3.
def save_each_monomer_as_pdb(list_of_monomers, output_dir, enc_type1 , enc_number1, enc_type2 = None ,enc_number2 = None ):
    """
    Save each monomer in a file (atom number for 1 file is limited to 99999 by PDBIO)
    Args:
        list_of_monomers: List of chain object (Monomers)
        output_dir
        enc_type: Type of encapsulin (MtEnc/TmEnc) have to be generalised
        enc_number: Number of the encapsulin file

    """
    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"aligned monomere a {len(list_of_monomers)} monomers")
    for i, monomer in enumerate(list_of_monomers, start=1):
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
def calculate_frustration(PDB_directory, output_directory, seqdist_flag, mode_flag):

    """
    this function call the R package FrustratometeR for the calculation of frustration.
    :param PDB_directory: the directory of the pdb file used for the frustration calculation
    :param output_directory: the directory to save the frustration results
    """

    # check if directories exists
    #if not os.path.isdir(PDB_directory):
    #    raise ValueError(f"PDB directory not found: {PDB_directory}")
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
    results <- calculate_frustration(PdbFile = pdb_file, Mode = "{mode}", SeqDist = {seqdist}, ResultsDir = output_dir, Graphics = FALSE)

    """

    # Pour chaque fichier PDB dans le répertoire
    # TO DO make the frustration calculation in order for the files
    if os.path.isdir(PDB_directory):
        for pdb_file in os.listdir(PDB_directory):
            if pdb_file.endswith('.pdb'):
                full_path = os.path.join(PDB_directory, pdb_file)

                # Personnaliser le script R pour ce fichier
                current_script = r_script.format(
                    pdb_file=full_path,
                    output_dir=output_directory,
                    seqdist=seqdist_flag,
                    mode=mode_flag
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
    elif os.path.isfile(PDB_directory):
        if PDB_directory.endswith('.pdb'):
                full_path = os.path.join(PDB_directory)

                # Personnaliser le script R pour ce fichier
                current_script = r_script.format(
                    pdb_file=full_path,
                    output_dir=output_directory,
                    seqdist=seqdist_flag,
                    mode=mode_flag
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
def dico_frustIndex(frustration_file, enc_type):
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
        #count = 1
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue

            parts = line.split()

            # Check we have enough columns (at least 8)
            if len(parts) >= 8:
                if enc_type =="TmEnc":
                    res_num = int(parts[0]) - 3
                else :
                    res_num = int(parts[0])
                #res_num = count  # Residue number
                aa = parts[3]  # Amino acid (M, E, F...)
                ##remplace the '?' with the correct AA f the data
                #if aa == '?':
                #    if count == 187 :
                #        aa = 'Tyr'
                #    if count == 50 :
                #        aa = 'Hsd'
                frst_index = float(parts[7])  # FrstIndex value

                # Create key like "M4" and add to dictionary
                key = f"{aa}{res_num}"
                frustration_dict[key] = frst_index
                #count += 1

    return frustration_dict

def dico_of_dico_frustIndex(frustration_directory, enc_type):
    """
    this function take the full path of the directory were the frustration results files are
    And for each result, buit the dico of the residues frsutrationIndex, and add this to a big dico
    dico = { (monomer)1: {'M4': 0.759, ...} , 2:... , ... }
    :param frustration_directory:
    :return:
    """
    dico_monomers = {}
    for directory in os.listdir(frustration_directory) :
        doc_path = os.path.join(frustration_directory, directory, "FrustrationData")
        for file in os.listdir(doc_path):
            if file.endswith("pdb_singleresidue"):
                full_path = os.path.join(doc_path, file)
                name_file = os.path.basename(full_path)
                #print("nom du fichier :",name_file)
                match = re.search(r'monomer(\d+)', name_file)
                number = int(match.group(1))
                dico_monomer = dico_frustIndex(full_path, enc_type)
                dico_monomers[number]= dico_monomer

    return dico_monomers

def dico_frustIndex_with_chains(frustration_file, enc_type):
    """
    Create a dictionary  with frustration indices from FrustratometeR output file.

    Args:
        frustration_file (str): Path to the frustration results file

    Returns:
        dict: {chainID:{residue_key: FrstIndex}} where residue_key is "AAresNum" (e.g. "M4")

    Example output:
        {'0':{'M4': 0.759, 'E5': -0.582, 'F6': 1.027, ...}, '1':{'M4': 0.87, 'E5': -0.465, 'F6': 1.754, ...}}
    """
    frustration_dict = {}

    with open(frustration_file, 'r') as f:
        # Skip header
        next(f)
        #count = 1
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) >= 8:
                if enc_type =="TmEnc":
                    res_num = int(parts[0]) - 3
                else :
                    res_num = int(parts[0])
                #res_num = count
                chain_id = parts[1]
                aa = parts[3]
                if aa == '?':
                #    if count == 187 :
                    aa = 'Y'
                #    if count == 50 :
                #        aa = 'Hsd'
                frst_index = float(parts[7])

                # Initialize chain dictionary if needed
                if chain_id not in frustration_dict:
                    frustration_dict[chain_id] = {}
                    #count = 1
                    #res_num = count
                # Add residue frustration
                residue_key = f"{aa}{res_num}"
                frustration_dict[chain_id][residue_key] = frst_index
                #count += 1
            #print(frustration_dict)
    return frustration_dict

def dico_of_dico_frustIndex_with_chains(frustration_directory, enc_type):
    """
    this function take the full path of the directory were the frustration results files are
    And for each result, buit the dico of the residues frsutrationIndex, and add this to a big dico
    dico = { (monomer)1: {'M4': 0.759, ...} , 2:... , ... }
    :param frustration_directory:
    :return:
    """

    for directory in os.listdir(frustration_directory) :
        doc_path = os.path.join(frustration_directory, directory, "FrustrationData")
        for file in os.listdir(doc_path):
            if file.endswith("pdb_singleresidue"):
                full_path = os.path.join(doc_path, file)
                dico_monomer = dico_frustIndex_with_chains(full_path, enc_type)
                #print(dico_monomer)
    return dico_monomer


def dico_percentage_frustration_types(dico_monomers, mode_flag):
    """
    Calculate the percentage of each frustration type for all monomers.

    Frustration types are defined as:


    - minimal: frustration > 0.78
    - high: frustration < -1
    - neutral: -1 ≤ frustration ≤ 0.78

    Args:
        dico_monomers: Dictionary with format {monomer_num: {'res1': frustration, ...}, ...}

    Returns:
        Dictionary with format {monomer_num: {'high': %, 'minimal': %, 'neutral': %}, ...}
    """
    dico_percentage = {}

    for monomer_num, residue_data in dico_monomers.items():
        # Initialize counters
        counts = {'high': 0, 'minimal': 0, 'neutral': 0}
        total_residues = len(residue_data)

        if total_residues == 0:
            dico_percentage[monomer_num] = {'high': 0, 'minimal': 0, 'neutral': 0}
            continue

        # Classify each residue
        for frustration in residue_data.values():
            if mode_flag == 'singleresidue':
                if frustration > 0.58:
                    counts['minimal'] += 1
                elif frustration < -1:
                    counts['high'] += 1
                else:
                    counts['neutral'] += 1
            elif mode_flag == 'configurational':
                if frustration > 0.78:
                    counts['minimal'] += 1
                elif frustration < -1:
                    counts['high'] += 1
                else:
                    counts['neutral'] += 1
        # Calculate percentages
        percentages = {
            'high': round((counts['high'] / total_residues) * 100, 2),
            'minimal': round((counts['minimal'] / total_residues) * 100, 2),
            'neutral': round((counts['neutral'] / total_residues) * 100, 2)
        }

        dico_percentage[monomer_num] = percentages

    return dico_percentage





def dico_mean_percentage_frustration_types(dico_percentage):
    """
    Calculate mean and standard deviation of frustration type percentages across all monomers.

    Args:
        dico_percentage: Dictionary with format
            {monomer_num: {'high': %, 'minimal': %, 'neutral': %}, ...}

    Returns:
        Dictionary with format
            {'minimal': {'mean': float, 'std': float},
            'high': {'mean': float, 'std': float},
            'neutral': {'mean': float, 'std': float}}
    """
    # Initialize lists to collect all percentages for each type
    high_percents = []
    minimal_percents = []
    neutral_percents = []

    # Collect all percentages
    for monomer_data in dico_percentage.values():
        minimal_percents.append(monomer_data['minimal'])
        high_percents.append(monomer_data['high'])
        neutral_percents.append(monomer_data['neutral'])

    # Calculate mean and std for each frustration type
    dico_mean_percentage = {
        'minimal': {
            'mean': round(float(np.mean(minimal_percents)), 2),
            'std': round(float(np.std(minimal_percents)), 2)
        },
        'high': {
            'mean': round(float(np.mean(high_percents)), 2),
            'std': round(float(np.std(high_percents)), 2)
        },
        'neutral': {
            'mean': round(float(np.mean(neutral_percents)), 2),
            'std': round(float(np.std(neutral_percents)), 2)
        }
    }

    return dico_mean_percentage

#7. plot of all the frustration variation, have to be improved
def plot_all_frustration_per_res (dico_monomers, enc_type , enc_number, plots_dir ) :
    """
    this function generate a plot of all the frustration lines for the 60 monomers of an encapsuline for a given frame.
    :param dico_monomers: dico with the frustration info
    :param enc_type: <TmEnc|MtEnc>
    :param enc_number: number of the frame
    :param plots_dir: the output directory
    """

    plt.figure(figsize=(12, 8))

    # Get all residues (sorted numerically)
    sample_monomer = next(iter(dico_monomers.values()))  # Prend un monomer quelconque
    residues = sorted(sample_monomer.keys(), key=lambda x: int(''.join(filter(str.isdigit, x))))

    # Create x-axis positions
    x = np.arange(len(residues))

    # Plot each monomer
    for monomer_id, frustration_data in dico_monomers.items():
        # Get frustration values in residue order
        y = [frustration_data[res] for res in residues]
        plt.plot(x, y,
            marker='o',
            markersize=0.5,
            linestyle='-',
            linewidth=0.2,
            alpha=0.7,
            label=f'Monomer {monomer_id}')

    # Customize plot
    plt.xticks(range(len(residues))[::3], residues[::3],
               rotation=45, fontsize=8, ha='right')
    #plt.xticks(x, residues, rotation=45, ha='right', fontsize=8)
    plt.xlabel('Residue ')
    plt.ylabel('Frustration Index')
    plt.title('Frustration Profiles Across Monomers')
    plt.grid(True, alpha=0.3)

    # Legend outside the plot
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    name_plot = f"frustration_per_res_{enc_type}_ALL_monomer_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Plot with standard deviation saved as {plot_path}")

#8. plot mean frustration + std per residue
#8.1. dico of mean frustration for each res
def dico_mean_frustration (dico_monomers) :
    """
    create a dictionary with residue IDs as keys and {'mean': float, 'std': float} as values.
                   Example: {'M4': {'mean': -0.135, 'std': 0.668}, ...}
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

def dico_list_frustration (dico_monomers) :
    """
    create a dictionary with residue IDs as keys and a list of the frustration indexes as values.
                   Example: {'M4': [0.51 , 0.67, ... ], ...}
    :param dico_monomers:  dico_monomers = { (monomer)1: {'M4': 0.759, ...} , 2:... , ... }
    :return:
    """
     # Initialize dictionary to store sum and count
    dico_list = {}

    # First pass: collect sum and count for each residue
    for monomer_data in dico_monomers.values():
        for residue, frustration in monomer_data.items():
            if residue not in dico_list:
                dico_list[residue] =  []
            dico_list[residue].append(frustration)

    return dico_list


def plot_frustration_per_res(dico, enc_type, enc_number, plots_dir, seqdist, isolate):
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

    # Find positions of residues ending with 187 and 188
    pos_187 = None
    pos_188 = None

    for i, res in enumerate(residues):
        if res.endswith('187'):
            pos_187 = i
        elif res.endswith('188'):
            pos_188 = i
        # Exit loop if both positions are found
        if pos_187 is not None and pos_188 is not None:
            break

    # Highlight the single zone from 187 to 188 if both residues are found
    if pos_187 is not None and pos_188 is not None:
        # Determine start and end positions (handle case where 188 comes before 187)
        start = min(pos_187, pos_188)
        end = max(pos_187, pos_188)
        print(start)
        print(end)

        plt.axvspan(
            start - 0.5, end + 0.5,
            facecolor='red', alpha=0.2,
            label=f'Residues {residues[start]} to {residues[end]}'
        )
    elif pos_187 is not None:  # Only 187 found
        plt.axvspan(
            pos_187 - 0.5, pos_187 + 0.5,
            facecolor='red', alpha=0.2,
            label=f'Residue {residues[pos_187]}'
        )
    elif pos_188 is not None:  # Only 188 found
        plt.axvspan(
            pos_188 - 0.5, pos_188 + 0.5,
            facecolor='red', alpha=0.2,
            label=f'Residue {residues[pos_188]}'
        )

    # Plot means with error bars for std and connect dots with lines
    plt.errorbar(
        range(len(residues)),  # Use numeric positions for better line connection
        means,
        yerr=stds,
        fmt='o',  # '-' connects dots, 'o' shows markers
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
    if isolate:
        plt.title(f'Mean Frustration seqdist {seqdist} per Residue ({enc_type}, frame {enc_number}, isolate)',
                  fontsize=10)
    else:
        plt.title(f'Mean Frustration seqdist {seqdist} per Residue ({enc_type}, frame {enc_number}, not isolate)',
                  fontsize=10)
    plt.xlabel('Residue', fontsize=9)
    plt.ylabel('Mean Frustration', fontsize=9)

    # Customize x-axis: Show every 2nd residue to avoid clutter
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
    name_plot = f"frustration_per_res_{enc_type}_monomer_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def plot_min_and_max_frustration_per_res(dico_monomers, enc_type, enc_number, plots_dir):
        """
        Generate a plot showing min and max frustration values per residue as points.

        :param dico_monomers: Dictionary with frustration info {monomer_num: {'M4': 0.759, ...}, ...}
        :param enc_type: Encounter type identifier
        :param enc_number: Encounter number
        :param plots_dir: Directory to save the plot
        """
        # Initialize dictionary to store min and max values
        dico_min_max = {}

        # First pass: find min and max for each residue across all monomers
        for monomer_data in dico_monomers.values():
            for residue, frustration in monomer_data.items():
                if residue not in dico_min_max:
                    dico_min_max[residue] = {'min': frustration, 'max': frustration}
                else:
                    if frustration > dico_min_max[residue]['max']:
                        dico_min_max[residue]['max'] = frustration
                    if frustration < dico_min_max[residue]['min']:
                        dico_min_max[residue]['min'] = frustration

        residues = list(dico_min_max.keys())
        min_vals = [dico_min_max[res]['min'] for res in residues]
        max_vals = [dico_min_max[res]['max'] for res in residues]

        # Create figure
        plt.figure(figsize=(12, 6))

        # Plot min values as points
        for i, val in enumerate(min_vals):
            color = 'green' if val < 0 else 'red'
            plt.plot(i, val, 'o', color=color, markersize=6, label='Min' if i == 0 else "")

        # Plot max values as points
        for i, val in enumerate(max_vals):
            color = 'green' if val < 0 else 'red'
            plt.plot(i, val, 's', color=color, markersize=6, label='Max' if i == 0 else "")

        # Add horizontal line at y=0
        plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

        # Customize plot
        plt.title(f'Min and Max Frustration Values per Residue ({enc_type}, Monomer {enc_number})')
        plt.xlabel('Residue')
        plt.ylabel('Frustration Value')

        # Customize x-axis: Show every 3rd residue to avoid clutter
        plt.xticks(
            range(len(residues))[::3],
            residues[::3],
            rotation=45,
            fontsize=8,
            ha='right'
        )

        # Create custom legend
        green_patch = mpatches.Patch(color='green', label='Negative frustration')
        red_patch = mpatches.Patch(color='red', label='Positive frustration')
        plt.legend(handles=[
            green_patch,
            red_patch,
            plt.Line2D([0], [0], marker='o', color='w', label='Min',
                       markersize=8, markerfacecolor='black'),
            plt.Line2D([0], [0], marker='s', color='w', label='Max',
                       markersize=8, markerfacecolor='black')
        ])

        plt.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()

        # Save and close
        name_plot = f"frustration_min_max_per_res_{enc_type}_monomer_{enc_number}.png"
        plot_path = os.path.join(plots_dir, name_plot)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()

    #print(dico_min_max)
    #return dico_min_max




def plot_min_max_frustration_bars(dico_monomers, enc_type, enc_number, plots_dir, seqdist_flag, mode_flag):
    """
    Generate a plot showing vertical bars from min to max frustration values per residue.

    :param dico_monomers: Dictionary with frustration info {monomer_num: {'M4': 0.759, ...}, ...}
    :param enc_type: Encapsuline type identifier
    :param enc_number: frame number
    :param plots_dir: Directory to save the plot
    """
    # Initialize dictionary to store min and max values
    dico_min_max = {}

    # Find min and max for each residue across all monomers
    for monomer_data in dico_monomers.values():
        for residue, frustration in monomer_data.items():
            if residue not in dico_min_max:
                dico_min_max[residue] = {'min': frustration, 'max': frustration}
            else:
                if frustration > dico_min_max[residue]['max']:
                    dico_min_max[residue]['max'] = frustration
                if frustration < dico_min_max[residue]['min']:
                    dico_min_max[residue]['min'] = frustration

    residues = list(dico_min_max.keys())
    min_vals = [dico_min_max[res]['min'] for res in residues]
    max_vals = [dico_min_max[res]['max'] for res in residues]

    # Create figure
    plt.figure(figsize=(16, 6))

    # Plot vertical bars for each residue
    if mode_flag == 'singleresidue' :
        for i, (res, vmin, vmax) in enumerate(zip(residues, min_vals, max_vals)):
            # Determine color based on values
            if vmax <= -1:
                color = 'red'  # Both min and max are negative
            elif vmin >= 0.58:
                color = 'green'  # Both min and max are positive
            elif vmax < 0.58 and vmin > -1 :
                color = 'grey'
            else:
                color = 'gold'  # Spanning both negative and positive
             # Draw vertical line from min to max
            plt.vlines(x=i, ymin=vmin, ymax=vmax,
                   colors=color, linewidth=2, alpha=0.7)

            # Add small markers for min and max points
            plt.plot(i, vmin, 'o', color='black', markersize=3)
            plt.plot(i, vmax, 'o', color='black', markersize=3)
    elif mode_flag == 'configurational' :
        for i, (res, vmin, vmax) in enumerate(zip(residues, min_vals, max_vals)):
            # Determine color based on values
            if vmax <= -1:
                color = 'red'  # Both min and max are negative
            elif vmin >= 0.78:
                color = 'green'  # Both min and max are positive
            elif vmax < 0.78 and vmin > -1:
                color = 'grey'
            else:
                color = 'gold'  # Spanning both negative and positive

            # Draw vertical line from min to max
            plt.vlines(x=i, ymin=vmin, ymax=vmax,
                    colors=color, linewidth=2, alpha=0.7)

            # Add small markers for min and max points
            plt.plot(i, vmin, 'o', color='black', markersize=3)
            plt.plot(i, vmax, 'o', color='black', markersize=3)

    # Add horizontal line at y=0
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

    # Customize plot
    plt.title(f'Frustration seqdist {seqdist_flag} Range per Residue ({enc_type}, frame {enc_number})')
    plt.xlabel('Residue')
    plt.ylabel('Frustration Value')

    # Customize x-axis
    plt.xticks(
        range(len(residues))[::3],
        residues[::3],
        rotation=45,
        fontsize=8,
        ha='right'
    )

    # Create custom legend
    green_patch = mpatches.Patch(color='green', label='minimal frustration', alpha=0.7)
    red_patch = mpatches.Patch(color='red', label='High frustration', alpha=0.7)
    grey_patch = mpatches.Patch(color='grey', label='neutral frustration', alpha=0.7)
    gold_patch = mpatches.Patch(color='gold', label='Mixed frustration', alpha=0.7)
    plt.legend(handles=[green_patch, grey_patch, red_patch, gold_patch])

    plt.grid(True, linestyle=':', alpha=0.4)
    plt.tight_layout()

    # Save plot
    name_plot = f"frustration_range_per_res_{enc_type}_frame_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import matplotlib.patches as mpatches

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

def plot_frustration_boxplots(dico_monomer, enc_type, enc_number, plots_dir, seqdist, isolate):
    """
    Generate boxplots of frustration values per residue from a dictionary of lists.
    Box color depends on the distribution:
      - Red: all values < -1
      - Green: all values > 0.58
      - Gray: all values between -1 and 0.58
      - Blue: mixed values
    Y-axis is inverted (negative up). Residues 187 and 188 are highlighted with red background.

    Parameters:
    - dico_monomer: Dict with residue IDs as keys and list of frustration values
    - enc_type: Encounter type
    - enc_number: Monomer number
    - plots_dir: Directory to save the plot
    - seqdist: Sequence distance (used in title)
    - isolate: Boolean flag for title customization
    """

    df = pd.DataFrame.from_dict(dico_monomer, orient='index').T
    df_melted = df.melt(var_name='Residue', value_name='Frustration')

    plt.figure(figsize=(16, 6))
    ax = sns.boxplot(
        data=df_melted,
        x='Residue',
        y='Frustration',
        showfliers=True,
        width=0.75,
        linewidth=1,
        fliersize=2,
        color='skyblue'
    )

    # Highlight residues ending with 187 or 188
    highlight_suffixes = ['189']
    residue_list = list(dico_monomer.keys())
    for i, res in enumerate(residue_list):
        if any(res.endswith(suffix) for suffix in highlight_suffixes):
            ax.axvspan(i - 0.5, i + 0.5, color='lightcoral', alpha=0.3, zorder=0)

    # Color each box
    for residue, patch in zip(residue_list, ax.patches):
        values = dico_monomer[residue]
        all_below = all(val < -1 for val in values)
        all_above = all(val > 0.58 for val in values)
        all_middle = all(-1 <= val <= 0.58 for val in values)

        if all_below:
            patch.set_facecolor('red')
            patch.set_edgecolor('darkred')
        elif all_above:
            patch.set_facecolor('green')
            patch.set_edgecolor('darkgreen')
        elif all_middle:
            patch.set_facecolor('gray')
            patch.set_edgecolor('black')
        else:
            patch.set_facecolor('skyblue')
            patch.set_edgecolor('steelblue')

    # Invert Y-axis
    plt.ylim(-3.5, 3.5)
    ax.invert_yaxis()


    # Horizontal reference lines
    plt.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    plt.axhline(-1, color='red', linestyle=':', linewidth=1, alpha=0.5)
    plt.axhline(0.58, color='green', linestyle=':', linewidth=1, alpha=0.5)

    # X-axis formatting
    xticks = ax.get_xticks()
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    ax.set_xticks(xticks[::3])
    ax.set_xticklabels(xlabels[::3], rotation=45, ha='right', fontsize=7)

    # Vertical grid lines
    for x in xticks:
        plt.axvline(x=x, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)

    # Title and labels
    if isolate:
        plt.title(f'Frustration Distribution per Residue ({enc_type}, frame {enc_number}, seqdist {seqdist}, isolate)',
                  fontsize=10)
    else:
        plt.title(f'Frustration Distribution per Residue ({enc_type}, frame {enc_number}, seqdist {seqdist}, not isolate)',
                  fontsize=10)

    plt.xlabel('Residue', labelpad=10)
    plt.ylabel('Frustration Index ', labelpad=10)
    plt.yticks(fontsize=10)
    plt.grid(axis='y', alpha=0.3)

    # Legend
    red_patch = mpatches.Patch(color='red', label='All values < -1')
    green_patch = mpatches.Patch(color='green', label='All values > 0.58')
    gray_patch = mpatches.Patch(color='gray', label='All between -1 and 0.58')
    blue_patch = mpatches.Patch(color='skyblue', label='Mixed values')
    plt.legend(handles=[red_patch, green_patch, gray_patch, blue_patch], fontsize=8)

    # Save plot
    plt.tight_layout()
    name_plot = f"frustration_boxplot_per_res_{enc_type}_monomer_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()



def plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_number1, enc_type2, enc_number2, plots_dir, seqdist_flag, mode_flag):
    """
    Create a scatter plot comparing mean frustration values between two conditions.
    Points are colored:
    - Green if both x and y are positive
    - Red if both x and y are negative
    - Gold otherwise

    Parameters:
    - dico_mean1: {'M4': {'mean': -0.135, 'std': 0.668}, ...} for first condition
    - dico_mean2: Same structure for second condition
    - enc_type1, enc_number1: Identifier for first condition
    - enc_type2, enc_number2: Identifier for second condition
    - plots_dir: Directory to save the plot
    """
    # Prepare data and determine colors
    colors = []
    x_vals = []
    y_vals = []

    # Get residue names list
    residues = list(dico_mean1.keys())
    residues2 = list(dico_mean2.keys())
    # Prepare data
    x_vals = [dico_mean2[res]['mean'] for res in dico_mean2.keys()]
    y_vals = [dico_mean1[res]['mean'] for res in dico_mean1.keys()]

    if mode_flag == 'singleresidue' :
        for k in range (len(x_vals)):
            if x_vals[k] >= 0.58 and y_vals[k] >= 0.58:
                colors.append('green')
            elif x_vals[k] <= -1 and y_vals[k] <= -1:
                colors.append('red')
            elif -1<x_vals[k]<0.58  and -1<y_vals[k]<0.58:
                colors.append('grey')
            elif x_vals[k]>y_vals[k]:
                colors.append('yellow')
            else :
                colors.append('goldenrod')
    if mode_flag == 'configurational' :
        for k in range (len(x_vals)):
            if x_vals[k] >= 0.78 and y_vals[k] >= 0.78:
                colors.append('green')
            elif x_vals[k] <= -1 and y_vals[k] <= -1:
                colors.append('red')
            elif -1<x_vals[k]<0.78  and -1<y_vals[k]<0.78:
                colors.append('grey')
            elif x_vals[k]>y_vals[k]:
                colors.append('yellow')
            else :
                colors.append('goldenrod')
    # Create figure
    plt.figure(figsize=(8, 8))

    # Create scatter plot with colored points
    plt.scatter(x_vals, y_vals, c=colors, alpha=0.7, s=60, edgecolors='black', linewidths=0.5)

    # Highlight points 187 to 188 with dark red outline and labels
    highlight_indices1 = range(186, 188)
    if highlight_indices1[-1] < len(residues):  # Check if indices are valid
        # Plot highlighted points with dark red outline
        plt.scatter(
            [x_vals[i] for i in highlight_indices1],
            [y_vals[i] for i in highlight_indices1],
            c=[colors[i] for i in highlight_indices1],
            s=80,  # Slightly larger size
            edgecolors='darkred',
            linewidths=1.5,
            alpha=0.9
        )

        # Add labels for highlighted points
        for i in highlight_indices1:
            plt.text(
                x_vals[i], y_vals[i],f"{residues[i]}/{residues2[i]}",
                fontsize=8,
                ha='center',
                va='bottom',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
            )

    # Highlight points 134 to 137 with dark red outline and labels
    #highlight_indices2 = range(133, 137)  # 134 to 137 (Python ranges are exclusive on end)
    #if highlight_indices2[-1] < len(residues):  # Check if indices are valid
        # Plot highlighted points with dark red outline
    #    plt.scatter(
    #        [x_vals[i] for i in highlight_indices2],
    #        [y_vals[i] for i in highlight_indices2],
    #        c=[colors[i] for i in highlight_indices2],
    #        s=80,  # Slightly larger size
    #        edgecolors='darkred',
    #        linewidths=1.5,
    #        alpha=0.9
    #    )

        # Add labels for highlighted points
    #    for i in highlight_indices2:
    #        plt.text(
    #            x_vals[i], y_vals[i],f"{residues[i]}/{residues2[i]}",
    #            fontsize=8,
    #            ha='center',
    #            va='bottom',
    #            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
    #        )

    # Add reference lines
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)

    # Add identity line
    max_val = max(max(np.abs(x_vals)), max(np.abs(y_vals))) * 1.1
    plt.plot([-max_val, max_val], [-max_val, max_val], 'k--', alpha=0.3)

    # Customize plot
    plt.title(
        f'Mean monomers frustration seqdist {seqdist_flag} comparison\n{enc_type1} frame {enc_number1} vs {enc_type2} frame {enc_number2}')
    plt.xlabel(f'Mean Frustration ({enc_type2} {enc_number2})')
    plt.ylabel(f'Mean Frustration ({enc_type1} {enc_number1})')
    plt.grid(True, linestyle=':', alpha=0.3)

    # Create legend
    import matplotlib.patches as mpatches
    legend_elements = [
        mpatches.Patch(color='green', label='minimal frustration conserved'),
        mpatches.Patch(color='red', label='High frustration conserved'),
        mpatches.Patch(color='grey', label='neutral frustration conserved'),
        mpatches.Patch(color='yellow', label='MtEnc more frustrated than TmEnc'),
        mpatches.Patch(color='goldenrod', label='TmEnc more frustrated than MtEnc')
    ]
    plt.legend(handles=legend_elements, loc='best')

    # Equal aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')

    # Save plot
    name_plot = f"scatter_frustration_per_res_{enc_type1}_frame_{enc_number1}_vs_{enc_type2}_frame_{enc_number2}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()
    return colors


def visualization_VMD_script(output_dir1, output_dir2, enc_type1, enc_type2, enc_number1, enc_number2, green_num,
                             red_num, grey_num, yellow_num, goldenrod_num, execute_vmd=False):
    """
    This function generates a tcl script for frustration by single residues visualization using VMD,
    which gives a specific color to each residue of the first monomer PDB files of the result dir.

    :param output_dir1: first directory where to create the tcl file
    :param output_dir2: second directory where to create the tcl file
    :param enc_type1: <TmEnc|MtEnc> for the first file
    :param enc_type2: <TmEnc|MtEnc> for the second file
    :param enc_number1: the number of the first frame studied
    :param enc_number2: the number of the second frame studied
    :param green_num: list of residues for group1 (green)
    :param red_num: list of residues for group2 (red)
    :param grey_num: list of residues for group3 (grey)
    :param yellow_num: list of residues for group4 (yellow)
    :param goldenrod_num: list of residues for group5 (goldenrod)
    :param execute_vmd: if True, automatically executes the script with VMD
    """

    def generate_script_content(output_dir, enc_type, enc_number):
        return f"""# Charger la molécule
mol new {output_dir}/{enc_type}{enc_number}_monomer1.pdb type pdb waitfor all

# Supprimer les représentations existantes
mol delrep 0 top

# Représentation de base pour toute la molécule
mol representation NewCartoon
mol color Name
mol selection "all"
mol material Opaque
mol addrep top

# Définir les groupes de résidus
set group1 {{{green_num}}}
set group2 {{{red_num}}}
set group3 {{{grey_num}}}
set group4 {{{yellow_num}}}
set group5 {{{goldenrod_num}}}

# Fonction pour créer une sélection et une représentation
proc add_residue_group {{resid_list color_id}} {{
    set selection_text ""
    foreach r $resid_list {{
        append selection_text "resid $r or "
    }}
    set selection_text [string range $selection_text 0 end-4]

    mol representation NewCartoon
    mol selection $selection_text
    mol color ColorID $color_id
    mol material Opaque
    mol addrep top
}}

# apply colors for each group
add_residue_group $group1 19  ;# green2
add_residue_group $group2 1   ;# red
add_residue_group $group3 6   ;# silver
add_residue_group $group4 4   ;# yellow
add_residue_group $group5 31   ;# orange
"""

    # Create output directories if they don't exist
    os.makedirs(output_dir1, exist_ok=True)
    os.makedirs(output_dir2, exist_ok=True)

    # Define output file paths
    output_file1 = os.path.join(output_dir1, "frustration_per_res_colors.tcl")
    output_file2 = os.path.join(output_dir2, "frustration_per_res_colors.tcl")

    # Generate and write script for output_file1
    script_content1 = generate_script_content(output_dir1, enc_type1, enc_number1)
    with open(output_file1, 'w') as f:
        f.write(script_content1)

    # Generate and write script for output_file2 (with different parameters)
    script_content2 = generate_script_content(output_dir2, enc_type2, enc_number2)
    with open(output_file2, 'w') as f:
        f.write(script_content2)

    # Execute VMD if requested
    if execute_vmd:
        try:
            subprocess.run(["vmd", "-e", output_file1], check=True)
            print(f"Successfully executed VMD with script: {output_file1}")
            subprocess.run(["vmd", "-e", output_file2], check=True)
            print(f"Successfully executed VMD with script: {output_file2}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing VMD: {e}")
        except FileNotFoundError:
            print("VMD not found. Please make sure VMD is installed and in your PATH.")



def plot_barplot_percentage_frustration_types(plots_dir, dico_mean_type_1 , dico_mean_type_2, enc_type1, enc_number1, enc_type2, enc_number2, seqdist_flag):
    """
    Crée un barplot comparant les types de frustration entre MtEnc et TmEnc
    avec barres d'erreur pour les écarts-types
    :param plots_dir: dossier où sauvegarder le graphique
    :param dico_mean_type_1: dictionnaire des stats pour MtEnc (format {'minimal': {'mean': x, 'std': y}, ...})
    :param dico_mean_type_2: dictionnaire des stats pour TmEnc (même format)
    """
    # Extraction des données depuis les dictionnaires
    categories = ['Minimal', 'High', 'Neutral']

    # Données pour MtEnc (dico_mean_type_1)
    mtenc_means = [
        dico_mean_type_1['minimal']['mean'],
        dico_mean_type_1['high']['mean'],
        dico_mean_type_1['neutral']['mean']
    ]
    mtenc_stds = [
        dico_mean_type_1['minimal']['std'],
        dico_mean_type_1['high']['std'],
        dico_mean_type_1['neutral']['std']
    ]

    # Données pour TmEnc (dico_mean_type_2)
    tmenc_means = [
        dico_mean_type_2['minimal']['mean'],
        dico_mean_type_2['high']['mean'],
        dico_mean_type_2['neutral']['mean']
    ]
    tmenc_stds = [
        dico_mean_type_2['minimal']['std'],
        dico_mean_type_2['high']['std'],
        dico_mean_type_2['neutral']['std']
    ]

    # Paramètres du graphique
    bar_width = 0.35
    index = np.arange(len(categories))

    # Création du plot
    plt.figure(figsize=(10, 6))

    # Barres pour MtEnc
    plt.bar(index - bar_width / 2, mtenc_means, bar_width,
            yerr=mtenc_stds, capsize=5,
            label='MtEnc', color='red')

    # Barres pour TmEnc
    plt.bar(index + bar_width / 2, tmenc_means, bar_width,
            yerr=tmenc_stds, capsize=5,
            label='TmEnc', color='goldenrod')

    # Personnalisation
    plt.title(f'Comparison of mean Frustration (seqdist {seqdist_flag}) types between MtEnc and TmEnc monomers', fontsize=14)
    plt.xlabel('Frustration Type', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.xticks(index, categories)
    plt.ylim(0, 70)

    # Ajout des valeurs sur les barres
    for i in range(len(categories)):
        plt.text(i - bar_width / 2, mtenc_means[i] + 2, f'{mtenc_means[i]:.2f}%',
            ha='center', va='bottom')
        plt.text(i + bar_width / 2, tmenc_means[i] + 2, f'{tmenc_means[i]:.2f}%',
            ha='center', va='bottom')

    # Légende et grille
    plt.legend(loc='upper right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Afficher le plot
    plt.tight_layout()
    # Save plot
    name_plot = f"frustration_barplot_per_res_{enc_type1}_frame_{enc_number1}_vs_{enc_type2}_frame_{enc_number2}.png"

    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()


def parse_arguments():
    """Parse command line arguments including optional -vmd and -frustration flags"""
    vmd_flag = False
    frustration_flag = False
    seqdist_flag = 12
    isolate_flag = True
    plots_flag = True
    mode_flag = 'singleresidue'
    pdb_files = []

    for arg in sys.argv[1:]:
        if arg.startswith('-vmd='):
            vmd_value = arg.split('=')[1].lower()
            if vmd_value == 'true':
                vmd_flag = True
            elif vmd_value == 'false':
                vmd_flag = False
            else :
                raise ValueError(
            f"the option -vmd expect a boolean (True or False), make sure is write in the correct format:\n"
            "-option=value")

        elif arg.startswith('-frustration='):
            frustration_value = arg.split('=')[1].lower()
            if frustration_value == 'true':
                frustration_flag = True
            elif frustration_value == 'false':
                frustration_flag = False
            else:
                 raise ValueError(
            f"the option -frustration expect a boolean (True or False), make sure is write in the correct format:\n"
            "-option=value")

        elif arg.startswith('-seqdist='):
            seqdist_flag = arg.split('=')[1]
            try:
                int(seqdist_flag)
            except :
                raise ValueError(
            f"the option -seqdist expect a integer, make sure is write in the correct format:\n"
            "-option=value")

        elif arg.startswith('-isolate='):
            isolate_value = arg.split('=')[1].lower()
            if isolate_value == 'false':
                isolate_flag = False
            elif isolate_value == 'true':
                isolate_flag = True
            else :
                raise ValueError(
            f"the option -isolate expect a boolean (True or False), make sure is write in the correct format:\n"
            "-option=value")
        elif arg.startswith('-plots='):
            plots_value = arg.split('=')[1].lower()
            if plots_value == 'false':
                plots_flag = False
            elif plots_value == 'true':
                plots_flag = True
            else :
                raise ValueError(
            f"the option -isolate expect a boolean (True or False), make sure is write in the correct format:\n"
            "-option=value")
        elif arg.startswith('-mode='):
            mode_flag = arg.split('=')[1].lower()
            print(mode_flag)
            if mode_flag != 'singleresidue' and mode_flag!='configurational' and mode_flag!='mutational':
                raise ValueError(
            f"the option -mode expect one of the following options : \n singleresidue \n configurational \n mutational \n , make sure is write in the correct format:\n"
            "-option=value")
        else:
            pdb_files.append(arg)

    return pdb_files, vmd_flag, frustration_flag, seqdist_flag, isolate_flag, plots_flag, mode_flag

#when only a file given
def main(pdb_file1, vmd_flag= False, frustration_flag=False, seqdist_flag = 12, isolate_flag=True, plots_flag=True, mode_flag='singleresidue' ):
    start_time = time.time()
    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filename
    enc_type, enc_number = extract_info_from_filename(pdb_file1)
    # 1.1. create the output directory path
    frustration_dir = f"FRUSTRATION_{enc_type.upper()}"
    capsids_dir = f"{enc_type.upper()}_CAPSIDS"
    monomer_dir = "FRUSTRATION_monomer_for_a_frame"

    # all path of the actual direcory (ls -a)
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # cd ..
    base_dir = os.path.dirname(current_dir)
    # cd results/
    results_dir = os.path.join(base_dir, "results")

    if isolate_flag:

        isolate_dir = "Isolate"

        results_pdb_dir = os.path.join(results_dir, frustration_dir, capsids_dir, monomer_dir,isolate_dir, f"{enc_type}{enc_number}_monomers")
        results_frustration_dir = os.path.join(results_dir, frustration_dir, capsids_dir, monomer_dir,isolate_dir, f"{enc_type}{enc_number}_frustration_seqdist_{seqdist_flag}_isolate_mode_{mode_flag}")

        plots_dir = os.path.join("../plots", enc_type, f"frustration_seqdist_{seqdist_flag}_isolate")

        #1.2 create the repository if it not exist
        os.makedirs(results_pdb_dir, exist_ok=True)
        os.makedirs(results_frustration_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)

        # 2. Charger les monomères
        monomers = load_monomers(pdb_file1)
        print(f"Loaded {len(monomers)} monomers")

        # 3. create PDB files
        save_each_monomer_as_pdb(monomers, results_pdb_dir, enc_type, enc_number)

        # 4. calculation of frustration
        if frustration_flag:
            # 4. calculation of frustration
            calculate_frustration(results_pdb_dir, results_frustration_dir, seqdist_flag, mode_flag)


        #6
        dico_monomers = dico_of_dico_frustIndex(results_frustration_dir, enc_type)
        #print(dico_monomers)


        #7 grafical views

        #plot_all_frustration_per_res(dico_monomers, enc_type, enc_number, plots_dir)

        #8.
        dico_mean = dico_mean_frustration(dico_monomers)
        #print(dico_mean)
        if plots_flag:
            #9
            plot_frustration_per_res(dico_mean, enc_type, enc_number, plots_dir, seqdist_flag, isolate_flag)

            #10
            #plot_min_and_max_frustration_per_res(dico_monomers, enc_type, enc_number, plots_dir)
            plot_min_max_frustration_bars(dico_monomers, enc_type, enc_number, plots_dir, seqdist_flag, mode_flag)

            dico_list = dico_list_frustration(dico_monomers)
            plot_frustration_boxplots(dico_list, enc_type, enc_number, plots_dir, seqdist_flag, isolate_flag)

        # % of frustration type
        dico_types = dico_percentage_frustration_types(dico_monomers, mode_flag)
        #print(dico_types)

        dico_mean_types = dico_mean_percentage_frustration_types(dico_types)
        #print(dico_mean_types)
        return dico_mean
    else:

        not_isolate_dir = "Not_Isolate"
        results_frustration_dir = os.path.join(results_dir, frustration_dir, capsids_dir, monomer_dir,not_isolate_dir,
                                               f"{enc_type}{enc_number}_frustration_seqdist_{seqdist_flag}_NOT_isolate_mode_{mode_flag}")

        plots_dir = os.path.join("../plots", enc_type, f"frustration_seqdist_{seqdist_flag}_NOT_isolate")

        # 1.2 create the repository if it not exist
        os.makedirs(results_frustration_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)

        # 4. calculation of frustration
        if frustration_flag:
            # 4. calculation of frustration
            calculate_frustration(pdb_file1, results_frustration_dir, seqdist_flag, mode_flag)

        dico_monomer = dico_of_dico_frustIndex_with_chains(results_frustration_dir, enc_type)
        print(dico_monomer)
        if plots_flag :
            dico_mean = dico_mean_frustration(dico_monomer)
            plot_frustration_per_res(dico_mean, enc_type, enc_number, plots_dir, seqdist_flag, isolate_flag)
            #print(dico_mean)
            dico_mean_isolate = main(pdb_file1, vmd_flag, frustration_flag, seqdist_flag, isolate_flag=True)
            #print(dico_mean_isolate)

            dico_list = dico_list_frustration(dico_monomer)
            plot_frustration_boxplots(dico_list, enc_type, enc_number, plots_dir, seqdist_flag, isolate_flag)
            # Créer un nouveau dictionnaire pour stocker les différences
            #diff_dict = {}
            # Parcourir les clés communes aux deux dictionnaires
            #for residue in dico_mean.keys() :
            #    diff = dico_mean_isolate[residue]['mean'] - dico_mean[residue]['mean']
            #    diff_dict[residue] = diff
            #print(diff_dict)
            # Préparation des données pour le plot
            #residues = list(diff_dict.keys())
            #values = list(diff_dict.values())
            #x_pos = np.arange(len(residues))  # positions des résidus sur l'axe x

            # Création du plot
            #plt.figure(figsize=(15, 6))

            # Barres ou points pour les valeurs de frustration
            #plt.scatter(x_pos, values, color='blue', label='Frustration résiduelle')
            # Ou utiliser des barres :
            # plt.bar(x_pos, values, color='blue', alpha=0.7, label='Frustration résiduelle')

            # Ligne horizontale à y=0
            #plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

            # Personnalisation du plot
            #plt.xticks(
            #range(len(residues))[::3],
            #residues[::3],
            #rotation=45,
            #fontsize=7,
            #ha='right'
            #)
            #plt.ylabel('delta Frustration ISOLATE - NOT ISOLATE ')
            #plt.title('delta frsutration per res')
            #plt.grid(axis='y', linestyle=':', alpha=0.5)
            #plt.tight_layout()  # Ajuste les marges

             # Save plot
            #name_plot = f"delta_frustration_per_res_{enc_type}_frame_{enc_number}.png"
            #plot_path = os.path.join(plots_dir, name_plot)
            #plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
            #plt.close()
        else:
            return dico_monomer


        print("Work in progress for isolate= False")

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")


#when 2 files given in parameters
def main2(pdb_file1, pdb_file2, vmd_flag = False, frustration_flag=False, seqdist_flag = 12, isolate_flag=True, mode_flag='singleresidue'):

    start_time = time.time()

    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filenames
    enc_type1, enc_number1 = extract_info_from_filename(pdb_file1)
    enc_type2, enc_number2 = extract_info_from_filename(pdb_file2)

    # 1.1. create the output directories paths

    # all path of the actual direcory (ls -a)
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # cd ..
    base_dir = os.path.dirname(current_dir)
    # cd results/
    results_dir = os.path.join(base_dir, "results")

    frustration_dir1 = f"FRUSTRATION_{enc_type1.upper()}"
    capsids_dir1 = f"{enc_type1.upper()}_CAPSIDS"

    frustration_dir2 = f"FRUSTRATION_{enc_type2.upper()}"
    capsids_dir2 = f"{enc_type2.upper()}_CAPSIDS"

    monomer_dir = "FRUSTRATION_monomer_for_a_frame"
    if isolate_flag:
        isolate_dir = "Isolate"

        results_pdb_dir1 = os.path.join(results_dir, frustration_dir1, capsids_dir1,monomer_dir, isolate_dir, f"{enc_type1}{enc_number1}_monomers")
        results_frustration_dir1 = os.path.join(results_dir, frustration_dir1, capsids_dir1,monomer_dir,isolate_dir, f"{enc_type1}{enc_number1}_frustration_seqdist_{seqdist_flag}_isolate_mode_{mode_flag}")

        results_pdb_dir2 = os.path.join(results_dir, frustration_dir2, capsids_dir2,monomer_dir,isolate_dir, f"{enc_type2}{enc_number2}_monomers")
        results_frustration_dir2 = os.path.join(results_dir, frustration_dir2, capsids_dir2,monomer_dir,isolate_dir, f"{enc_type2}{enc_number2}_frustration_seqdist_{seqdist_flag}_isolate_mode_{mode_flag}")

        plots_dir = os.path.join("../plots", "COMMON", f"frustration_seqdist_{seqdist_flag}_isolate")

        #1.2 create the repository if it not exist
        os.makedirs(results_pdb_dir1, exist_ok=True)
        os.makedirs(results_frustration_dir1, exist_ok=True)
        os.makedirs(results_pdb_dir2, exist_ok=True)
        os.makedirs(results_frustration_dir2, exist_ok=True)

        os.makedirs(plots_dir, exist_ok=True)

        # 2. Charger les monomères
        monomers1 = load_monomers(pdb_file1)
        print(f"Loaded {len(monomers1)} monomers")

        monomers2 = load_monomers(pdb_file2)
        print(f"Loaded {len(monomers2)} monomers")


         # 3. create PDB files
        save_each_monomer_as_pdb(monomers1, results_pdb_dir1, enc_type1, enc_number1)
        save_each_monomer_as_pdb(monomers2, results_pdb_dir2, enc_type2, enc_number2)

        # 4. calculation of frustration
        if frustration_flag :
            calculate_frustration(results_pdb_dir1, results_frustration_dir1, seqdist_flag, mode_flag)
            calculate_frustration(results_pdb_dir2, results_frustration_dir2, seqdist_flag, mode_flag)

        #6
        dico_monomers1 = dico_of_dico_frustIndex(results_frustration_dir1, enc_type1)
        dico_monomers2 = dico_of_dico_frustIndex(results_frustration_dir2, enc_type2)

        #8.
        dico_mean1 = dico_mean_frustration(dico_monomers1)
        dico_mean2 = dico_mean_frustration(dico_monomers2)


       # graphical view, scatter plot
        colors = plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_number1, enc_type2, enc_number2, plots_dir, seqdist_flag, mode_flag)
        #print(colors)
        green_list = []
        red_list = []
        grey_list = []
        yellow_list = []
        goldenrod_list = []
        for k in range(len(colors)):
            if colors[k] == "green" :
                green_list.append(k+1)
            elif colors[k] == "red" :
                red_list.append(k+1)
            elif colors[k] == "grey" :
                grey_list.append(k+1)
            elif colors[k] == "yellow" :
                yellow_list.append(k+1)
            elif colors[k] == "goldenrod" :
                goldenrod_list.append(k+1)
        green_num = ""
        for num in green_list :
            green_num+=str(num)+" "
        #print( green_num)

        red_num = ""
        for num in red_list :
            red_num+=str(num)+" "
        #print( red_num)

        grey_num = ""
        for num in grey_list :
            grey_num+=str(num)+" "
        #print( grey_num)

        yellow_num = ""
        for num in yellow_list :
            yellow_num+=str(num)+" "
        #print( yellow_num)

        goldenrod_num = ""
        for num in goldenrod_list :
            goldenrod_num+=str(num)+" "
        #print( goldenrod_num)

        #print(green_list)
        #print(red_list)
        #print(grey_list)
        #print(yellow_list)
        #print(goldenrod_list)
        #coloring structure
        #color_structure_by_index(colors,  results_pdb_dir1, results_pdb_dir2)

        #plot visualisation
        visualization_VMD_script(results_pdb_dir1, results_pdb_dir2, enc_type1, enc_type2, enc_number1, enc_number2, green_num, red_num, grey_num, yellow_num, goldenrod_num, vmd_flag)

        # % of frustration type
        dico_types_1 = dico_percentage_frustration_types(dico_monomers1, mode_flag)
        #print(dico_types_1)
        dico_types_2 = dico_percentage_frustration_types(dico_monomers2, mode_flag)
        #print(dico_types_2)

        # % of each type of frustration barplot

        dico_mean_types_1 = dico_mean_percentage_frustration_types(dico_types_1)
        print(dico_mean_types_1)
        dico_mean_types_2 = dico_mean_percentage_frustration_types(dico_types_2)
        print(dico_mean_types_2)

        plot_barplot_percentage_frustration_types(plots_dir, dico_mean_types_1, dico_mean_types_2, enc_type1, enc_number1, enc_type2, enc_number2, seqdist_flag )
    else :
        not_isolate_dir = "Not_Isolate"
        results_frustration_dir1 = os.path.join(results_dir, frustration_dir1, capsids_dir1, monomer_dir, not_isolate_dir,
                                                f"{enc_type1}{enc_number1}_frustration_seqdist_{seqdist_flag}_NOT_isolate_mode_{mode_flag}")

        results_frustration_dir2 = os.path.join(results_dir, frustration_dir2, capsids_dir2, monomer_dir, not_isolate_dir,
                                                f"{enc_type2}{enc_number2}_frustration_seqdist_{seqdist_flag}_NOT_isolate_mode_{mode_flag}")

        plots_dir = os.path.join("../plots", "COMMON", f"frustration_seqdist_{seqdist_flag}_NOT_isolate")

        # 1.2 create the repository if it not exist

        os.makedirs(results_frustration_dir1, exist_ok=True)

        os.makedirs(results_frustration_dir2, exist_ok=True)

        os.makedirs(plots_dir, exist_ok=True)

         # 4. calculation of frustration
        if frustration_flag :
            calculate_frustration(pdb_file1, results_frustration_dir1, seqdist_flag, mode_flag)
            calculate_frustration(pdb_file2, results_frustration_dir2, seqdist_flag, mode_flag)

        #6
        dico_monomers1 = dico_of_dico_frustIndex_with_chains(results_frustration_dir1, enc_type1)
        dico_monomers2 = dico_of_dico_frustIndex_with_chains(results_frustration_dir2, enc_type2)

        #8.
        dico_mean1 = dico_mean_frustration(dico_monomers1)
        dico_mean2 = dico_mean_frustration(dico_monomers2)

        # graphical view, scatter plot
        colors = plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_number1, enc_type2, enc_number2, plots_dir, seqdist_flag, mode_flag)
        # print(colors)
        green_list = []
        red_list = []
        grey_list = []
        yellow_list = []
        goldenrod_list = []
        for k in range(len(colors)):
            if colors[k] == "green":
                green_list.append(k + 1)
            elif colors[k] == "red":
                red_list.append(k + 1)
            elif colors[k] == "grey":
                grey_list.append(k + 1)
            elif colors[k] == "yellow":
                yellow_list.append(k + 1)
            elif colors[k] == "goldenrod":
                goldenrod_list.append(k + 1)
        green_num = ""
        for num in green_list:
            green_num += str(num) + " "
        # print( green_num)

        red_num = ""
        for num in red_list:
            red_num += str(num) + " "
        # print( red_num)

        grey_num = ""
        for num in grey_list:
            grey_num += str(num) + " "
        # print( grey_num)

        yellow_num = ""
        for num in yellow_list:
            yellow_num += str(num) + " "
        # print( yellow_num)

        goldenrod_num = ""
        for num in goldenrod_list:
            goldenrod_num += str(num) + " "
        # print( goldenrod_num)

        # print(green_list)
        # print(red_list)
        # print(grey_list)
        # print(yellow_list)
        # print(goldenrod_list)
        # coloring structure
        # color_structure_by_index(colors,  results_pdb_dir1, results_pdb_dir2)

        # plot visualisation
        #visualization_VMD_script(results_pdb_dir1, results_pdb_dir2, enc_type1, enc_type2, enc_number1, enc_number2,
        #                         green_num, red_num, grey_num, yellow_num, goldenrod_num, vmd_flag)

        # % of frustration type
        dico_types_1 = dico_percentage_frustration_types(dico_monomers1, mode_flag)
        # print(dico_types_1)
        dico_types_2 = dico_percentage_frustration_types(dico_monomers2, mode_flag)
        # print(dico_types_2)

        # % of each type of frustration barplot

        dico_mean_types_1 = dico_mean_percentage_frustration_types(dico_types_1)
        print(dico_mean_types_1)
        dico_mean_types_2 = dico_mean_percentage_frustration_types(dico_types_2)
        print(dico_mean_types_2)

        plot_barplot_percentage_frustration_types(plots_dir, dico_mean_types_1, dico_mean_types_2, enc_type1,
                                                  enc_number1, enc_type2, enc_number2, seqdist_flag)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

if __name__ == "__main__":
    pdb_files, vmd_flag, frustration_flag, seqdist_flag, isolate_flag, plots_flag, mode_flag = parse_arguments()

    if len(pdb_files) == 1:
        try:
            main(pdb_files[0], vmd_flag=vmd_flag, frustration_flag=frustration_flag, seqdist_flag= seqdist_flag, isolate_flag= isolate_flag, plots_flag= plots_flag, mode_flag=mode_flag)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)
    elif len(pdb_files) == 2:
        try:
            main2(pdb_files[0], pdb_files[1], vmd_flag=vmd_flag, frustration_flag=frustration_flag, seqdist_flag=seqdist_flag, isolate_flag= isolate_flag, mode_flag=mode_flag)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)
    else:
        print("Usage: python3 structural_align.py path/to/encapsulin.pdb [-vmd=True] [-frustration=True]")
        print("Or: python3 structural_align.py path/to/encapsulin1.pdb path/to/encapsulin2.pdb [-vmd=True] [-frustration=True]")
        sys.exit(1)