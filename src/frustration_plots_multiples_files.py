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
import matplotlib.patches as mpatches



'''
This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a wanted monomer, and make a frustration study of the monomers in each file (frames) using FrustratometeR 
the results files are added to results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_frames_for_a_monomer 
A plot is add to the plots/frustration directory, with name of type rmsf_with_std_per_frame_<TmEnc|MmEnc>_monomer_<(i)>.png
Usage:

python3 frustration_plots_multiples_files.py path/to/directory/ number_of_the_wanted_monomer

If 2 file given, frustration will be calculated for the monomers of the 2 files and a scatter-plot of the frsutration means for each residue will be made

python3 frustration_plots_multiples_files.py path/to/file1.pdb path/to/file2.pdb 

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
def calculate_frustration(PDB_directory, output_directory,seqdist_flag):

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
    results <- calculate_frustration(PdbFile = pdb_file, Mode = "singleresidue", SeqDist = {seqdist}, ResultsDir = output_dir, Graphics = FALSE)

    """

    # Pour chaque fichier PDB dans le répertoire
    # TO DO make the frustration calculation in order for the files
    for pdb_file in os.listdir(PDB_directory):
        if pdb_file.endswith('.pdb'):
            full_path = os.path.join(PDB_directory, pdb_file)

            # Personnaliser le script R pour ce fichier
            current_script = r_script.format(
                pdb_file=full_path,
                output_dir=output_directory,
                seqdist=seqdist_flag
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

def dico_percentage_frustration_types(dico_frames):
    """
    Calculate the percentage of each frustration type for all frames.

    Frustration types are defined as:
    - minimal: frustration > 0.78
    - high: frustration < -1
    - neutral: -1 ≤ frustration ≤ 0.78

    Args:
        dico_monomers: Dictionary with format {frame_num: {'res1': frustration, ...}, ...}

    Returns:
        Dictionary with format {monomer_num: {'high': %, 'minimal': %, 'neutral': %}, ...}
    """
    dico_percentage = {}

    for monomer_num, residue_data in dico_frames.items():
        # Initialize counters
        counts = {'high': 0, 'minimal': 0, 'neutral': 0}
        total_residues = len(residue_data)

        if total_residues == 0:
            dico_percentage[monomer_num] = {'high': 0, 'minimal': 0, 'neutral': 0}
            continue

        # Classify each residue
        for frustration in residue_data.values():
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



def plot_frustration_per_res(dico, enc_type, monomer_number, plots_dir, file_dic, seqdist_flag):
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
    plt.title(f'Mean Frustration seqdist {seqdist_flag} per Residue ({enc_type}, Monomer {monomer_number}) for all the {len(file_dic)} frames', fontsize=10)
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

def plot_min_max_frustration_bars(dico_monomers, enc_type, monomer_number, plots_dir, file_dic, seqdist_flag):
    """
    Generate a plot showing vertical bars from min to max frustration values per residue.

    :param dico_monomers: Dictionary with frustration info {monomer_num: {'M4': 0.759, ...}, ...}
    :param enc_type: Encapsuline type identifier
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
    for i, (res, vmin, vmax) in enumerate(zip(residues, min_vals, max_vals)):
         # Determine color based on values
        if vmax <= -1:
            color = 'red'  # Both min and max are negative
        elif vmin >= 0.78:
            color = 'green'  # Both min and max are positive
        elif vmax < 0.78 and vmin > -1 :
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
    plt.title(f'Frustration seqdist {seqdist_flag} Range per Residue ({enc_type}, Monomer {monomer_number}) for all the {len(file_dic)} frames')
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
    green_patch = mpatches.Patch(color='green', label='Minimal frustration', alpha=0.7)
    red_patch = mpatches.Patch(color='red', label='High frustration', alpha=0.7)
    grey_patch = mpatches.Patch(color='grey', label='neutral frustration', alpha=0.7)
    gold_patch = mpatches.Patch(color='gold', label='Mixed frustration', alpha=0.7)
    plt.legend(handles=[green_patch, grey_patch, red_patch, gold_patch])

    plt.grid(True, linestyle=':', alpha=0.4)
    plt.tight_layout()

    # Save plot
    name_plot = f"frustration_range_per_frame_{enc_type}_monomer_{monomer_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_type2, monomer_number, plots_dir, file_dict, seqdist):
    """
    Create a scatter plot comparing mean frustration values between two conditions.
    Points are colored:
    - Green if both x and y are positive
    - Red if both x and y are negative
    - Gold otherwise

    Parameters:
    - dico_mean1: {'M4': {'mean': -0.135, 'std': 0.668}, ...} for first condition
    - dico_mean2: Same structure for second condition
    - enc_type1, enc_type2: Identifier for conditions
    - monomer_number: Monomer number being analyzed
    - plots_dir: Directory to save the plot
    - file_dict: Dictionary containing frame information to display count in title
    """
    # Prepare data
    x_vals = [dico_mean2[res]['mean'] for res in dico_mean2.keys()]
    y_vals = [dico_mean1[res]['mean'] for res in dico_mean1.keys()]

    # Determine colors
    colors = []
    for x, y in zip(x_vals, y_vals):
        if x >= 0.78 and y >= 0.78:
            colors.append('green')
        elif x <= -1 and y <= -1:
            colors.append('red')
        elif -1<x<0.78  and -1<y<0.78:
            colors.append('grey')
        elif x>y:
            colors.append('yellow')
        else :
            colors.append('goldenrod')
    # Create figure
    plt.figure(figsize=(8, 8))

    # Create scatter plot with colored points
    plt.scatter(x_vals, y_vals, c=colors, alpha=0.7, s=60, edgecolors='black', linewidths=0.5)

    # Add reference lines
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)

    # Add identity line
    max_val = max(max(np.abs(x_vals)), max(np.abs(y_vals))) * 1.1
    plt.plot([-max_val, max_val], [-max_val, max_val], 'k--', alpha=0.3)

    # Customize plot with frame count in title
    plt.title(
        f'Mean monomers frustration seqdist {seqdist} comparison\n{enc_type1} vs {enc_type2} for monomer {monomer_number}\n({len(file_dict)} frames analyzed)',
        pad=20
    )
    plt.xlabel(f'Mean Frustration ({enc_type2} monomer {monomer_number})')
    plt.ylabel(f'Mean Frustration ({enc_type1} monomer {monomer_number})')
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
    name_plot = f"scatter_frustration_per_frame_{enc_type1}_monomer_{monomer_number}_vs_{enc_type2}_monomer_{monomer_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return colors

def visualization_VMD_script(output_dir1, output_dir2, enc_type1, enc_type2, monomer_number, green_num,
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
mol new {output_dir}/{enc_type}0_monomer{monomer_number}.pdb type pdb waitfor all

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
    output_file1 = os.path.join(output_dir1, "frustration_per_frame_colors.tcl")
    output_file2 = os.path.join(output_dir2, "frustration_per_frame_colors.tcl")

    # Generate and write script for output_file1
    script_content1 = generate_script_content(output_dir1, enc_type1, monomer_number)
    with open(output_file1, 'w') as f:
        f.write(script_content1)

    # Generate and write script for output_file2 (with different parameters)
    script_content2 = generate_script_content(output_dir2, enc_type2, monomer_number)
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


def plot_barplot_percentage_frustration_types(plots_dir,dico_mean_type_1 , dico_mean_type_2, enc_type1, enc_type2, monomer_number, seqdist, file_dict):
        """
        Crée un barplot comparant les types de frustration entre MtEnc et TmEnc
        avec barres d'erreur pour les écarts-types.
        :param plots_dir: dossier où sauvegarder le graphique
        :param dico_mean_type_1: dictionnaire des stats pour MtEnc (format {'minimal': {'mean': x, 'std': y}, ...})
        :param dico_mean_type_2: dictionnaire des stats pour TmEnc (même format)
        """

        # Organized data
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
        plt.title(f'Comparison of mean Frustration (secdict {seqdist}) type between MtEnc and TmEnc monomer {monomer_number} for {len(file_dict)} frames ', fontsize=14)
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
        name_plot = f"frustration_barplot_per_frame_{enc_type1}_monomer_{monomer_number}_vs_{enc_type2}_monomer_{monomer_number}.png"

        plot_path = os.path.join(plots_dir, name_plot)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()

def parse_arguments():
    """Parse command line arguments including optional -vmd and -frustration flags"""
    vmd_flag = False
    frustration_flag = False
    seqdist_flag = 12
    pdb_files = []

    monomer_number = None

    for arg in sys.argv[1:]:
        if arg.startswith('-vmd='):
            vmd_value = arg.split('=')[1].lower()
            if vmd_value == 'true':
                vmd_flag = True
        elif arg.startswith('-frustration='):
            frustration_value = arg.split('=')[1].lower()
            if frustration_value == 'true':
                frustration_flag = True
        elif arg.isdigit():
            monomer_number = arg
        elif arg.startswith('-seqdist='):
            seqdist_flag = arg.split('=')[1]
        else:
            pdb_files.append(arg)

    if not monomer_number:
        raise ValueError("Monomer number must be specified as a numeric argument")

    return pdb_files, monomer_number, vmd_flag, frustration_flag, seqdist_flag

def main(pdb_directory, monomer_number, vmd_flag=False, frustration_flag=False, seqdist_flag = 12 ):
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
    results_frustration_dir = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir, capsids_dir,frames_dir, f"{enc_type}_monomer_{monomer_number}_frustration_seqdist_{seqdist_flag}")

    plots_dir = os.path.join("../plots", enc_type, f"frustration_seqdist_{seqdist_flag}")

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
    if frustration_flag:
        calculate_frustration(results_pdb_dir, results_frustration_dir, seqdist_flag)


    #6
    dico_monomers = dico_of_dico_frustIndex(results_frustration_dir)
    #print(dico_monomers)


    #7 grafical views

    #8.
    dico_mean = dico_mean_frustration(dico_monomers)
    #print(dico_mean)

    #9
    plot_frustration_per_res(dico_mean, enc_type, monomer_number, plots_dir, files_dict, seqdist_flag)

    plot_min_max_frustration_bars(dico_monomers, enc_type, monomer_number, plots_dir, files_dict, seqdist_flag)

    # % of frustration type
    dico_types = dico_percentage_frustration_types(dico_monomers)
    #print(dico_types)

    dico_mean_types = dico_mean_percentage_frustration_types(dico_types)
    print(dico_mean_types)


    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

#when 2 files given in parameters
def main2(pdb_directory1, pdb_directory2, monomer_number,vmd_flag=False, frustration_flag=False, seqdist_flag = 12 ):

    start_time = time.time()

    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the directories
    enc_type1, files_dict1, chain_dict1 = extract_info_from_directory(pdb_directory1)
    print(f"Fist encapsulin type : {enc_type1}")
    print(f"Number of Monomers by structure : {len(chain_dict1)}")
    print(f"Number of frames : {len(files_dict1)} ")

    enc_type2, files_dict2, chain_dict2 = extract_info_from_directory(pdb_directory2)
    print(f"Secon encapsulin type : {enc_type2}")
    print(f"Number of Monomers by structure : {len(chain_dict2)}")
    print(f"Number of frames : {len(files_dict2)} ")


    #1.1. create the output directories paths

    frustration_dir1 = f"FRUSTRATION_{enc_type1.upper()}"
    capsids_dir1 = f"{enc_type1.upper()}_CAPSIDS"
    frames_dir = "FRUSTRATION_frames_for_a_monomer"
    results_pdb_dir1 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir1, capsids_dir1,frames_dir, f"{enc_type1}_monomer_{monomer_number}_monomers")
    results_frustration_dir1 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir1, capsids_dir1,frames_dir, f"{enc_type1}_monomer_{monomer_number}_frustration_seqdist_{seqdist_flag}")

    frustration_dir2 = f"FRUSTRATION_{enc_type2.upper()}"
    capsids_dir2 = f"{enc_type2.upper()}_CAPSIDS"
    frames_dir = "FRUSTRATION_frames_for_a_monomer"
    results_pdb_dir2 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir2, capsids_dir2,frames_dir, f"{enc_type2}_monomer_{monomer_number}_monomers")
    results_frustration_dir2 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir2, capsids_dir2,frames_dir, f"{enc_type2}_monomer_{monomer_number}_frustration_seqdist_{seqdist_flag}")


    plots_dir = os.path.join("../plots", "COMMON", f"frustration_seqdist_{seqdist_flag}")

    #1.2 create the repository if it not exist
    os.makedirs(results_pdb_dir1, exist_ok=True)
    os.makedirs(results_frustration_dir1, exist_ok=True)
    os.makedirs(results_pdb_dir2, exist_ok=True)
    os.makedirs(results_frustration_dir2, exist_ok=True)

    os.makedirs(plots_dir, exist_ok=True)

    # 2. Charger les monomères
    monomers1 = load_monomers(files_dict1, chain_dict1, monomer_number)
    print(f"Loaded {len(monomers1)} monomers")

    monomers2 = load_monomers(files_dict2, chain_dict2, monomer_number)
    print(f"Loaded {len(monomers2)} monomers")

     # 3. create PDB files
    save_each_monomer_as_pdb(monomers1, results_pdb_dir1, enc_type1, files_dict1, monomer_number)
    save_each_monomer_as_pdb(monomers2, results_pdb_dir2, enc_type2, files_dict2, monomer_number)

    # 4. calculation of frustration
    if frustration_flag:
        calculate_frustration(results_pdb_dir1, results_frustration_dir1, seqdist_flag )
        calculate_frustration(results_pdb_dir2, results_frustration_dir2, seqdist_flag )


    #6
    dico_monomers1 = dico_of_dico_frustIndex(results_frustration_dir1)
    dico_monomers2 = dico_of_dico_frustIndex(results_frustration_dir2)

    #8.
    dico_mean1 = dico_mean_frustration(dico_monomers1)
    dico_mean2 = dico_mean_frustration(dico_monomers2)
    #print(dico_mean1)
    #print(dico_mean2)

   # graphical view, scatter plot
    colors = plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_type2, monomer_number, plots_dir, files_dict1, seqdist_flag)
    print(colors)
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
    #print( print_num)

    red_num = ""
    for num in red_list :
        red_num+=str(num)+" "
    #print( print_num)

    grey_num = ""
    for num in grey_list :
        grey_num+=str(num)+" "
    #print( print_num)

    yellow_num = ""
    for num in yellow_list :
        yellow_num+=str(num)+" "
    #print( print_num)

    goldenrod_num = ""
    for num in goldenrod_list :
        goldenrod_num+=str(num)+" "
    #print( print_num)

    #print(green_list)
    #print(red_list)
    #print(grey_list)
    #print(yellow_list)
    #print(goldenrod_list)

    # visualise mean frustration for each residue
    visualization_VMD_script(results_pdb_dir1, results_pdb_dir2, enc_type1, enc_type2 ,monomer_number,green_num, red_num, grey_num, yellow_num, goldenrod_num, vmd_flag)

    # % of frustration type
    dico_types_1 = dico_percentage_frustration_types(dico_monomers1)
    #print(dico_types_1)
    dico_types_2 = dico_percentage_frustration_types(dico_monomers2)
    #print(dico_types_2)

    # % of each type of frustration barplot

    dico_mean_types_1 = dico_mean_percentage_frustration_types(dico_types_1)
    print(dico_mean_types_1)
    dico_mean_types_2 = dico_mean_percentage_frustration_types(dico_types_2)
    print(dico_mean_types_2)

    plot_barplot_percentage_frustration_types(plots_dir, dico_mean_types_1, dico_mean_types_2, enc_type1, enc_type2, monomer_number, seqdist_flag, files_dict1 )

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")


if __name__ == "__main__":
    try:
        pdb_files, monomer_number, vmd_flag, frustration_flag, seqdist_flag = parse_arguments()

        if len(pdb_files) == 1:
            try:
                main(pdb_files[0], monomer_number, vmd_flag=vmd_flag, frustration_flag=frustration_flag,seqdist_flag= seqdist_flag )
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)
        elif len(pdb_files) == 2:
            try:
                main2(pdb_files[0], pdb_files[1], monomer_number, vmd_flag=vmd_flag, frustration_flag=frustration_flag,seqdist_flag= seqdist_flag )
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)
        else:
            print(
                "Usage: python3 frustration_plots_multiples_files.py path/to/directory/ number_of_wanted_monomer [-vmd=True] [-frustration=True]")
            print(
                "Or: python3 frustration_plots_multiples_files.py path/to/directory1 path/to/directory2 number_of_wanted_monomer [-vmd=True] [-frustration=True]")
            sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        print(
            "Usage: python3 frustration_plots_multiples_files.py path/to/directory/ number_of_wanted_monomer [-vmd=True] [-frustration=True]")
        print(
            "Or: python3 frustration_plots_multiples_files.py path/to/directory1 path/to/directory2 number_of_wanted_monomer [-vmd=True] [-frustration=True]")
        sys.exit(1)