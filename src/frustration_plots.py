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

#from FrustraMotion.src.parser_pdb import return_sequence_3_letter_format

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

If 2 file given, frustration will be calculated for the monomers of the 2 files and a scatter-plot of the frsutration means for each residue will be made

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
#6
def dico_of_dico_frustIndex(frustration_directory):
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
                dico_monomer = dico_frustIndex(full_path)
                dico_monomers[number]= dico_monomer

    return dico_monomers

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

def plot_frustration_per_res(dico, enc_type, enc_number, plots_dir):
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
    plt.title(f'Mean Frustration per Residue ({enc_type}, frame {enc_number})', fontsize=10)
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




def plot_min_max_frustration_bars(dico_monomers, enc_type, enc_number, plots_dir):
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
    for i, (res, vmin, vmax) in enumerate(zip(residues, min_vals, max_vals)):
        # Determine color based on values
        if vmax <= 0:
            color = 'green'  # Both min and max are negative
        elif vmin >= 0:
            color = 'red'  # Both min and max are positive
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
    plt.title(f'Frustration Range per Residue ({enc_type}, frame {enc_number})')
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
    green_patch = mpatches.Patch(color='green', label='Negative frustration', alpha=0.7)
    red_patch = mpatches.Patch(color='red', label='Positive frustration', alpha=0.7)
    gold_patch = mpatches.Patch(color='gold', label='Mixed frustration', alpha=0.7)
    plt.legend(handles=[green_patch, red_patch, gold_patch])

    plt.grid(True, linestyle=':', alpha=0.4)
    plt.tight_layout()

    # Save plot
    name_plot = f"frustration_range_per_res_{enc_type}_frame_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()




def plot_frustration_boxplots(dico_list, enc_type, enc_number, plots_dir):
    """
    Generate boxplots of frustration values per residue from a dictionary of lists.
    Shows only boxplots without individual data points.

    Parameters:
    - dico_list: Dictionary with residue IDs as keys and lists of frustration values
                Example: {'M4': [0.51, 0.67, ...], ...}
    - enc_type: Encounter type identifier
    - enc_number: Encounter number
    - plots_dir: Directory to save the plot
    """
    # Convert dictionary to DataFrame for easier plotting
    df = pd.DataFrame.from_dict(dico_list, orient='index').T

    # Melt the DataFrame for seaborn compatibility
    df_melted = df.melt(var_name='Residue', value_name='Frustration')

    # Create figure
    plt.figure(figsize=(14, 7))

    # Create boxplot with Seaborn (showing only boxes)
    sns.boxplot(
        data=df_melted,
        x='Residue',
        y='Frustration',
        #palette='vlag',  # Red-white-blue diverging palette
        showfliers=True,  # Hide outliers (since we don't want individual points)
        width=0.6,
        linewidth=1.5,
        fliersize=0  # Ensure no outlier markers are shown
    )

    # Reference line at y=0
    plt.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    # Customize plot
    plt.title(f'Frustration Distribution per Residue ({enc_type}, Monomer {enc_number})', pad=20, fontsize=14)
    plt.xlabel('Residue', labelpad=10)
    plt.ylabel('Frustration Index', labelpad=10)

    # Rotate x-axis labels and adjust font
    plt.xticks(
        rotation=45,
        ha='right',
        fontsize=10
    )

    # Improve y-axis
    plt.yticks(fontsize=10)
    plt.grid(axis='y', alpha=0.3)

    # Adjust layout
    plt.tight_layout()

    # Save plot
    name_plot = f"frustration_boxplot_per_res_{enc_type}_monomer_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_number1, enc_type2, enc_number2, plots_dir):
    """
    Create a scatter plot comparing mean frustration values between two conditions.

    Parameters:
    - dico_mean1: {'M4': {'mean': -0.135, 'std': 0.668}, ...} for first condition
    - dico_mean2: Same structure for second condition
    - enc_type1, enc_number1: Identifier for first condition
    - enc_type2, enc_number2: Identifier for second condition
    - plots_dir: Directory to save the plot
    """
    # Extract common residues
    #common_residues = sorted(set(dico_mean1.keys()) & set(dico_mean2.keys()))
    #print(dico_mean1.keys())
    #print(dico_mean2.keys())
    #print(common_residues)
    # Prepare data
    x_vals = [dico_mean2[res]['mean'] for res in dico_mean2.keys()]
    y_vals = [dico_mean1[res]['mean'] for res in dico_mean1.keys()]
    #print(x_vals)
    #print(y_vals)

    # Create figure
    plt.figure(figsize=(8, 8))

    # Create scatter plot
    plt.scatter(x_vals, y_vals, color='blue', alpha=0.6, s=50)

    # Add reference lines
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)

    # Add identity line
    max_val = max(max(np.abs(x_vals)), max(np.abs(y_vals))) * 1.1
    plt.plot([-max_val, max_val], [-max_val, max_val], 'r--', alpha=0.5)

    # Customize plot
    plt.title(f'Mean monomers frustration comparison\n{enc_type1} frame {enc_number1} vs {enc_type2} frame {enc_number2}')
    plt.xlabel(f'Mean Frustration ({enc_type2} {enc_number2})')
    plt.ylabel(f'Mean Frustration ({enc_type1} {enc_number1})')
    plt.grid(True, linestyle='--', alpha=0.3)

    # Equal aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')

    # Save plot
    name_plot = f"scatter_frustration_per_res_{enc_type1}_frame_{enc_number1}_vs_{enc_type2}_frame_{enc_number2}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()



#when only a file given
def main(pdb_file1):
    start_time = time.time()
    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filename
    enc_type, enc_number = extract_info_from_filename(pdb_file1)

    #1.1. create the output directory path
    frustration_dir = f"FRUSTRATION_{enc_type.upper()}"
    capsids_dir = f"{enc_type.upper()}_CAPSIDS"
    monomer_dir = "FRUSTRATION_monomer_for_a_frame"
    results_pdb_dir = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir, capsids_dir,monomer_dir, f"{enc_type}{enc_number}_monomers")
    results_frustration_dir = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir, capsids_dir,monomer_dir, f"{enc_type}{enc_number}_frustration")

    plots_dir = os.path.join("../plots", enc_type, "frustration")

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
    calculate_frustration(results_pdb_dir, results_frustration_dir)

    #6
    dico_monomers = dico_of_dico_frustIndex(results_frustration_dir)
    #print(dico_monomers)


    #7 grafical views

    plot_all_frustration_per_res(dico_monomers, enc_type, enc_number, plots_dir)

    #8.
    dico_mean = dico_mean_frustration(dico_monomers)
    #print(dico_mean)

    #9
    plot_frustration_per_res(dico_mean, enc_type, enc_number, plots_dir)

    #10
    #plot_min_and_max_frustration_per_res(dico_monomers, enc_type, enc_number, plots_dir)
    plot_min_max_frustration_bars(dico_monomers, enc_type, enc_number, plots_dir)

    dico_list = dico_list_frustration(dico_monomers)
    plot_frustration_boxplots(dico_list, enc_type, enc_number, plots_dir)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")


#when 2 files given in parameters
def main2(pdb_file1, pdb_file2):

    start_time = time.time()

    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filenames
    enc_type1, enc_number1 = extract_info_from_filename(pdb_file1)
    enc_type2, enc_number2 = extract_info_from_filename(pdb_file2)

    #1.1. create the output directories paths

    frustration_dir1 = f"FRUSTRATION_{enc_type1.upper()}"
    capsids_dir1 = f"{enc_type1.upper()}_CAPSIDS"
    monomer_dir = "FRUSTRATION_monomer_for_a_frame"
    results_pdb_dir1 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir1, capsids_dir1,monomer_dir, f"{enc_type1}{enc_number1}_monomers")
    results_frustration_dir1 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir1, capsids_dir1,monomer_dir, f"{enc_type1}{enc_number1}_frustration")

    frustration_dir2 = f"FRUSTRATION_{enc_type2.upper()}"
    capsids_dir2 = f"{enc_type2.upper()}_CAPSIDS"
    monomer_dir = "FRUSTRATION_monomer_for_a_frame"
    results_pdb_dir2 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir2, capsids_dir2,monomer_dir, f"{enc_type2}{enc_number2}_monomers")
    results_frustration_dir2 = os.path.join("/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results", frustration_dir2, capsids_dir2,monomer_dir, f"{enc_type2}{enc_number2}_frustration")

    plots_dir = os.path.join("../plots", "COMMON", "frustration")

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
    calculate_frustration(results_pdb_dir1, results_frustration_dir1)
    calculate_frustration(results_pdb_dir2, results_frustration_dir2)

    #6
    dico_monomers1 = dico_of_dico_frustIndex(results_frustration_dir1)
    dico_monomers2 = dico_of_dico_frustIndex(results_frustration_dir2)

    #8.
    dico_mean1 = dico_mean_frustration(dico_monomers1)
    dico_mean2 = dico_mean_frustration(dico_monomers2)


   # graphical view, scatter plot
    plot_scatter_frustration_mean(dico_mean1, dico_mean2, enc_type1, enc_number1, enc_type2, enc_number2, plots_dir)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

if __name__ == "__main__":
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("Usage: python3 frustration_plots.py path/to/encapsulin.pdb \n or python3 frustration_plots.py path/to/encapsulin1.pdb path/to/encapsulin2.pdb")
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
    elif len(sys.argv) ==3:
        try:
            main2(sys.argv[1], sys.argv[2])
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)