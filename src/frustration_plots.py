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


'''
TO DO : 
make the frustration calculation in order for the files


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
    os.makedirs(plots_dir, exist_ok=True)

    # 2. Charger les monomères
    monomers = load_monomers(pdb_file1)
    print(f"Loaded {len(monomers)} monomers")

     # 3. create PDB files
    save_each_monomer_as_pdb(monomers, results_pdb_dir, enc_type, enc_number)

    # 4. calculation of frustration
    #calculate_frustration(results_pdb_dir, results_frustration_dir)

    #6
    dico_monomers = dico_of_dico_frustIndex(results_frustration_dir)
    #print(dico_monomers)


    #7 grafical view
    #monomers_id = list(dico_monomers.keys())
    #residues = list(dico_monomers[1].keys())

    plt.figure(figsize=(25, 8))

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
            markersize=2,
            linestyle='-',
            linewidth=0.01,
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

    name_plot = f"frustration_per_res_{enc_type}_monomer_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Plot with standard deviation saved as {plot_path}")


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
