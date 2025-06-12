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
from collections import OrderedDict



'''
This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a all the data of all the residues and monomers and make a frustration study of the residues of the monomers in each file (frames) using FrustratometeR results 
the results files are added to results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_frames_for_a_residue 
A plot is add to the plots/<MTENC|TMENC>/frustration_MN frustration directory, with name of type frustration_monomer_1.png
Usage:

python3 frustration_plots_MN_data.py path/to/directory/ number_of_the_wanted_monomer


'''
# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)

#parse every directory name contain in the given directory and make a list of absolut paths of all the directories having as first digit the given digit in parameters
def list_of_wanted_monomer(frustration_directory, monomer_number):
    """
    parse every directory name contain in the given directory and make a list of absolut paths of all the directories having as first digit the given digit in parameters

    :param frustration_directory: directory where alle the frsutration directory results are
    :param monomer_number: the monomer wanted
    :return: a list of the absolute path to the file containing the results of frustration
    """
    #dico of number to letter equivalency in the monomer names
    monomer_id = {'0':1, '1':2, '2':3, '3':4, '4':5, '5':6, '6':7, '7':8, 'A':9, 'B':10,
                  'C':11,'D':12,'E':13,'F':14,'G':15,'H':16,'I':17,'J':18,'K':19,'L':20,
                  'M':21,'N':22,'O':23,'P':24,'Q':25,'R':26,'S':27,'T':28,'U':29,'V':30,
                  'W':31,'X':32,'Y':33,'Z':34,'a':35,'b':36,'c':37,'d':38,'e':39,'f':40,
                  'g':41,'h':42,'i':43,'j':44,'k':45,'l':46,'m':47,'n':48,'o':49,'p':50,
                  'q':51,'r':52,'s':53,'t':54,'u':55,'v':56,'w':57,'x':58,'y':59,'z':60 }

    wanted_dirs = []
    # Ensure the frustration_directory exists
    if not os.path.isdir(frustration_directory):
        print(f"Directory not found: {frustration_directory}")
        return wanted_dirs

    for entry in os.listdir(frustration_directory):
        full_path = os.path.join(frustration_directory, entry)
        if os.path.isdir(full_path):
            # Split the directory name to extract the monomer character
            parts = entry.split('_')
            if len(parts) < 3:
                continue
            monomer_char = parts[-1].split('.')[0]
            #print(f'la valeur dans le nom du directory est {monomer_char}')
            # Check if monomer_char is in the dictionary
            if monomer_char in monomer_id:
                #print(f'{monomer_char} est dans monomer id')
                # Convert dictionary value to 0-based index
                current_monomer_index = monomer_id[monomer_char]
                #print(f'sa valeur est {current_monomer_index}')
                if str(current_monomer_index) == monomer_number:
                    #print(f'{current_monomer_index} est pareil que {monomer_number}, on ajoute le path a la liste')
                    full_path = os.path.join(full_path, "FrustrationData")
                    for file in os.listdir(full_path):
                        if file.endswith("singleresidue"):
                            full_path= os.path.join(full_path, file)
                            wanted_dirs.append(full_path)

    return wanted_dirs

#take the list of result files, and for a given residue, take the frustration value to form a dico
def dico_frustration_of_a_residue_by_frame(results_files, residue_flag):
    """
    take one result file, and for a given residue, take the frsutration value to form a dico
    :param results_files:
    :param residue_flag:
    :return:  residue_key , dico_frustration_by_frame = {(frame)0:1.5(frustration), 10:2.5 , 20:2.4 , ....}
    """

    frustration_dict = {}
    for file in results_files:
        with open(file, 'r') as f:
            # Skip header line if present
            header = f.readline()
            for line in f:
                # Skip empty lines
                if not line.strip():
                    continue

                parts = line.split()

                # Check we have enough columns (at least 8)
                if len(parts) >= 8:
                    if os.path.basename(file).startswith('TmEnc'):
                        res_num = int(parts[0]) - 3
                        type = "TmEnc"
                    else :
                        res_num = int(parts[0])
                        type = "MtEnc"
                    if res_num == int(residue_flag) :
                        aa = parts[3]  # Amino acid (M, E, F...)
                        ##remplace the '?' with the correct AA f the data
                        if aa == '?':
                            aa = 'H'
                        frst_index = float(parts[7])  # FrstIndex value
                        residue_key = f"{aa}{res_num}"
                        frame = os.path.basename(file).split("_")[1]
                        frustration_dict[int(frame)] = frst_index

    return type,residue_key, frustration_dict


def plot_frustration(dico_frustration, type, residue, monomer):
    """
    Génère un plot des valeurs de frustration par frame avec le même code couleur que plot_all_frustrations.
    Points rouges si < -1, verts si > -0.58, gris sinon.
    Lignes horizontales à y=0, y=-1 (rouge) et y=-0.58 (vert).

    :param dico_frustration: OrderedDict contenant les paires (frame, frustration_value)
    :param type: Type de structure (MTENC/TMENC)
    :param residue: Numéro du résidu étudié
    :param monomer: Numéro du monomère
    """
    # Extraire les frames et les valeurs de frustration
    frames = list(dico_frustration.keys())
    frustrations = list(dico_frustration.values())

    # Créer le plot
    plt.figure(figsize=(12, 6))

    # Tracer les points avec le code couleur
    for x, y in zip(frames, frustrations):
        if y < -1:
            color = 'r'  # Rouge
        elif y > -0.58:
            color = 'g'  # Vert
        else:
            color = 'gray'  # Gris
        plt.plot(x, y, 'o', color=color, markersize=4)

    # Connecter les points avec une ligne bleue fine
    plt.plot(frames, frustrations, linestyle='-', color='b', alpha=0.3, linewidth=0.5)

    # Lignes horizontales
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)  # Ligne noire en pointillés à y=0
    plt.axhline(y=-1, color='r', linestyle=':', alpha=0.5)  # Ligne rouge en pointillés à y=-1
    plt.axhline(y=-0.58, color='g', linestyle=':', alpha=0.5)  # Ligne verte en pointillés à y=-0.58

    # Inverser l'axe des ordonnées
    plt.gca().invert_yaxis()

    # Labels et titre
    plt.xlabel('Frame')
    plt.ylabel('Frustration Index (inverted)')
    plt.title(f'Frustration Index by Frame of residue {residue} ({type}), monomer {monomer}')

    # Gestion des ticks
    if len(frames) > 20:
        step = max(1, len(frames) // 10)  # Évite step=0 pour les petites listes
        plt.xticks(frames[::step], rotation=45)
    else:
        plt.xticks(rotation=45)

    # Ajouter une légende pour les couleurs
    red_patch = mpatches.Patch(color='red', label='Highly frustrated (< -1)')
    green_patch = mpatches.Patch(color='green', label='Minimally frustrated (> -0.58)')
    gray_patch = mpatches.Patch(color='gray', label='Neutral frustration')
    plt.legend(handles=[red_patch, green_patch, gray_patch], loc='best')

    # Grille et affichage
    plt.grid(True, linestyle=':', alpha=0.3)
    plt.tight_layout()
    plt.show()

import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy as np


def plot_all_frustrations(all_data, type, residue):
    """Crée une figure avec 60 subplots (12x5) selon l'ordre spécifié"""
    # Création de la figure globale
    fig, axs = plt.subplots(12, 5, figsize=(24, 30))
    fig.suptitle(f'Frustration Index by Frame for residue {residue} ({type})', y=1.02, fontsize=16)

    # Aplatir le tableau d'axes
    axs = axs.flatten()

    # Ordre spécifique des monomères (12 lignes x 5 colonnes)
    custom_order = [
        '0', '5', 'l', 'q', 'v',
        'J', 'F', 'G', 'H', 'I',
        '1', '6', 'm', 'r', 'w',
        '2', '7', 'n', 's', 'x',
        '3', 'j', 'o', 't', 'y',
        '4', 'k', 'p', 'u', 'z',
        'A', 'B', 'C', 'D', 'E',
        'K', 'P', 'U', 'Z', 'e',
        'L', 'Q', 'V', 'a', 'f',
        'M', 'R', 'W', 'b', 'g',
        'N', 'S', 'X', 'c', 'h',
        'O', 'T', 'Y', 'd', 'i'
    ]

    # Mapping inverse : label -> numéro de monomère
    monomer_labels = {1: '0', 2: '1', 3: '2', 4: '3', 5: '4', 6: '5', 7: '6', 8: '7',
                      9: 'A', 10: 'B', 11: 'C', 12: 'D', 13: 'E', 14: 'F', 15: 'G', 16: 'H',
                      17: 'I', 18: 'J', 19: 'K', 20: 'L', 21: 'M', 22: 'N', 23: 'O', 24: 'P',
                      25: 'Q', 26: 'R', 27: 'S', 28: 'T', 29: 'U', 30: 'V', 31: 'W', 32: 'X',
                      33: 'Y', 34: 'Z', 35: 'a', 36: 'b', 37: 'c', 38: 'd', 39: 'e', 40: 'f',
                      41: 'g', 42: 'h', 43: 'i', 44: 'j', 45: 'k', 46: 'l', 47: 'm', 48: 'n',
                      49: 'o', 50: 'p', 51: 'q', 52: 'r', 53: 's', 54: 't', 55: 'u', 56: 'v',
                      57: 'w', 58: 'x', 59: 'y', 60: 'z'}

    # Inverser le mapping pour avoir label->numéro
    label_to_num = {v: k for k, v in monomer_labels.items()}

    # Parcourir dans l'ordre personnalisé
    for i, label in enumerate(custom_order):
        ax = axs[i]
        monomer_num = label_to_num[label]

        if monomer_num in all_data:
            frames = list(all_data[monomer_num].keys())
            frustrations = list(all_data[monomer_num].values())

            try:
                frames_num = [int(f) for f in frames]
            except:
                frames_num = range(len(frames))

            # Plot des points avec couleurs conditionnelles
            for x, y in zip(frames_num, frustrations):
                if y < -1:
                    color = 'r'  # Rouge
                elif y > -0.58:
                    color = 'g'  # Vert
                else:
                    color = 'gray'  # Gris

                ax.plot(x, y, 'o', color=color, markersize=1.5)

            # Lignes horizontales
            ax.axhline(0, color='k', linestyle=':', alpha=0.3)
            ax.axhline(-1, color='r', linestyle=':', alpha=0.3, linewidth=0.5)
            ax.axhline(-0.58, color='g', linestyle=':', alpha=0.3, linewidth=0.5)

            ax.invert_yaxis()
            ax.set_title(f'M{monomer_num}', fontsize=8, pad=2)
        else:
            ax.axis('off')

        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.grid(True, linestyle=':', alpha=0.1)

        # Gestion des ticks sur l'axe X
        if monomer_num in all_data and len(frames_num) > 0:
            min_frame = min(frames_num)
            max_frame = max(frames_num)
            step = max(1, (max_frame - min_frame) // 4)
            ax.set_xticks([min_frame, min_frame + step,max_frame - step , max_frame])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    plt.show()


def parse_arguments():
    """Parse command line arguments with flexible options"""
    monomer_flag = None
    residue_flag = None
    pdb_file = []

    for arg in sys.argv[1:]:
        if arg.startswith('-monomer='):
            monomer_flag = arg.split('=')[1]
            try:
                int(monomer_flag)
            except:
                raise ValueError(
                    "The option -monomer expects an integer, use format: -monomer=value")
        elif arg.startswith('-residue='):
            residue_flag = arg.split('=')[1]
            try:
                int(residue_flag)
            except:
                raise ValueError(
                    "The option -residue expects an integer, use format: -residue=value")
        elif not arg.startswith('-'):
            pdb_file.append(arg)
        else:
            raise ValueError(f"Unknown option: {arg}")

    return pdb_file, monomer_flag, residue_flag


def main(frustration_dir, monomer_flag=None, residue_flag=None):
    # Create output directory if it doesn't exist
    output_dir = "../results/MN_data/SR/"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "TM_data_single_residue.txt")

    # Mode 1: No parameters - Generate complete data file
    if monomer_flag is None and residue_flag is None:
        with open(output_file, 'w') as f_out:
            f_out.write("residue\tmonomer\tframe\tfrustration\n")
            for residue_num in range(1, 265):
                for monomer in range(1, 61):
                    list_of_result_files = list_of_wanted_monomer(frustration_dir, str(monomer))
                    local_type, residue_key, dico_frst = dico_frustration_of_a_residue_by_frame(
                        list_of_result_files, residue_num)
                    dico_ordonne = OrderedDict(sorted(dico_frst.items()))

                    for frame, frustration in dico_ordonne.items():
                        f_out.write(f"{residue_num}\t{monomer}\t{frame}\t{frustration}\n")
                print(f'DATA OF RESIDUE {residue_num} ADDED')
        print(f"All data saved to {output_file}")
        return

    # Mode 2: Only residue provided - Show plot with all 60 monomers
    elif residue_flag is not None and monomer_flag is None:
        all_data = {}
        for monomer in range(1, 61):
            list_of_result_files = list_of_wanted_monomer(frustration_dir, str(monomer))
            local_type, residue_key, dico_frst = dico_frustration_of_a_residue_by_frame(
                list_of_result_files, residue_flag)
            all_data[monomer] = OrderedDict(sorted(dico_frst.items()))
            print(f'DATA OF MONOMER {monomer} ADDED')

        plot_all_frustrations(all_data, "TM", residue_flag)
        return

    # Mode 3: Both residue and monomer provided - Show single plot
    elif residue_flag is not None and monomer_flag is not None:
        list_of_result_files = list_of_wanted_monomer(frustration_dir, monomer_flag)
        local_type, residue_key, dico_frst = dico_frustration_of_a_residue_by_frame(
            list_of_result_files, residue_flag)
        dico_ordonne = OrderedDict(sorted(dico_frst.items()))

        plot_frustration(dico_ordonne, "TM", residue_key, monomer_flag)
        return


if __name__ == "__main__":
    try:
        frustration_dir, monomer_flag, residue_flag = parse_arguments()

        if len(frustration_dir) == 1:
            try:
                main(frustration_dir[0], monomer_flag, residue_flag)
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)
        else:
            print("Usage: python3 frustration_plots_MN_data.py path/to/directory/ [-residue=value] [-monomer=value]")
            sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        print("Usage options:")
        print("1. No parameters: Generate complete data file")
        print("2. -residue=value: Show plot for all 60 monomers of specified residue")
        print("3. -residue=value -monomer=value: Show single plot for specified residue and monomer")
        sys.exit(1)