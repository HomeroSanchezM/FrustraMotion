import sys
import os
import matplotlib.patches as mpatches
from collections import OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import argparse
import os
import umap  # Import UMAP library




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
                if current_monomer_index == monomer_number:
                    #print(f'{current_monomer_index} est pareil que {monomer_number}, on ajoute le path a la liste')
                    full_path = os.path.join(full_path, "FrustrationData")
                    for file in os.listdir(full_path):
                        if file.endswith("singleresidue"):
                            #print(file)
                            full_path= os.path.join(full_path, file)
                            wanted_dirs.append(full_path)
    print(f"For the monomer {monomer_number} we parse {len(wanted_dirs)} files")
    return wanted_dirs

#take a list of result_files and make write a line of the data frame
def line_of_matrix(results_files, file_directory):
    # Mode 1: No parameters - Generate complete data file
    line_dico =  {}
    with open(file_directory, 'a') as f_out:
        for file in results_files :
            with open(file, 'r') as f:
                #Skip header line if present
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

                        else :
                            res_num = int(parts[0])
                        frst_index = float(parts[7])  # FrstIndex value
                        frame = os.path.basename(file).split("_")[1]
                        if res_num in line_dico :
                            #print(f"Ading the frustration of the frame {frame}, of the residue {res_num}")
                            line_dico[res_num][frame]=frst_index
                        else:
                            #print(f"Start parsing residue {res_num}")
                            line_dico[res_num]= {}
        nb_val=0
        for residue,sous_dico in sorted(line_dico.items()):
            for frame, frust in sorted(sous_dico.items()):
                f_out.write(f"{frust}\t")
                nb_val +=1
        f_out.write("\n")
        #print(f"the line have {nb_val} values of frustration")
    f_out.close()






def parse_arguments():
    """Parse command line arguments with flexible options"""

    pdb_file = []

    for arg in sys.argv[1:]:

        if not arg.startswith('-'):
            pdb_file.append(arg)
        else:
            raise ValueError(f"Unknown option: {arg}")

    return pdb_file


def main(frustration_dir):
    # Create output directory if it doesn't exist
    output_dir = "../results/MN_data/SR/new_script/"
    os.makedirs(output_dir, exist_ok=True)
    name = os.path.basename(frustration_dir)
    if name == "" :
        raise ValueError("the directory name must not end with /")
    output_file = os.path.join(output_dir, f"{str(name)}_data_single_residue.txt")
    for monomer in range(1, 61):
        result_files = list_of_wanted_monomer(frustration_dir, monomer)
        line_of_matrix(result_files, output_file)

    # 1. Charger les données
    #df = pd.read_csv(output_file, sep='\t')
    #print("Aperçu des données chargées:")
    #print(df.head())
    #print(len(df))



if __name__ == "__main__":
    try:
        frustration_dir = parse_arguments()

        if len(frustration_dir) == 1:
            try:
                main(frustration_dir[0])
            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)
        else:
            print("Usage: python3 frustration_plots_MN_data.py path/to/directory/ ")
            sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)