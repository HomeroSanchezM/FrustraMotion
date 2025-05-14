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



'''
TO DO : 


or isolate the funtions that make the structural alignment
or organise all the plot generation with options of clculate only a part of it of alls of it

'''
'''
This scripts take one PDB File containing multiples monomers (chains) and make a structural alignment of the monomers and mesure the  mean RMSF for all the monomers per #res
A plot is add to the plot directory, with name of type rmsf_with_std_per_res_<TmEnc|MtEnc>_monomer_<(t)>.png
The aligned PDB files are added to the result directory. 
Usage:

python3 structural_align.py chemin/vers/fichier.pdb

If 2 files given as parameter, the script will align the 2 files structures entirely (not considering monmers)

Usage:

python3 structural_align.py chemin/vers/1erfichier.pdb parser_pdb.py chemin/vers/2emefichier.pdb
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





#5.
def coord_of_atom (aligned_monomers) :
    """return a dico with the list of cords for all the atoms of each residue"""
    # dico_of_atoms = {MET: {CA:[[x1,y1,z1],    , ... },...}
    #                            [x2,y2,z2],
    #                            ...
    #                            [x60,y60,z60]]
    dico_of_atoms = {}

    for k in range (len(list(list(aligned_monomers)[0]))) : #for obtain the number of residues = 264 for TmEnc
        for monomer in aligned_monomers :
            name_res = list(monomer)[k].get_resname() +" "+  str(k+1)
            if name_res in dico_of_atoms :
                atoms = list(monomer)[k].get_atoms()
                for atom in atoms:
                    if atom.id in dico_of_atoms[name_res]:
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)
                    else:
                        dico_of_atoms[name_res][atom.id] = []
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)

            else :
                dico_of_atoms[name_res]={}
                atoms = list(monomer)[k].get_atoms()
                for atom in atoms:
                    if atom.id in dico_of_atoms[name_res]:
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)
                    else:
                        dico_of_atoms[name_res][atom.id] = []
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)
    return  dico_of_atoms

#6.1. Calculate standard deviation for atom positions
def calculate_residue_std(list_of_RMSF):
    """Take a list of atomical RMSF (the RMSF of all the atoms of a residue) and calculate standard deviation for the RMSF of the residue"""
    return np.std(list_of_RMSF)
#6.2.
def calculate_atom_RMSF(list_of_coord):
    """Take a list of 3D coordinates and calculate the corresponding RMSF."""
    # Calcul de la position moyenne
    mean_pos = np.mean(list_of_coord, axis=0)
    #print(mean_pos)
    # Fluctuations (vecteurs de déplacement)
    displacements = list_of_coord - mean_pos
    # Norme au carré de chaque vecteur
    squared_displacements = np.sum(displacements ** 2, axis=1)
    # Moyenne des normes au carré
    mean_squared_disp = np.mean(squared_displacements)
     # RMSF
    rmsf = np.sqrt(mean_squared_disp)
    return rmsf


#7.1 Calculate RMSF and standard deviation for atom positions
def RMSF_std_of_atom(dico_of_atoms):
    """Calculate RMSF and standard deviation for all atoms"""
    dico_of_atom_stats = {}
    for residue in dico_of_atoms:
        dico_of_atom_stats[residue] = {}
        for atom in dico_of_atoms[residue]:
            coords = dico_of_atoms[residue][atom]
            dico_of_atom_stats[residue][atom] = {
                'rmsf': calculate_atom_RMSF(coords),
                #'std': calculate_atom_std(coords)
            }
    return dico_of_atom_stats
#7.2 Calculate RMSF and standard deviation for residue
def RMSF_std_of_Residue(dico_of_atom_RMSF):
    """Calculate average RMSF and std per residue"""
    dico_of_residue_RMSF_std = {}
    for residue in dico_of_atom_RMSF:
        total_rmsf = 0
        total_std = 0
        count = 0
        list_of_RMSF = []
        for atom in dico_of_atom_RMSF[residue]:
            total_rmsf += dico_of_atom_RMSF[residue][atom]['rmsf']
            list_of_RMSF.append(dico_of_atom_RMSF[residue][atom]['rmsf'])
            #total_std += dico_of_atom_RMSF_std[residue][atom]['std']
            count += 1
        dico_of_residue_RMSF_std[residue] = {
            'rmsf': total_rmsf / count,
            'std': calculate_residue_std(list_of_RMSF)
        }
    return dico_of_residue_RMSF_std

#when only a file given
def main(pdb_file1):
    start_time = time.time()
    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filename
    enc_type, enc_number = extract_info_from_filename(pdb_file1)

    #1.1. create the output directory path
    frustration_dir = f"FRUSTRATION_{enc_type.upper()}"
    capsids_dir = f"{enc_type.upper()}_CAPSIDS"
    monomer_dir = "aligned_monomer_for_a_frame"
    results_dir = os.path.join("../results", frustration_dir, capsids_dir,monomer_dir, f"{enc_type}{enc_number}_aligned_monomers")

    plots_dir = os.path.join("../plots", enc_type)

    #1.2 create the repository if it not exist
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    # 2. Charger les monomères
    monomers = load_monomers(pdb_file1)
    print(f"Loaded {len(monomers)} monomers")

    # 3. Alignement structural
    aligned = structural_alignment(monomers)
    print("Structural alignment completed")

    # 4. create aligned PDB files
    save_each_monomer_as_pdb(aligned, results_dir, enc_type, enc_number)

    # 5. create dico of coord of each atom
    dico_of_coords = coord_of_atom(aligned)
    print("atoms coordinates parsed")

    #7.1 Calculate RMSF and standard deviation for atom positions
    dico_atom_RMSF_STD = RMSF_std_of_atom(dico_of_coords)
    print("atoms RMSF calculated")

    #7.2 Calculate RMSF and standard deviation for residue
    dico_res_RMSF_STD = RMSF_std_of_Residue(dico_atom_RMSF_STD)
    print("residues RMSF and std calculated")
    #print(dico_res_RMSF_STD["MET 1"])

    #8 grafical view
    residues = list(dico_res_RMSF_STD.keys())
    rmsf_values = []
    std_values = []
    for res in dico_res_RMSF_STD :
        rmsf_values.append(dico_res_RMSF_STD[res]['rmsf'] )
        std_values.append(dico_res_RMSF_STD[res]['std'])

    plt.figure(figsize=(15, 6))

    # Plot with error bars
    plt.errorbar(range(len(residues)), rmsf_values, yerr=std_values,
                 fmt='-o', markersize=3, linewidth=1,
                 color='blue', alpha=0.7,
                 ecolor='red', elinewidth=0.5, capsize=2,
                 label='RMSF ± std')

    # Customize x-axis
    plt.xticks(range(len(residues))[::3], residues[::3],
               rotation=45, fontsize=8, ha='right')

    plt.xlabel('Residue Number', fontsize=10)
    plt.ylabel('RMSF (Å)', fontsize=10)
    plt.title(f"Residue Flexibility Analysis of {enc_type} at frame {enc_number} ", fontsize=12, pad=20)
    plt.grid(True, alpha=0.3)
    #plt.legend(fontsize=9)
    plt.legend([f'RMSF (#res={len(residues)} , #monomeres={len(monomers)} )'], fontsize=9)
    plt.tight_layout()

    name_plot = f"rmsf_with_std_per_res_{enc_type}_frame_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Plot with standard deviation saved as {plot_path}")
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

#when 2 files given in parameters
def main2(pdb_file1, pdb_file2):
    start_time = time.time()
    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filenames
    enc_type1, enc_number1 = extract_info_from_filename(pdb_file1)
    enc_type2, enc_number2 = extract_info_from_filename(pdb_file2)

    #1.1. create the output directory path
    common_dir = "COMMON" #could not be common but has to be generalised
    type_dir = "CHAIN_A" #or MONOMER or CHAIN_A has to be automated, CAPSID would not work cause have more than 99999 atoms
    results_dir = os.path.join("../results", common_dir, type_dir, f"{enc_type1}{enc_number1}_{enc_type2}{enc_number2}_structural_align")#f"{enc_type1}{enc_number1}_{enc_type2}{enc_number2}_structural_align"

    #1.2 create the repository if it not exist
    os.makedirs(results_dir, exist_ok=True)

    # 2. Charger les pdb
    structure_list = []
    pdb1 = load_pdb(pdb_file1)
    pdb2 = load_pdb(pdb_file2)
    structure_list.append(pdb1)
    structure_list.append(pdb2)
    print(f"Loaded two structures")

    # 3. Alignement structural
    aligned = structural_alignment(structure_list)
    print("Structural alignment completed")

     # 4. create aligned PDB files
    save_each_monomer_as_pdb(aligned, results_dir, enc_type1, enc_number1 , enc_type2, enc_number2)
    print("finish adding files")

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

if __name__ == "__main__":
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("Usage: python3 analyze_encapsulin.py path/to/encapsulin.pdb \n or python3 structural_align.py chemin/vers/1erfichier.pdb parser_pdb.py chemin/vers/2emefichier.pdb")
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