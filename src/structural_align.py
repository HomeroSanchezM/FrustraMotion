from Bio.PDB import *
import numpy as np
import warnings
from Bio import BiopythonWarning
import sys
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain
import os
import re


'''
TO DO : 
take the <(t)> number and the type <MtEnc|TmEnc> of the file name for adding to the name of the plot automaticly
to do it have to change :
l 111 : f"<MtEnc|TmEnc><0|...|1000|...>_monomer{i}.pdb"
l 243 : "../results/<FRUSTRATION_MTENC|FRUSTRATION_TMENC>/<MTENC_CAPSIDS|TMENC_CAPSIDS>/<MtEnc|TmEnc><0|...|1000|...>_aligned_monomers")
l 293 :  "rmsf_with_std_<MtEnc|TmEnc>_capsids_<0|...|1000|...>.png"
l 234 : f"../plots/<MtEnc|TmEnc>/{name_plot}"
AND
or isolate the funtions that make the structural alignment
or organise all the plot generation with options of clculate only a part of it of alls of it

'''
'''
This scripts take one pdf File containing multiples monomers (chains) and make a structural alignment and mesure the  mean RMSF for all the monomers per #res
A plot is add to the plot directory, with name of type rmsf_per_res_TmEnc_capsids_<(t)>.png
The aligned PDB files are added to the result directory. 
Usage:

python3 structural_align.py chemin/vers/fichier.pdb
'''

# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)
#1.
def extract_info_from_filename(pdb_file):
    """Extract the type (MtEnc/TmEnc) and number (t) from the filename"""
    filename = os.path.basename(pdb_file)
    # Pattern search of a name followed by a number, have to be generalised of others types of encapsulins / proteins
    match = re.search(r'(MtEnc|TmEnc)_(\d+)', filename)
    if match:
        enc_type = match.group(1)
        enc_number = match.group(2)
        return enc_type, enc_number
    else:
        # error if the name of the given file not follow this structure, have to be generalised
        raise ValueError(
            f"Filename '{filename}' doesn't match expected pattern. "
            "Expected format: '.../MtEnc<number>.pdb' or '.../TmEnc<number>.pdb' "
            "where <number> is one or more digits."
        )

#2.
def load_monomers(pdb_file):
    """Load the n monomers, corresponding to the chains"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)

    chains = []
    for model in structure :
        for chain in model:
            chains.append(chain)
    return chains

#3.
def structural_alignment(monomers):
    """Align all the monomers to the first one"""
    super_imposer = Superimposer()
    reference = monomers[0]
    aligned = [reference]

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

    return aligned

#4.
def save_each_monomer_as_pdb(aligned_monomers, output_dir, enc_type , enc_number):
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

    for i, monomer in enumerate(aligned_monomers, start=1):
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

        # name of the output file
        output_file = os.path.join(output_dir, f"{enc_type}{enc_number}_monomer{i}.pdb")

        # create de pdb
        io = PDBIO()
        io.set_structure(structure)
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
def calculate_atom_std(list_of_coord):
    """Take a list of 3D coordinates and calculate standard deviation for atom positions"""
    mean_pos = np.mean(list_of_coord, axis=0)
    displacements = list_of_coord - mean_pos
    distances = np.sqrt(np.sum(displacements**2, axis=1))
    return np.std(distances)

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
                'std': calculate_atom_std(coords)
            }
    return dico_of_atom_stats
#7.2 Calculate RMSF and standard deviation for residue
def RMSF_std_of_Residue(dico_of_atom_RMSF_std):
    """Calculate average RMSF and std per residue"""
    dico_of_residue_RMSF_std = {}
    for residue in dico_of_atom_RMSF_std:
        total_rmsf = 0
        total_std = 0
        count = 0
        for atom in dico_of_atom_RMSF_std[residue]:
            total_rmsf += dico_of_atom_RMSF_std[residue][atom]['rmsf']
            total_std += dico_of_atom_RMSF_std[residue][atom]['std']
            count += 1
        dico_of_residue_RMSF_std[residue] = {
            'rmsf': total_rmsf / count,
            'std': total_std / count
        }
    return dico_of_residue_RMSF_std


def main(pdb_file):

    # 1. Extract the type (MtEnc/TmEnc) and number (t) from the filename
    enc_type, enc_number = extract_info_from_filename(pdb_file)

    #1.1. create the output directory path
    frustration_dir = f"FRUSTRATION_{enc_type.upper()}"
    capsids_dir = f"{enc_type.upper()}_CAPSIDS"
    results_dir = os.path.join("../results", frustration_dir, capsids_dir, f"{enc_type}{enc_number}_aligned_monomers")

    plots_dir = os.path.join("../plots", enc_type)

    #1.2 create the repository if it not exist
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    # 2. Charger les monomères
    monomers = load_monomers(pdb_file)
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
    print("atoms RMSF and std calculated")

    #7.2 Calculate RMSF and standard deviation for residue
    dico_res_RMSF_STD = RMSF_std_of_Residue(dico_atom_RMSF_STD)
    print("residues RMSF and std calculated")

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
    plt.title('Residue Flexibility Analysis with Standard Deviation', fontsize=12, pad=20)
    plt.grid(True, alpha=0.3)
    #plt.legend(fontsize=9)
    plt.legend([f'RMSF (#res={len(residues)} , #monomeres={len(monomers)} )'], fontsize=9)
    plt.tight_layout()

    name_plot = f"rmsf_with_std_{enc_type}_capsids_{enc_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Plot with standard deviation saved as {plot_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_encapsulin.py path/to/encapsulin.pdb")
        sys.exit(1)
    try:
        main(sys.argv[1])
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
