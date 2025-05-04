from Bio.PDB import *
import numpy as np
import warnings
from Bio import BiopythonWarning
import sys
import os
import matplotlib.pyplot as plt

# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)


def load_monomers(pdb_file):
    """Load the n monomers, corresponding to the chains"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)

    chains = []
    for model in structure :
        for chain in model:
            chains.append(chain)
    return chains

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

def coord_of_atom (aligned_monomers) :
    """return a dico with the list of cords for all the atoms of each residue"""
    # dico_of_atoms = {MET: {CA:[[x1,y1,z1],    , ... },...}
    #                            [x2,y2,z2],
    #                            ...
    #                            [x60,y60,z60]]
    dico_of_atoms = {}

    for k in range (len(list(aligned_monomers))) :
        #print(f"on a k = {k}")
        for monomer in aligned_monomers :
            #print(f"on a monomer = {monomer}")
            name_res = list(monomer)[k].get_resname() +" "+  str(k+1)
            #print(f"on a name_res = {name_res}")
            if name_res in dico_of_atoms :
                #print(f"name_res = {name_res} est dans dico_of_atoms")
                atoms = list(monomer)[k].get_atoms()
                #print(list(atoms))
                for atom in atoms:
                    if atom.id in dico_of_atoms[name_res]:
                        #print(f"on a atom.id = {atom.id} qui est dans dico_of_atoms[name_res] = {name_res} ")
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)
                    else:
                        #print(f"PAS atom.id = {atom.id} qui est dans dico_of_atoms[name_res] = {name_res} , on cree la liste")
                        dico_of_atoms[name_res][atom.id] = []
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)

            else :
                #print(f"PAS name_res = {name_res} est dans dico_of_atoms, on cree le dico")
                dico_of_atoms[name_res]={}
                atoms = list(monomer)[k].get_atoms()
                # print(list(atoms))
                for atom in atoms:
                    if atom.id in dico_of_atoms[name_res]:
                        # print(f"on a atom.id = {atom.id} qui est dans dico_of_atoms[name_res] = {name_res} ")
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)
                    else:
                        #print(f"PAS atom.id = {atom.id} qui est dans dico_of_atoms[name_res] = {name_res} , on cree la liste")
                        dico_of_atoms[name_res][atom.id] = []
                        coord = atom.get_coord()
                        dico_of_atoms[name_res][atom.id].append(coord)


    #print(dico_of_atoms['MET1'])
    return  dico_of_atoms


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

def RMSF_of_atom (dico_of_atoms) :
    """return a dico with the RMSF for all the atoms of each residue"""
    # RMSF_of_atom = {MET: {CA:3.21,... },...}

    dico_of_atom_RMSF = {}
    for residue in dico_of_atoms :
        if residue in dico_of_atom_RMSF :
            for atom in dico_of_atoms[residue]:
                dico_of_atom_RMSF[residue][atom] = calculate_atom_RMSF(dico_of_atoms[residue][atom])
        else:
            dico_of_atom_RMSF[residue]= {}
            for atom in dico_of_atoms[residue]:
                dico_of_atom_RMSF[residue][atom] = calculate_atom_RMSF(dico_of_atoms[residue][atom])

    return dico_of_atom_RMSF

def calculate_residue_RMSF(dico_of_atoms_RMSF) :
    residue_RMSF = 0
    for atom in dico_of_atoms_RMSF :
        residue_RMSF += dico_of_atoms_RMSF[atom]
    residue_RMSF = residue_RMSF / len(dico_of_atoms_RMSF)
    return residue_RMSF

def RMSF_of_Residue (dico_of_atom_RMSF) :
    """return a dico with the RMSF for all the residues"""
    # RMSF_of_atom = {MET: 3.14 ,...}
    dico_of_residue_RMSF = {}
    for residue in dico_of_atom_RMSF:
            dico_of_residue_RMSF[residue] = calculate_residue_RMSF(dico_of_atom_RMSF[residue])

    return dico_of_residue_RMSF


def main(pdb_file):
    # 1. Charger les monomères
    monomers = load_monomers(pdb_file)
    print(f"Loaded {len(monomers)} monomers")

    #test of the funtion calculate_atom_RMSF
    #test_list = [[5,6,2],[4,2,1],[7,2,3],[4,1,5]]
    #test_res =calculate_atom_RMSF(test_list)
    #print(test_res)

    # 2. Alignement structural
    aligned = structural_alignment(monomers)
    print("Structural alignment completed")
    # print(aligned)

    # 3. create dico of coord of each atom
    dico_of_coords = coord_of_atom(aligned)
    print("coordinates parsed")


    # 4. create dico of RMSF of each atom
    dico_of_atom_RMSF= RMSF_of_atom(dico_of_coords)
    print("atoms RMSF calculated")

    # test of the funtion calculate_residue_RMSF
    #test_dico = {'atom1' : 3.21 , 'atom2' : 3.58, 'atom3' : 7.21}
    #test_res =calculate_residue_RMSF(test_dico)
    #print(test_res)

    # 4. create dico of RMSF of each residue
    dico_of_residue_RMSF = RMSF_of_Residue(dico_of_atom_RMSF)
    #print(dico_of_residue_RMSF)
    print("residue RMSF calculated")

    # 5. grafical view
    # Extract data
    residues = list(dico_of_residue_RMSF.keys())
    rmsf_values = list(dico_of_residue_RMSF.values())

    # Tracer le graphique
    plt.figure(figsize=(10, 5))
    plt.plot(residues, rmsf_values, marker='o', linestyle='-', color='blue', label='RMSF')
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)

    # Ajouter titres et légendes
    plt.xlabel('Residues')
    plt.ylabel('RMSF')
    plt.title('RMSF for monomers residues')
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    name_plot = "rmsf_per_res_TmEnc_capsids_1000.png"
    plt.savefig(f"../plots/{name_plot}", dpi=300)
    print(f"plot {name_plot} created in folder plot/")




if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_encapsulin.py path/to/encapsulin.pdb")
        sys.exit(1)

    main(sys.argv[1])