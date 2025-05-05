from Bio.PDB import *
import numpy as np
import warnings
from Bio import BiopythonWarning
import sys
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain
import os



# delete non essential warning
warnings.simplefilter('ignore', BiopythonWarning)

#1
def load_monomers(pdb_file):
    """Load the n monomers, corresponding to the chains"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("encapsulin", pdb_file)

    chains = []
    for model in structure :
        for chain in model:
            chains.append(chain)
            #print(list(list(chain)[0])[0].get_serial_number())
    return chains

#2.
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
#test 2.1


def save_each_monomer_as_pdb(aligned_monomers, output_dir="aligned_monomers"):
    """
    Enregistre chaque monomère dans un fichier PDB séparé

    Args:
        aligned_monomers: Liste d'objets Chain (monomères alignés)
        output_dir: Répertoire de sortie (sera créé si inexistant)
    """
    # Créer le répertoire de sortie
    os.makedirs(output_dir, exist_ok=True)

    for i, monomer in enumerate(aligned_monomers, start=1):
        # Créer une structure minimale pour ce monomère
        structure = Structure.Structure(f"MONOMER_{i}")
        model = Model.Model(0)
        structure.add(model)

        # Créer une nouvelle chaîne (ID 'A' par défaut)
        new_chain = Chain.Chain("A")

        # Copier tous les résidus du monomère
        for residue in monomer:
            new_chain.add(residue)

        model.add(new_chain)

        # Renumérotation propre des atomes
        atom_number = 1
        for atom in structure.get_atoms():
            atom.serial_number = atom_number
            atom_number += 1

        # Nom du fichier de sortie
        output_file = os.path.join(output_dir, f"TmEnc0_monomer{i}.pdb")

        # Sauvegarde
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)

        print(f"Monomère {i} sauvegardé dans {output_file}")





#3.
def coord_of_atom (aligned_monomers) :
    """return a dico with the list of cords for all the atoms of each residue"""
    # dico_of_atoms = {MET: {CA:[[x1,y1,z1],    , ... },...}
    #                            [x2,y2,z2],
    #                            ...
    #                            [x60,y60,z60]]
    dico_of_atoms = {}

    for k in range (len(list(list(aligned_monomers)[0]))) : #for obtain the number of residues = 264 for TmEnc
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

#4.1. Calculate standard deviation for atom positions
def calculate_atom_std(list_of_coord):
    """Take a list of 3D coordinates and calculate standard deviation for atom positions"""
    mean_pos = np.mean(list_of_coord, axis=0)
    displacements = list_of_coord - mean_pos
    distances = np.sqrt(np.sum(displacements**2, axis=1))
    return np.std(distances)

#4.2.
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

#5.1 Calculate RMSF and standard deviation for atom positions
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
#5.2
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

    #test 2.1
    # Exemple d'utilisation
    save_each_monomer_as_pdb(aligned, "../results/FRUSTRATION_TMENC/TMENC_CAPSIDS/TmEnc0_aligned_monomers/")

    # 3. create dico of coord of each atom
    dico_of_coords = coord_of_atom(aligned)
    print("coordinates parsed")

    # test of the funtion calculate_atom_std
    #test_list = [[5, 6, 2], [4, 2, 1], [7, 2, 3], [4, 1, 5]]
    #test_res =calculate_atom_std(test_list)
    #print(test_res)

    #5.1
    dico_atom_RMSF_STD = RMSF_std_of_atom(dico_of_coords)
    print("atoms RMSF and std calculated")

    #5.2
    dico_res_RMSF_STD = RMSF_std_of_Residue(dico_atom_RMSF_STD)
    #print(dico_res_RMSF_STD)
    print("residues RMSF and std calculated")

    #6 grafical view
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

    name_plot = "rmsf_with_std_TmEnc_capsids_0.png"
    plt.savefig(f"../plots/{name_plot}", dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Plot with standard deviation saved as {name_plot}")



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_encapsulin.py path/to/encapsulin.pdb")
        sys.exit(1)

    main(sys.argv[1])