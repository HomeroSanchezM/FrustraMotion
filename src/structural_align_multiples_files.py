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


'''
This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a wanted monomer, and make a structural alignment of this monomer and mesure the  mean RMSF for this monomers in each file (frames) 
A plot is add to the plot directory, with name of type rmsf_with_std_per_frame_<TmEnc|MmEnc>_monomer_<(i)>.png
The aligned PDB files are added to the result directory. 
Usage:

python3 structural_align_multiples_files.py path/to/directory/ number_of_the_wanted_monomer

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

def visualization_VMD_script(output_dir, enc_type, monomer_number, execute_vmd = False):
    """
    This function generates a tcl script for structural alignment visualization using VMD,
    which loads multiple monomer PDB files and applies different colors to each.

    :param output_dir: where to create the tcl file
    :param enc_type1: <TmEnc|MtEnc>
    :param enc_number1: the number of the frame studied
    """

    script_content = f"""# Load 491 monomers {enc_type}0_monomer{monomer_number}.pdb ... {enc_type}4900_monomer{monomer_number}.pdb

set num_monomers 491
set base_filename "{output_dir}/{enc_type}"
set monomer_suffix "_monomer{monomer_number}.pdb"

# Colors - using VMD's default color IDs
set available_colors [list 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]

# Calculate frame numbers (0, 10, 20, ..., 4900)
set frame_numbers {{}}
for {{set i 0}} {{$i < $num_monomers}} {{incr i}} {{
    lappend frame_numbers [expr {{$i * 10}}]
}}

# Load each structure with proper naming
set color_index 0
foreach frame $frame_numbers {{
    set filename "${{base_filename}}${{frame}}${{monomer_suffix}}"
    puts "Loading $filename..."
    
    # Load PDB file
    mol new $filename type pdb waitfor all
    
    # Apply color (cycle through available colors)
    mol modcolor 0 [molinfo top] "ColorID" [lindex $available_colors $color_index]
    set color_index [expr {{($color_index + 1) % [llength $available_colors]}}]
    
    # Apply representation (optional)
    mol modstyle 0 [molinfo top] "NewCartoon"
}}

puts "Successfully loaded $num_monomers structures."

"""

    # Create the output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path
    output_file = os.path.join(output_dir, "load_data.tcl")

    # Write the script to file
    with open(output_file, 'w') as f:
        f.write(script_content)
    # Execute VMD if requested
    if execute_vmd:
        try:

            subprocess.run(["vmd", "-e", output_file], check=True)
            print(f"Successfully executed VMD with script: {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing VMD: {e}")
        except FileNotFoundError:
            print("VMD not found. Please make sure VMD is installed and in your PATH.")


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



def main(pdb_directory, monomer_number):
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
    frames_dir = "aligned_frames_for_a_monomer"
    results_dir = os.path.join("../results", frustration_dir, capsids_dir,frames_dir, f"{enc_type}_monomer_{monomer_number}_aligned_frames")

    plots_dir = os.path.join("../plots", enc_type)

    #1.2 create the repository if it not exist
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    #2
    monomers = load_monomers(files_dict,chain_dict, monomer_number)
    print(f"Loaded {len(monomers)} monomers")

     # 3. Alignement structural
    aligned = structural_alignment(monomers)
    print("Structural alignment completed")

    # 4. create aligned PDB files
    save_each_monomer_as_pdb(aligned, results_dir, enc_type, files_dict, monomer_number)

    #4.2 alignment visualization
    visualization_VMD_script(results_dir, enc_type, monomer_number, True)
    # 5. create dico of coord of each atom
    dico_of_coords = coord_of_atom(aligned)
    print("atoms coordinates parsed")

    #7.1 Calculate RMSF and standard deviation for atom positions
    dico_atom_RMSF_STD = RMSF_std_of_atom(dico_of_coords)
    print("atoms RMSF and std calculated")

    #7.2 Calculate RMSF and standard deviation for residue
    dico_res_RMSF_STD = RMSF_std_of_Residue(dico_atom_RMSF_STD)
    print("residues RMSF and std calculated")
    print(dico_res_RMSF_STD["MET 1"])

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
    plt.title(f"Residue Flexibility Analysis of  monomer {monomer_number} of {enc_type} for {len(files_dict)} frames ", fontsize=12, pad=20)
    plt.grid(True, alpha=0.3)
    #plt.legend(fontsize=9)
    plt.legend([f'RMSF (#res={len(residues)} , #frames={len(monomers)} )'], fontsize=9)
    plt.tight_layout()

    name_plot = f"rmsf_with_std_per_frame_{enc_type}_monomer_{monomer_number}.png"
    plot_path = os.path.join(plots_dir, name_plot)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Plot with standard deviation saved as {plot_path}")
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution: {execution_time:.2f} secondes")

if __name__ == "__main__":
    if len(sys.argv) != 3 :
        print("Usage: python analyze_encapsulin.py path/to/encapsulin.pdb number_of_wanted_monomer")
        sys.exit(1)
    elif len(sys.argv) ==3:
        try:
            main(sys.argv[1], sys.argv[2])
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)
