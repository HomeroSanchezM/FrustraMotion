from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning


file_directory = "../data/FRUSTRATION_TMENC/FRUSTRATION_ANALYSIS_TMENC_CAPSIDS/MONOMER_CHAINA/TmEnc_monomer_0.pdb"


def return_sequence (file_directory) :
    # Suppression des warnings non critiques
    warnings.simplefilter('ignore', BiopythonWarning)
    list = []
    # Création d'un objet PDBParser
    parser = PDBParser(QUIET=True)  # QUIET=True supprime certains warnings
    # Création d'un objet PDBParser
    # Chargement d'un fichier PDB local
    structure = parser.get_structure("struct", file_directory)
    # iterate over all residues in a structure
    for residue in structure.get_residues():
        list.append(residue.get_resname())
    return list


print(return_sequence(file_directory ))