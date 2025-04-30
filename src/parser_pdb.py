from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning


file_directory = "../data/FRUSTRATION_TMENC/FRUSTRATION_ANALYSIS_TMENC_CAPSIDS/MONOMER_CHAINA/TmEnc_monomer_0.pdb"


'''
controles a faire:
le fichier PDB contient quelque chose qui n'est pas un AA, faut indiquer dans quelle ligne
'''

#fonction thar take return the residues sequence on a given PDB
def return_sequence_3_letter_format (file_directory) :
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

#fonction that converte a list of 3 letter code AminoAcids to a string on 1 letter code AminoAcids
def three_to_one (list_of_tree):
    dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H','HSD': 'H',
           'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
           'TYR': 'Y', 'VAL': 'V'}

    seq = ''
    count = 0
    for residue in list_of_tree:

        if residue in dic:
            seq+=dic[residue]
            count += 1
        else:
            warnings.warn(f"the residue number {count} is in a invalid format")
            seq+='_'

    return seq


print(three_to_one(return_sequence_3_letter_format(file_directory )))