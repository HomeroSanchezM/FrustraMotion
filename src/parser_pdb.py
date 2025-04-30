from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning
import sys
import os

file_directory = "../data/FRUSTRATION_TMENC/FRUSTRATION_ANALYSIS_TMENC_CAPSIDS/MONOMER_CHAINA/TmEnc_monomer_0.pdb"

'''
Usage:
python parser_pdb.py chemin/vers/fichier.pdb
'''

#fonction thar take return the residues sequence on a given PDB
def return_sequence_3_letter_format (file_directory) :
    """Extract residue sequence from PDB file in 3-letter code format"""
    # Suppression des warnings non critiques
    warnings.simplefilter('ignore', BiopythonWarning)
    list = []

    try:
        # Création d'un objet PDBParser
        parser = PDBParser(QUIET=True)  # QUIET=True supprime certains warnings
        # Chargement d'un fichier PDB local
        structure = parser.get_structure("struct", file_directory)
        # iterate over all residues in a structure
        for residue in structure.get_residues():
            list.append(residue.get_resname())
    except Exception as e:
        print(f"Error parsing PDB file: {e}", file=sys.stderr)
        sys.exit(1)

    return list

#fonction that converte a list of 3 letter code AminoAcids to a string on 1 letter code AminoAcids
def three_to_one (list_of_tree):
    """Convert 3-letter AA codes to 1-letter format"""
    dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H','HSD': 'H','HSE': 'H',
           'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
           'TYR': 'Y', 'VAL': 'V'}

    seq = ''
    count = 0
    for residue in list_of_tree:

        if residue in dic:
            seq+=dic[residue]
            count += 1
        else:
            warnings.warn(f"the residue number {count} is in a invalid format ({residue})")
            seq+='_'

    return seq

def main():
    if len(sys.argv) != 2:
        print("Usage: python parser_pdb.py path/to/your_file.pdb", file=sys.stderr)
        sys.exit(1)

    pdb_file = sys.argv[1]

    if not os.path.isfile(pdb_file):
        print(f"Error: File not found - {pdb_file}", file=sys.stderr)
        sys.exit(1)

    three_letter_seq = return_sequence_3_letter_format(pdb_file)
    one_letter_seq = three_to_one(three_letter_seq)
    print("\nResults:")
    print(f"Input PDB file: {pdb_file}")
    print(f"Residue count: {len(three_letter_seq)}")
    print(f"3-letter sequence: {' '.join(three_letter_seq)}")
    print(f"1-letter sequence: {one_letter_seq}")

if __name__ == "__main__":
    main()