from Bio.PDB import PDBParser
#for ignoring warning after parsing
import warnings
from Bio import BiopythonWarning
#for using parameters during execution
import sys
import os
# for protein alignment
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


'''
Usage:
get the sequence of one PDB :python3 parser_pdb.py chemin/vers/fichier.pdb
Get the seqeunce for two PBD and compare it : python3 parser_pdb.py chemin/vers/1erfichier.pdb parser_pdb.py chemin/vers/2emefichier.pdb
'''

#fonction thar take a PDB file directory and return the residues sequence on a given PDB
def return_sequence_3_letter_format (file_directory) :
    """Extract residue sequence from PDB file in 3-letter code format"""
    # Ignore non essential warnings
    warnings.simplefilter('ignore', BiopythonWarning)
    residues = []

    try:
        # Création d'un objet PDBParser
        parser = PDBParser(QUIET=True)  # QUIET=True supprime certains warnings
        # Chargement d'un fichier PDB local
        structure = parser.get_structure("struct", file_directory)
        # iterate over all residues in a structure
        for residue in structure.get_residues():
            residues.append(residue.get_resname())
    except Exception as e:
        print(f"Error parsing PDB file: {e}", file=sys.stderr)
        sys.exit(1)

    return residues

#fonction that convert a list of 3 letter code AminoAcids to a string on 1 letter code AminoAcids
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

#fonction that take 2 sequences and make an alignement.
def align_seqs(SEQ1 , SEQ2):
    seq1 = Seq(SEQ1)
    seq2 = Seq(SEQ2)
    # Init align
    ALIGN = PairwiseAligner()
    alignments = ALIGN.align(seq1, seq2) #by defaut is a global alignement

    # give the best align
    best_alignment = alignments[0]
    return best_alignment.score , best_alignment

def main():
    if len(sys.argv) != 2 and len(sys.argv) != 3 :
        print("Usage: python parser_pdb.py path/to/your_file.pdb \n or python3 parser_pdb.py chemin/vers/1erfichier.pdb parser_pdb.py chemin/vers/2emefichier.pdb", file=sys.stderr)
        sys.exit(1)

    elif len(sys.argv) ==2:
        pdb_file = sys.argv[1]

        if not os.path.isfile(pdb_file):
            print(f"Error: File not found - {pdb_file}", file=sys.stderr)
            sys.exit(1)
        else :
            option = input("Choose between the next options: \n1 : display the sequence on the given PDB \n2 : make a structural alignment of the monomers of the given PDB \n3 : align with a custom sequence \n")
            if option == "1":
                three_letter_seq = return_sequence_3_letter_format(pdb_file)
                one_letter_seq = three_to_one(three_letter_seq)
                print("\nResults:")
                print(f"Input PDB file: {pdb_file}")
                print(f"Residue count: {len(three_letter_seq)}")
                #print(f"3-letter sequence: {' '.join(three_letter_seq)}")
                print(f"1-letter sequence: {one_letter_seq}")
            elif option == "2":
                print('structural_align.py')
            elif option == "3":
                pdb_file1 = sys.argv[1]
                if not os.path.isfile(pdb_file1):
                    print(f"Error: File not found - {pdb_file1}", file=sys.stderr)
                    sys.exit(1)
                else:
                    three_letter_seq1 = return_sequence_3_letter_format(pdb_file1)
                    one_letter_seq1 = three_to_one(three_letter_seq1)
                    custom_seq = input('input the custom sequence :\n')
                    score, alignment = align_seqs(one_letter_seq1, custom_seq)
                    print("\nResults:")
                    print(f"Input PDB files: \n {pdb_file1} and a custom sequence")
                    print(f"Residue count: \n {len(three_letter_seq1)} \n {len(custom_seq)}")
                    print(f"Score: {score}")
                    print(f"alignment : \n{alignment}")
            else :
                print("Error: Invalid option, have to write 1 or 2 ", file=sys.stderr)
                sys.exit(1)

    elif len(sys.argv) ==3:
        pdb_file1 = sys.argv[1]
        pdb_file2 = sys.argv[2]
        if not os.path.isfile(pdb_file1):
            print(f"Error: File not found - {pdb_file1}", file=sys.stderr)
            sys.exit(1)
        elif not os.path.isfile(pdb_file2):
            print(f"Error: File not found - {pdb_file2}", file=sys.stderr)
            sys.exit(1)
        else :
            three_letter_seq1 = return_sequence_3_letter_format(pdb_file1)
            one_letter_seq1 = three_to_one(three_letter_seq1)
            three_letter_seq2 = return_sequence_3_letter_format(pdb_file2)
            one_letter_seq2 = three_to_one(three_letter_seq2)
            score, alignment = align_seqs(one_letter_seq1, one_letter_seq2)
            print("\nResults:")
            print(f"Input PDB files: \n {pdb_file1} \n {pdb_file2}")
            print(f"Residue count: \n {len(three_letter_seq1)} \n {len(three_letter_seq2)}")
            print(f"Score: {score}")
            print(f"alignment : \n{alignment}")



if __name__ == "__main__":
    main()