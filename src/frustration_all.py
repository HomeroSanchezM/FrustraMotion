import sys
import os
from frustration_plots import main
from frustration_plots import parse_arguments


def main_all(pdb_directory, vmd_flag=False, frustration_flag=False, seqdist_flag=12, isolate_flag=True,
             plots_flag=False, output_file="../results/tab_results.txt"):
    """
    Process all PDB files in the given directory by calling main() for each file.
    Results are saved in a tabulated text file with columns:
    pdb_file, chain, residue, frustration

    Args:
        pdb_directory: Directory containing PDB files to process
        vmd_flag: Whether to execute VMD visualization
        frustration_flag: Whether to calculate frustration
        seqdist_flag: Sequence distance parameter
        isolate_flag: Whether to process as isolated monomers
        plots_flag: Whether to generate plots
        output_file: Name of the output text file
    """
    # Verify the directory exists
    if not os.path.isdir(pdb_directory):
        raise ValueError(f"Directory not found: {pdb_directory}")

    # Create and open the output file
    with open(output_file, 'w') as f:
        # Write header
        header = "pdb_file\tchain\tresidue\tfrustration\n"
        f.write(header)

        # Process each PDB file in the directory
        for filename in os.listdir(pdb_directory):
            if filename.endswith('.pdb'):
                full_path = os.path.join(pdb_directory, filename)
                print(f"\nProcessing file: {filename}")

                try:
                    # Call main() with the current PDB file
                    result = main(full_path,
                                  vmd_flag=vmd_flag,
                                  frustration_flag=frustration_flag,
                                  seqdist_flag=seqdist_flag,
                                  isolate_flag=isolate_flag,
                                  plots_flag=plots_flag)

                    # Write results to file
                    for chain_id, residues in result.items():
                        for residue, frustration in residues.items():
                            line = f"{filename}\t{chain_id}\t{residue}\t{frustration:.3f}\n"
                            f.write(line)

                    print(f"Results for {filename} written to {output_file}")

                except Exception as e:
                    print(f"Error processing {filename}: {e}")
                    continue

    print(f"\nAll results saved to {output_file}")


if __name__ == "__main__":
    # Parse command line arguments
    pdb_files, vmd_flag, frustration_flag, seqdist_flag, isolate_flag, plots_flag = parse_arguments()

    if len(pdb_files) == 1:
        # Check if the argument is a directory
        if os.path.isdir(pdb_files[0]):
            try:
                # Define output file name based on directory
                dir_name = os.path.basename(os.path.normpath(pdb_files[0]))
                output_file = f"frustration_results_{dir_name}.txt"

                main_all(pdb_files[0],
                         vmd_flag=vmd_flag,
                         frustration_flag=frustration_flag,
                         seqdist_flag=seqdist_flag,
                         isolate_flag=isolate_flag,
                         plots_flag=plots_flag,
                         output_file=output_file)

            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)
        else:
            # If it's a single file, process it and save to specific output file
            try:
                base_name = os.path.splitext(os.path.basename(pdb_files[0]))[0]
                output_file = f"frustration_results_{base_name}.txt"

                # Open output file
                with open(output_file, 'w') as f:
                    # Write header
                    header = "pdb_file\tchain\tresidue\tfrustration\n"
                    f.write(header)

                    # Process file
                    result = main(pdb_files[0],
                                  vmd_flag=vmd_flag,
                                  frustration_flag=frustration_flag,
                                  seqdist_flag=seqdist_flag,
                                  isolate_flag=isolate_flag,
                                  plots_flag=plots_flag)

                    # Write results
                    for chain_id, residues in result.items():
                        for residue, frustration in residues.items():
                            line = f"{os.path.basename(pdb_files[0])}\t{chain_id}\t{residue}\t{frustration:.3f}\n"
                            f.write(line)

                print(f"Results saved to {output_file}")

            except ValueError as e:
                print(f"Error: {e}")
                sys.exit(1)
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                sys.exit(1)
    else:
        print("Usage: python3 frustration_all.py path/to/pdb_directory [-vmd=True] [-frustration=True]")
        print("Or: python3 frustration_all.py path/to/single.pdb [-vmd=True] [-frustration=True]")
        sys.exit(1)