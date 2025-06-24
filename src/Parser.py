import sys
import os
import matplotlib.patches as mpatches
from collections import OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

"""
This script takes as a parameter the directory containing the frustration data and creates various dataframes 
with the structure: frame × residues, while keeping the chain name, the frame number, and the residue ID.

Options:
- isolate True/False: indicates whether the frustration was calculated by separating the chains of the structure or not (False by default).
    If isolate == True, the folders containing the frustration data follow this format: ProteinName_FrameNumber_chainName.done
                        the data frames will be separated by chainName 
    If isolate == False, the folders containing the frustration data follow this format: ProteinName_FrameNumber.done
                        if the are more than one chain, the data frames will be separated by chain, 
                        else there will be only one dataframe with all the data

More features in progress...
"""


def parse_arguments():
    """
    this function parse the given argument when executing the script, store the directory containing the frustration data int dir and use argparse for the rest
    :return: dir and isolates (boolean)
    """
    parser = argparse.ArgumentParser(description='Process frustration.')
    parser.add_argument('dir', type=str, help='Directory containing the frustration data')
    parser.add_argument('--isolate', action='store_true', default=False,
                        help='Whether the frustration was calculated by separating the chains')

    args = parser.parse_args()

    # Verify that the directory exists
    if not os.path.isdir(args.dir):
        print(f"Error: Directory {args.dir} does not exist")
        sys.exit(1)

    return args.dir, args.isolate


def parse_frustration_results(dir, isolate):
    """
    this function take the directory containing the frustration data (dir) and a boolean parameter (isolate)
    if isolate == True, it will search for directories in dir having this format : ProteinName_FrameNumber_chainName.done
    after counting how many different chainNames there are it will print the number
    for each chainName it will make a list with the path to all the directories with the same chainName.
    this list is given to function Data_frame(list_of_frustration_data_files) and we keep the return of this function in the variable df

    if isolate == False, it will search for directories in dir having this format : ProteinName_FrameNumber.done
    for each directory with this format it will make a list of the path and give it to function Data_frame(list_of_frustration_data_files)
    and we keep the return of this function in the variable df

    :return: df
    """
    if isolate:
        df = {}
        chain_files = {}
        wrong_format_dirs = []

        # Find all .done directories and group by chain name
        for dirname in os.listdir(dir):
            if dirname.endswith('.done'):
                parts = dirname.split('_')
                if len(parts) >= 3:  # Expected format: Protein_Frame_chain.done
                    chain_name = parts[-1].split('.')[0]  # Last part before .done is chain name
                    full_path = os.path.join(dir, dirname)

                    # Look for FrustrationData/singleresidue file
                    frustration_data_path = os.path.join(full_path, "FrustrationData")
                    if os.path.exists(frustration_data_path):
                        for file in os.listdir(frustration_data_path):
                            if file.endswith("singleresidue"):
                                single_res_path = os.path.join(frustration_data_path, file)

                                if chain_name not in chain_files:
                                    chain_files[chain_name] = []
                                chain_files[chain_name].append(single_res_path)
                else:
                    wrong_format_dirs.append(dirname)

        # Check if any directories don't match the expected format
        if wrong_format_dirs:
            print("\nError: Found directories with wrong format in isolate mode:")
            for dirname in wrong_format_dirs:
                print(f"- {dirname}")
            print("Expected format: ProteinName_FrameNumber_chainName.done")
            sys.exit(1)

        if not chain_files:
            print("\nError: No valid directories found with expected format in isolate mode")
            print("Expected format: ProteinName_FrameNumber_chainName.done")
            sys.exit(1)

        print(f"Found {len(chain_files)} different chains: {', '.join(chain_files.keys())}")

        # Create a dataframe for each chain
        for chain_name, files in chain_files.items():
            df[chain_name] = Data_frame(files, isolate)

    else:
        # List all .done directories (format: Protein_Frame.done)
        files = []
        wrong_format_dirs = []

        for dirname in os.listdir(dir):
            if dirname.endswith('.done'):
                parts = dirname.split('_')
                if len(parts) == 2:  # Expected format: Protein_Frame.done
                    full_path = os.path.join(dir, dirname)
                    frustration_data_path = os.path.join(full_path, "FrustrationData")
                    if os.path.exists(frustration_data_path):
                        for file in os.listdir(frustration_data_path):
                            if file.endswith("singleresidue"):
                                single_res_path = os.path.join(frustration_data_path, file)
                                files.append(single_res_path)
                else:
                    wrong_format_dirs.append(dirname)

        # Check if any directories don't match the expected format
        if wrong_format_dirs:
            print("\nError: Found directories with wrong format in non-isolate mode:")
            for dirname in wrong_format_dirs:
                print(f"- {dirname}")
            print("Expected format: ProteinName_FrameNumber.done")
            sys.exit(1)

        if not files:
            print("\nError: No valid directories found with expected format in non-isolate mode")
            print("Expected format: ProteinName_FrameNumber.done")
            sys.exit(1)

        # Create a single dataframe
        df = Data_frame(files, isolate)

    return df


def Data_frame(list_of_frustration_data_files, isolate):
    """
    this function take a list of files path,
    Using all the files we will create a dataframe with headers with this format :
    chain   residue   frame0   frame10   frame20   ... frame4000

    the files have this format :
Res ChainRes DensityRes AA NativeEnergy DecoyEnergy SDEnergy FrstIndex
1 0 2.969 M -1.273 -0.155 0.931 1.201
2 0 0.000 N -0.658 -0.651 0.147 0.049
...

    if isolate == true,
    the files would have only one chain and it will be the same for all the files,it can be optained from the 2nd column ChainRes
    the residue will be the 3rd colum AA concatenated with the first chain Res . If the fist line of res si not 1, then we make a subtraction to all the line for starting with one (if it is 4 5 6 7, it will be 1 2 3 4)
    the frame is given in the name of the file having the format ProteinName_FrameNumber_chainName.pdb_singleresidue , so fisrt we have to put the list of file in order by the frames number order,
    then each file will complete a column of the dataframe and the column will be the content of the 7th column of the file FrstIndex
    if isolate == False,
    the files could have multiples chains ,it can be optained from the 2nd column ChainRes, a data frame have to be created for each chain of the files
    the residue will be the 3rd colum AA concatenated with the first chain Res . If the fist line of res si not 1, then we make a subtraction to all the line for starting with one (if it is 4 5 6 7, it will be 1 2 3 4)
    the frame is given in the name of the file having the format ProteinName_FrameNumber.pdb_singleresidue , so fisrt we have to put the list of file in order by the frames number order,
    then each file will complete a column of the dataframe and the column will be the content of the 7th column of the file FrstIndex
    :return: dataframe or dictionary of dataframes
    """
    # Step 1: Extract frame numbers and sort files
    file_frame_pairs = []

    for file_path in list_of_frustration_data_files:
        # Get the grandparent directory (the .done directory)
        frustration_dir = os.path.dirname(file_path)
        done_dir = os.path.dirname(frustration_dir)
        done_dir_name = os.path.basename(done_dir)

        # Remove .done extension
        if done_dir_name.endswith('.done'):
            base_name = done_dir_name[:-5]
        else:
            base_name = done_dir_name

        # Split to get frame number (middle part)
        parts = base_name.split('_')
        if isolate:
            # Format: Protein_Frame_chain
            frame_num_str = parts[-2]  # Second to last part is frame number
        else:
            # Format: Protein_Frame
            frame_num_str = parts[-1]  # Last part is frame number

        try:
            frame_num = int(frame_num_str)
        except ValueError:
            frame_num = float('inf')  # Place non-integers at the end

        file_frame_pairs.append((file_path, frame_num, frame_num_str))

    # Sort files by frame number
    file_frame_pairs.sort(key=lambda x: x[1])
    sorted_files = [x[0] for x in file_frame_pairs]
    frame_labels = [f"frame{x[2]}" for x in file_frame_pairs]

    # Step 2: Initialize data structure
    all_chains_data = {}
    n_frames = len(sorted_files)

    # Step 3: Process each file
    for idx, file_path in enumerate(sorted_files):
        try:
            # Read the frustration file (modification ici)
            df_file = pd.read_csv(
                file_path,
                sep='\s+',  # Remplace delim_whitespace=True
                header=0,
                names=[
                    'Res', 'ChainRes', 'DensityRes', 'AA',
                    'NativeEnergy', 'DecoyEnergy', 'SDEnergy', 'FrstIndex'
                ],
                dtype={'Res': str, 'ChainRes': str, 'AA': str}
            )
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            continue

        # Group by chain
        grouped = df_file.groupby('ChainRes')

        for chain_id, group in grouped:
            # Initialize chain data if not present
            if chain_id not in all_chains_data:
                all_chains_data[chain_id] = {}

            # Get residue numbers and find min residue
            residues = group['Res'].tolist()
            try:
                min_res = int(residues[0])
            except ValueError:
                min_res = 0

            # Process each residue in the chain
            for _, row in group.iterrows():
                try:
                    # Calculate normalized residue number
                    res_num = int(row['Res'])
                    norm_res = res_num - min_res + 1
                except ValueError:
                    # If residue is not integer, use original
                    norm_res = row['Res']

                # Create residue label (AA + normalized residue number)
                residue_label = f"{row['AA']}{norm_res}"

                # Initialize residue data if not present
                if residue_label not in all_chains_data[chain_id]:
                    all_chains_data[chain_id][residue_label] = [np.nan] * n_frames

                # Store frustration value
                all_chains_data[chain_id][residue_label][idx] = row['FrstIndex']

    # Step 4: Create dataframes
    result_dfs = {}

    for chain_id, residue_data in all_chains_data.items():
        # Create dataframe for the chain
        df_chain = pd.DataFrame.from_dict(
            residue_data,
            orient='index',
            columns=frame_labels
        )

        # Add chain and residue info
        df_chain = df_chain.reset_index().rename(columns={'index': 'residue'})
        df_chain['chain'] = chain_id

        # Reorder columns: chain, residue, then frames
        cols = ['chain', 'residue'] + frame_labels
        df_chain = df_chain[cols]

        result_dfs[chain_id] = df_chain

    # Return appropriate structure based on number of chains
    if len(result_dfs) == 1:
        return list(result_dfs.values())[0]
    else:
        return result_dfs


def main():
    # Parse command line arguments
    dir_path, isolate = parse_arguments()
    print(f"\nProcessing directory: {dir_path}")
    print(f"Isolate mode: {isolate}\n")

    # Parse frustration results
    print("Starting to parse frustration results...")
    frustration_data = parse_frustration_results(dir_path, isolate)
    print("\nFinished parsing frustration results.")

    # Display information about the parsed data
    if isinstance(frustration_data, dict):
        print(f"\nFound {len(frustration_data)} chains:")
        for chain_id, df in frustration_data.items():
            print(f"Chain {chain_id}: DataFrame with shape {df.shape}")
            print(f"First few rows:\n{df.head()}\n")
    else:
        print("\nSingle DataFrame result:")
        print(f"DataFrame shape: {frustration_data.shape}")
        print(f"First few rows:\n{frustration_data.head()}\n")

    # Additional analysis example
    if isinstance(frustration_data, dict):
        for chain_id, df in frustration_data.items():
            print(f"\nBasic statistics for chain {chain_id}:")
            print(df.describe())
    else:
        print("\nBasic statistics for single chain:")
        print(frustration_data.describe())

    print("\nProcessing complete.")


if __name__ == "__main__":
    main()