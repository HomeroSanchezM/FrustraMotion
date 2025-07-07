import os
import sys
import pandas as pd
import argparse
from collections import defaultdict


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process contact frustration data.')
    parser.add_argument('dir', type=str, help='Directory containing the frustration data')
    parser.add_argument('--isolate', action='store_true', default=False,
                        help='Whether the frustration was calculated by separating the chains')

    args = parser.parse_args()

    if not os.path.isdir(args.dir):
        print(f"Error: Directory {args.dir} does not exist")
        sys.exit(1)

    return args.dir, args.isolate


def get_protein_name(dir_path, isolate):
    """Extract protein name from directory structure"""
    for dirname in os.listdir(dir_path):
        if dirname.endswith('.done'):
            parts = dirname.split('_')
            if isolate and len(parts) >= 3:
                return parts[0]
            elif not isolate and len(parts) >= 2:
                return parts[0]

    print("Error: Could not determine protein name")
    sys.exit(1)


def save_contact_data(contact_data, protein_name, isolate_mode):
    """Save contact data to TSV files"""
    output_dir = os.path.join('../contact_data', 'Isolated' if isolate_mode else 'Not_isolated', protein_name)
    os.makedirs(output_dir, exist_ok=True)

    if isinstance(contact_data, dict):
        # Multiple chains - save each chain separately
        for chain_id, df in contact_data.items():
            output_path = os.path.join(output_dir, f'contacts_chain_{chain_id}.tsv')
            df.to_csv(output_path, sep='\t', index=False)
            print(f"Saved contact data for chain {chain_id} to {output_path}")
    else:
        # Single dataframe
        output_path = os.path.join(output_dir, 'contacts_data.tsv')
        contact_data.to_csv(output_path, sep='\t', index=False)
        print(f"Saved contact data to {output_path}")


def parse_contact_file(file_path):
    """Parse a single contact frustration file"""
    try:
        df = pd.read_csv(
            file_path,
            sep='\s+',
            header=0,
            names=[
                'Res1', 'Res2', 'ChainRes1', 'ChainRes2',
                'DensityRes1', 'DensityRes2', 'AA1', 'AA2',
                'NativeEnergy', 'DecoyEnergy', 'SDEnergy',
                'FrstIndex', 'Welltype', 'FrstState'
            ],
            usecols=[
                'Res1', 'Res2', 'ChainRes1', 'ChainRes2',
                'AA1', 'AA2', 'FrstIndex', 'Welltype', 'FrstState'
            ]
        )

        # Convert residue numbers to integers and normalize if needed
        try:
            # Check first residue of Res1
            first_res = int(df['Res1'].iloc[0])
            if first_res != 1:
                offset = first_res - 1
                df['Res1'] = (df['Res1'].astype(int) - offset).astype(str)
                df['Res2'] = (df['Res2'].astype(int) - offset).astype(str)
        except (ValueError, IndexError):
            # If conversion fails, keep original values
            pass

        return df
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None


def _extract_frame_number(file_path, isolate_mode):
    """Extract frame number from the file path"""
    # Get the parent directory name (the .done directory)
    dir_name = os.path.basename(os.path.dirname(os.path.dirname(file_path)))

    # Remove .done extension
    base_name = dir_name[:-5] if dir_name.endswith('.done') else dir_name

    parts = base_name.split('_')
    if isolate_mode:
        # Format: Protein_Frame_Chain.done
        return parts[-2]  # Frame number is second to last
    else:
        # Format: Protein_Frame.done
        return parts[-1]  # Frame number is last


def process_contact_data(files, isolate_mode):
    """Process all contact files and create consolidated data"""
    all_data = defaultdict(list)

    for file_path in files:
        df = parse_contact_file(file_path)
        if df is not None:
            # Extract frame number from file path
            frame = _extract_frame_number(file_path, isolate_mode)

            # Add frame column
            df['Frame'] = frame

            # Convert chain columns to string
            df['ChainRes1'] = df['ChainRes1'].astype(str)
            df['ChainRes2'] = df['ChainRes2'].astype(str)

            # Create residue IDs (now using normalized residue numbers)
            df['ResID1'] = (
                df['ChainRes1'] + ':' +
                df['AA1'] +
                df['Res1'].astype(str)
            )
            df['ResID2'] = (
                df['ChainRes2'] + ':' +
                df['AA2'] +
                df['Res2'].astype(str)
            )

            if isolate_mode:
                # In isolate mode, all contacts are within the same chain
                chain = df['ChainRes1'].iloc[0]
                all_data[chain].append(df)
            else:
                # In non-isolate mode, we need to separate by chain pairs
                for _, row in df.iterrows():
                    chain_pair = (row['ChainRes1'], row['ChainRes2'])
                    all_data[chain_pair].append(pd.DataFrame([row]))

    # Combine all frames for each chain/pair
    combined = {}
    for key, dfs in all_data.items():
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_key = key[0] if isolate_mode and isinstance(key, tuple) else key
        combined[combined_key] = combined_df

    return combined if len(combined) > 1 else combined.popitem()[1]

def _parse_isolated_contacts(dir_path):
    """Parse contacts in isolate mode (Protein_Frame_Chain.done)"""
    chain_files = defaultdict(list)

    for dirname in os.listdir(dir_path):
        if dirname.endswith('.done'):
            parts = dirname.split('_')
            if len(parts) >= 3:
                chain_name = parts[-1].split('.')[0]
                # Find the actual contact file (name may vary)
                for f in os.listdir(os.path.join(dir_path, dirname, "FrustrationData")):
                    if f.endswith("configurational"):
                        chain_files[chain_name].append(os.path.join(dir_path, dirname, "FrustrationData", f))
                        break

    if not chain_files:
        print("No valid contact files found in isolate mode")
        sys.exit(1)

    print(f"Found {len(chain_files)} chains: {', '.join(chain_files.keys())}")
    return {chain: process_contact_data(files, True) for chain, files in chain_files.items()}


def _parse_non_isolated_contacts(dir_path):
    """Parse contacts in non-isolate mode (Protein_Frame.done)"""
    contact_files = []

    for dirname in os.listdir(dir_path):
        #print(dirname)
        if dirname.endswith('.done'):

            parts = dirname.split('_')
            if len(parts) >= 2:

                # Find the actual contact file

                #print( os.listdir(os.path.join(dir_path, dirname, "FrustrationData")))
                #print( os.path.join(dir_path, dirname, "FrustrationData"))
                for f in os.listdir(os.path.join(dir_path, dirname, "FrustrationData")):
                    #print(f)
                    if f.endswith("configurational"):

                        contact_files.append(os.path.join(dir_path, dirname, "FrustrationData", f))
                        break

    if not contact_files:
        print("No valid contact files found in non-isolate mode")
        sys.exit(1)

    return process_contact_data(contact_files, False)


def parse_contact_results(dir_path, isolate_mode):
    """Main function to parse contact frustration results"""
    if isolate_mode:
        return _parse_isolated_contacts(dir_path)
    else:
        return _parse_non_isolated_contacts(dir_path)


def main():
    dir_path, isolate_mode = parse_arguments()
    protein_name = get_protein_name(dir_path, isolate_mode)

    print(f"\nProcessing contact frustration data for {protein_name}")
    print(f"Isolate mode: {isolate_mode}\n")

    contact_data = parse_contact_results(dir_path, isolate_mode)
    save_contact_data(contact_data, protein_name, isolate_mode)

    # Display summary
    if isinstance(contact_data, dict):
        print(f"\nFound {len(contact_data)} chain(s):")
        for key, df in contact_data.items():
            print(f"{key}: {len(df)} contacts from {df['Frame'].nunique()} frames")
    else:
        print(f"\nFound {len(contact_data)} contacts from {contact_data['Frame'].nunique()} frames")


if __name__ == "__main__":
    main()