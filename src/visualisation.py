import pandas as pd
import sys
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches


def parse_arguments():
    """
    Parse command line arguments including the chain selection and TCL script generation
    """
    parser = argparse.ArgumentParser(description='Visualize frustration data.')
    parser.add_argument('dataframes_dir', type=str,
                      help='Directory containing the frustration dataframes')
    parser.add_argument('--chain', type=str, default=None,
                      help='Specific chain to analyze (e.g., "0", "A"). If not specified, uses first chain.')
    parser.add_argument('--generate-tcl', action='store_true', default=False,
                      help='Generate TCL script for visualization in VMD')
    parser.add_argument('--pdb-dir', type=str, default=None,
                      help='Directory containing PDB files')
    parser.add_argument('--frame-step', type=int, default=10,
                      help='Step between frames in the animation (default: 10)')
    parser.add_argument('--max-frames', type=int, default=None,
                      help='Maximum number of frames to include in the TCL script (default: all frames)')
    parser.add_argument('--residue', type=int, default=None,
                      help='Specific residue number to analyse')
    parser.add_argument('--boxplot', action='store_true', default=False,
                      help='Generate boxplot of frustration values per residue')
    parser.add_argument('--variability', action='store_true', default=False,
                      help='Show variability across chains')
    parser.add_argument('--frame', type=int, default=None,
                      help='Specific frame number to analyze for boxplot/variability')
    parser.add_argument('--dynamic', action='store_true', default=False,
                        help='Generate dynamic boxplot of frustration values across all frames for a chain')


    args = parser.parse_args()

    # Verify that the directory exists
    if not os.path.isdir(args.dataframes_dir):
        print(f"Error: Directory {args.dataframes_dir} does not exist")
        sys.exit(1)

    return (args.dataframes_dir, args.chain, args.generate_tcl, args.pdb_dir,
            args.frame_step, args.max_frames, args.residue, args.boxplot,
            args.variability, args.frame, args.dynamic)


def load_dataframes(dataframes_dir):
    """
    Load frustration dataframes from the specified directory
    Returns a dictionary of dataframes (one per chain) or a single dataframe
    """
    frustration_data = {}

    # Check if there's a single frustration_data.csv file
    single_df_path = os.path.join(dataframes_dir, 'frustration_data.csv')
    if os.path.exists(single_df_path):
        return pd.read_csv(single_df_path)

    # Otherwise look for chain-specific files
    for filename in os.listdir(dataframes_dir):
        if filename.startswith('frustration_chain_') and filename.endswith('.csv'):
            chain_id = filename.split('_')[-1].split('.')[0]
            df_path = os.path.join(dataframes_dir, filename)
            frustration_data[chain_id] = pd.read_csv(df_path)

    if not frustration_data:
        print(f"Error: No frustration dataframes found in {dataframes_dir}")
        sys.exit(1)

    return frustration_data


def classify_residues(df_chain):
    """
    Classify residues for each frame into three categories:
    - minimally frustrated (> 0.58)
    - highly frustrated (< -1.0)
    - neutral (between -1.0 and 0.58)
    """
    frustration_data = {}

    # Get all frame columns (assuming they start with 'frame')
    frame_columns = [col for col in df_chain.columns if col.startswith('frame')]

    for frame in frame_columns:
        frame_data = {
            'minimally_frustrated': [],
            'highly_frustrated': [],
            'neutral': []
        }

        for _, row in df_chain.iterrows():
            residue = row['residue']
            value = row[frame]

            if value > 0.58:
                frame_data['minimally_frustrated'].append(residue)
            elif value < -1.0:
                frame_data['highly_frustrated'].append(residue)
            else:
                frame_data['neutral'].append(residue)

        frustration_data[frame] = frame_data

    return frustration_data


def select_chain(frustration_data, chain_id):
    """
    Select a specific chain from the frustration data
    """
    if isinstance(frustration_data, dict):
        if chain_id in frustration_data:
            return frustration_data[chain_id]
        else:
            available_chains = list(frustration_data.keys())
            print(f"\nError: Chain '{chain_id}' not found. Available chains: {available_chains}")
            sys.exit(1)
    else:
        if chain_id is not None and frustration_data['chain'].iloc[0] != chain_id:
            print(f"\nError: Only chain '{frustration_data['chain'].iloc[0]}' available")
            sys.exit(1)
        return frustration_data


def generate_tcl_script(classified_residues, pdb_dir, frame_step, max_frames=None,
                        output_path="frustration_visualization.tcl"):
    """
    Generate a TCL script to visualize frustration in VMD with each frame as a separate molecule
    colored according to its specific frustration data
    """
    with open(output_path, 'w') as tcl_file:
        # Write header
        tcl_file.write("# Frustration Visualization Script\n")
        tcl_file.write("# Generated automatically from frustration data\n\n")
        tcl_file.write("color change rgb grey 0.8 0.8 0.8\n")
        tcl_file.write("color Display Background white\n")
        tcl_file.write("display projection Orthographic\n")
        tcl_file.write("axes location off\n")
        tcl_file.write("rotate x by -90\n")
        tcl_file.write("rotate y by -90\n")
        tcl_file.write("display resetview\n\n")

        # Function to color a single molecule
        tcl_file.write("# Function to color a molecule based on frustration data\n")
        tcl_file.write("proc color_molecule {molid min_frustrated high_frustrated neutral} {\n")
        tcl_file.write("    # Clear all existing representations\n")
        tcl_file.write("    set num_reps [molinfo $molid get numreps]\n")
        tcl_file.write("    for {set i [expr $num_reps-1]} {$i >= 0} {incr i -1} {\n")
        tcl_file.write("        mol delrep $i $molid\n")
        tcl_file.write("    }\n\n")
        tcl_file.write("    # Base representation for the whole protein (light blue)\n")
        tcl_file.write("    mol representation NewCartoon\n")
        tcl_file.write("    mol color ColorID 7\n")
        tcl_file.write("    mol selection \"all\"\n")
        tcl_file.write("    mol material Opaque\n")
        tcl_file.write("    mol addrep $molid\n\n")
        tcl_file.write("    # Create selection text for each group\n")
        tcl_file.write("    set sel_text_min \"\"\n")
        tcl_file.write("    foreach r $min_frustrated {\n")
        tcl_file.write("        append sel_text_min \"resid $r or \"\n")
        tcl_file.write("    }\n")
        tcl_file.write("    set sel_text_min [string range $sel_text_min 0 end-4]\n\n")
        tcl_file.write("    set sel_text_high \"\"\n")
        tcl_file.write("    foreach r $high_frustrated {\n")
        tcl_file.write("        append sel_text_high \"resid $r or \"\n")
        tcl_file.write("    }\n")
        tcl_file.write("    set sel_text_high [string range $sel_text_high 0 end-4]\n\n")
        tcl_file.write("    set sel_text_neutral \"\"\n")
        tcl_file.write("    foreach r $neutral {\n")
        tcl_file.write("        append sel_text_neutral \"resid $r or \"\n")
        tcl_file.write("    }\n")
        tcl_file.write("    set sel_text_neutral [string range $sel_text_neutral 0 end-4]\n\n")
        tcl_file.write("    # Add representations for each group\n")
        tcl_file.write("    if {$sel_text_min ne \"\"} {\n")
        tcl_file.write("        mol representation NewCartoon\n")
        tcl_file.write("        mol selection $sel_text_min\n")
        tcl_file.write("        mol color ColorID 19  ;# Green\n")
        tcl_file.write("        mol material Opaque\n")
        tcl_file.write("        mol addrep $molid\n")
        tcl_file.write("    }\n\n")
        tcl_file.write("    if {$sel_text_high ne \"\"} {\n")
        tcl_file.write("        mol representation NewCartoon\n")
        tcl_file.write("        mol selection $sel_text_high\n")
        tcl_file.write("        mol color ColorID 1   ;# Red\n")
        tcl_file.write("        mol material Opaque\n")
        tcl_file.write("        mol addrep $molid\n")
        tcl_file.write("    }\n\n")
        tcl_file.write("    if {$sel_text_neutral ne \"\"} {\n")
        tcl_file.write("        mol representation NewCartoon\n")
        tcl_file.write("        mol selection $sel_text_neutral\n")
        tcl_file.write("        mol color ColorID 6   ;# Grey\n")
        tcl_file.write("        mol material Opaque\n")
        tcl_file.write("        mol addrep $molid\n")
        tcl_file.write("    }\n")
        tcl_file.write("}\n\n")

        # Process each frame and load as separate molecule
        frames_processed = 0
        for frame, data in classified_residues.items():
            if max_frames is not None and frames_processed >= max_frames:
                break

            frame_num = frame[5:]  # Remove 'frame' prefix
            pdb_file = f"frame_{frame_num}.pdb"
            pdb_path = os.path.join(pdb_dir, pdb_file)

            if not os.path.exists(pdb_path):
                print(f"Warning: PDB file not found for frame {frame}: {pdb_path}")
                continue

            # Convert residue names to numbers (M1 -> 1, K4 -> 4, etc.)
            min_frustrated = [int(res[1:]) + 3 for res in data['minimally_frustrated']]
            high_frustrated = [int(res[1:]) + 3 for res in data['highly_frustrated']]
            neutral = [int(res[1:]) + 3 for res in data['neutral']]

            # Load and color this molecule
            tcl_file.write(f"# Loading and coloring frame {frame_num}\n")
            tcl_file.write(f"mol new {pdb_file} type pdb waitfor all\n")
            tcl_file.write(f"mol rename top \"Frame {frame_num}\"\n")
            tcl_file.write(f"color_molecule top [list {' '.join(map(str, min_frustrated))}] \\\n")
            tcl_file.write(f"    [list {' '.join(map(str, high_frustrated))}] \\\n")
            tcl_file.write(f"    [list {' '.join(map(str, neutral))}]\n\n")

            frames_processed += 1

    print(f"\nTCL visualization script generated at: {output_path}")
    print(f"Included {frames_processed} frames" +
          (f" (limited to {max_frames} frames)" if max_frames is not None else ""))
    print("Run in VMD with: vmd -e frustration_visualization.tcl")


def plot_frustration_vs_frames(frust_frames_data, frames, residue, chain_id, dataframe_dir):
    """
    Plot frustration values across frames with specified coloring.

    Args:
        frust_frames_data: List of frustration values (as numpy.float64 or float)
        frames: List of frame names (e.g., ['frame0', 'frame10', ...])
    """
    #extract Protein_name and type
    protein_name = dataframe_dir.split('/')[-2]
    Type = dataframe_dir.split('/')[-3]
    # Extract numerical values from frame names
    frame_numbers = [int(f[5:]) for f in frames]

    # Convert frustration data to simple floats if they're numpy types
    frustration_values = [float(val) for val in frust_frames_data]

    # Create figure
    plt.figure(figsize=(24, 4))

    # Plot lines connecting all points (thin line)
    plt.plot(frame_numbers, frustration_values,
             color='gray', linewidth=0.5, linestyle='-', zorder=1)

    # Lignes horizontales
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)  # Ligne noire en pointillés à y=0
    plt.axhline(y=-1, color='r', linestyle=':', alpha=0.5)  # Ligne rouge en pointillés à y=-1
    plt.axhline(y=0.58, color='g', linestyle=':', alpha=0.5)  # Ligne verte en pointillés à y=-0.58

    # Plot each point with appropriate color
    for x, y in zip(frame_numbers, frustration_values):
        if y > 0.58:
            color = 'green'  # Minimally frustrated
        elif y < -1.0:
            color = 'red'  # Highly frustrated
        else:
            color = 'gray'  # Neutral

        plt.scatter(x, y, color=color, s=15, zorder=2)

    # Set Y-axis limits
    plt.ylim(-4, 2)

    # Customize plot
    plt.gca().invert_yaxis()  # Invert y-axis as requested
    plt.xlabel('Frame Number')
    plt.ylabel('Frustration Value')
    plt.title(f'Frustration Value of residue {residue} Across Frames ({protein_name}, chain {chain_id}, {Type})')
    plt.grid(True, linestyle='--', alpha=0.6)

    # Add legend
    plt.scatter([], [], color='green', label='Minimally frustrated (>0.58)')
    plt.scatter([], [], color='red', label='Highly frustrated (<-1.0)')
    plt.scatter([], [], color='gray', label='Neutral')
    plt.legend()

    plt.tight_layout()
    plt.show()


def plot_all_chains_frustration(frustration_data, residue_num):
    """
    Plot frustration values for a specific residue across all chains in a specific order
    with chain labels converted to monomer numbers. All subplots have the same Y-axis scale (-3 to 1).
    """
    # Define the custom order of chains
    custom_order = [
        '0', '5', 'l', 'q', 'v',
        'J', 'F', 'G', 'H', 'I',
        '1', '6', 'm', 'r', 'w',
        '2', '7', 'n', 's', 'x',
        '3', 'j', 'o', 't', 'y',
        '4', 'k', 'p', 'u', 'z',
        'A', 'B', 'C', 'D', 'E',
        'K', 'P', 'U', 'Z', 'e',
        'L', 'Q', 'V', 'a', 'f',
        'M', 'R', 'W', 'b', 'g',
        'N', 'S', 'X', 'c', 'h',
        'O', 'T', 'Y', 'd', 'i'
    ]

    # Define the monomer labels mapping
    monomer_labels = {
        1: '0', 2: '1', 3: '2', 4: '3', 5: '4', 6: '5', 7: '6', 8: '7',
        9: 'A', 10: 'B', 11: 'C', 12: 'D', 13: 'E', 14: 'F', 15: 'G', 16: 'H',
        17: 'I', 18: 'J', 19: 'K', 20: 'L', 21: 'M', 22: 'N', 23: 'O', 24: 'P',
        25: 'Q', 26: 'R', 27: 'S', 28: 'T', 29: 'U', 30: 'V', 31: 'W', 32: 'X',
        33: 'Y', 34: 'Z', 35: 'a', 36: 'b', 37: 'c', 38: 'd', 39: 'e', 40: 'f',
        41: 'g', 42: 'h', 43: 'i', 44: 'j', 45: 'k', 46: 'l', 47: 'm', 48: 'n',
        49: 'o', 50: 'p', 51: 'q', 52: 'r', 53: 's', 54: 't', 55: 'u', 56: 'v',
        57: 'w', 58: 'x', 59: 'y', 60: 'z'
    }

    # Create reverse mapping (label -> number)
    label_to_num = {v: k for k, v in monomer_labels.items()}

    # Determine available chains
    if isinstance(frustration_data, dict):
        available_chains = list(frustration_data.keys())
    else:
        available_chains = [frustration_data['chain'].iloc[0]]

    # Filter and sort chains according to custom_order
    chains = [chain for chain in custom_order if chain in available_chains]

    # Add any remaining chains not in custom_order (sorted alphabetically)
    remaining_chains = sorted(set(available_chains) - set(custom_order))
    chains.extend(remaining_chains)

    # Create figure with subplots (12 rows x 5 columns)
    rows = 12
    cols = 5
    fig, axs = plt.subplots(rows, cols, figsize=(25, 30))
    fig.suptitle(f'Frustration Index for residue {residue_num} across all chains', y=1.02)

    # Set common Y-axis limits for all subplots
    y_min, y_max = -4, 2

    # Plot each chain in custom order
    for i, chain_label in enumerate(custom_order):
        if chain_label not in available_chains:
            continue  # Skip chains that don't exist in our data

        # Get the monomer number
        monomer_num = label_to_num.get(chain_label, '?')

        # Get the dataframe for this chain
        if isinstance(frustration_data, dict):
            df_chain = frustration_data[chain_label]
        else:
            df_chain = frustration_data

        # Get the row index (residue_num - 1)
        row_idx = residue_num - 1
        if row_idx >= len(df_chain):
            print(f"Warning: Residue {residue_num} not found in chain {chain_label}")
            continue

        # Get data for this residue
        row_data = df_chain.iloc[row_idx]
        frust_values = row_data.values[2:]  # Skip chain and residue columns
        frames = row_data.index[2:]  # Skip chain and residue columns

        # Extract numerical frame numbers
        frame_numbers = [int(f[5:]) for f in frames]

        # Get current axis
        ax = axs[i // cols, i % cols]

        # Plot each point with appropriate color
        for x, y in zip(frame_numbers, frust_values):
            if y > 0.58:
                color = 'green'  # Minimally frustrated
            elif y < -1.0:
                color = 'red'  # Highly frustrated
            else:
                color = 'gray'  # Neutral

            ax.plot(x, y, 'o', color=color, markersize=1.5)

        # Customize plot
        ax.set_ylim(y_min, y_max)
        ax.invert_yaxis()
        ax.set_title(f'M{monomer_num}', fontsize=8, pad=2)
        ax.grid(True, linestyle=':', alpha=0.1)

        # Add horizontal lines
        ax.axhline(y=0, color='k', linestyle=':', alpha=0.3)
        ax.axhline(y=-1, color='r', linestyle=':', alpha=0.3, linewidth=0.5)
        ax.axhline(y=0.58, color='g', linestyle=':', alpha=0.3, linewidth=0.5)

        # Remove y-axis label
        ax.set_ylabel('')
        ax.set_yticklabels([])

        # Set x labels only for bottom row
        if i < (rows - 1) * cols:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Frame')

    # Hide unused axes
    for j in range(len(custom_order), rows * cols):
        axs[j // cols, j % cols].axis('off')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    plt.show()


def plot_frustration_boxplots(frustration_data, frame_num, output_dir, chain_id , dataframe_dir):
    """
    Generate boxplots of frustration values per residue across all chains for a specific frame.
    Box color depends on the distribution:
      - Red: all values < -1
      - Green: all values > 0.58
      - Gray: all values between -1 and 0.58
      - Blue: mixed values
    Y-axis is inverted (negative up).

    Parameters:
    - frustration_data: Dict of DataFrames (one per chain) or single DataFrame
    - frame_num: Frame number to analyze
    - output_dir: Directory to save the plot (if None, shows plot)
    """


    # extract Protein_name and type
    protein_name = dataframe_dir.split('/')[-2]
    Type = dataframe_dir.split('/')[-3]
    # Prepare data
    frame_col = f'frame{frame_num}'
    all_data = []

    if isinstance(frustration_data, dict):
        # Multiple chains case
        for chain_id, df in frustration_data.items():
            temp_df = df[['residue', frame_col]].copy()
            temp_df['chain'] = chain_id
            all_data.append(temp_df)
    else:
        # Single chain case
        temp_df = frustration_data[['residue', frame_col]].copy()
        temp_df['chain'] = frustration_data['chain'].iloc[0]
        all_data.append(temp_df)

    combined_df = pd.concat(all_data)

    plt.figure(figsize=(16, 6))
    ax = sns.boxplot(
        data=combined_df,
        x='residue',
        y=frame_col,
        showfliers=True,
        width=0.75,
        linewidth=1,
        fliersize=2,
        color='skyblue'
    )

    # Color each box based on values
    residue_list = combined_df['residue'].unique()
    for i, residue in enumerate(residue_list):
        values = combined_df[combined_df['residue'] == residue][frame_col].values
        all_below = all(val < -1 for val in values)
        all_above = all(val > 0.58 for val in values)
        all_middle = all(-1 <= val <= 0.58 for val in values)

        if all_below:
            ax.patches[i].set_facecolor('red')
            ax.patches[i].set_edgecolor('darkred')
        elif all_above:
            ax.patches[i].set_facecolor('green')
            ax.patches[i].set_edgecolor('darkgreen')
        elif all_middle:
            ax.patches[i].set_facecolor('gray')
            ax.patches[i].set_edgecolor('black')
        else:
            ax.patches[i].set_facecolor('skyblue')
            ax.patches[i].set_edgecolor('steelblue')

    # Invert Y-axis and set limits
    plt.ylim(-3.5, 3.5)
    ax.invert_yaxis()

    # Horizontal reference lines
    plt.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    plt.axhline(-1, color='red', linestyle=':', linewidth=1, alpha=0.5)
    plt.axhline(0.58, color='green', linestyle=':', linewidth=1, alpha=0.5)

    # X-axis formatting
    xticks = ax.get_xticks()
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    ax.set_xticks(xticks[::3])
    ax.set_xticklabels(xlabels[::3], rotation=45, ha='right', fontsize=7)

    # Vertical grid lines
    for x in xticks:
        plt.axvline(x=x, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)

    # Title and labels
    plt.title(f'Frustration Distribution per Residue (Frame {frame_num}, all chains, {protein_name}, {Type})', fontsize=12)
    plt.xlabel('Residue', labelpad=10)
    plt.ylabel('Frustration Index', labelpad=10)
    plt.yticks(fontsize=10)
    plt.grid(axis='y', alpha=0.3)

    # Legend
    red_patch = mpatches.Patch(color='red', label='All values < -1')
    green_patch = mpatches.Patch(color='green', label='All values > 0.58')
    gray_patch = mpatches.Patch(color='gray', label='All between -1 and 0.58')
    blue_patch = mpatches.Patch(color='skyblue', label='Mixed values')
    plt.legend(handles=[red_patch, green_patch, gray_patch, blue_patch], fontsize=8)

    plt.tight_layout()
    plt.show()
    #if output_dir:
    #    os.makedirs(output_dir, exist_ok=True)
    #    plot_path = os.path.join(output_dir, f"frustration_boxplot_frame_{frame_num}.png")
    #    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    #    plt.close()
    #    print(f"Boxplot saved to {plot_path}")
    #else:
    #    plt.show()

def plot_dynamic_frustration_boxplots(frustration_data, output_dir, chain_id, dataframe_dir, high_residue):
    """
    Generate boxplots of frustration values per residue across all frames for a specific chain.
    Box color depends on the distribution:
      - Red: all values < -1
      - Green: all values > 0.58
      - Gray: all values between -1 and 0.58
      - Blue: mixed values
    Y-axis is inverted (negative up). Residue  is highlighted with red background.

    Parameters:
    - frustration_data: Dict of DataFrames (one per chain) or single DataFrame
    - chain_id: Specific chain to analyze (if None, uses first available chain)
    - output_dir: Directory to save the plot (if None, shows plot)
    """


    # extract Protein_name and type
    protein_name = dataframe_dir.split('/')[-2]
    Type = dataframe_dir.split('/')[-3]
    # Select the chain to analyze
    if isinstance(frustration_data, dict):
        if chain_id is None:
            chain_id = list(frustration_data.keys())[0]
        df_chain = frustration_data[chain_id]
    else:
        df_chain = frustration_data

    # Get all frame columns
    frame_columns = [col for col in df_chain.columns if col.startswith('frame')]

    # Prepare data for plotting
    plot_data = []
    for _, row in df_chain.iterrows():
        residue = row['residue']
        for frame in frame_columns:
            plot_data.append({
                'residue': residue,
                'frustration': row[frame],
                'frame': int(frame[5:])  # Extract frame number
            })

    plot_df = pd.DataFrame(plot_data)

    plt.figure(figsize=(16, 6))
    ax = sns.boxplot(
        data=plot_df,
        x='residue',
        y='frustration',
        showfliers=True,
        width=0.75,
        linewidth=1,
        fliersize=2,
        color='skyblue'
    )

    # Highlight residue
    highlight_suffixes = [str(high_residue)]

    residue_list = plot_df['residue'].unique()
    for i, res in enumerate(residue_list):
        if any(res[1:]==suffix for suffix in highlight_suffixes):
            high_res = res
            ax.axvspan(i - 0.5, i + 0.5, color='lightcoral', alpha=0.3, zorder=0)

    # Color each box based on values
    for i, residue in enumerate(residue_list):
        values = plot_df[plot_df['residue'] == residue]['frustration'].values
        all_below = all(val < -1 for val in values)
        all_above = all(val > 0.58 for val in values)
        all_middle = all(-1 <= val <= 0.58 for val in values)

        if all_below:
            ax.patches[i].set_facecolor('red')
            ax.patches[i].set_edgecolor('darkred')
        elif all_above:
            ax.patches[i].set_facecolor('green')
            ax.patches[i].set_edgecolor('darkgreen')
        elif all_middle:
            ax.patches[i].set_facecolor('gray')
            ax.patches[i].set_edgecolor('black')
        else:
            ax.patches[i].set_facecolor('skyblue')
            ax.patches[i].set_edgecolor('steelblue')

    # Invert Y-axis and set limits
    plt.ylim(-3.5, 3.5)
    ax.invert_yaxis()

    # Horizontal reference lines
    plt.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    plt.axhline(-1, color='red', linestyle=':', linewidth=1, alpha=0.5)
    plt.axhline(0.58, color='green', linestyle=':', linewidth=1, alpha=0.5)

    # X-axis formatting
    xticks = ax.get_xticks()
    xlabels = [label.get_text() for label in ax.get_xticklabels()]
    ax.set_xticks(xticks[::3])
    ax.set_xticklabels(xlabels[::3], rotation=45, ha='right', fontsize=7)

    # Vertical grid lines
    for x in xticks:
        plt.axvline(x=x, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)

    # Title and labels
    plt.title(f'Dynamic Frustration Distribution per Residue, highlight residue {high_res} (All Frames, Chain {chain_id}, {protein_name}, {Type})', fontsize=12)
    plt.xlabel('Residue', labelpad=10)
    plt.ylabel('Frustration Index', labelpad=10)
    plt.yticks(fontsize=10)
    plt.grid(axis='y', alpha=0.3)

    # Legend
    red_patch = mpatches.Patch(color='red', label='All values < -1')
    green_patch = mpatches.Patch(color='green', label='All values > 0.58')
    gray_patch = mpatches.Patch(color='gray', label='All between -1 and 0.58')
    blue_patch = mpatches.Patch(color='skyblue', label='Mixed values')
    plt.legend(handles=[red_patch, green_patch, gray_patch, blue_patch], fontsize=8)

    plt.tight_layout()
    plt.show()
    #if output_dir:
    #    os.makedirs(output_dir, exist_ok=True)
    #    plot_path = os.path.join(output_dir, f"dynamic_frustration_boxplot_chain_{chain_id}.png")
    #    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    #    plt.close()
    #    print(f"Dynamic boxplot saved to {plot_path}")
    #else:
    #    plt.show()

def main():
    # Parse command line arguments with new parameters
    (dataframes_dir, chain_id, generate_tcl, pdb_dir, frame_step,
     max_frames, residue, boxplot, variability, frame_num, dynamic) = parse_arguments()

    print(f"\nLoading dataframes from: {dataframes_dir}")
    if chain_id is not None:
        print(f"Selected chain: {chain_id}")
    else:
        print("No specific chain selected - using first available chain")
    if max_frames is not None:
        print(f"Maximum frames in TCL script: {max_frames}")

    # Load frustration data from dataframes
    frustration_data = load_dataframes(dataframes_dir)

    # Select the chain to analyze
    if chain_id is not None:
        df_chain = select_chain(frustration_data, chain_id)
    else:
        if isinstance(frustration_data, dict):
            df_chain = list(frustration_data.values())[0]
        else:
            df_chain = frustration_data

    print(f"\nAnalyzing chain {df_chain['chain'].iloc[0]} with shape {df_chain.shape}")

    # Classify residues for each frame
    classified_residues = classify_residues(df_chain)

    # Print some statistics
    print("\nClassification results:")
    frames_to_show = min(5, max_frames) if max_frames is not None else 5
    for frame, data in list(classified_residues.items())[:frames_to_show]:
        print(f"\nFrame {frame}:")
        print(f"Minimally frustrated residues (>0.58): {len(data['minimally_frustrated'])}")
        print(f"Highly frustrated residues (<-1.0): {len(data['highly_frustrated'])}")
        print(f"Neutral residues: {len(data['neutral'])}")

    # Handle boxplot/variability analysis
    if boxplot :

        if variability:
            if frame_num is None:
                print("\nError: --frame parameter must be specified with --boxplot and --variability")
                sys.exit(1)
            plot_frustration_boxplots(frustration_data, frame_num, "plots/boxplots", chain_id, dataframes_dir)
        if dynamic:
            if chain_id is None:
                print("\nError: --chain parameter must be specified with --boxplot and --dynamic")
                sys.exit(1)
            plot_dynamic_frustration_boxplots(frustration_data, "plots/dynamic_boxplots", chain_id, dataframes_dir, residue )

    # Generate TCL script if requested
    if generate_tcl:
        if pdb_dir is None:
            print("\nError: PDB directory must be specified with --pdb-dir when generating TCL script")
            sys.exit(1)

        if not os.path.isdir(pdb_dir):
            print(f"\nError: PDB directory {pdb_dir} does not exist")
            sys.exit(1)

        # Generate TCL script
        generate_tcl_script(classified_residues, pdb_dir, frame_step, max_frames)

    # Plot frustration for specific residue if requested
    if residue is not None:
        if chain_id is None and isinstance(frustration_data, dict):
            # Plot for all chains if no specific chain is selected and we have multiple chains
            plot_all_chains_frustration(frustration_data, residue)
        else:
            # Plot for specific chain (or single chain case)
            if isinstance(frustration_data, dict):
                if chain_id is None:
                    chain_id = list(frustration_data.keys())[0]
                df_chain = frustration_data[chain_id]

            print(f"\nData for residue {residue}:")
            print(df_chain.iloc[residue - 1])

            frust_frames_data = df_chain.iloc[residue - 1].values[2:]  # Skip chain and residue columns
            frames = df_chain.columns.tolist()[2:]  # Skip chain and residue columns

            plot_frustration_vs_frames(frust_frames_data, frames, residue, chain_id, dataframes_dir )

if __name__ == "__main__":
    main()