# Not the most obvious io but this write to df allowing to write csv file

import pandas as pd
import numpy as np
import csv
import os, sys
import random # Added: required for random.uniform and random.choice
from . import default_config # Assuming default_config contains necessary constants
from .geometry.tip import generate_multiple_tip_lengths_in_memory

from pathlib import Path


# Global constants
REQUIRED_COLUMNS = [
    "DoubletNumber", "X", "Y", "Z", 
    "Idx_A", "Idx_B", "Angle", 
    "A_Shift", "B_Shift"
]

def load_template_data(template_filename: str):
    # 1. Get the path to the current Python file (where this function is defined)
    #    This ensures we know exactly where 'main_script.py' is.
    script_path = Path(__file__).resolve()

    # 2. Navigate UP one directory (from src/ to CiliaBuilder/)
    #    and then DOWN into the 'template' directory to find the file.
    template_file_path = script_path.parent / 'geometry' / template_filename
    
    # Alternatively, if the template file is right next to the script:
    # template_file_path = script_path.parent / template_filename 

    if not template_file_path.exists():
        raise FileNotFoundError(f"Template file not found at: {template_file_path}")
    
    # 3. Load the data using the absolute path
    df = read_3d_csv(template_file_path)
    
    return df

# Example of how to call the function:
try:
    df_template = load_template_data('initial_template.csv')
    print(f"Successfully loaded template from {df_template.shape[0]} rows.")
except FileNotFoundError as e:
    print(e)
    

def read_3d_csv(file_path):
    """
    Reads a microtubule geometry CSV file and validates the required columns.

    Args:
        file_path (str): The path to the input CSV file.

    Returns:
        pandas.DataFrame: The DataFrame containing the validated data.
    """

    # --- 1. Check if file exists ---
    if not os.path.exists(file_path):
        print(f"Error: File not found at path '{file_path}'", file=sys.stderr)
        return None

    # --- 2. Read the CSV ---
    try:
        df = pd.read_csv(file_path)
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{file_path}' is empty.", file=sys.stderr)
        return None
    except pd.errors.ParserError:
        print(f"Error: Could not parse '{file_path}'. Check if it is a valid CSV format.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}", file=sys.stderr)
        return None

    # --- 3. Validate Columns ---
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if missing_columns:
        print("\n" + "="*50, file=sys.stderr)
        print("Error: The CSV file is missing critical columns.", file=sys.stderr)
        print(f"Missing columns: {missing_columns}", file=sys.stderr)
        print("="*50 + "\n", file=sys.stderr)
        return None
    
    return df[REQUIRED_COLUMNS] # Return only the validated columns
    
def write_3d_csv(df: pd.DataFrame, file_path: str):
    """
    Writes a DataFrame to a CSV file, ensuring only the required 3D geometry 
    columns are included. It prints a warning if the input DataFrame contains 
    extra, unused columns.

    Args:
        df (pd.DataFrame): The DataFrame containing the microtubule geometry data.
        file_path (str): The full path where the CSV file should be saved.
    """
    
    # 1. Check for missing required columns (for robustness)
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_columns:
        error_msg = (
            f"Cannot write CSV: Input DataFrame is missing required columns: "
            f"{missing_columns}"
        )
        # Printing error to stderr
        print(error_msg, file=sys.stderr)
        return

    # 2. Check for extra columns and print a warning
    extra_columns = [col for col in df.columns if col not in REQUIRED_COLUMNS]
    if extra_columns:
        # Print warning to stderr so it's clearly visible as a non-fatal error
        print("="*50, file=sys.stderr)
        print("WARNING: Extra columns found in the input DataFrame.", file=sys.stderr)
        print(f"Ignored columns: {extra_columns}", file=sys.stderr)
        print(f"Only columns: {REQUIRED_COLUMNS} will be written.", file=sys.stderr)
        print("="*50, file=sys.stderr)

    # 3. Select only the required columns and write to CSV
    try:
        # Select columns in the required order
        df_output = df[REQUIRED_COLUMNS]
        df_output['X'] = df_output['X'].map('{:.2f}'.format)
        df_output['Y'] = df_output['Y'].map('{:.2f}'.format)
        df_output['Z'] = df_output['Z'].map('{:.2f}'.format)
        df_output.to_csv(file_path, index=False)
        print(f"\nSuccessfully wrote data to '{file_path}'.")
        print(f"File includes only the required columns.")
        
    except Exception as e:
        # Printing error to stderr
        print(f"Error occurred while writing to CSV: {e}", file=sys.stderr)

# --- 1. Constant Definitions (Readability Improvement) ---
# Use constants from default_config, falling back to assumed values if needed.
# NOTE: The constants INITIAL_LENGTH, TRANSITION_LENGTH, and CILIA_DOUBLET_SHIFT
# used in save_curves_to_csv were not imported but seemed to come from globals/config.
# I am assuming they are meant to be accessed via the default_config object.

INITIAL_LENGTH = getattr(default_config, 'TIP_INITIAL_LENGTH', 300)
TRANSITION_LENGTH = getattr(default_config, 'TIP_TRANSITION_LENGTH', 2000)
CILIA_DOUBLET_SHIFT = getattr(default_config, 'CILIA_DOUBLET_SHIFT', 70)
Z_MAX_IDX_B_DEFAULT = INITIAL_LENGTH + TRANSITION_LENGTH
MAX_INTERVAL = getattr(default_config, 'MAX_INTERVAL', 20)


# --- 2. Helper Function for CSV Row Generation ---
def _generate_doublet_rows(curve_data, doublet_number, z_max_idx_b, doublet_shift, num_doublets=None):
    """
    Internal helper function to process a single curve and generate the CSV rows.
    Separates the data processing logic for better reuse and clarity.
    """
    curve = curve_data['curve']
    all_data = []

    # Calculate Angle: Use 360/num_doublets if provided, otherwise use original logic.
    if num_doublets:
        # Logic from create_mixed_csv_from_memory
        angle = 90 + (360 / num_doublets) * (doublet_number - 1)
    else:
        # Logic from save_curves_to_csv (assuming 9 doublets in that context)
        angle = 90 + 40 * (doublet_number - 1)

    # Define constant values
    IDX_A = 1
    A_SHIFT = -doublet_shift
    B_SHIFT = doublet_shift
    
    # 1. Randomly select the transition Z-coordinate for this specific curve
    idx_b_transition_z = random.uniform(0, z_max_idx_b)
    
    # 2. Iterate over all points and assign Idx_B conditionally
    for x, y, z in curve:
        # Idx_B is 1 if Z is within the randomized region, 0 otherwise
        idx_b = 1 if z <= idx_b_transition_z else 0
        
        # Append data in the required column order
        all_data.append([
            doublet_number, x, y, z, 
            IDX_A, idx_b, angle, 
            A_SHIFT, B_SHIFT
        ])
        
    return all_data

# --- 4. Refactored create_mixed_csv_from_memory Function ---
def create_mixed_csv_from_memory(all_curves_data, output_filename=None, 
                                 initial_length=INITIAL_LENGTH, 
                                 transition_length=TRANSITION_LENGTH, 
                                 num_doublets=default_config.CILIA_NUM_DOUBLETS, 
                                 doublet_shift=CILIA_DOUBLET_SHIFT):
    """
    Create mixed CSV file by randomly selecting curves from different tip lengths.
    Ensures no tip length is chosen more than twice.
    """
    mixed_data = []
    Z_MAX_IDX_B = initial_length + transition_length
    USE_LIMIT = 1
    
    # Create a list to track how many times each tip length has been used
    # Format: [(data_index, count), ...]
    usage_count = {i: 0 for i in range(len(all_curves_data))}
    
    for doublet_number in range(1, num_doublets + 1):
        # Filter available tip lengths that haven't been used twice yet
        available_indices = [idx for idx, count in usage_count.items() if count < USE_LIMIT]
        
        if not available_indices:
            # Shouldn't happen with 9 doublets and 10 tip lengths, but handle edge case
            print(f"Warning: All tip lengths used twice, resetting counts for doublet {doublet_number}")
            usage_count = {i: 0 for i in range(len(all_curves_data))}
            available_indices = list(usage_count.keys())
        
        # Randomly select from available tip lengths
        selected_idx = random.choice(available_indices)
        usage_count[selected_idx] += 1
        
        selected_data = all_curves_data[selected_idx]
        selected_tip_length = selected_data['tip_length']
        selected_curves = selected_data['curves']
        
        print(f"Doublet {doublet_number}: Selected from tip_length={selected_tip_length} (used {usage_count[selected_idx]}/{USE_LIMIT} time)")
        
        curve_data = selected_curves[doublet_number - 1] 
        doublet_rows = _generate_doublet_rows(
            curve_data, doublet_number, Z_MAX_IDX_B, doublet_shift, 
            num_doublets=num_doublets
        )
        mixed_data.extend(doublet_rows)
    
    # Print usage summary
    print("\nTip length usage summary:")
    for idx, count in usage_count.items():
        if count > 0:
            print(f"  Tip length {all_curves_data[idx]['tip_length']:.1f}: used {count} time(s)")
    
    columns = ['DoubletNumber', 'X', 'Y', 'Z', 'Idx_A', 'Idx_B', 'Angle', 'A_Shift', 'B_Shift']
    mixed_df = pd.DataFrame(mixed_data, columns=columns)
    
    if output_filename:
        mixed_df.to_csv(output_filename, index=False)
    
    return mixed_df


# --- 5. generate_tip_csv Function (Cleaning up hardcoded constants) ---
def generate_tip_csv(
    cilia_radius=default_config.CILIA_RADIUS, 
    tip_length_end=default_config.TIP_LENGTH, 
    transition_radius=default_config.TIP_TRANSITION_RADIUS, 
    final_radius=default_config.TIP_FINAL_RADIUS, 
    initial_length=INITIAL_LENGTH, 
    transition_length=TRANSITION_LENGTH, 
    num_doublets=default_config.CILIA_NUM_DOUBLETS, 
    doublet_shift=CILIA_DOUBLET_SHIFT
):
    
    # NOTE: The function generate_multiple_tip_lengths_in_memory is missing, 
    # but the logic flow is preserved.
    tip_length_start = initial_length + transition_length
    
    # This call is assumed to be defined elsewhere (e.g., geometry/tip.py)
    # The 'number_of_steps=10' is a hardcoded constant and could be a parameter
    all_curves_data = generate_multiple_tip_lengths_in_memory( 
        cilia_radius, transition_radius, final_radius,
        tip_length_start, tip_length_end, number_of_steps=10
    )
            
    # Create mixed CSV directly from in-memory data
    mixed_df = create_mixed_csv_from_memory(
        all_curves_data, 
        initial_length=initial_length, 
        transition_length=transition_length, 
        num_doublets=num_doublets, 
        doublet_shift=doublet_shift
    )
    
    return mixed_df
    

def generate_cilia_with_tip(
    cilia_length=10000,
    cilia_radius=default_config.CILIA_RADIUS,
    tip_length=default_config.TIP_LENGTH,
    transition_radius=default_config.TIP_TRANSITION_RADIUS,
    final_radius=default_config.TIP_FINAL_RADIUS,
    initial_length=INITIAL_LENGTH,
    transition_length=TRANSITION_LENGTH,
    num_doublets=default_config.CILIA_NUM_DOUBLETS,
    doublet_shift=CILIA_DOUBLET_SHIFT,
    cp_doublet_length_diff=default_config.CILIA_CP_DOUBLET_LENGTH_DIFF,
    cp_shift=default_config.CILIA_CP_SHIFT,
    membrane_radius=default_config.CILIA_MEMBRANE_RADIUS,
    membrane_fraction=default_config.CILIA_MEMBRANE_FRACTION,
    max_interval=MAX_INTERVAL
):
    """
    Generate complete cilia structure with straight base + tapered tip.
    
    Parameters:
    -----------
    cilia_length : float
        Length of the straight base section (default: 10000)
    cilia_radius : float
        Radius of the cilia at the base (default: from config)
    tip_length_end : float
        End Z-coordinate for the tip section (default: from config)
    transition_radius : float
        Radius at the transition point in the tip (default: from config)
    final_radius : float
        Final radius at the tip end (default: from config)
    output_filename : str or None
        If provided, save the DataFrame to this CSV file
    initial_length : float
        Initial length parameter for tip geometry (default: 300)
    transition_length : float
        Transition length parameter for tip geometry (default: 2000)
    num_doublets : int
        Number of doublets (default: 9)
    doublet_shift : float
        Shift distance for A and B tubules (default: 70)
    max_interval : float
        Maximum interval between points along Z-axis (default: 20)
        
    Returns:
    --------
    pd.DataFrame
        Complete cilia structure data
    """
    
    print("=" * 60)
    print("GENERATING CILIA WITH TIP")
    print("=" * 60)
    print(f"Cilia base length: {cilia_length} Å")
    print(f"Tip length: {tip_length} Å")
    print(f"Total length: {cilia_length + tip_length} Å")
    print("=" * 60)
    
    tip_length_end = tip_length - cp_doublet_length_diff
    
    # Step 1: Generate tip data
    print("\n1. Generating tip geometry...")
    tip_df = generate_tip_csv(
        cilia_radius=cilia_radius,
        tip_length_end=tip_length_end,
        transition_radius=transition_radius,
        final_radius=final_radius,
        initial_length=initial_length,
        transition_length=transition_length,
        num_doublets=num_doublets,
        doublet_shift=doublet_shift
    )
    
    # Step 2: Shift all tip Z-coordinates by cilia_length
    print(f"\n2. Shifting tip Z-coordinates by {cilia_length} Å...")
    tip_df['Z'] = tip_df['Z'] + cilia_length - initial_length
    
    # Step 3: Generate straight base section for each doublet
    print("\n3. Generating straight base sections...")
    base_data = []
    
    for doublet_num in range(1, num_doublets + 1):
        # Get the first point of this doublet from the tip data
        doublet_tip_data = tip_df[tip_df['DoubletNumber'] == doublet_num]
        
        if len(doublet_tip_data) == 0:
            print(f"Warning: No tip data for doublet {doublet_num}, skipping")
            continue
        
        # Get the first point (which is now at Z = cilia_length after shift)
        first_point = doublet_tip_data.iloc[0]
        x1, y1, z1 = first_point['X'], first_point['Y'], first_point['Z']
        angle = first_point['Angle']
        a_shift = first_point['A_Shift']
        b_shift = first_point['B_Shift']
        
        # Calculate number of points needed for the straight section
        n_points = int(np.ceil(cilia_length / max_interval)) + 1
        
        # Generate Z coordinates from 0 to cilia_length (which connects to z1)
        z_coords = np.linspace(0, cilia_length - initial_length, n_points)
        
        # Since it's a straight line, X and Y remain constant
        for i, z in enumerate(z_coords):
            # First point (Z=0): Idx_A=1, Idx_B=0
            # All other points: Idx_A=1, Idx_B=1
            if i == 0:
                idx_a, idx_b = 1, 0
            else:
                idx_a, idx_b = 1, 1
            
            base_data.append([
                doublet_num, x1, y1, z,
                idx_a, idx_b, angle,
                a_shift, b_shift
            ])
        
        print(f"  Doublet {doublet_num}: Generated {n_points} base points from Z=0 to Z={cilia_length}")
    
    # Step 4: Create base DataFrame
    columns = ['DoubletNumber', 'X', 'Y', 'Z', 'Idx_A', 'Idx_B', 'Angle', 'A_Shift', 'B_Shift']
    base_df = pd.DataFrame(base_data, columns=columns)
    
    # Step 5: Combine base and tip (remove duplicate point at Z=cilia_length from tip)
    print("\n4. Combining base and tip sections...")
    # Remove the first point of each doublet in tip_df (which is at Z=cilia_length, duplicated)
    tip_df_filtered = tip_df.groupby('DoubletNumber').apply(
        lambda group: group.iloc[1:]  # Skip first row of each group
    ).reset_index(drop=True)
    
    # Concatenate base and filtered tip
    complete_data = pd.concat([base_df, tip_df_filtered], ignore_index=True)
        
    # Step 6: Add Central Pair (DoubletNumber = -1)
    print("\n5. Generating central pair...")
    total_length = cilia_length + tip_length
    n_cp_points = int(np.ceil(total_length / max_interval)) + 1
    cp_z_coords = np.linspace(0, total_length, n_cp_points)
    
    cp_data = []
    for z in cp_z_coords:
        cp_data.append([
            -1, 0, 0, z,  # DoubletNumber=-1, X=0, Y=0, Z
            1, 1, 0,       # Idx_A=1, Idx_B=1, Angle=0
            cp_shift, -cp_shift  # A_Shift=cp_shift, B_Shift=-cp_shift
        ])
        

    # Create Cap DataFrame and combine with complete data (DoubletNumber = -2)
    cp_data.append([
            -2, 0, 0, cp_z_coords[-1],  # DoubletNumber=-2, X=0, Y=0, Z
            1, 1, 0,       # Idx_A=1, Idx_B=1, Angle=0
            0, 0  # A_Shift=cp_shift, B_Shift=-cp_shift
        ])
    
    # Create CP DataFrame and combine with complete data
    cp_df = pd.DataFrame(cp_data, columns=columns)
    
    # Step 7: Add Membrane (DoubletNumber = 0)
    total_membrane_length = membrane_fraction*cilia_length
    n_membrane_points = int(np.ceil(total_membrane_length / max_interval)) + 1
    membrane_z_coords = np.linspace(0, total_membrane_length, n_membrane_points)
    
    membrane_data = []        
    for z in membrane_z_coords:
            membrane_data.append([
                0, 0, 0, z, # DoubletNumber=0, X=0, Y=0, Z
                1, 1, 0, # Idx_A=1, Idx_B=1, Angle=0
                0, 0
            ])
            
    # Create Membrane DataFrame and combine with complete data
    membrane_df = pd.DataFrame(membrane_data, columns=columns)
    complete_df = pd.concat([complete_data, cp_df, membrane_df], ignore_index=True)
    complete_df = complete_df.sort_values(['DoubletNumber', 'Z']).reset_index(drop=True)

        
    return complete_df



# Example usage (uncomment and fix dependency for testing)
# generate_tip_csv(cilia_radius=875, tip_length_end=5000, transition_radius=656, final_radius=460, output_filename="test_tip.csv")