# geometry/primarycilia.py


import csv
import pandas as pd
import numpy as np
import sys, os
import random
import copy

from scipy.interpolate import UnivariateSpline

# Define thresholds for Idx_A and Idx_B for each doublet (1-9)
IDX_A_THRESHOLDS = [0.95, 0.92, 0.4, 0.94, 1.0, 0.55, 0.80, 0.65, 0.75]
IDX_B_THRESHOLDS = [0.20, 0.20, 0.35, 0.20, 0.20, 0.22, 0.23, 0.45, 0.20]
MAX_INTERVAL = 20
PC_MEMBRANE_FRACTION = 0.18
PRIMARY_CILIA_LENGTH = 10000
MEMBRANE_RADIUS = 1100
REQUIRED_COLUMNS = [
    "DoubletNumber", "X", "Y", "Z", 
    "Idx_A", "Idx_B", "Angle", 
    "A_Shift", "B_Shift"
]
    
def process_single_doublet(df_doublet: pd.DataFrame, new_length: float, interval: float) -> pd.DataFrame:
    """
    Processes points for a single microtubule doublet:
    1. Scales the Z-values to the new_length.
    2. Calculates the spline curve length parameter (s) based on X, Y, Z.
    3. Fits a smoothing spline to (X, Y, Angle, A_Shift, B_Shift) vs s.
    4. Recalculates points at fixed 'interval' along the new length.
    5. Recalculates Idx_A and Idx_B.
    
    Returns:
        pd.DataFrame: A new DataFrame with resampled and scaled points.
    """
    
    # --- 1. Scale Z to the new_length ---
    max_z_old = df_doublet['Z'].max()
    if max_z_old == 0:
        # Avoid division by zero if all Z are zero
        df_doublet['Z_scaled'] = 0.0
    else:
        scale_factor = new_length / max_z_old
        df_doublet['Z_scaled'] = df_doublet['Z'] * scale_factor
        
    # --- 2. Calculate Cumulative Arc Length (s) ---
    # The true length is the arc length in 3D (X, Y, Z_scaled)
    coords = df_doublet[['X', 'Y', 'Z_scaled']].values
    
    # Calculate segment lengths (distance between consecutive points)
    segment_lengths = np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))
    
    # Cumulative arc length 's' starts at 0
    s_values = np.insert(np.cumsum(segment_lengths), 0, 0.0)
    
    # --- 3. Define New Sample Points (New s_values) ---
    # Calculate the number of new points needed
    num_points = int(np.ceil(new_length / interval))
    
    # The new points will span [0, new_length] at the specified interval
    s_new = np.linspace(0, new_length, num_points + 1)
    
    # --- 4. Fit Splines and Interpolate New Values ---
    new_data = {}
    
    # Columns to be interpolated
    interp_cols = ['X', 'Y', 'Angle', 'A_Shift', 'B_Shift']
    
    # Use the original Z as the input for the Z spline for a smoother fit
    # (Since Z_scaled is derived linearly from Z)
    interp_cols.append('Z_scaled') 
    
    for col in interp_cols:
        # Ensure we have enough unique points for spline fitting
        if len(s_values) > 3 and not df_doublet[col].nunique() == 1:
            # Fit a smoothing spline (k=3 for cubic spline)
            spline = UnivariateSpline(s_values, df_doublet[col], k=3, s=0) 
            new_data[col] = spline(s_new)
        else:
            # If data is constant or too few points, just take the first value
            new_data[col] = np.full_like(s_new, df_doublet[col].iloc[0])

    # --- 5. Assemble and Finalize New DataFrame ---
    df_new = pd.DataFrame(new_data)
    
    # Recalculate Index columns
    df_new['Idx_A'] = np.arange(len(df_new))
    df_new['Idx_B'] = np.arange(len(df_new))
    
    # Constant columns
    df_new['DoubletNumber'] = df_doublet['DoubletNumber'].iloc[0]
    
    # Rename 'Z_scaled' back to 'Z' for the final output
    df_new = df_new.rename(columns={'Z_scaled': 'Z'})
    
    return df_new[REQUIRED_COLUMNS]
    
    
        
def randomize_list_order(input_list: list) -> list:
    """
    Takes a list and returns a new list with the elements shuffled into 
    a random order.

    Args:
        input_list (list): The list of values to be randomized.

    Returns:
        list: A new list with the elements in a random order.
    """
    # Create a shallow copy of the list to avoid modifying the original list
    shuffled_list = copy.copy(input_list)
    
    # Use random.shuffle() to randomize the order of the elements in-place
    random.shuffle(shuffled_list)
    
    return shuffled_list
    
def process_cilia_csv(df, idx_a_thresholds=IDX_A_THRESHOLDS, idx_b_thresholds=IDX_B_THRESHOLDS):
    """
    Process primary cilia CSV file by setting Angle, Idx_A, and Idx_B values.
    
    Parameters:
    input_file (str): Path to input CSV file
    output_file (str): Path to output CSV file
    """
    
    
    # Set Angle values for each DoubletNumber
    # Doublet1 = 90, Doublet2 = 130, Doublet3 = 170, etc. (increment by 40)
    for doublet_num in range(1, 10):
        angle_value = 90 + (doublet_num - 1) * 40
        df.loc[df['DoubletNumber'] == doublet_num, 'Angle'] = angle_value
    
    # Zero out Idx_A and Idx_B columns
    df['Idx_A'] = 0
    df['Idx_B'] = 0
    
    # Find maximum Z value
    max_z = df['Z'].max()
    
    # Process each DoubletNumber
    for doublet_num in range(1, 10):
        idx = doublet_num - 1  # Array index (0-8)
        
        # Calculate Z thresholds
        z_threshold_a = idx_a_thresholds[idx] * max_z
        z_threshold_b = idx_b_thresholds[idx] * max_z
        
        # Set Idx_A to 1 for Z from 0 to z_threshold_a
        mask_a = (df['DoubletNumber'] == doublet_num) & (df['Z'] >= 0) & (df['Z'] <= z_threshold_a)
        df.loc[mask_a, 'Idx_A'] = 1
        
        # Set Idx_B to 1 for Z from 20 to z_threshold_b
        mask_b = (df['DoubletNumber'] == doublet_num) & (df['Z'] >= 20) & (df['Z'] <= z_threshold_b)
        df.loc[mask_b, 'Idx_B'] = 1
        
    # Remove all rows where both Idx_A=0 and Idx_B=0
    df_filtered = df[~((df['Idx_A'] == 0) & (df['Idx_B'] == 0))]
    return df_filtered


def scale_cilia_df(df_in: pd.DataFrame, cilia_length: float, interval: float) -> pd.DataFrame:
    """
    Scales, smooths, and resamples the coordinates for all microtubule doublets 
    in the input DataFrame.

    Args:
        df_in (pd.DataFrame): The input DataFrame with geometry points.
        cilia_length (float): The target maximum Z-value (final length).
        interval (float): The desired spacing between new points.

    Returns:
        pd.DataFrame: A new DataFrame with resampled, scaled, and smoothed points.
    """
    if 'DoubletNumber' not in df_in.columns:
        raise ValueError("Input DataFrame must contain 'DoubletNumber' column for grouping.")

    if df_in['Z'].max() == 0:
        print("Warning: Max Z is zero. Scaling will likely result in a line of zero length.")

    print(f"Starting scaling and smoothing for {df_in['DoubletNumber'].nunique()} doublets...")
    
    # Group by DoubletNumber and apply the processing function to each group
    scaled_df = (
        df_in.groupby('DoubletNumber', group_keys=False)
             .apply(process_single_doublet, new_length=cilia_length, interval=interval)
             .reset_index(drop=True) # Reset indices after groupby/apply
    )
    
    print("Processing complete.")
    return scaled_df

def generate_primary_cilia(
			df_template,
            cilia_length=PRIMARY_CILIA_LENGTH,
            membrane_radius=MEMBRANE_RADIUS,
            membrane_fraction=PC_MEMBRANE_FRACTION,
            interval=MAX_INTERVAL
    ):
    
    idx_a_list = randomize_list_order(IDX_A_THRESHOLDS)
    idx_b_list = randomize_list_order(IDX_B_THRESHOLDS)
    
    df_scaled = scale_cilia_df(df_template, cilia_length, interval)
    df_filtered = process_cilia_csv(df_scaled, idx_a_list, idx_b_list)
    
        # Step 7: Add Membrane (DoubletNumber = 0)
    total_membrane_length = membrane_fraction*cilia_length
    n_membrane_points = int(np.ceil(total_membrane_length / interval)) + 1
    membrane_z_coords = np.linspace(0, total_membrane_length, n_membrane_points)
    
    membrane_data = []        
    for z in membrane_z_coords:
            membrane_data.append([
                0, 0, 0, z, # DoubletNumber=0, X=0, Y=0, Z
                1, 1, 0, # Idx_A=1, Idx_B=1, Angle=0
                0, 0
            ])
            
    # Create Membrane DataFrame and combine with complete data
    membrane_df = pd.DataFrame(membrane_data, columns=REQUIRED_COLUMNS)
    complete_df = pd.concat([df_filtered, membrane_df], ignore_index=True)
    complete_df = complete_df.sort_values(['DoubletNumber', 'Z']).reset_index(drop=True)

    return complete_df