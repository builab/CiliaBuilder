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

def process_cilia_csv(df_scaled, idx_a_list, idx_b_list):
    """
    Applies Idx_A and Idx_B length filtering based on thresholds.
    """
    
    # Step 5: Filter points based on randomized length thresholds
    df_filtered = []
    doublet_numbers = sorted(df_scaled['DoubletNumber'].unique())
    
    for i, doublet_num in enumerate(doublet_numbers):
        df_doublet = df_scaled[df_scaled['DoubletNumber'] == doublet_num].copy()
        
        # Calculate maximum Z to establish length percentage
        max_z = df_doublet['Z'].max()
        
        # Angle
        angle_value = 90 + (doublet_num - 1) * 40
        df_doublet['Angle'] = angle_value
        
        # Determine Z-thresholds for Idx_A and Idx_B
        idx_a_threshold_z = max_z * idx_a_list[i % len(idx_a_list)]
        idx_b_threshold_z = max_z * idx_b_list[i % len(idx_b_list)]
        
        # Apply filtering
        df_doublet.loc[df_doublet['Z'] < idx_a_threshold_z, 'Idx_A'] = 1
        df_doublet.loc[df_doublet['Z'] >= idx_a_threshold_z, 'Idx_A'] = 0
        
        df_doublet.loc[df_doublet['Z'] < idx_b_threshold_z, 'Idx_B'] = 1
        df_doublet.loc[df_doublet['Z'] >= idx_b_threshold_z, 'Idx_B'] = 0
        
        df_filtered.append(df_doublet)

    return pd.concat(df_filtered, ignore_index=True) if df_filtered else pd.DataFrame(columns=REQUIRED_COLUMNS)


def generate_primary_cilia(
			df_template,
            cilia_length=PRIMARY_CILIA_LENGTH,
            membrane_radius=MEMBRANE_RADIUS,
            membrane_fraction=PC_MEMBRANE_FRACTION,
            interval=MAX_INTERVAL
    ):
    
    idx_a_list = randomize_list_order(IDX_A_THRESHOLDS)
    idx_b_list = randomize_list_order(IDX_B_THRESHOLDS)
    
    # --- FIX: Normalize Z-coordinates of the TEMPLATE to ensure it starts exactly at Z=0 ---
    if not df_template.empty and 'Z' in df_template.columns:
        # Find the minimum Z of the template data
        min_z_template = df_template['Z'].min()
        
        if min_z_template != 0:
            print(f"Normalizing TEMPLATE Z-coordinates by subtracting minimum Z: {min_z_template:.2f}")
            # Use .copy() to ensure we modify the Z in the current scope without
            # potential SettingWithCopyWarning if df_template was a slice.
            df_template = df_template.copy() 
            df_template['Z'] -= min_z_template
    # --- END FIX ---
    
    # Step 1-4: Scale and Process Doublet Data (which now starts from a Z=0 base)
    df_scaled = scale_cilia_df(df_template, cilia_length, interval)
    df_filtered = process_cilia_csv(df_scaled, idx_a_list, idx_b_list)

    # The normalization is now done on the template, so the output of df_filtered 
    # should have min(Z) = 0, assuming proportional scaling.
    
    # Step 7: Add Membrane (DoubletNumber = 0)
    total_membrane_length = membrane_fraction*cilia_length
    n_membrane_points = int(np.ceil(total_membrane_length / interval)) + 1
    # Ensure membrane Z starts at 0, which is guaranteed by np.linspace(0, ...)
    membrane_z_coords = np.linspace(0, total_membrane_length, n_membrane_points)
    
    membrane_data = []        
    for z in membrane_z_coords:
            membrane_data.append([
                0, 0, 0, z, # DoubletNumber=0, X=0, Y=0, Z
                1, 1, 0, # Idx_A=1, Idx_B=1, Angle=0
                0, 0 # A_Shift=0, B_Shift=0 (unused for membrane)
            ])
            
    df_membrane = pd.DataFrame(membrane_data, columns=REQUIRED_COLUMNS)
    
    # Step 8: Combine all data and return
    df_final = pd.concat([df_filtered, df_membrane], ignore_index=True)
    
    print("Processing complete.")
    return df_final