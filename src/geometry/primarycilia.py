# geometry/primarycilia.py

import pandas as pd
import numpy as np
import random
from scipy.interpolate import UnivariateSpline

# --- Constants ---
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

    
def process_single_doublet(df_doublet, new_length, interval):
    """
    Process points for a single microtubule doublet.
    
    Performs the following operations:
    1. Scales Z-values to match new_length
    2. Calculates cumulative arc length (s) based on 3D distance
    3. Fits smoothing splines to X, Y, Angle, A_Shift, B_Shift vs arc length
    4. Resamples points at fixed intervals along the curve
    5. Recalculates Idx_A and Idx_B indices
    
    Parameters:
    -----------
    df_doublet : pd.DataFrame
        DataFrame containing points for a single doublet
    new_length : float
        Target length for the scaled doublet
    interval : float
        Desired spacing between resampled points
    
    Returns:
    --------
    pd.DataFrame
        Resampled and scaled DataFrame with standard columns
    """
    
    # Scale Z to new length
    max_z_old = df_doublet['Z'].max()
    if max_z_old == 0:
        df_doublet['Z_scaled'] = 0.0
    else:
        scale_factor = new_length / max_z_old
        df_doublet['Z_scaled'] = df_doublet['Z'] * scale_factor
        
    # Calculate cumulative arc length in 3D
    coords = df_doublet[['X', 'Y', 'Z_scaled']].values
    segment_lengths = np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))
    s_values = np.insert(np.cumsum(segment_lengths), 0, 0.0)
    
    # Define new sample points along arc length
    num_points = int(np.ceil(new_length / interval))
    s_new = np.linspace(0, new_length, num_points + 1)
    
    # Fit splines and interpolate
    new_data = {}
    interp_cols = ['X', 'Y', 'Angle', 'A_Shift', 'B_Shift', 'Z_scaled']
    
    for col in interp_cols:
        if len(s_values) > 3 and df_doublet[col].nunique() > 1:
            # Fit cubic smoothing spline
            spline = UnivariateSpline(s_values, df_doublet[col], k=3, s=0)
            new_data[col] = spline(s_new)
        else:
            # Handle constant or insufficient data
            new_data[col] = np.full_like(s_new, df_doublet[col].iloc[0])

    # Assemble new DataFrame
    df_new = pd.DataFrame(new_data)
    df_new['Idx_A'] = np.arange(len(df_new))
    df_new['Idx_B'] = np.arange(len(df_new))
    df_new['DoubletNumber'] = df_doublet['DoubletNumber'].iloc[0]
    df_new = df_new.rename(columns={'Z_scaled': 'Z'})
    
    return df_new[REQUIRED_COLUMNS]


def randomize_list_order(input_list):
    """
    Shuffle a list into random order.
    
    Parameters:
    -----------
    input_list : list
        List to be randomized
    
    Returns:
    --------
    list
        New list with elements in random order
    """
    shuffled_list = input_list.copy()
    random.shuffle(shuffled_list)
    return shuffled_list


def scale_cilia_df(df_in, cilia_length, interval):
    """
    Scale, smooth, and resample coordinates for all microtubule doublets.
    
    Processes each doublet independently by:
    - Scaling to target length
    - Fitting smoothing splines
    - Resampling at regular intervals
    
    Parameters:
    -----------
    df_in : pd.DataFrame
        Input DataFrame with doublet geometry
    cilia_length : float
        Target maximum Z-value (final length)
    interval : float
        Desired spacing between new points
    
    Returns:
    --------
    pd.DataFrame
        Resampled, scaled, and smoothed DataFrame
    """
    if 'DoubletNumber' not in df_in.columns:
        raise ValueError("Input DataFrame must contain 'DoubletNumber' column")

    if df_in['Z'].max() == 0:
        print("Warning: Max Z is zero. Scaling will result in zero-length structure")

    print(f"Scaling and smoothing {df_in['DoubletNumber'].nunique()} doublets...")
    
    # Process each doublet separately
    scaled_df = (
        df_in.groupby('DoubletNumber', group_keys=False)
             .apply(process_single_doublet, new_length=cilia_length, interval=interval)
             .reset_index(drop=True)
    )
    
    print("Processing complete.")
    return scaled_df


def process_cilia_csv(df_scaled, idx_a_list, idx_b_list):
    """
    Apply length-based filtering for Idx_A and Idx_B based on thresholds.
    
    For each doublet:
    - Sets Idx_A=1 for points below threshold percentage of length
    - Sets Idx_B=1 for points below threshold percentage of length
    - Updates Angle based on doublet number
    
    Parameters:
    -----------
    df_scaled : pd.DataFrame
        Scaled doublet data
    idx_a_list : list
        List of Idx_A threshold fractions (0-1) for each doublet
    idx_b_list : list
        List of Idx_B threshold fractions (0-1) for each doublet
    
    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame with updated Idx_A, Idx_B, and Angle values
    """
    df_filtered = []
    doublet_numbers = sorted(df_scaled['DoubletNumber'].unique())
    
    for i, doublet_num in enumerate(doublet_numbers):
        df_doublet = df_scaled[df_scaled['DoubletNumber'] == doublet_num].copy()
        
        max_z = df_doublet['Z'].max()
        
        # Set angle
        angle_value = 90 + (doublet_num - 1) * 40
        df_doublet['Angle'] = angle_value
        
        # Calculate Z thresholds
        idx_a_threshold_z = max_z * idx_a_list[i % len(idx_a_list)]
        idx_b_threshold_z = max_z * idx_b_list[i % len(idx_b_list)]
        
        # Apply filtering
        df_doublet['Idx_A'] = (df_doublet['Z'] < idx_a_threshold_z).astype(int)
        df_doublet['Idx_B'] = (df_doublet['Z'] < idx_b_threshold_z).astype(int)
        df_doublet['Idx_B'].iloc[0] = 0 # Make plotting better

        df_filtered.append(df_doublet)

    return pd.concat(df_filtered, ignore_index=True) if df_filtered else pd.DataFrame(columns=REQUIRED_COLUMNS)


def generate_primary_cilia(
    df_template,
    cilia_length=PRIMARY_CILIA_LENGTH,
    membrane_radius=MEMBRANE_RADIUS,
    membrane_fraction=PC_MEMBRANE_FRACTION,
    interval=MAX_INTERVAL
):
    """
    Generate primary cilia structure from a template DataFrame.
    
    Creates a primary cilia model by:
    1. Normalizing template to start at Z=0
    2. Scaling doublets to target length with smooth interpolation
    3. Randomizing Idx_A and Idx_B threshold assignments
    4. Adding membrane coverage
    
    Primary cilia characteristics:
    - 9 doublet microtubules with variable A and B tubule lengths
    - No central pair
    - Membrane covering basal region
    
    Parameters:
    -----------
    df_template : pd.DataFrame
        Template DataFrame with doublet geometry (typically from motile cilia)
    cilia_length : float
        Target length of primary cilia (Å)
    membrane_radius : float
        Radius of membrane (not used in geometry, kept for compatibility)
    membrane_fraction : float
        Fraction of length covered by membrane (0-1)
    interval : float
        Maximum spacing between points
    
    Returns:
    --------
    pd.DataFrame
        Primary cilia structure with columns:
        DoubletNumber, X, Y, Z, Idx_A, Idx_B, Angle, A_Shift, B_Shift
        
        Special DoubletNumber values:
        - Positive (1-9): Doublet microtubules with variable lengths
        - 0: Membrane
    """
    
    # Randomize threshold assignments for variability
    idx_a_list = randomize_list_order(IDX_A_THRESHOLDS)
    idx_b_list = randomize_list_order(IDX_B_THRESHOLDS)
    
    # Normalize template to start at Z=0
    if not df_template.empty and 'Z' in df_template.columns:
        min_z_template = df_template['Z'].min()
        
        if min_z_template != 0:
            print(f"Normalizing template Z-coordinates (subtracting {min_z_template:.2f})")
            df_template = df_template.copy()
            df_template['Z'] -= min_z_template
    
    # Scale and process doublets
    df_scaled = scale_cilia_df(df_template, cilia_length, interval)
    df_filtered = process_cilia_csv(df_scaled, idx_a_list, idx_b_list)
    
    # Add membrane
    print(f"Generating membrane (covering {membrane_fraction*100:.1f}% of base)...")
    total_membrane_length = membrane_fraction * cilia_length
    n_membrane_points = int(np.ceil(total_membrane_length / interval)) + 1
    membrane_z_coords = np.linspace(0, total_membrane_length, n_membrane_points)
    
    membrane_data = []
    for z in membrane_z_coords:
        membrane_data.append([
            0, 0, 0, z,  # DoubletNumber=0, X=0, Y=0, Z
            1, 1, 0,     # Idx_A=1, Idx_B=1, Angle=0
            0, 0         # A_Shift=0, B_Shift=0
        ])
    
    df_membrane = pd.DataFrame(membrane_data, columns=REQUIRED_COLUMNS)
    
    # Combine all components
    df_final = pd.concat([df_filtered, df_membrane], ignore_index=True)
    
    print(f"Primary cilia generation complete ({len(df_final)} total points)")
    return df_final