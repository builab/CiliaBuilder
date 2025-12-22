# geometry/tip.py

import numpy as np
import pandas as pd
import random
from .. import default_config

# --- Constants ---
INITIAL_LENGTH = getattr(default_config, 'TIP_INITIAL_LENGTH', 300)
TRANSITION_LENGTH = getattr(default_config, 'TIP_TRANSITION_LENGTH', 2000)
CILIA_DOUBLET_SHIFT = getattr(default_config, 'CILIA_DOUBLET_SHIFT', 70)
Z_MAX_IDX_B_DEFAULT = INITIAL_LENGTH + TRANSITION_LENGTH
MAX_INTERVAL = getattr(default_config, 'MAX_INTERVAL', 20)


def generate_tip_curves(tip_length, cilia_radius, final_radius, interval=MAX_INTERVAL):
    """
    Generate 9 cilia curves centered around the Z axis with tapered tips.
    Uses single cosine interpolation for smooth, monotonic radius decrease.
    
    Parameters:
    -----------
    tip_length : float
        Total length of the tip section
    cilia_radius : float
        Radius at Z=0 (starting radius)
    final_radius : float
        Radius at Z=tip_length (final radius)
    interval : float
        Maximum spacing between points (default: MAX_INTERVAL)
    
    Returns:
    --------
    list of dict
        List of 9 curves, each containing:
        - 'curve': numpy array of shape (n_points, 3) with (x, y, z) coordinates
        - 'control_points': numpy array of 4 control points
    """
    num_lines = 9
    curves = []
    
    # Segment 1: p1 to p2 (straight section, length = INITIAL_LENGTH)
    n_points_linear = int(np.ceil(INITIAL_LENGTH / interval)) + 1
    n_points_linear = max(n_points_linear, 5)

    # Segment 2: p2 to p4 (tapered section, length = tip_length)
    n_points_spline = int(np.ceil(tip_length / interval)) + 1
    n_points_spline = max(n_points_spline, 5)
    
    # Z-coordinates for the tapered segment
    z_spline = np.linspace(INITIAL_LENGTH, tip_length + INITIAL_LENGTH, n_points_spline)

    for i in range(num_lines):
        angle = (i / num_lines) * 2 * np.pi + np.pi  # Add π to match curve.py
        
        # Define 4 control points
        p1 = np.array([cilia_radius * np.cos(angle), cilia_radius * np.sin(angle), 0])
        p2 = np.array([p1[0], p1[1], INITIAL_LENGTH])
        p3 = np.array([final_radius * np.cos(angle), final_radius * np.sin(angle), tip_length + INITIAL_LENGTH])
        control_points = np.array([p1, p2, p3])
        
        # Linear segment from p1 to p2 (straight)
        t_linear = np.linspace(0, 1, n_points_linear)
        linear_segment = np.array([p1 + t * (p2 - p1) for t in t_linear])
        
        # Smooth tapered segment using single cosine interpolation
        t = (z_spline - INITIAL_LENGTH) / tip_length
        t_smooth = (1 - np.cos(t * np.pi)) / 2
        radii = cilia_radius + t_smooth * (final_radius - cilia_radius)
        
        # Calculate X and Y based on radius and angle
        x_new = radii * np.cos(angle)
        y_new = radii * np.sin(angle)
        
        # Combine X, Y, Z to form the smooth segment
        smooth_segment = np.column_stack([x_new, y_new, z_spline])
        
        # Combine segments (skip first point of smooth segment to avoid duplication)
        curve = np.vstack([linear_segment, smooth_segment[1:]])
        
        curves.append({
            'curve': curve,
            'control_points': control_points
        })
    
    return curves


def generate_multiple_tip_lengths_in_memory(
    cilia_radius=default_config.CILIA_RADIUS,
    final_radius=default_config.TIP_FINAL_RADIUS,
    tip_length_start=default_config.TIP_INITIAL_LENGTH,
    tip_length_end=default_config.TIP_LENGTH + default_config.TIP_INITIAL_LENGTH,
    number_of_steps=10
):
    """
    Generate multiple tip length variations in memory for randomized tip assembly.
    
    Parameters:
    -----------
    cilia_radius : float
        Starting radius of cilia
    final_radius : float
        Final radius at tip end
    tip_length_start : float
        Minimum tip length
    tip_length_end : float
        Maximum tip length
    number_of_steps : int
        Number of different tip lengths to generate
    
    Returns:
    --------
    list of dict
        List of tip length variations, each containing:
        - 'tip_length': float
        - 'curves': list of curve dictionaries
    """
    all_curves_data = []
    tip_lengths = np.linspace(tip_length_start, tip_length_end, max(1, number_of_steps))

    print("GENERATING MULTIPLE TIP LENGTH VARIATIONS (IN MEMORY)")

    for tip_length in tip_lengths:
        print(f"Generating curves for tip_length = {tip_length}...")
        
        curves = generate_tip_curves(
            tip_length=tip_length, 
            cilia_radius=cilia_radius, 
            final_radius=final_radius
        )
        
        all_curves_data.append({
            'tip_length': tip_length,
            'curves': curves
        })

    print(f"Generated {len(all_curves_data)} tip length variations")

    return all_curves_data


def _generate_doublet_rows(curve_data, doublet_number, z_max_idx_b, doublet_shift, num_doublets=None):
    """
    Generate CSV row data for a single doublet curve.
    
    Parameters:
    -----------
    curve_data : dict
        Curve data containing 'curve' array
    doublet_number : int
        Doublet number (1-9)
    z_max_idx_b : float
        Maximum Z for Idx_B=1 region
    doublet_shift : float
        Shift distance for A and B tubules
    num_doublets : int, optional
        Total number of doublets (for angle calculation)
    
    Returns:
    --------
    list
        List of row data [DoubletNumber, X, Y, Z, Idx_A, Idx_B, Angle, A_Shift, B_Shift]
    """
    curve = curve_data['curve']
    all_data = []

    # Calculate angle
    if num_doublets:
        angle = 90 + (360 / num_doublets) * (doublet_number - 1)
    else:
        angle = 90 + 40 * (doublet_number - 1)

    # Define constant values
    IDX_A = 1
    A_SHIFT = -doublet_shift
    B_SHIFT = doublet_shift
    
    # Randomly select the transition Z-coordinate for this curve
    idx_b_transition_z = random.uniform(0, z_max_idx_b)
    
    # Generate row data for each point
    for x, y, z in curve:
        idx_b = 1 if z <= idx_b_transition_z else 0
        
        all_data.append([
            doublet_number, x, y, z, 
            IDX_A, idx_b, angle, 
            A_SHIFT, B_SHIFT
        ])
        
    return all_data


def create_mixed_tip_from_memory(all_curves_data, 
                                 initial_length=INITIAL_LENGTH, 
                                 transition_length=TRANSITION_LENGTH, 
                                 num_doublets=default_config.CILIA_NUM_DOUBLETS, 
                                 doublet_shift=CILIA_DOUBLET_SHIFT):
    """
    Create mixed tip df by randomly selecting curves from different tip lengths.
    Ensures each tip length is used at most once for variety.
    
    Parameters:
    -----------
    all_curves_data : list
        List of tip length variations from generate_multiple_tip_lengths_in_memory
    initial_length : float
        Initial straight section length
    transition_length : float
        Transition region length
    num_doublets : int
        Number of doublets
    doublet_shift : float
        Shift distance for A and B tubules
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with mixed tip geometry
    """
    mixed_data = []
    Z_MAX_IDX_B = initial_length + transition_length
    USE_LIMIT = 1  # Each tip length used at most once
    
    # Track usage count for each tip length
    usage_count = {i: 0 for i in range(len(all_curves_data))}
    
    for doublet_number in range(1, num_doublets + 1):
        # Filter available tip lengths
        available_indices = [idx for idx, count in usage_count.items() if count < USE_LIMIT]
        
        if not available_indices:
            # Reset if all options exhausted
            print(f"Warning: All tip lengths used, resetting counts for doublet {doublet_number}")
            usage_count = {i: 0 for i in range(len(all_curves_data))}
            available_indices = list(usage_count.keys())
        
        # Randomly select from available tip lengths
        selected_idx = random.choice(available_indices)
        usage_count[selected_idx] += 1
        
        selected_data = all_curves_data[selected_idx]
        selected_tip_length = selected_data['tip_length']
        selected_curves = selected_data['curves']
        
        print(f"Doublet {doublet_number}: Selected tip_length={selected_tip_length:.1f} (used {usage_count[selected_idx]}/{USE_LIMIT})")
        
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

    return mixed_df


def generate_tip_csv(
    cilia_radius=default_config.CILIA_RADIUS, 
    tip_length_end=default_config.TIP_LENGTH, 
    final_radius=default_config.TIP_FINAL_RADIUS, 
    initial_length=INITIAL_LENGTH, 
    transition_length=TRANSITION_LENGTH, 
    num_doublets=default_config.CILIA_NUM_DOUBLETS, 
    doublet_shift=CILIA_DOUBLET_SHIFT
):
    """
    Generate tip geometry with randomized tip lengths for each doublet.
    
    Parameters:
    -----------
    cilia_radius : float
        Starting radius
    tip_length_end : float
        Maximum tip length
    final_radius : float
        Final radius at tip end
    initial_length : float
        Initial straight section length
    transition_length : float
        Transition region length
    num_doublets : int
        Number of doublets
    doublet_shift : float
        Shift distance for A and B tubules
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with tip geometry for all doublets
    """
    tip_length_start = initial_length + transition_length
    
    # Generate multiple tip length variations
    all_curves_data = generate_multiple_tip_lengths_in_memory( 
        cilia_radius, final_radius,
        tip_length_start, tip_length_end, number_of_steps=10
    )
            
    # Create mixed CSV with random tip lengths per doublet
    mixed_df = create_mixed_tip_from_memory(
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
    Generate complete cilia structure with straight base section and tapered tip.
    
    Creates a cilia model with:
    - Straight base section from Z=0 to Z=cilia_length
    - Tapered tip section from Z=cilia_length to Z=cilia_length+tip_length
    - Central pair (DoubletNumber=-1) spanning full length
    - Cap at tip end (DoubletNumber=-2)
    - Membrane (DoubletNumber=0) covering fraction of base
    
    Parameters:
    -----------
    cilia_length : float
        Length of straight base section (Å)
    cilia_radius : float
        Radius at base
    tip_length : float
        Length of tapered tip section (includes CP extension)
    final_radius : float
        Final radius at tip end
    initial_length : float
        Initial straight section before taper starts
    transition_length : float
        Transition region length
    num_doublets : int
        Number of doublets (typically 9)
    doublet_shift : float
        Radial shift for A and B tubules
    cp_doublet_length_diff : float
        How much longer CP extends beyond doublets
    cp_shift : float
        Radial shift for central pair tubules
    membrane_radius : float
        Membrane radius (not currently used in geometry)
    membrane_fraction : float
        Fraction of base covered by membrane (0-1)
    max_interval : float
        Maximum spacing between points
        
    Returns:
    --------
    pd.DataFrame
        Complete cilia structure with columns:
        DoubletNumber, X, Y, Z, Idx_A, Idx_B, Angle, A_Shift, B_Shift
        
        Special DoubletNumber values:
        - Positive (1-9): Doublet microtubules
        - -1: Central pair
        - -2: Cap at tip end
        - 0: Membrane
    """
    
    print("=" * 60)
    print("GENERATING CILIA WITH TIP")
    print("=" * 60)
    print(f"Base length: {cilia_length} Å")
    print(f"Tip length: {tip_length} Å")
    print(f"Total length: {cilia_length + tip_length} Å")
    print("=" * 60)
    
    # Adjust tip length to account for CP extension
    tip_length_end = tip_length - cp_doublet_length_diff
    
    # Step 1: Generate tip geometry
    print("\n1. Generating tip geometry...")
    tip_df = generate_tip_csv(
        cilia_radius=cilia_radius,
        tip_length_end=tip_length_end,
        final_radius=final_radius,
        initial_length=initial_length,
        transition_length=transition_length,
        num_doublets=num_doublets,
        doublet_shift=doublet_shift
    )
    
    # Step 2: Shift tip Z-coordinates to connect to base
    print(f"\n2. Shifting tip Z-coordinates by {cilia_length - initial_length} Å...")
    tip_df['Z'] = tip_df['Z'] + cilia_length - initial_length
    
    # Step 3: Generate straight base section
    print("\n3. Generating straight base sections...")
    base_data = []
    
    for doublet_num in range(1, num_doublets + 1):
        doublet_tip_data = tip_df[tip_df['DoubletNumber'] == doublet_num]
        
        if len(doublet_tip_data) == 0:
            print(f"Warning: No tip data for doublet {doublet_num}, skipping")
            continue
        
        # Get connection point from tip
        first_point = doublet_tip_data.iloc[0]
        x1, y1 = first_point['X'], first_point['Y']
        angle = first_point['Angle']
        a_shift = first_point['A_Shift']
        b_shift = first_point['B_Shift']
        
        # Generate base points
        n_points = int(np.ceil((cilia_length - initial_length) / max_interval)) + 1
        z_coords = np.linspace(0, cilia_length - initial_length, n_points)
        
        for i, z in enumerate(z_coords):
            # First point: Idx_B=0 (B-tubule hasn't started yet)
            idx_a, idx_b = (1, 0) if i == 0 else (1, 1)
            
            base_data.append([
                doublet_num, x1, y1, z,
                idx_a, idx_b, angle,
                a_shift, b_shift
            ])
        
        print(f"  Doublet {doublet_num}: {n_points} base points from Z=0 to Z={cilia_length - initial_length:.1f}")
    
    # Step 4: Combine base and tip
    print("\n4. Combining base and tip...")
    columns = ['DoubletNumber', 'X', 'Y', 'Z', 'Idx_A', 'Idx_B', 'Angle', 'A_Shift', 'B_Shift']
    base_df = pd.DataFrame(base_data, columns=columns)
    
    # Remove duplicate points at junction
    tip_df_filtered = tip_df.groupby('DoubletNumber').apply(
        lambda group: group.iloc[1:]
    ).reset_index(drop=True)
    
    complete_data = pd.concat([base_df, tip_df_filtered], ignore_index=True)
    
    # Step 5: Add central pair
    print("\n5. Generating central pair...")
    total_length = cilia_length + tip_length
    n_cp_points = int(np.ceil(total_length / max_interval)) + 1
    cp_z_coords = np.linspace(0, total_length, n_cp_points)
    
    cp_data = []
    for z in cp_z_coords:
        cp_data.append([
            -1, 0, 0, z,
            1, 1, 0,
            cp_shift, -cp_shift
        ])
    
    print(f"  Generated {n_cp_points} central pair points")
    
    # Step 6: Add cap at tip end
    cp_data.append([
        -2, 0, 0, cp_z_coords[-1],
        1, 1, 0,
        0, 0
    ])
    
    cp_df = pd.DataFrame(cp_data, columns=columns)
    
    # Step 7: Add membrane
    print("\n6. Generating membrane...")
    total_membrane_length = membrane_fraction * cilia_length
    n_membrane_points = int(np.ceil(total_membrane_length / max_interval)) + 1
    membrane_z_coords = np.linspace(0, total_membrane_length, n_membrane_points)
    
    membrane_data = []
    for z in membrane_z_coords:
        membrane_data.append([
            0, 0, 0, z,
            1, 1, 0,
            0, 0
        ])
    
    print(f"  Generated {n_membrane_points} membrane points (covering {membrane_fraction*100:.1f}% of base)")
    
    membrane_df = pd.DataFrame(membrane_data, columns=columns)
    
    # Step 8: Combine all components
    complete_df = pd.concat([complete_data, cp_df, membrane_df], ignore_index=True)
    complete_df = complete_df.sort_values(['DoubletNumber', 'Z']).reset_index(drop=True)
    
    print("\n7. Complete!")
    print(f"   Total points: {len(complete_df)}")
    print(f"   Z range: {complete_df['Z'].min():.1f} to {complete_df['Z'].max():.1f} Å")
    print("=" * 60)
    
    return complete_df