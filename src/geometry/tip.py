import numpy as np
from .. import default_config

# --- CONSTANTS (Global variables for easy modification) ---
MAX_INTERVAL = default_config.MAX_INTERVAL # Default interval for dynamic point calculation
INITIAL_LENGTH = default_config.TIP_INITIAL_LENGTH     # Z-coordinate for control point p2
TRANSITION_LENGTH = default_config.TIP_TRANSITION_LENGTH # Z-coordinate for control point p3

def generate_tip_curves(tip_length, cilia_radius, transition_radius, final_radius, interval=MAX_INTERVAL):
    """
    Generate 9 cilia curves centered around the Z axis.
    Uses single cosine interpolation for smooth, monotonic radius decrease.
    
    Parameters:
    -----------
    tip_length : float
        Total length of the tip
    cilia_radius : float
        Radius of the circle at Z=0 (starting points)
    transition_radius : float
        Radius at Z=TRANSITION_LENGTH (not used in single-drop version, kept for compatibility)
    final_radius : float
        Radius of the circle at Z=tip_length (final points)
    
    Returns:
    --------
    curves : list of dicts
        List of 9 curves, each containing 'curve' (array of shape (n_points, 3)) 
        and 'control_points'.
    """
    num_lines = 9
    curves = []
    
    # Segment 1: p1 to p2 (Length = INITIAL_LENGTH)
    n_points_linear = int(np.ceil(INITIAL_LENGTH / interval)) + 1
    n_points_linear = max(n_points_linear, 5)

    # Segment 2: p2 to p4 (Length = tip_length)
    spline_sampling_length = tip_length
    n_points_spline = int(np.ceil(spline_sampling_length / interval)) + 1
    n_points_spline = max(n_points_spline, 5)
    
    # Z-coordinates for the transition segment
    z_spline = np.linspace(INITIAL_LENGTH, tip_length + INITIAL_LENGTH, n_points_spline)

    for i in range(num_lines):
        #angle = (i / num_lines) * 2 * np.pi  # 
        angle = (i / num_lines) * 2 * np.pi  + np.pi # Add pi to match curve.py

        # Define 4 control points
        p1 = np.array([cilia_radius * np.cos(angle), cilia_radius * np.sin(angle), 0])
        p2 = np.array([p1[0], p1[1], INITIAL_LENGTH])
        p3 = np.array([transition_radius * np.cos(angle), transition_radius * np.sin(angle), TRANSITION_LENGTH])
        p4 = np.array([final_radius * np.cos(angle), final_radius * np.sin(angle), tip_length + INITIAL_LENGTH])
        control_points = np.array([p1, p2, p3, p4])
        
        # --- Linear segment from p1 to p2 ---
        t_linear = np.linspace(0, 1, n_points_linear)
        linear_segment = np.array([
            p1 + t * (p2 - p1) for t in t_linear
        ])
        
        # --- Single smooth radius transition using cosine interpolation ---
        # Interpolate directly from cilia_radius to final_radius
        # Calculate normalized position along the tip length
        t = (z_spline - INITIAL_LENGTH) / (tip_length + INITIAL_LENGTH - INITIAL_LENGTH)
        # Cosine interpolation: smooth transition from 0 to 1
        t_smooth = (1 - np.cos(t * np.pi)) / 2
        # Calculate radius at each Z position
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

# --- Main Function: generate_multiple_tip_lengths_in_memory ---
def generate_multiple_tip_lengths_in_memory(
    cilia_radius=default_config.CILIA_RADIUS,
    transition_radius=default_config.TIP_TRANSITION_RADIUS,
    final_radius=default_config.TIP_FINAL_RADIUS,
    tip_length_start=default_config.TIP_INITIAL_LENGTH,
    tip_length_end=default_config.TIP_LENGTH+default_config.TIP_INITIAL_LENGTH,
    number_of_steps=10
):
    """
    Generate multiple tip length variations and return all curve data in memory.
    """
    all_curves_data = []
    tip_lengths = np.linspace(tip_length_start, tip_length_end, max(1, number_of_steps))

    print("GENERATING MULTIPLE TIP LENGTH VARIATIONS (IN MEMORY)")

    for tip_length in tip_lengths:
        print(f"\nGenerating curves for tip_length = {tip_length}...")
        
        # Use keyword arguments for clarity
        curves = generate_tip_curves(
            tip_length=tip_length, 
            cilia_radius=cilia_radius, 
            transition_radius=transition_radius, 
            final_radius=final_radius
        )
        
        all_curves_data.append({
            'tip_length': tip_length,
            'curves': curves
        })

    print(f"Generated {len(all_curves_data)} tip length variations")

    return all_curves_data
    
