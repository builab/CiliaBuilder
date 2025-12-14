"""
Calculate points for cilia structure with curved centerlines
"""
import numpy as np


def generate_centerline_points(length=10.0, num_points=100, 
                               centerline_type='straight',
                               curve_radius=5.0, 
                               # curve_angle removed
                               sine_frequency=2.0, sine_amplitude=2.0):
    """
    Generate centerline points for straight, curved, or sinusoidal paths.
    
    Parameters:
    -----------
    length : float
        Length along the center line path (arc length for 'curve') (default: 10.0)
    num_points : int
        Number of points along the centerline (default: 100). Determines smoothness.
    centerline_type : str
        Type of center line: 'straight', 'curve', or 'sinusoidal' (default: 'straight')
    curve_radius : float
        Radius of curvature for 'curve' type (default: 5.0)
    # curve_angle removed: total angle is derived from length / curve_radius
    sine_frequency : float
        Frequency of sinusoidal oscillation (default: 2.0)
    sine_amplitude : float
        Amplitude of sinusoidal oscillation (default: 2.0)
    
    Returns:
    --------
    points : numpy.ndarray
        Array of shape (num_points, 3) containing (x, y, z) coordinates
    """
    
    # Parameter t goes from 0 to length
    t = np.linspace(0, length, num_points)
    
    # Generate center line coordinates based on type
    if centerline_type == 'straight':
        # Straight line along z-axis
        x_center = np.zeros(num_points)
        y_center = np.zeros(num_points)
        z_center = t
        
    elif centerline_type == 'curve':
        # Curved arc in the x-z plane
        # Calculate total angle (in radians) from arc length (length) and radius
        # s = R * theta => theta = s / R
        if curve_radius == 0:
             raise ValueError("curve_radius cannot be zero for 'curve' type.")
        total_angle_rad = length / curve_radius
        
        # theta parameter goes from 0 to total_angle_rad
        # This defines the angle of the arc in the XZ plane.
        theta = (t / length) * total_angle_rad
        
        x_center = curve_radius * (1 - np.cos(theta))  # Displacement in X
        y_center = np.zeros(num_points)
        z_center = curve_radius * np.sin(theta)
        
    elif centerline_type == 'sinusoidal':
        # Planar sinusoidal curve in the x-z plane
        x_center = sine_amplitude * np.sin(sine_frequency * 2 * np.pi * t / length)
        y_center = np.zeros(num_points)
        z_center = t
        
    else:
        raise ValueError("centerline_type must be 'straight', 'curve', or 'sinusoidal'")
    
    # Combine into array
    points = np.column_stack([x_center, y_center, z_center]).astype(np.float32)
    
    return points


def calculate_doublet_positions(centerline_points, doublet_index, 
                                total_doublets=9, cilia_radius=190.0):
# ... (function remains the same)
    # Calculate angle for this doublet (evenly distributed around 360°)
    angle = (360.0 / total_doublets) * doublet_index
    
    return angle, cilia_radius


def generate_cilia_structure(length=5000.0, 
                             # num_points argument removed from here
                             centerline_type='straight',
                             curve_radius=5000.0, 
                             # curve_angle removed
                             sine_frequency=2.0, sine_amplitude=500.0,
                             num_doublets=9, cilia_radius=190.0):
    """
    Generate complete cilia structure with centerline and doublet positions.
    
    Parameters:
    -----------
    length : float
        Length of the cilia (arc length for 'curve') (default: 5000.0 Angstroms)
    # num_points is now calculated internally for optimal smoothness
    centerline_type : str
        Type of centerline: 'straight', 'curve', or 'sinusoidal' (default: 'straight')
    curve_radius : float
        Radius of curvature for 'curve' type (default: 5000.0)
    # curve_angle removed
    sine_frequency : float
        Frequency of sinusoidal oscillation (default: 2.0)
    sine_amplitude : float
        Amplitude of sinusoidal oscillation (default: 500.0)
    num_doublets : int
        Number of doublet microtubules around perimeter (default: 9)
    cilia_radius : float
        Radius from center to doublets (default: 190.0 Angstroms)
    
    Returns:
    --------
    structure : dict
        Dictionary containing:
        - 'centerline': centerline points for central pair
        - 'doublets': list of (angle, shift_distance) for each doublet
        - 'num_doublets': number of doublets
        - 'centerline_type': type of centerline used
    """
    
    # === CRUCIAL CHANGE FOR SMOOTHNESS ===
    # Set a maximum interval between points (e.g., 10 Angstroms) 
    # to ensure high-density sampling for smooth curves.
    MAX_INTERVAL = 10.0 
    num_points = int(length / MAX_INTERVAL) + 1
    # ====================================

    # Generate centerline for central pair
    centerline = generate_centerline_points(
        length=length,
        num_points=num_points, # Use the high-density number
        centerline_type=centerline_type,
        curve_radius=curve_radius,
        # curve_angle removed
        sine_frequency=sine_frequency,
        sine_amplitude=sine_amplitude
    )
    
    # ... (Rest of the function remains the same)
    
    # Calculate positions for each doublet
    doublets = []
    for i in range(num_doublets):
        angle, shift_dist = calculate_doublet_positions(
            centerline, i, num_doublets, cilia_radius
        )
        doublets.append({
            'index': i,
            'angle': angle,
            'shift_distance': shift_dist,
            'name': f"doublet_{i+1}"
        })
    
    structure = {
        'centerline': centerline,
        'doublets': doublets,
        'num_doublets': num_doublets,
        'centerline_type': centerline_type,
        'length': length,
        'cilia_radius': cilia_radius
    }
    
    return structure


def get_doublet_centerline(cilia_centerline, angle, shift_distance):
# ... (function remains the same as high-density points should resolve tangent issues)
    """
    Calculate the centerline for a doublet microtubule by shifting the cilia centerline.
    
    Parameters:
    -----------
    cilia_centerline : numpy.ndarray
        Main cilia centerline points (N, 3)
    angle : float
        Angle in degrees for the doublet position
    shift_distance : float
        Radial distance from the cilia centerline
    
    Returns:
    --------
    doublet_centerline : numpy.ndarray
        Shifted centerline points for the doublet (N, 3)
    """
    n_points = len(cilia_centerline)
    doublet_centerline = np.zeros_like(cilia_centerline)
    
    # Convert angle to radians
    angle_rad = np.radians(angle)
    
    for i in range(n_points):
        # Calculate tangent vector
        if i == 0:
            tangent = cilia_centerline[i+1] - cilia_centerline[i]
        elif i == n_points - 1:
            tangent = cilia_centerline[i] - cilia_centerline[i-1]
        else:
            tangent = cilia_centerline[i+1] - cilia_centerline[i-1]
        
        tangent = tangent / np.linalg.norm(tangent)
        
        # --- Kink Fix: Use a stable reference vector (Y-axis) ---
        # Since the curve is planar in XZ, the Y-axis (0, 1, 0) is a stable perpendicular reference.
        # This prevents the abrupt flip when the Z-component of the tangent changes.
        up = np.array([0.0, 1.0, 0.0]) 
        
        # Ensure 'up' is not parallel to the tangent (only happens if tangent is also (0,1,0))
        # If the tangent is close to (0,1,0), choose a different reference, e.g., (1,0,0)
        if np.linalg.norm(np.cross(tangent, up)) < 1e-6:
             up = np.array([1.0, 0.0, 0.0])
        # ---------------------------------------------------------
        
        normal = np.cross(tangent, up)
        normal = normal / np.linalg.norm(normal)
        binormal = np.cross(tangent, normal)
        
        # Calculate shift direction based on angle
        shift_vector = shift_distance * (np.cos(angle_rad) * normal + np.sin(angle_rad) * binormal)
        
        # Apply shift
        doublet_centerline[i] = cilia_centerline[i] + shift_vector
    
    return doublet_centerline


# Example usage
if __name__ == "__main__":
    print("=" * 60)
    print("Cilia Structure Point Calculator")
    print("=" * 60)
    
    # Generate straight cilia structure
    print("\n1. Generating STRAIGHT cilia structure...")
    structure = generate_cilia_structure(
        length=5000.0,
        centerline_type='straight',
        num_doublets=9,
        cilia_radius=190.0
    )
    
    print(f"Centerline type: {structure['centerline_type']}")
    print(f"Centerline points shape: {structure['centerline'].shape}")
    print(f"Number of doublets: {structure['num_doublets']}")
    
    # Generate curved cilia structure (using length as arc length)
    print("\n2. Generating CURVED cilia structure...")
    structure_curved = generate_cilia_structure(
        length=5000.0, # Arc length = 5000
        centerline_type='curve',
        curve_radius=10000.0, # Radius = 10000
        num_doublets=9,
        cilia_radius=190.0
    )
    
    print(f"Centerline type: {structure_curved['centerline_type']}")
    print(f"Curve radius: 10000.0 Å")
    print(f"Calculated Total Angle (approx): {np.degrees(structure_curved['length'] / 10000.0):.2f}°")
    
    # Generate sinusoidal cilia structure
    print("\n3. Generating SINUSOIDAL cilia structure...")
    structure_sine = generate_cilia_structure(
        length=5000.0,
        centerline_type='sinusoidal',
        sine_frequency=3.0,
        sine_amplitude=500.0,
        num_doublets=9,
        cilia_radius=190.0
    )
    
    print(f"Centerline type: {structure_sine['centerline_type']}")
    
    print("\n" + "=" * 60)