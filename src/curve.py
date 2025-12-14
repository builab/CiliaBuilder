"""
Calculate points for cilia structure with curved centerlines
"""
import numpy as np


def generate_centerline_points(length=10.0, num_points=100, 
                               centerline_type='straight',
                               curve_radius=5.0, 
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
        # Calculate total angle (in radians) from arc length (length) and radius: s = R * theta
        if curve_radius == 0:
             raise ValueError("curve_radius cannot be zero for 'curve' type.")
        total_angle_rad = length / curve_radius
        
        # theta parameter goes from 0 to total_angle_rad
        theta = (t / length) * total_angle_rad
        
        x_center = curve_radius * np.sin(theta)
        y_center = np.zeros(num_points)
        z_center = curve_radius * (1 - np.cos(theta))
        
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
    """
    Calculate the position offset for a doublet microtubule around the cilia perimeter.
    """
    # Calculate angle for this doublet (evenly distributed around 360°)
    angle = (360.0 / total_doublets) * doublet_index
    
    return angle, cilia_radius


def generate_cilia_structure(length=5000.0,
                             centerline_type='straight',
                             curve_radius=5000.0, 
                             sine_frequency=2.0, sine_amplitude=500.0,
                             num_doublets=9, cilia_radius=190.0):
    """
    Generate complete cilia structure with centerline and doublet positions.
    """
    
    # === CRUCIAL CHANGE FOR SMOOTHNESS ===
    # Set a maximum interval between points (e.g., 10 Angstroms) 
    # to ensure high-density sampling for smooth curves and stable tangents.
    MAX_INTERVAL = 10.0 
    num_points = int(length / MAX_INTERVAL) + 1
    # ====================================

    # Generate centerline for central pair
    centerline = generate_centerline_points(
        length=length,
        num_points=num_points,
        centerline_type=centerline_type,
        curve_radius=curve_radius,
        sine_frequency=sine_frequency,
        sine_amplitude=sine_amplitude
    )
    
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
    """
    Calculate the centerline for a doublet microtubule by shifting the cilia centerline.
    
    This function uses a stable reference vector (Y-axis) to calculate the Normal/Binormal 
    frame, which fixes the sudden kinks/twists in the geometry.
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
        up = np.array([0.0, 1.0, 0.0]) 
        
        # If the tangent is near parallel to 'up' (e.g., tangent is along Y-axis), choose X-axis instead
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
    print(f"Cilia radius: {structure['cilia_radius']}")
    print("\nDoublet positions:")
    for doublet in structure['doublets']:
        print(f"  {doublet['name']}: angle={doublet['angle']:.1f}°, "
              f"shift={doublet['shift_distance']:.1f}Å")
    
    # Generate curved cilia structure
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
    # Total angle in radians = 5000 / 10000 = 0.5 rad
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
    print(f"Sine frequency: 3.0")
    print(f"Sine amplitude: 500.0 Å")
    
    # Example: Get specific doublet centerline
    print("\n4. Calculating specific doublet centerline...")
    doublet_0 = structure['doublets'][0]
    doublet_centerline = get_doublet_centerline(
        structure['centerline'],
        doublet_0['angle'],
        doublet_0['shift_distance']
    )
    print(f"Doublet 0 centerline shape: {doublet_centerline.shape}")
    print(f"First point: {doublet_centerline[0]}")
    print(f"Last point: {doublet_centerline[-1]}")
    
    print("\n" + "=" * 60)