# geometry/centerline.py

import numpy as np
from scipy.interpolate import UnivariateSpline

from .base import calculate_radial_position

def generate_centerline_points(length=10000.0, num_points=100, 
                               centerline_type='straight',
                               curve_radius=10000.0, 
                               sine_frequency=1.0, sine_amplitude=1000.0,
                               template_file='template.csv'):
    """
    Generate centerline points for straight, curved, or sinusoidal paths.
    
    Parameters:
    -----------
    length : float
        Length along the center line path (arc length for 'curve') (default: 10.0)
    num_points : int
        Number of points along the centerline (default: 100). Determines smoothness.
    centerline_type : str
        Type of center line: 'straight', 'curve', 'sinusoidal', or 'template' (default: 'straight')
    curve_radius : float
        Radius of curvature for 'curve' type (default: 10000.0)
    sine_frequency : float
        Frequency of sinusoidal oscillation (default: 1.0)
    sine_amplitude : float
        Amplitude of sinusoidal oscillation (default: 1000.0)
    template_file : str
        Path to CSV file containing template points (default: 'template.csv')
        Required for centerline_type='template'. Must have at least 100 points.
    
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
        # Calculate Z-offset using the provided formula with cilia_radius = 1000
        cilia_radius = 1000.0
        angle_rad = np.arctan2(sine_amplitude, length / (sine_frequency * 4))
        z_offset = np.cos(np.pi/2 - angle_rad) * (2 * cilia_radius)
        
        # Generate upper part: sinusoidal curve shifted by z_offset
        x_upper = sine_amplitude * np.sin(sine_frequency * 2 * np.pi * t / length)
        y_upper = np.zeros(num_points)
        z_upper = t + z_offset
        
        # Generate lower part: straight segment from 0 to z_offset/2
        straight_length = z_offset / 2
        
        # Find indices for different regions
        # Lower region: straight from 0 to z_offset/2
        lower_mask = t <= straight_length
        # Upper region: original sinusoidal (no interpolation needed)
        upper_mask = t >= z_offset
        # Interpolation region: between z_offset/2 and z_offset
        interp_mask = ~lower_mask & ~upper_mask
        
        # Initialize arrays
        x_center = np.zeros(num_points)
        y_center = np.zeros(num_points)
        z_center = np.zeros(num_points)
        
        # Lower part: straight segment
        x_center[lower_mask] = 0
        y_center[lower_mask] = 0
        z_center[lower_mask] = t[lower_mask]
        
        # Upper part: sinusoidal with offset
        x_center[upper_mask] = x_upper[upper_mask]
        y_center[upper_mask] = y_upper[upper_mask]
        z_center[upper_mask] = z_upper[upper_mask]
        
        # Interpolation region: smooth transition
        if np.any(interp_mask):
            t_interp = t[interp_mask]
            # Interpolation parameter from 0 to 1
            alpha = (t_interp - straight_length) / (z_offset - straight_length)
            
            # Lower endpoint values (end of straight segment)
            x_lower_end = 0
            z_lower_end = straight_length
            
            # Upper endpoint values (start of sinusoidal curve)
            x_upper_start = x_upper[interp_mask]
            z_upper_start = z_upper[interp_mask]
            
            # Linear interpolation
            x_center[interp_mask] = (1 - alpha) * x_lower_end + alpha * x_upper_start
            y_center[interp_mask] = 0
            z_center[interp_mask] = (1 - alpha) * z_lower_end + alpha * z_upper_start
    
    elif centerline_type == 'template':
        # Load template from CSV file
        try:
            template_data = np.loadtxt(template_file, delimiter=',', skiprows=0)
            
            # Check if first row is a header (contains non-numeric data)
            # Try to detect if first row should be skipped
            first_row = template_data[0] if template_data.ndim > 1 else template_data
            
            # If the file might have headers, try loading with skiprows=1
            try:
                # Check if we can convert first row to float
                _ = float(first_row[0]) if hasattr(first_row, '__getitem__') else float(first_row)
            except (ValueError, TypeError):
                # First row is not numeric, skip it
                template_data = np.loadtxt(template_file, delimiter=',', skiprows=1)
                
        except Exception as e:
            raise ValueError(f"Could not load template file '{template_file}': {e}")
        
        # Validate template has at least 100 points
        if len(template_data) < 100:
            raise ValueError(f"Template file must have at least 100 points, got {len(template_data)}")
        
        # Ensure template has 2 columns (x, y) which will be converted to (x, 0, z)
        if template_data.ndim == 1:
            raise ValueError("Template file must have 2 columns (x, y)")
        if template_data.shape[1] != 2:
            raise ValueError(f"Template file must have 2 columns (x, y), got {template_data.shape[1]}")
        
        # Extract x, y from template and convert to x, 0, z (planar in XZ plane)
        x_template = template_data[:, 0]
        y_template = np.zeros(len(template_data))  # Y is always 0 (planar in XZ)
        z_template = template_data[:, 1]  # Second column becomes Z
        
        # Calculate cumulative arc length of the template
        dx = np.diff(x_template)
        dy = np.diff(y_template)
        dz = np.diff(z_template)
        segment_lengths = np.sqrt(dx**2 + dy**2 + dz**2)
        template_arc_length = np.concatenate([[0], np.cumsum(segment_lengths)])
        template_total_length = template_arc_length[-1]
        
        # Scale factor to match desired length
        scale_factor = length / template_total_length
        
        # Scale the arc length parameter for interpolation
        scaled_arc_length = template_arc_length * scale_factor
        
        # Create smooth splines for each coordinate with scaled coordinates
        # Use UnivariateSpline with smoothing factor s=0 for exact interpolation
        spline_x = UnivariateSpline(scaled_arc_length, scale_factor * x_template, s=0, k=3)
        spline_y = UnivariateSpline(scaled_arc_length, scale_factor * y_template, s=0, k=3)
        spline_z = UnivariateSpline(scaled_arc_length, scale_factor * z_template, s=0, k=3)
        
        # Generate new arc length parameter for desired number of points
        new_arc_length = np.linspace(0, length, num_points)
        
        # Evaluate splines at new points
        x_center = spline_x(new_arc_length)
        y_center = spline_y(new_arc_length)
        z_center = spline_z(new_arc_length)
        
    else:
        raise ValueError("centerline_type must be 'straight', 'curve', 'sinusoidal', or 'template'")
    
    # Combine into array
    points = np.column_stack([x_center, y_center, z_center]).astype(np.float32)
    
    return points


def generate_cilia_structure(length=15000.0, 
                             centerline_type='straight',
                             curve_radius=10000.0, 
                             sine_frequency=1.0, sine_amplitude=1000.0,
                             template_file='template.csv',
                             num_doublets=9, cilia_radius=875.0):
    """
    Generate complete cilia structure with centerline and doublet positions.
    
    Parameters:
    -----------
    length : float
        Length of the cilia (arc length for 'curve') (default: 5000.0 Angstroms)
    centerline_type : str
        Type of centerline: 'straight', 'curve', or 'sinusoidal' (default: 'straight')
    curve_radius : float
        Radius of curvature for 'curve' type (default: 5000.0)
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
    
    # Set a maximum interval between points (e.g., 10 Angstroms) 
    # to ensure high-density sampling for smooth curves.
    MAX_INTERVAL = 20.0 
    num_points = int(length / MAX_INTERVAL) + 1

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
        angle, shift_dist = calculate_radial_position(
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
        
        # Use a stable reference vector (Y-axis)
        up = np.array([0.0, 1.0, 0.0]) 
        
        # Ensure 'up' is not parallel to the tangent
        if np.linalg.norm(np.cross(tangent, up)) < 1e-6:
             up = np.array([1.0, 0.0, 0.0])
        
        normal = np.cross(tangent, up)
        normal = normal / np.linalg.norm(normal)
        binormal = np.cross(tangent, normal)
        
        # Calculate shift direction based on angle
        shift_vector = shift_distance * (np.cos(angle_rad) * normal + np.sin(angle_rad) * binormal)
        
        # Apply shift
        doublet_centerline[i] = cilia_centerline[i] + shift_vector
    
    return doublet_centerline


