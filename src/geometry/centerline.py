# geometry/centerline.py

import numpy as np
from scipy.interpolate import UnivariateSpline
from ..io import read_2d_csv


def _calculate_radial_position(centerline_points, doublet_index, 
                                total_doublets=9, cilia_radius=875.0):
    """
    Calculate angle for a doublet evenly distributed around 360°.
    
    Parameters:
    -----------
    centerline_points : np.ndarray
        Centerline points (not used, kept for compatibility)
    doublet_index : int
        Index of the doublet (0-based)
    total_doublets : int
        Total number of doublets
    cilia_radius : float
        Radial distance from centerline
    
    Returns:
    --------
    tuple
        (angle in degrees, radius)
    """
    angle = (360.0 / total_doublets) * doublet_index
    return angle, cilia_radius


def _generate_straight_centerline(t):
    """
    Generate points for a straight centerline along the Z-axis.
    
    Parameters:
    -----------
    t : np.ndarray
        Arc length parameter values
    
    Returns:
    --------
    tuple
        (x, y, z) coordinate arrays
    """
    x_center = np.zeros_like(t)
    y_center = np.zeros_like(t)
    z_center = t
    return x_center, y_center, z_center


def _generate_curved_centerline(t, length, curve_radius):
    """
    Generate points for a curved centerline (arc) in the X-Z plane.
    
    Parameters:
    -----------
    t : np.ndarray
        Arc length parameter values
    length : float
        Total arc length
    curve_radius : float
        Radius of curvature
    
    Returns:
    --------
    tuple
        (x, y, z) coordinate arrays
    """
    if curve_radius == 0:
        raise ValueError("curve_radius cannot be zero for 'curve' type")
    
    # Arc length s = R * θ => θ = s / R
    total_angle_rad = length / curve_radius
    theta = (t / length) * total_angle_rad
    
    # Curve starts at origin and extends into positive X-Z quadrant
    x_center = curve_radius * (1 - np.cos(theta))
    y_center = np.zeros_like(t)
    z_center = curve_radius * np.sin(theta)
    return x_center, y_center, z_center


def _generate_sinusoidal_centerline(t, length, sine_frequency, sine_amplitude):
    """
    Generate points for a sinusoidal centerline with smooth lower transition.
    
    Creates a centerline with:
    - Straight section at the base
    - Smooth transition region
    - Sinusoidal upper section
    
    Parameters:
    -----------
    t : np.ndarray
        Arc length parameter values
    length : float
        Total length
    sine_frequency : float
        Frequency of oscillation
    sine_amplitude : float
        Amplitude of oscillation
    
    Returns:
    --------
    tuple
        (x, y, z) coordinate arrays
    """
    num_points = len(t)
    CILIA_RADIUS_REF = 1000.0
    
    # Calculate Z-offset to ensure smooth start at Z=0
    angle_rad = np.arctan2(sine_amplitude, length / (sine_frequency * 4))
    z_offset = np.cos(np.pi/2 - angle_rad) * (2 * CILIA_RADIUS_REF)
    straight_length = z_offset / 2
    
    # Initialize coordinate arrays
    x_center = np.zeros(num_points)
    y_center = np.zeros(num_points)
    z_center = np.zeros(num_points)
    
    # Base sinusoidal curve
    x_upper_base = sine_amplitude * np.sin(sine_frequency * 2 * np.pi * t / length)
    z_upper_shifted = t + z_offset
    
    # Define three regions
    lower_mask = t <= straight_length
    upper_mask = t >= z_offset
    interp_mask = ~(lower_mask | upper_mask)
    
    # Lower part: straight
    z_center[lower_mask] = t[lower_mask]
    
    # Upper part: sinusoidal
    x_center[upper_mask] = x_upper_base[upper_mask]
    z_center[upper_mask] = z_upper_shifted[upper_mask]
    
    # Interpolation region: smooth transition
    if np.any(interp_mask):
        t_interp = t[interp_mask]
        alpha = (t_interp - straight_length) / (z_offset - straight_length)
        
        x_upper_start = x_upper_base[interp_mask]
        z_upper_start = z_upper_shifted[interp_mask]
        
        x_center[interp_mask] = alpha * x_upper_start
        z_center[interp_mask] = (1 - alpha) * straight_length + alpha * z_upper_start
    
    return x_center, y_center, z_center


def _generate_template_centerline(t, length, template_data):
    """
    Generate points by scaling and interpolating a 2D template file.
    
    Parameters:
    -----------
    t : np.ndarray
        Arc length parameter values
    length : float
        Target length
    template_data : np.ndarray
        Template data with shape (n, 2) containing (x, z) coordinates
    
    Returns:
    --------
    tuple
        (x, y, z) coordinate arrays
    """
    # Extract template coordinates (XZ plane)
    x_template = template_data[:, 0]
    z_template = template_data[:, 1]
    
    # Calculate arc length of template
    dx = np.diff(x_template)
    dz = np.diff(z_template)
    segment_lengths = np.sqrt(dx**2 + dz**2)
    template_arc_length = np.concatenate([[0], np.cumsum(segment_lengths)])
    template_total_length = template_arc_length[-1]
    
    # Scale template to match desired length
    scale_factor = length / template_total_length
    scaled_arc_length = template_arc_length * scale_factor
    scaled_x = scale_factor * x_template
    scaled_z = scale_factor * z_template
    
    # Create cubic splines for interpolation
    spline_x = UnivariateSpline(scaled_arc_length, scaled_x, s=0, k=3)
    spline_z = UnivariateSpline(scaled_arc_length, scaled_z, s=0, k=3)
    
    # Interpolate at desired points
    x_center = spline_x(t)
    z_center = spline_z(t)
    y_center = np.zeros_like(t)
    
    return x_center, y_center, z_center


def generate_centerline_points(length=10000.0, num_points=501, 
                               centerline_type='straight',
                               curve_radius=10000.0, 
                               sine_frequency=1.0, sine_amplitude=1000.0,
                               template_file='template.csv'):
    """
    Generate centerline points for various path types.
    
    Parameters:
    -----------
    length : float
        Total length of centerline (Å)
    num_points : int
        Number of points to generate
    centerline_type : str
        Type of centerline: 'straight', 'curve', 'sinusoidal', or '2Dtemplate'
    curve_radius : float
        Radius of curvature for 'curve' type
    sine_frequency : float
        Frequency for 'sinusoidal' type
    sine_amplitude : float
        Amplitude for 'sinusoidal' type
    template_file : str
        Path to CSV file for '2Dtemplate' type
    
    Returns:
    --------
    np.ndarray
        Array of shape (num_points, 3) with (x, y, z) coordinates
    """
    # Arc length parameter from 0 to length
    t = np.linspace(0, length, num_points)
    
    # Generate centerline based on type
    if centerline_type == 'straight':
        x_center, y_center, z_center = _generate_straight_centerline(t)
    elif centerline_type == 'curve':
        x_center, y_center, z_center = _generate_curved_centerline(t, length, curve_radius)
    elif centerline_type == 'sinusoidal':
        x_center, y_center, z_center = _generate_sinusoidal_centerline(t, length, sine_frequency, sine_amplitude)
    elif centerline_type == '2Dtemplate':
        template_data = read_2d_csv(template_file)
        x_center, y_center, z_center = _generate_template_centerline(t, length, template_data)
    else:
        raise ValueError(f"Invalid centerline_type: {centerline_type}. "
                        "Must be 'straight', 'curve', 'sinusoidal', or '2Dtemplate'")
    
    # Combine into single array
    points = np.column_stack([x_center, y_center, z_center]).astype(np.float32)
    
    return points


def generate_cilia_structure(
    length=15000.0, 
    centerline_type='straight',
    curve_radius=10000.0, 
    sine_frequency=1.0, 
    sine_amplitude=1000.0,
    template_file='template.csv',
    num_doublets=9, 
    cilia_radius=875.0,
    max_interval=20.0
):
    """
    Generate complete cilia structure with centerline and doublet positions.
    
    Parameters:
    -----------
    length : float
        Total length of cilia (Å)
    centerline_type : str
        Type of centerline path
    curve_radius : float
        Radius for curved centerline
    sine_frequency : float
        Frequency for sinusoidal centerline
    sine_amplitude : float
        Amplitude for sinusoidal centerline
    template_file : str
        Path to template file for 2D template centerline
    num_doublets : int
        Number of doublet microtubules
    cilia_radius : float
        Radial distance from centerline to doublets
    max_interval : float
        Maximum spacing between centerline points
    
    Returns:
    --------
    dict
        Structure containing:
        - 'centerline': np.ndarray of centerline points
        - 'doublets': list of doublet info dicts
        - 'num_doublets': int
        - 'centerline_type': str
        - 'length': float
        - 'cilia_radius': float
    """
    num_points = int(length / max_interval) + 1

    # Generate centerline
    centerline = generate_centerline_points(
        length=length,
        num_points=num_points,
        centerline_type=centerline_type,
        curve_radius=curve_radius,
        sine_frequency=sine_frequency,
        sine_amplitude=sine_amplitude,
        template_file=template_file
    )
    
    # Calculate doublet positions
    doublets = []
    for i in range(num_doublets):
        angle, shift_dist = _calculate_radial_position(
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
    Calculate doublet centerline by radially shifting the cilia centerline.
    
    Uses the Frenet-Serret frame (Tangent, Normal, Binormal) to ensure
    the shift is perpendicular to the centerline at every point.
    
    Parameters:
    -----------
    cilia_centerline : np.ndarray
        Main cilia centerline points, shape (n, 3)
    angle : float
        Angle in degrees for radial position
    shift_distance : float
        Radial distance from centerline
    
    Returns:
    --------
    np.ndarray
        Shifted centerline points for the doublet, shape (n, 3)
    """
    n_points = len(cilia_centerline)
    doublet_centerline = np.zeros_like(cilia_centerline)
    angle_rad = np.radians(angle)
    
    # Calculate tangent vectors at all points
    tangents = np.gradient(cilia_centerline, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1)[:, np.newaxis]

    # Reference vector for normal calculation
    UP_VECTOR = np.array([0.0, 1.0, 0.0])
    
    for i in range(n_points):
        tangent = tangents[i]
        
        # Calculate normal vector (perpendicular to tangent)
        normal_est = np.cross(tangent, UP_VECTOR)
        
        # Handle case where tangent is parallel to UP_VECTOR
        if np.linalg.norm(normal_est) < 1e-6:
            normal_est = np.cross(tangent, np.array([1.0, 0.0, 0.0]))
        
        normal = normal_est / np.linalg.norm(normal_est)
        
        # Calculate binormal (perpendicular to both tangent and normal)
        binormal = np.cross(tangent, normal)
        
        # Calculate shift vector in the normal-binormal plane
        shift_vector = shift_distance * (
            np.cos(angle_rad) * normal + np.sin(angle_rad) * binormal
        )
        
        # Apply shift
        doublet_centerline[i] = cilia_centerline[i] + shift_vector
    
    return doublet_centerline