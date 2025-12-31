# geometry/centerline.py

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import splprep, splev
from ..io import read_2d_csv

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
    return np.zeros_like(t), np.zeros_like(t), t


def _generate_curved_centerline(t, length, curve_radius):
    """
    Generate points for a curved centerline (arc) in the X-Z plane.
    Uses true arc-length parameterization.
    Arc length s = R * theta => theta = s / R.
    
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
        raise ValueError("curve_radius cannot be zero")
    
    # theta = s / R. Since t is our arc length parameter [0, length]:
    theta = t / curve_radius
    
    x_center = curve_radius * (1 - np.cos(theta))
    y_center = np.zeros_like(t)
    z_center = curve_radius * np.sin(theta)
    return x_center, y_center, z_center


def _generate_sinusoidal_centerline(t, length, sine_frequency, sine_amplitude):
    """
    Improved Sinusoid:
    1. Quintic taper for 15% straight base (C2 smoothness).
    2. Numerical integration to ensure TOTAL path length = length.
    
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
    # 1. Create a high-resolution temporary z-axis
    # We use a larger range initially because the sine wave 'eats up' arc length
    # High-resolution temp axis for length calculation
    num_fine = max(2000, len(t))
    z_temp = np.linspace(0, length, num_fine)
    
    # 1. Define shape with 15% straight transition
    phase = sine_frequency * 2 * np.pi * z_temp / length
    t_norm = np.clip(z_temp / (length * 0.15), 0, 1)
    taper = 6*t_norm**5 - 15*t_norm**4 + 10*t_norm**3 # Quintic ramp
    
    x_temp = (sine_amplitude * np.sin(phase)) * taper
    
    # 2. Calculate true cumulative arc length (s)
    dx = np.gradient(x_temp)
    dz = np.gradient(z_temp)
    ds = np.sqrt(dx**2 + dz**2)
    s_accumulated = np.cumsum(ds)
    s_accumulated -= s_accumulated[0] # Zero the start
    
    # 3. Scale structure so the final s matches 'length'
    scale = length / s_accumulated[-1]
    x_final = x_temp * scale
    z_final = z_temp * scale
    s_final = s_accumulated * scale
    
    # 4. Interpolate to uniform t points (0 to length)
    x_center = np.interp(t, s_final, x_final)
    y_center = np.zeros_like(t)
    z_center = np.interp(t, s_final, z_final)
    
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
    interval=20.0
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
    interval : float
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
    num_points = int(length / interval) + 1

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
        angle = -360 / num_doublets * i # Correct doublet order by using -
        doublets.append({
            'index': i,
            'angle': angle,
            'shift_distance': cilia_radius,
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
    # Fit a B-Spline for analytical derivative calculation
    tck, u = splprep([cilia_centerline[:,0], cilia_centerline[:,1], cilia_centerline[:,2]], s=0)
    
    # Sample at 5x density to ensure the shift vector is calculated smoothly
    u_fine = np.linspace(0, 1, n_points * 5)
    points_fine = np.array(splev(u_fine, tck)).T
    derivs = np.array(splev(u_fine, tck, der=1)).T
    
    # Generate the Frenet-Serret Frame (Tangent, Normal, Binormal)
    tangents = derivs / np.linalg.norm(derivs, axis=1)[:, np.newaxis]
    up = np.array([0, 1, 0])
    normals = np.cross(tangents, up)
    
    # Safety check for segments aligned with the up-vector
    norms = np.linalg.norm(normals, axis=1)
    mask = norms < 1e-6
    if np.any(mask): 
        normals[mask] = np.cross(tangents[mask], [1, 0, 0])
        norms[mask] = np.linalg.norm(normals[mask], axis=1)
    normals /= norms[:, np.newaxis]
    binormals = np.cross(tangents, normals)
    
    # Apply the radial shift
    angle_rad = np.radians(angle)
    shift_vec = shift_distance * (np.cos(angle_rad) * normals + np.sin(angle_rad) * binormals)
    shifted_points = points_fine + shift_vec
    
    # Resample by the arc length of the SHIFTED curve to prevent jagged bunching
    diffs = np.diff(shifted_points, axis=0)
    dist = np.concatenate(([0], np.cumsum(np.linalg.norm(diffs, axis=1))))
    
    # Final smoothing (s=0.1) filters out numerical jitter on high-compression inner curves
    tck_final, _ = splprep([shifted_points[:,0], shifted_points[:,1], shifted_points[:,2]], u=dist, s=0.1)
    u_final = np.linspace(0, dist[-1], n_points)
    
    return np.array(splev(u_final, tck_final)).T