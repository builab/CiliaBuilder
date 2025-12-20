# geometry/centerline.py

import numpy as np
from scipy.interpolate import UnivariateSpline

# --- Helper Functions for Centerline Types ---

def calculate_radial_position(centerline_points, doublet_index, 
                                total_doublets=9, cilia_radius=875.0):
    """Calculate angle for this doublet (evenly distributed around 360°)"""
    angle = (360.0 / total_doublets) * doublet_index
    
    return angle, cilia_radius

def _generate_straight_centerline(t):
    """Generate points for a straight centerline along the Z-axis."""
    x_center = np.zeros_like(t)
    y_center = np.zeros_like(t)
    z_center = t
    return x_center, y_center, z_center

def _generate_curved_centerline(t, length, curve_radius):
    """Generate points for a curved centerline (arc) in the X-Z plane."""
    if curve_radius == 0:
         raise ValueError("curve_radius cannot be zero for 'curve' type.")
    
    # s = R * theta => theta = s / R
    total_angle_rad = length / curve_radius
    theta = (t / length) * total_angle_rad
    
    # Curve starts at (0, 0, 0) and curves into the positive X-Z quadrant
    x_center = curve_radius * (1 - np.cos(theta))
    y_center = np.zeros_like(t)
    z_center = curve_radius * np.sin(theta)
    return x_center, y_center, z_center

def _generate_sinusoidal_centerline(t, length, sine_frequency, sine_amplitude):
    """Generate points for a sinusoidal centerline with a smooth lower transition."""
    num_points = len(t)
    
    # Constants used in the original calculation
    CILIA_RADIUS_REF = 1000.0
    
    # Calculate Z-offset (used to shift the curve so its start is at Z=0)
    angle_rad = np.arctan2(sine_amplitude, length / (sine_frequency * 4))
    z_offset = np.cos(np.pi/2 - angle_rad) * (2 * CILIA_RADIUS_REF)
    straight_length = z_offset / 2 # End of pure straight section
    
    # Initialize arrays
    x_center = np.zeros(num_points)
    y_center = np.zeros(num_points)
    z_center = np.zeros(num_points)
    
    # 1. Sinusoidal base curve (x_upper) and shifted z (z_upper)
    x_upper_base = sine_amplitude * np.sin(sine_frequency * 2 * np.pi * t / length)
    z_upper_shifted = t + z_offset
    
    # 2. Define Masks for three regions: Lower, Upper, Interpolation
    lower_mask = t <= straight_length
    upper_mask = t >= z_offset
    interp_mask = ~(lower_mask | upper_mask) # The middle transition region
    
    # Lower part: Straight segment (0, 0, t)
    z_center[lower_mask] = t[lower_mask]
    # x_center/y_center remain 0

    # Upper part: Sinusoidal
    x_center[upper_mask] = x_upper_base[upper_mask]
    z_center[upper_mask] = z_upper_shifted[upper_mask]
    
    # Interpolation region: Smooth transition
    if np.any(interp_mask):
        t_interp = t[interp_mask]
        
        # Interpolation parameter (0 at straight_length, 1 at z_offset)
        alpha = (t_interp - straight_length) / (z_offset - straight_length)
        
        # Upper endpoint values (start of pure sinusoidal region)
        x_upper_start = x_upper_base[interp_mask]
        z_upper_start = z_upper_shifted[interp_mask]
        
        # Linear interpolation of X and Z coordinates
        x_center[interp_mask] = alpha * x_upper_start # X starts at 0 (alpha=0)
        z_center[interp_mask] = (1 - alpha) * straight_length + alpha * z_upper_start
    
    return x_center, y_center, z_center

def _generate_template_centerline(t, length, template_file):
    """Generate points by scaling and interpolating a 2D template file (x, z)."""
    
    # --- 1. Load and Validate Template Data ---
    try:
        # Load with a potential header skip
        template_data = np.loadtxt(template_file, delimiter=',', skiprows=1)
    except Exception:
        # If loading with skiprows=1 fails (e.g., no header), try loading without skipping
        try:
            template_data = np.loadtxt(template_file, delimiter=',', skiprows=0)
        except Exception as e:
             raise ValueError(f"Could not load template file '{template_file}'. Check format.") from e

    # Validation checks
    if len(template_data) < 100:
        raise ValueError(f"Template file must have at least 100 points, got {len(template_data)}")
    if template_data.ndim == 1 or template_data.shape[1] != 2:
        raise ValueError(f"Template file must have 2 columns (x, y/z), got {template_data.shape[1]}")

    # Extract template coordinates (assuming template is in the XZ plane)
    x_template = template_data[:, 0]
    z_template = template_data[:, 1]
    y_template = np.zeros_like(x_template)

    # --- 2. Calculate Arc Length and Scaling ---
    # Calculate cumulative arc length of the template
    dx = np.diff(x_template)
    dz = np.diff(z_template)
    segment_lengths = np.sqrt(dx**2 + dz**2)
    template_arc_length = np.concatenate([[0], np.cumsum(segment_lengths)])
    template_total_length = template_arc_length[-1]
    
    # Scale factor to match desired length
    scale_factor = length / template_total_length

    # Scale the template coordinates and the arc length parameter for interpolation
    scaled_arc_length = template_arc_length * scale_factor
    scaled_x = scale_factor * x_template
    scaled_z = scale_factor * z_template

    # --- 3. Interpolation ---
    # Create smooth splines for each coordinate (using arc length as parameter)
    # s=0 ensures exact interpolation through control points
    spline_x = UnivariateSpline(scaled_arc_length, scaled_x, s=0, k=3)
    spline_z = UnivariateSpline(scaled_arc_length, scaled_z, s=0, k=3)
    
    # Interpolate for the desired number of points (t is the new arc length)
    x_center = spline_x(t)
    z_center = spline_z(t)
    y_center = np.zeros_like(t)

    return x_center, y_center, z_center


# --- Main Function: generate_centerline_points ---
def generate_centerline_points(length=10000.0, num_points=501, 
                               centerline_type='straight',
                               curve_radius=10000.0, 
                               sine_frequency=1.0, sine_amplitude=1000.0,
                               template_file='template.csv'):
    """
    Generate centerline points for straight, curved, sinusoidal, or template paths.
    
    Returns:
    --------
    points : numpy.ndarray
        Array of shape (num_points, 3) containing (x, y, z) coordinates
    """
    
    # t represents the normalized arc length parameter from 0 to length
    t = np.linspace(0, length, num_points)
    
    # Dispatch based on centerline type
    if centerline_type == 'straight':
        x_center, y_center, z_center = _generate_straight_centerline(t)
    elif centerline_type == 'curve':
        x_center, y_center, z_center = _generate_curved_centerline(t, length, curve_radius)
    elif centerline_type == 'sinusoidal':
        x_center, y_center, z_center = _generate_sinusoidal_centerline(t, length, sine_frequency, sine_amplitude)
    elif centerline_type == 'template':
        x_center, y_center, z_center = _generate_template_centerline(t, length, template_file)
    else:
        raise ValueError("centerline_type must be 'straight', 'curve', 'sinusoidal', or 'template'")
    
    # Combine into array
    points = np.column_stack([x_center, y_center, z_center]).astype(np.float32)
    
    return points


# --- Main Function: generate_cilia_structure ---
def generate_cilia_structure(
    length=15000.0, 
    centerline_type='straight',
    curve_radius=10000.0, 
    sine_frequency=1.0, 
    sine_amplitude=1000.0,
    template_file='template.csv',
    num_doublets=9, 
    cilia_radius=875.0,
    max_interval=20.0 # Use a descriptive default parameter
):
    """
    Generate complete cilia structure with centerline and doublet position information.
    """
    num_points = int(length / max_interval) + 1

    # 1. Generate centerline for central pair
    centerline = generate_centerline_points(
        length=length,
        num_points=num_points,
        centerline_type=centerline_type,
        curve_radius=curve_radius,
        sine_frequency=sine_frequency,
        sine_amplitude=sine_amplitude,
        template_file=template_file
    )
    
    # 2. Calculate radial position for each doublet
    doublets = []
    for i in range(num_doublets):
        # NOTE: calculate_radial_position is assumed to be defined in geometry/base.py
        # and returns angle in DEGREES and shift_distance.
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


# --- Function: get_doublet_centerline (Improved Tangent Calculation) ---
def get_doublet_centerline(cilia_centerline, angle, shift_distance):
    """
    Calculate the centerline for a doublet microtubule by shifting the cilia centerline.
    Uses the Frenet-Serret frame method (Tangent, Normal, Binormal) to ensure
    the shift is radial and perpendicular to the centerline at every point.
    """
    n_points = len(cilia_centerline)
    doublet_centerline = np.zeros_like(cilia_centerline)
    angle_rad = np.radians(angle)
    
    # Pre-calculate tangent vectors for all points
    tangents = np.gradient(cilia_centerline, axis=0)
    tangents /= np.linalg.norm(tangents, axis=1)[:, np.newaxis]

    # Use a global reference vector (Y-axis) for the initial cross product (Normal estimate)
    UP_VECTOR = np.array([0.0, 1.0, 0.0])
    
    for i in range(n_points):
        tangent = tangents[i]
        
        # 1. Find a stable Normal vector
        # Cross product of T and UP_VECTOR gives the Binormal (or a vector perpendicular to T)
        normal_est = np.cross(tangent, UP_VECTOR)
        
        # If tangent is parallel to UP_VECTOR (e.g., vertical section)
        if np.linalg.norm(normal_est) < 1e-6:
             # Use a different reference vector (X-axis)
             normal_est = np.cross(tangent, np.array([1.0, 0.0, 0.0]))

        normal = normal_est / np.linalg.norm(normal_est)
        
        # 2. Calculate the Binormal (perpendicular to T and N)
        binormal = np.cross(tangent, normal)
        # Note: binormal is already normalized since T and N are orthogonal unit vectors
        
        # 3. Calculate shift vector: combination of Normal and Binormal, rotated by angle
        # The (N, B) plane is the cross-sectional plane of the cilia
        shift_vector = shift_distance * (np.cos(angle_rad) * normal + np.sin(angle_rad) * binormal)
        
        # 4. Apply shift
        doublet_centerline[i] = cilia_centerline[i] + shift_vector
    
    return doublet_centerline