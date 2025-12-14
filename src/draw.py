"""
General function to draw multiple tubules (doublet or triplet)
"""
import numpy as np
from chimerax.core.models import Surface
from chimerax.surface import calculate_vertex_normals

def create_tube_geometry(path_points, radius=10.0, segments=16, capped=True):
    """
    Create tube geometry from a path of points.
    
    Parameters:
    -----------
    path_points : numpy.ndarray
        Array of shape (N, 3) containing the path points
    radius : float
        Radius of the tube (default: 10.0)
    segments : int
        Number of segments around the tube circumference (default: 16)
    capped : bool
        Whether to add caps at the ends (default: True)
    
    Returns:
    --------
    vertices : numpy.ndarray
        Array of vertex positions
    triangles : numpy.ndarray
        Array of triangle indices
    """
    n_points = len(path_points)
    
    # Create circle points around the path
    theta = np.linspace(0, 2*np.pi, segments, endpoint=False)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)
    
    vertices = []
    
    # For each point along the path
    for i in range(n_points):
        if i == 0:
            # Use next point for direction
            tangent = path_points[i+1] - path_points[i]
        elif i == n_points - 1:
            # Use previous point for direction
            tangent = path_points[i] - path_points[i-1]
        else:
            # Use average of forward and backward
            tangent = path_points[i+1] - path_points[i-1]
        
        tangent = tangent / np.linalg.norm(tangent)
        
        # Create perpendicular vectors
        if abs(tangent[2]) < 0.9:
            up = np.array([0, 0, 1])
        else:
            up = np.array([1, 0, 0])
        
        normal = np.cross(tangent, up)
        normal = normal / np.linalg.norm(normal)
        binormal = np.cross(tangent, normal)
        
        # Create circle of vertices around this point
        for j in range(segments):
            vertex = path_points[i] + circle_x[j] * normal + circle_y[j] * binormal
            vertices.append(vertex)
    
    vertices = np.array(vertices, dtype=np.float32)
    
    # Create triangles connecting the rings
    triangles = []
    for i in range(n_points - 1):
        for j in range(segments):
            # Current ring
            v0 = i * segments + j
            v1 = i * segments + (j + 1) % segments
            # Next ring
            v2 = (i + 1) * segments + j
            v3 = (i + 1) * segments + (j + 1) % segments
            
            # Two triangles per quad
            triangles.append([v0, v2, v1])
            triangles.append([v1, v2, v3])
    
    # Add end caps if requested
    if capped:
        # Store current vertex count
        base_vertex_count = len(vertices)
        
        # Add center vertices for caps
        start_center_idx = base_vertex_count
        end_center_idx = base_vertex_count + 1
        
        vertices = np.vstack([vertices, [path_points[0]], [path_points[-1]]])
        
        # Start cap (pointing inward, so reverse winding)
        for j in range(segments):
            v0 = start_center_idx
            v1 = j
            v2 = (j + 1) % segments
            triangles.append([v0, v2, v1])
        
        # End cap (pointing outward)
        last_ring_start = (n_points - 1) * segments
        for j in range(segments):
            v0 = end_center_idx
            v1 = last_ring_start + j
            v2 = last_ring_start + (j + 1) % segments
            triangles.append([v0, v1, v2])
    
    triangles = np.array(triangles, dtype=np.int32)
    
    return vertices, triangles


def generate_tube_surface(session, path_points, radius=10.0, segments=16, 
                         color=(255, 255, 0, 255), name="tube", capped=True, add_to_session=True):
    """
    Generate a tube surface model in ChimeraX.
    
    Parameters:
    -----------
    session : chimerax.core.session.Session
        ChimeraX session
    path_points : numpy.ndarray
        Array of shape (N, 3) containing the path points
    radius : float
        Radius of the tube
    segments : int
        Number of segments around circumference
    color : tuple
        RGBA color (0-255)
    name : str
        Name of the surface model
    capped : bool
        Whether to add caps at the ends (default: True)
    add_to_session : bool
        Whether to add the surface to session immediately (default: True)
    """
    # Create geometry
    vertices, triangles = create_tube_geometry(path_points, radius, segments, capped)
    
    # Calculate normals
    normals = calculate_vertex_normals(vertices, triangles)
    
    # Create surface model
    surf = Surface(name, session)
    surf.set_geometry(vertices, normals, triangles)
    
    # Set color (convert to 0-255 uint8 array)
    color_array = np.array(color, dtype=np.uint8)
    surf.color = color_array
    
    # Add to session if requested
    if add_to_session:
        session.models.add([surf])
    
    return surf


def generate_centerline_points(length=1500.0, interval=160, start_point=(0, 0, 0)):
    """
    Generate points along a straight center line in the Z direction.
    
    Parameters:
    -----------
    length : float
        Total length of the center line (default: 1500.0)
    interval : float
        Distance between consecutive points (default: 160)
    start_point : tuple
        Starting point (x, y, z) coordinates (default: (0, 0, 0))
    
    Returns:
    --------
    points : numpy.ndarray
        Array of shape (N, 3) containing the (x, y, z) coordinates
    """
    # Calculate number of points
    num_points = int(length / interval) + 1
    
    # Generate z coordinates
    z_coords = np.linspace(start_point[2], start_point[2] + length, num_points)
    
    # x and y remain constant (straight line in z direction)
    x_coords = np.full(num_points, start_point[0])
    y_coords = np.full(num_points, start_point[1])
    
    # Combine into a single array
    points = np.column_stack([x_coords, y_coords, z_coords])
    
    return points


def shift_path_perpendicular(path_points, shift_distance, angle_degrees):
    """
    Shift a path perpendicular to its axis by a given distance and angle.
    
    Parameters:
    -----------
    path_points : numpy.ndarray
        Array of shape (N, 3) containing the original path points
    shift_distance : float
        Distance to shift perpendicular to the path
    angle_degrees : float
        Angle in degrees to determine shift direction (0-360)
    
    Returns:
    --------
    shifted_points : numpy.ndarray
        Array of shifted path points
    """
    n_points = len(path_points)
    shifted_points = np.zeros_like(path_points)
    
    # Convert angle to radians
    angle_rad = np.radians(angle_degrees)
    
    for i in range(n_points):
        # Calculate tangent vector
        if i == 0:
            tangent = path_points[i+1] - path_points[i]
        elif i == n_points - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        
        tangent = tangent / np.linalg.norm(tangent)
        
        # Create perpendicular vectors
        if abs(tangent[2]) < 0.9:
            up = np.array([0, 0, 1])
        else:
            up = np.array([1, 0, 0])
        
        normal = np.cross(tangent, up)
        normal = normal / np.linalg.norm(normal)
        binormal = np.cross(tangent, normal)
        
        # Calculate shift direction based on angle
        shift_vector = shift_distance * (np.cos(angle_rad) * normal + np.sin(angle_rad) * binormal)
        
        # Apply shift
        shifted_points[i] = path_points[i] + shift_vector
    
    return shifted_points


COLOR_A_TUBULE = (255, 100, 100, 255)  # Reddish
COLOR_B_TUBULE = (100, 100, 255, 255)  # Blueish
COLOR_C_TUBULE = (100, 255, 100, 255)  # Greenish
COLOR_CP_TUBULE = (100, 255, 255, 255)  # Cyanish


def draw_tubules(session, 
                 length=1500, 
                 interval=80, 
                 angle=0,
                 radii=None,
                 shift_distances=None,
                 length_diffs=None,
                 tubule_names=None,
                 colors=None,
                 group_name="tubules"):
    """
    Draw multiple tubules (doublet, triplet, etc.) shifted perpendicular to centerline.
    
    Parameters:
    -----------
    session : chimerax.core.session.Session
        ChimeraX session
    length : float
        Base length for all tubules (default: 1500)
    interval : float
        Distance between path points (default: 80)
    angle : float
        Direction angle in degrees for shifting all tubules (default: 0)
        This defines the perpendicular direction from the centerline
    radii : list of float
        Radius for each tubule. If None, uses [125, 130] for doublet or [125, 130, 135] for triplet
    shift_distances : list of float
        Distance from centerline for each tubule along the angle direction
        Positive = shift in +angle direction, Negative = shift in opposite direction
        (default: [+70, -70] for doublet meaning opposite sides)
    length_diffs : list of float
        Length adjustment for each tubule (+/- from base length). 
        E.g., [0, -5, +10] means tubule1=length, tubule2=length-5, tubule3=length+10
        (default: [0, -5] for doublet or [0, -5, -10] for triplet)
    tubule_names : list of str
        Name for each tubule (default: ["A", "B"] or ["A", "B", "C"])
    colors : list of tuple
        RGBA color for each tubule (default: preset colors)
    group_name : str
        Name for the group containing all tubules (default: "tubules")
    
    Returns:
    --------
    list : List of tube models
    """
    
    # Determine number of tubules from the first provided list parameter
    n_tubules = None
    for param in [radii, shift_distances, length_diffs, tubule_names, colors]:
        if param is not None:
            n_tubules = len(param)
            break
    
    # Default to doublet (2 tubules) if nothing specified
    if n_tubules is None:
        n_tubules = 2
    
    # Set defaults for all parameters
    if radii is None:
        if n_tubules == 2:
            radii = [125, 130]
        elif n_tubules == 3:
            radii = [125, 130, 135]
        else:
            radii = [125 + i*5 for i in range(n_tubules)]
    
    if shift_distances is None:
        if n_tubules == 2:
            shift_distances = [70, -70]  # Opposite sides
        elif n_tubules == 3:
            shift_distances = [70, 0, -70]  # Spread across
        else:
            # Distribute symmetrically around center
            shift_distances = [70 - (140 / (n_tubules - 1)) * i for i in range(n_tubules)]
    
    if length_diffs is None:
        if n_tubules == 2:
            length_diffs = [0, -5]
        elif n_tubules == 3:
            length_diffs = [0, -5, -10]
        else:
            length_diffs = [-i*5 for i in range(n_tubules)]
    
    if tubule_names is None:
        default_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
        tubule_names = default_names[:n_tubules]
    
    if colors is None:
        default_colors = [
            COLOR_A_TUBULE,   # Red
            COLOR_B_TUBULE,   # Blue
            COLOR_C_TUBULE,   # Green
            (255, 255, 100, 255),  # Yellow
            (255, 100, 255, 255),  # Magenta
            (100, 255, 255, 255),  # Cyan
        ]
        colors = default_colors[:n_tubules]
    
    # Validate all lists have the same length
    if not all(len(lst) == n_tubules for lst in [radii, shift_distances, 
                                                   length_diffs, tubule_names, colors]):
        raise ValueError(f"All parameter lists must have length {n_tubules}")
    
    # Create tubules
    surfs = []
    
    for i in range(n_tubules):
        # Calculate actual length for this tubule
        tubule_length = length + length_diffs[i]
        
        # Generate center line path for this tubule
        center_path = generate_centerline_points(tubule_length, interval, start_point=(0, 0, 0))
        
        # Shift the path perpendicular to axis using the single angle
        path = shift_path_perpendicular(center_path, shift_distances[i], angle)
        
        # Create tubule surface
        tube = generate_tube_surface(session, path,
                                     radius=radii[i],
                                     segments=32,
                                     color=colors[i],
                                     name=tubule_names[i],
                                     capped=True,
                                     add_to_session=False)
        surfs.append(tube)
    
    # Add all as a group
    session.models.add_group(surfs, name=group_name)
    
    # Log information
    info_lines = [f"Created {n_tubules}-tubule group \"{group_name}\" at angle {angle}°:"]
    for i in range(n_tubules):
        info_lines.append(f"  {tubule_names[i]}: length={length + length_diffs[i]}, "
                         f"radius={radii[i]}, shift={shift_distances[i]:+.1f}")
    session.logger.info("\n".join(info_lines))
    
    return surfs


def draw_mt(session, length=1500, interval=80, radius=125, name="singlet", color=COLOR_A_TUBULE):
    """
    Draw a single microtubule (singlet) along the centerline.
    """
    # Calculate path points along centerline (no shift)
    path_points = generate_centerline_points(length, interval, start_point=(0, 0, 0))
    
    # Create tube
    tube = generate_tube_surface(session, path_points,
                                 radius=radius,
                                 segments=32,
                                 color=color,
                                 name=name,
                                 capped=True)
    session.logger.info(f"Created singlet \"{name}\" (length={length}, radius={radius})")
    return tube


def draw_doublet(session, length=3500, interval=80, angle=0, 
                radius_a=125, radius_b=130, shift_distance=70,
                length_diff=5,
                name="doublet", color_a=COLOR_A_TUBULE, color_b=COLOR_B_TUBULE):
    """
    Draw a microtubule doublet (convenience wrapper for draw_tubules).
    """
    return draw_tubules(
        session=session,
        length=length,
        interval=interval,
        angle=angle,
        radii=[radius_a, radius_b],
        shift_distances=[shift_distance, -shift_distance],  # Opposite sides
        length_diffs=[0, -length_diff],
        tubule_names=["A", "B"],
        colors=[color_a, color_b],
        group_name=name
    )


def draw_cp(session, length=3500, interval=80, angle=0, 
            radius_a=125, radius_b=125, shift_distance=160,
            length_diff=0,
            name="central_pair", color_a=COLOR_CP_TUBULE, color_b=COLOR_CP_TUBULE):
    """
    Draw central pair (convenience wrapper for draw_tubules).
    """
    return draw_tubules(
        session=session,
        length=length,
        interval=interval,
        angle=angle,
        radii=[radius_a, radius_b],
        shift_distances=[shift_distance, -shift_distance],  # Opposite sides
        length_diffs=[0, -length_diff],
        tubule_names=["C1", "C2"],
        colors=[color_a, color_b],
        group_name=name
    )


def draw_triplet(session, length=3000, interval=80, angle=0,
                radii=None, shift_distances=None,
                length_diffs=None,
                name="triplet", colors=None):
    """
    Draw a microtubule triplet (A, B, C tubules).
    
    Default configuration:
    - A-tubule: length=length, radius=125, shift=+70
    - B-tubule: length=length-5, radius=130, shift=0 (center)
    - C-tubule: length=length-10, radius=135, shift=-70
    """
    if radii is None:
        radii = [125, 130, 130]
    if shift_distances is None:
        shift_distances = [140, 0, -140]
    if length_diffs is None:
        length_diffs = [0, -5, -1000]
    if colors is None:
        colors = [COLOR_A_TUBULE, COLOR_B_TUBULE, COLOR_C_TUBULE]
    
    return draw_tubules(
        session=session,
        length=length,
        interval=interval,
        angle=angle,
        radii=radii,
        shift_distances=shift_distances,
        length_diffs=length_diffs,
        tubule_names=["A", "B", "C"],
        colors=colors,
        group_name=name
    )


# Example usage:
"""
# Singlet (single microtubule along centerline)
draw_mt(session, length=1500, radius=125, name="center_mt")

# Doublet at 0° angle (default X direction)
draw_doublet(session, length=1500, angle=0)

# Doublet at 45° angle
draw_doublet(session, length=1500, angle=45)

# Central pair at 90° angle
draw_cp(session, length=1500, angle=90)

# Triplet at 30° angle
draw_triplet(session, length=1500, angle=30)

# Custom configuration - 3 tubules all same length
draw_tubules(session, 
            length=2000,
            angle=0,
            radii=[125, 130, 130],
            shift_distances=[80, 0, -80],  # One side, center, other side
            length_diffs=[0, 0, 0],  # All same length
            tubule_names=["MT1", "MT2", "MT3"],
            colors=[(255,0,0,255), (0,255,0,255), (0,0,255,255)],
            group_name="custom_triplet")

# Or draw a singlet using draw_tubules with single-element lists
draw_tubules(session,
            length=1500,
            radii=[125],
            shift_distances=[0],  # No shift = centerline
            length_diffs=[0],
            tubule_names=["Center"],
            colors=[(255,255,0,255)],
            group_name="singlet")

# The angle defines the perpendicular direction:
# angle=0° means shift along X-axis
# angle=90° means shift along Y-axis
# angle=45° means shift along XY diagonal
# Positive shift_distance = in the +angle direction
# Negative shift_distance = in the -angle direction (opposite side)
"""