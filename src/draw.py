"""
General functions to draw tubules, membranes, and geometric primitives for cilia/centriole visualization.
"""

import numpy as np
from chimerax.core.models import Surface
from chimerax.surface import calculate_vertex_normals


def _calculate_local_frame(tangent):
    """
    Calculate local coordinate frame (normal, binormal) for a tangent vector.
    
    Parameters:
    -----------
    tangent : np.ndarray
        Tangent vector (will be normalized)
    
    Returns:
    --------
    tuple
        (normal, binormal) orthonormal vectors
    """
    tangent = tangent / np.linalg.norm(tangent)
    
    up = np.array([0.0, 1.0, 0.0])
    
    # If tangent is parallel to up, use X-axis instead
    if np.linalg.norm(np.cross(tangent, up)) < 1e-6:
        up = np.array([1.0, 0.0, 0.0])
    
    normal = np.cross(tangent, up)
    normal = normal / np.linalg.norm(normal)
    binormal = np.cross(tangent, normal)
    
    return normal, binormal


def _create_ellipsoid_geometry(center, radii, segments_u, segments_v):
    """
    Create vertices and triangles for an ellipsoid.
    
    Parameters:
    -----------
    center : ndarray
        Center point (x, y, z)
    radii : ndarray
        Semi-axes lengths (rx, ry, rz)
    segments_u : int
        Number of longitudinal segments
    segments_v : int
        Number of latitudinal segments
    
    Returns:
    --------
    vertices : ndarray
        Vertex positions (N, 3)
    triangles : ndarray
        Triangle indices (M, 3)
    """
    vertices = []
    
    # Generate vertices using spherical coordinates, scaled by radii
    for i in range(segments_v + 1):
        theta = np.pi * i / segments_v  # Latitude angle (0 to π)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        
        for j in range(segments_u):
            phi = 2 * np.pi * j / segments_u  # Longitude angle (0 to 2π)
            sin_phi = np.sin(phi)
            cos_phi = np.cos(phi)
            
            # Parametric equation for ellipsoid
            x = center[0] + radii[0] * sin_theta * cos_phi
            y = center[1] + radii[1] * sin_theta * sin_phi
            z = center[2] + radii[2] * cos_theta
            
            vertices.append([x, y, z])
    
    vertices = np.array(vertices, dtype=np.float32)
    
    # Generate triangles
    triangles = []
    for i in range(segments_v):
        for j in range(segments_u):
            # Current vertex indices
            current = i * segments_u + j
            next_u = i * segments_u + (j + 1) % segments_u
            next_v = (i + 1) * segments_u + j
            next_both = (i + 1) * segments_u + (j + 1) % segments_u
            
            # Two triangles per quad
            triangles.append([current, next_v, next_u])
            triangles.append([next_u, next_v, next_both])
    
    triangles = np.array(triangles, dtype=np.int32)
    
    return vertices, triangles


def _create_tube_geometry(path_points, radius=10.0, segments=16, capped=True):
    """
    Create tube geometry from a path of points.
    
    Parameters:
    -----------
    path_points : np.ndarray
        Array of shape (n, 3) defining the tube path
    radius : float
        Tube radius
    segments : int
        Number of segments around circumference
    capped : bool
        Whether to add end caps
    
    Returns:
    --------
    tuple
        (vertices, triangles) as numpy arrays
    """
    n_points = len(path_points)
    
    if n_points < 2:
        return np.empty((0, 3), dtype=np.float32), np.empty((0, 3), dtype=np.int32)

    # Create circle points around the path
    theta = np.linspace(0, 2*np.pi, segments, endpoint=False)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)
    
    vertices = []
    
    # Generate vertices along path
    for i in range(n_points):
        # Calculate tangent vector
        if i == 0:
            tangent = path_points[i+1] - path_points[i]
        elif i == n_points - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        
        # Calculate local frame
        normal, binormal = _calculate_local_frame(tangent)
        
        # Create ring of vertices
        for j in range(segments):
            vertex = path_points[i] + circle_x[j] * normal + circle_y[j] * binormal
            vertices.append(vertex)
    
    vertices = np.array(vertices, dtype=np.float32)
    
    # Create triangles connecting the rings
    triangles = []
    for i in range(n_points - 1):
        for j in range(segments):
            v0 = i * segments + j
            v1 = i * segments + (j + 1) % segments
            v2 = (i + 1) * segments + j
            v3 = (i + 1) * segments + (j + 1) % segments
            
            triangles.append([v0, v2, v1])
            triangles.append([v1, v2, v3])
    
    # Add end caps if requested
    if capped:
        base_vertex_count = len(vertices)
        start_center_idx = base_vertex_count
        end_center_idx = base_vertex_count + 1
        
        vertices = np.vstack([vertices, [path_points[0]], [path_points[-1]]])
        
        # Start cap (reverse winding for inward normal)
        for j in range(segments):
            v0 = start_center_idx
            v1 = j
            v2 = (j + 1) % segments
            triangles.append([v0, v2, v1])
        
        # End cap (outward normal)
        last_ring_start = (n_points - 1) * segments
        for j in range(segments):
            v0 = end_center_idx
            v1 = last_ring_start + j
            v2 = last_ring_start + (j + 1) % segments
            triangles.append([v0, v1, v2])
    
    triangles = np.array(triangles, dtype=np.int32)
    
    return vertices, triangles
    

def _create_ellipsoid_cap_geometry(center, radius_xy, radius_z, segments_u, segments_v, hemisphere='lower'):
    """
    Create vertices and triangles for an ellipsoid hemisphere cap.
    
    Parameters:
    -----------
    center : ndarray
        Center point of the cap (where it attaches to cylinder)
    radius_xy : float
        Radius in X and Y directions (matches cylinder radius)
    radius_z : float
        Radius in Z direction (e.g., radius/2 for flattened cap)
    segments_u : int
        Number of longitudinal segments
    segments_v : int
        Number of latitudinal segments
    hemisphere : str
        'lower' for bottom cap (negative Z), 'upper' for top cap (positive Z)
    
    Returns:
    --------
    vertices : ndarray
        Vertex positions (N, 3)
    triangles : ndarray
        Triangle indices (M, 3)
    """
    vertices = []
    
    # Determine theta range based on hemisphere
    if hemisphere == 'lower':
        theta_start = np.pi / 2  # Equator
        theta_end = np.pi         # South pole
        z_offset = -radius_z      # Offset center down by radius_z
    else:  # 'upper'
        theta_start = 0           # North pole
        theta_end = np.pi / 2     # Equator
        z_offset = radius_z       # Offset center up by radius_z
    
   
    # Generate vertices using spherical coordinates, scaled by radii
    for i in range(segments_v + 1):
        # Map i to theta range
        t = i / segments_v
        theta = theta_start + t * (theta_end - theta_start)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        
        for j in range(segments_u):
            phi = 2 * np.pi * j / segments_u
            sin_phi = np.sin(phi)
            cos_phi = np.cos(phi)
            
            # Parametric equation for ellipsoid
            x = center[0] + radius_xy * sin_theta * cos_phi
            y = center[1] + radius_xy * sin_theta * sin_phi
            z = center[2] + radius_z * cos_theta
            
            vertices.append([x, y, z])
    
    vertices = np.array(vertices, dtype=np.float32)
    
    # Generate triangles
    triangles = []
    for i in range(segments_v):
        for j in range(segments_u):
            current = i * segments_u + j
            next_u = i * segments_u + (j + 1) % segments_u
            next_v = (i + 1) * segments_u + j
            next_both = (i + 1) * segments_u + (j + 1) % segments_u
            
            triangles.append([current, next_v, next_u])
            triangles.append([next_u, next_v, next_both])
    
    triangles = np.array(triangles, dtype=np.int32)
    
    return vertices, triangles
    
    
def _create_ring_cap_geometry(center, direction, inner_radius, outer_radius, segments=32):
    """
    Create cylindrical ring (annulus) cap geometry.
    
    Parameters:
    -----------
    center : np.ndarray
        Center point of the ring
    direction : np.ndarray
        Direction vector for ring normal
    inner_radius : float
        Inner radius
    outer_radius : float
        Outer radius
    segments : int
        Number of segments around circumference
    
    Returns:
    --------
    tuple
        (vertices, triangles) as numpy arrays
    """
    # Normalize direction vector
    direction = direction / np.linalg.norm(direction)
    
    # Create orthonormal basis
    up = np.array([0.0, 1.0, 0.0])
    if np.linalg.norm(np.cross(direction, up)) < 1e-6:
        up = np.array([1.0, 0.0, 0.0])
    
    normal = np.cross(direction, up)
    normal = normal / np.linalg.norm(normal)
    binormal = np.cross(direction, normal)
    
    # Generate vertices
    theta = np.linspace(0, 2*np.pi, segments, endpoint=False)
    vertices = []
    
    # Inner ring
    for angle in theta:
        x = inner_radius * np.cos(angle)
        y = inner_radius * np.sin(angle)
        vertex = center + x * normal + y * binormal
        vertices.append(vertex)
    
    # Outer ring
    for angle in theta:
        x = outer_radius * np.cos(angle)
        y = outer_radius * np.sin(angle)
        vertex = center + x * normal + y * binormal
        vertices.append(vertex)
    
    vertices = np.array(vertices, dtype=np.float32)
    
    # Create triangles
    triangles = []
    for i in range(segments):
        next_i = (i + 1) % segments
        
        v0 = i
        v1 = next_i
        v2 = segments + i
        v3 = segments + next_i
        
        triangles.append([v0, v2, v1])
        triangles.append([v1, v2, v3])
    
    triangles = np.array(triangles, dtype=np.int32)
    
    return vertices, triangles
    

def _combine_geometries(geom_list):
    """
    Combine multiple (vertices, triangles) tuples into a single geometry.
    
    Parameters:
    -----------
    geom_list : list
        List of (vertices, triangles) tuples
    
    Returns:
    --------
    tuple
        Combined (vertices, triangles)
    """
    all_vertices = []
    all_triangles = []
    vertex_count = 0
    
    for vertices, triangles in geom_list:
        if len(vertices) > 0:
            all_vertices.append(vertices)
            all_triangles.append(triangles + vertex_count)
            vertex_count += len(vertices)
            
    if len(all_vertices) == 0:
        return np.empty((0, 3), dtype=np.float32), np.empty((0, 3), dtype=np.int32)
        
    return np.vstack(all_vertices), np.vstack(all_triangles)
    
    
def generate_centerline_points(length=1500.0, interval=160, start_point=(0, 0, 0)):
    """
    Generate points along a straight centerline in the Z direction.
    
    Parameters:
    -----------
    length : float
        Total length
    interval : float
        Spacing between points
    start_point : tuple
        Starting (x, y, z) coordinates
    
    Returns:
    --------
    np.ndarray
        Array of shape (n, 3) with centerline points
    """
    num_points = int(length / interval) + 1
    z_coords = np.linspace(start_point[2], start_point[2] + length, num_points)
    x_coords = np.full(num_points, start_point[0])
    y_coords = np.full(num_points, start_point[1])
    
    points = np.column_stack([x_coords, y_coords, z_coords])
    return points


def shift_path_perpendicular(path_points, shift_distance, angle_degrees):
    """
    Shift a path perpendicular to its axis.
    
    Uses stable reference frame to prevent kinks/twists in the shifted path.
    
    Parameters:
    -----------
    path_points : np.ndarray
        Array of path points, shape (n, 3)
    shift_distance : float
        Radial shift distance
    angle_degrees : float
        Angle for shift direction
    
    Returns:
    --------
    np.ndarray
        Shifted path points
    """
    n_points = len(path_points)
    shifted_points = np.zeros_like(path_points)
    angle_rad = np.radians(angle_degrees)
    
    for i in range(n_points):
        # Calculate tangent
        if i == 0:
            tangent = path_points[i+1] - path_points[i]
        elif i == n_points - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        
        tangent = tangent / np.linalg.norm(tangent)
        
        # Use stable reference vector
        up = np.array([0.0, 1.0, 0.0])
        if np.linalg.norm(np.cross(tangent, up)) < 1e-6:
            up = np.array([1.0, 0.0, 0.0])

        normal = np.cross(tangent, up)
        normal = normal / np.linalg.norm(normal)
        binormal = np.cross(tangent, normal)
        
        # Calculate and apply shift
        shift_vector = shift_distance * (np.cos(angle_rad) * normal + np.sin(angle_rad) * binormal)
        shifted_points[i] = path_points[i] + shift_vector
    
    return shifted_points


def generate_tube_surface(session, path_points, radius=10.0, segments=16, 
                         color=(255, 255, 0, 255), name="tube", capped=True, add_to_session=True):
    """
    Generate a tube surface model in ChimeraX.
    
    Parameters:
    -----------
    session : ChimeraX session
        Session object
    path_points : np.ndarray
        Path points for tube centerline
    radius : float
        Tube radius
    segments : int
        Number of circumferential segments
    color : tuple
        RGBA color (0-255)
    name : str
        Surface name
    capped : bool
        Whether to cap tube ends
    add_to_session : bool
        Whether to add to session immediately
    
    Returns:
    --------
    Surface or None
        Created surface model
    """
    vertices, triangles = _create_tube_geometry(path_points, radius, segments, capped)
    
    if len(vertices) == 0:
        session.logger.warning(f"Skipping '{name}': Path must have at least 2 points")
        return None
    
    normals = calculate_vertex_normals(vertices, triangles)
    
    surf = Surface(name, session)
    surf.set_geometry(vertices, normals, triangles)
    surf.color = np.array(color, dtype=np.uint8)
    
    if add_to_session:
        session.models.add([surf])
    
    return surf
       

def generate_ellipsoid_surface(session, center=(0.0, 0.0, 0.0), radii=(10.0, 8.0, 6.0), segments_u=32, segments_v=16,
                               color=(255, 255, 0, 255), name="ellipsoid", add_to_session=True):
    """
    Generate an ellipsoid surface model for a visualization session.
    
    Parameters:
    -----------
    session : ChimeraX session
        Session object
    center : tuple
        Center of the ellipsoid (x, y, z)
    radii : tuple
        Semi-axes lengths (rx, ry, rz) along x, y, z axes
    segments_u : int
        Number of longitudinal segments
    segments_v : int
        Number of latitudinal segments
    color : tuple
        RGBA color (0-255)
    name : str
        Surface name
    add_to_session : bool
        Whether to add to session immediately
    
    Returns:
    --------
    Surface or None
        Created surface model
    """
    center = np.array(center)
    radii = np.array(radii)
    
    # Create geometry
    vertices, triangles = _create_ellipsoid_geometry(center, radii, segments_u, segments_v)
    
    # Check for valid geometry
    if len(vertices) == 0:
        session.logger.warning(f"Skipping surface creation for '{name}': No vertices generated.")
        return None
    
    # Calculate normals for ellipsoid
    # For an ellipsoid, the normal at point (x,y,z) is proportional to (x/rx^2, y/ry^2, z/rz^2)
    vectors_from_center = vertices - center
    # Scale by inverse square of radii
    normals_unnormalized = vectors_from_center / (radii ** 2)
    # Normalize
    normals = normals_unnormalized / np.linalg.norm(normals_unnormalized, axis=1)[:, np.newaxis]
    
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

        
def generate_capsule_surface(session, center=(0.0, 0.0, 0.0), capsule_length=25.0, radius=10.0, 
                             segments=32, color=(255, 255, 0, 255), name="capsule", 
                             flat_end_z='low', add_to_session=True):
    """
    Generate a capsule surface model (cylinder with two hemispherical caps) for cap complex
    in ChimeraX, centered around the Z-axis.
    
    Parameters:
    -----------
    session : ChimeraX session
        Session object
    center : tuple
        Center of the capsule (Z axis)
    capsule_length : float
        Length of the capsule in total
    radius : float
        Sphere radius
    segments : int
        Number of circumferential segments
    color : tuple
        RGBA color (0-255)
    name : str
        Surface name
    flat_end_z : str or None
        Which end(s) to flatten: 'low', 'high', 'both', or None
        - None: Both ends are spherical hemispheres (rx=ry=rz=radius)
        - 'low'/'high'/'both': Specified end(s) are flattened ellipsoids (rx=ry=radius, rz=radius/2)
    add_to_session : bool
        Whether to add to session immediately
    
    Returns:
    --------
    Surface or None
        Created surface model
    """
    
    center = np.array(center)
    
    # --- 1. Validate capsule length based on end cap types ---
    
    # Determine minimum length based on cap configuration
    if flat_end_z is None:
        # Both ends spherical: need 2*radius (full sphere height)
        min_length = 2 * radius
        cap_type_desc = "spherical hemispheres"
    elif flat_end_z == 'both':
        # Both ends flattened: need 2*(radius/2) = radius
        min_length = radius
        cap_type_desc = "flattened ellipsoids"
    elif flat_end_z in ['low', 'high']:
        # One spherical (radius), one flattened (radius/2): need 1.5*radius
        min_length = 1.5 * radius
        cap_type_desc = "mixed (one spherical, one flattened)"
    else:
        session.logger.error(f"Invalid flat_end_z value: '{flat_end_z}'. Must be None, 'low', 'high', or 'both'.")
        return None
    
    if capsule_length < min_length:
        session.logger.error(
            f"Capsule length ({capsule_length:.2f}) is less than minimum required length "
            f"({min_length:.2f}) for {cap_type_desc} caps. Cannot create capsule."
        )
        return None
    
    # --- 2. Calculate dimensions and centers ---
    
    # Determine cap heights
    if flat_end_z is None:
        cap_height_low = radius
        cap_height_high = radius
    elif flat_end_z == 'both':
        cap_height_low = radius / 2.0
        cap_height_high = radius / 2.0
    elif flat_end_z == 'low':
        cap_height_low = radius / 2.0
        cap_height_high = radius
    else:  # flat_end_z == 'high'
        cap_height_low = radius
        cap_height_high = radius / 2.0
    
    # Calculate cylinder length
    cylinder_centerline_length = capsule_length - cap_height_low - cap_height_high
    
    # Cylinder endpoints (cap attachment points)
    cyl_start_z = center[2] - cylinder_centerline_length / 2.0
    cyl_end_z = center[2] + cylinder_centerline_length / 2.0
    cap_center_low_z = np.array([center[0], center[1], cyl_start_z])
    cap_center_high_z = np.array([center[0], center[1], cyl_end_z])
    
    # Generate the centerline path for the cylinder
    path_points = generate_centerline_points(
        length=cylinder_centerline_length, 
        interval=cylinder_centerline_length, 
        start_point=cap_center_low_z
    )
    
    # --- 3. Generate Geometries ---
    
    geom_list = []
    segments_v = max(4, segments // 2)
    
    # Low Z End Cap (A)
    radius_z_low = radius / 2.0 if (flat_end_z == 'low' or flat_end_z == 'both') else radius
    
    vertices_a, triangles_a = _create_ellipsoid_cap_geometry(
        cap_center_low_z, 
        radius, 
        radius_z_low,
        segments_u=segments, 
        segments_v=segments_v,
        hemisphere='lower'
    )
    geom_list.append((vertices_a, triangles_a))
    
    # High Z End Cap (B)
    radius_z_high = radius / 2.0 if (flat_end_z == 'high' or flat_end_z == 'both') else radius
    
    vertices_b, triangles_b = _create_ellipsoid_cap_geometry(
        cap_center_high_z, 
        radius, 
        radius_z_high,
        segments_u=segments, 
        segments_v=segments_v,
        hemisphere='upper'
    )
    geom_list.append((vertices_b, triangles_b))
    
    # Generate Cylinder Body
    if cylinder_centerline_length > 0.0:
        vertices_cyl, triangles_cyl = _create_tube_geometry(path_points, radius, segments, capped=False)
        geom_list.append((vertices_cyl, triangles_cyl))
    
    # --- 4. Combine Geometries and Create Surface ---
    
    vertices, triangles = _combine_geometries(geom_list)
    
    if len(vertices) == 0:
        session.logger.warning(f"Skipping surface creation for '{name}': No geometry generated.")
        return None
    
    # Calculate normals for the combined geometry
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
        
    session.logger.info(
        f"Created capsule \"{name}\" (Total Length: {capsule_length:.2f}, Radius: {radius:.2f}, "
        f"Cylinder Length: {cylinder_centerline_length:.2f}, Flat End: {flat_end_z})"
    )
    
    return surf


def draw_tubules(session, 
                 length=None,
                 centerline_points=None,
                 interval=80, 
                 angle=0,
                 radii=None,
                 shift_distances=None,
                 length_diffs=None,
                 tubule_names=None,
                 colors=None,
                 group_name="tubules",
                 add_to_session=False):
    """
    Draw multiple tubules (doublet, triplet, etc.) shifted perpendicular to centerline.
    
    Parameters:
    -----------
    session : ChimeraX session
        Session object
    length : float, optional
        Base length for tubules (if centerline_points not provided)
    centerline_points : np.ndarray, optional
        Centerline path points
    interval : float
        Spacing between points
    angle : float
        Rotation angle in degrees
    radii : list, optional
        Radius for each tubule
    shift_distances : list, optional
        Radial shift for each tubule
    length_diffs : list, optional
        Length differences from base length
    tubule_names : list, optional
        Names for each tubule
    colors : list, optional
        RGBA colors for each tubule
    group_name : str
        Name for the tubule group
    add_to_session : bool
        Whether to add to session immediately
    
    Returns:
    --------
    list
        List of created Surface models
    """
    if length is None and centerline_points is None:
        raise ValueError("Must provide either 'length' or 'centerline_points'")
    
    # Determine number of tubules
    n_tubules = None
    for param in [radii, shift_distances, length_diffs, tubule_names, colors]:
        if param is not None:
            n_tubules = len(param)
            break
    
    if n_tubules is None:
        n_tubules = 2
    
    # Set defaults
    if radii is None:
        radii = [125 + i*5 for i in range(n_tubules)]
    
    if shift_distances is None:
        if n_tubules == 2:
            shift_distances = [70, -70]
        elif n_tubules == 3:
            shift_distances = [70, 0, -70]
        else:
            shift_distances = [70 - (140 / (n_tubules - 1)) * i for i in range(n_tubules)]
    
    if length_diffs is None:
        length_diffs = [-i*5 for i in range(n_tubules)]
    
    if tubule_names is None:
        default_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
        tubule_names = default_names[:n_tubules]
    
    if colors is None:
        default_colors = [
            (100, 100, 255, 255),  # Blue
            (100, 100, 255, 255),  # Blue
            (179, 179, 255, 255),  # Light blue
            (255, 255, 100, 255),  # Yellow
            (255, 100, 255, 255),  # Magenta
            (100, 255, 255, 255),  # Cyan
        ]
        colors = default_colors[:n_tubules]
    
    # Validate list lengths
    if not all(len(lst) == n_tubules for lst in [radii, shift_distances, length_diffs, tubule_names, colors]):
        raise ValueError(f"All parameter lists must have length {n_tubules}")
    
    # Create tubules
    surfs = []
    
    for i in range(n_tubules):
        # Determine centerline for this tubule
        if centerline_points is not None:
            center_path = centerline_points
            
            # Truncate based on length_diffs if length is provided
            if length is not None and length_diffs is not None:
                tubule_length = length + length_diffs[i]
    
                if len(center_path) >= 2:
                    segment_lengths = np.linalg.norm(np.diff(center_path, axis=0), axis=1)
                    cumulative_length = np.concatenate(([0], np.cumsum(segment_lengths)))
        
                    n_keep = np.searchsorted(cumulative_length, tubule_length, side='right')
                    n_keep = max(2, min(n_keep, len(center_path)))
                else:
                    n_keep = len(center_path)

                center_path_to_shift = center_path[:n_keep]
            else:
                center_path_to_shift = center_path
        else:
            tubule_length = length + length_diffs[i]
            center_path_to_shift = generate_centerline_points(tubule_length, interval, start_point=(0, 0, 0))
        
        # Shift path perpendicular
        path = shift_path_perpendicular(center_path_to_shift, shift_distances[i], angle)
        
        # Create tubule surface
        tube = generate_tube_surface(
            session, path,
            radius=radii[i],
            segments=32,
            color=colors[i],
            name=tubule_names[i],
            capped=True,
            add_to_session=False
        )
        if tube is not None:
            surfs.append(tube)
    
    # Add as group if requested
    if surfs and add_to_session:
        session.models.add_group(surfs, name=group_name)
    
    # Log information
    info_source = "from centerline" if centerline_points is not None else ""
    session.logger.info(f"Created {len(surfs)}-tubule group \"{group_name}\" {info_source} at angle {angle}°")
    
    return surfs


def draw_membrane(session, centerline_points, radius=1100.0, segments=32, 
                  color=(105, 105, 105, 255), name="membrane", capped=True,
                  membrane_thickness=40.0):
    """
    Draw double membrane (inner and outer) with optional ring caps.
    
    Parameters:
    -----------
    session : ChimeraX session
        Session object
    centerline_points : np.ndarray
        Membrane centerline path
    radius : float
        Inner membrane radius
    segments : int
        Number of circumferential segments
    color : tuple
        RGBA color (0-255)
    name : str
        Base name for membrane surfaces
    capped : bool
        Whether to add cylindrical ring caps
    membrane_thickness : float
        Thickness between inner and outer membranes
    
    Returns:
    --------
    list
        List of created membrane surfaces [inner, outer, start_cap, end_cap]
    """
    session.logger.info(f"Drawing membrane with {len(centerline_points)} points, "
                       f"Z range: {centerline_points[0][2]:.1f} to {centerline_points[-1][2]:.1f}")

    surfaces = []
    
    # Inner membrane
    inner_membrane = generate_tube_surface(
        session, centerline_points,
        radius=radius,
        segments=segments,
        color=color,
        name=f"{name}_inner",
        capped=False,
        add_to_session=False
    )
    
    if inner_membrane is not None:
        surfaces.append(inner_membrane)
        session.logger.info(f"Created inner membrane (radius={radius})")
    
    # Outer membrane
    outer_radius = radius + membrane_thickness
    outer_membrane = generate_tube_surface(
        session, centerline_points,
        radius=outer_radius,
        segments=segments,
        color=color,
        name=f"{name}_outer",
        capped=False,
        add_to_session=False
    )
    
    if outer_membrane is not None:
        surfaces.append(outer_membrane)
        session.logger.info(f"Created outer membrane (radius={outer_radius})")
    
    # Add ring caps if requested
    if capped and len(centerline_points) >= 2:
        # Start cap
        start_cap_vertices, start_cap_triangles = _create_ring_cap_geometry(
            center=centerline_points[0],
            direction=centerline_points[1] - centerline_points[0],
            inner_radius=radius,
            outer_radius=outer_radius,
            segments=segments
        )
        
        if len(start_cap_vertices) > 0:
            start_cap_normals = calculate_vertex_normals(start_cap_vertices, start_cap_triangles)
            start_cap = Surface(f"{name}_start_cap", session)
            start_cap.set_geometry(start_cap_vertices, start_cap_normals, start_cap_triangles)
            start_cap.color = np.array(color, dtype=np.uint8)
            surfaces.append(start_cap)
            session.logger.info(f"Created start cap")
        
        # End cap
        end_cap_vertices, end_cap_triangles = _create_ring_cap_geometry(
            center=centerline_points[-1],
            direction=centerline_points[-1] - centerline_points[-2],
            inner_radius=radius,
            outer_radius=outer_radius,
            segments=segments
        )
        
        if len(end_cap_vertices) > 0:
            end_cap_normals = calculate_vertex_normals(end_cap_vertices, end_cap_triangles)
            end_cap = Surface(f"{name}_end_cap", session)
            end_cap.set_geometry(end_cap_vertices, end_cap_normals, end_cap_triangles)
            end_cap.color = np.array(color, dtype=np.uint8)
            surfaces.append(end_cap)
            session.logger.info(f"Created end cap")
    
    return surfaces


def draw_ladders(session, centerline_points=None, angle=0.0, periodicity=320.0,
                shift_distances=None, radius=20.0, color=None, segments=16,
                name="LadderStructure", add_to_session=False):
    """
    Generate ladder-like structure connecting two shifted centerlines at regular intervals.
    
    Parameters:
    -----------
    session : ChimeraX session
        Session object
    centerline_points : np.ndarray
        Centerline coordinates
    angle : float
        Rotation angle in degrees for shift plane
    periodicity : float
        Distance between rungs along centerline
    shift_distances : list
        Must contain at least two values for ladder sides
    radius : float
        Radius of cylindrical rungs
    color : tuple
        RGBA color (0-255)
    segments : int
        Number of segments for rungs
    name : str
        Name for resulting surface
    add_to_session : bool
        Whether to add to session immediately
    
    Returns:
    --------
    Surface or None
        Combined ladder surface model
    """
    # Validate inputs
    if centerline_points is None or len(centerline_points) < 2:
        session.logger.error("Centerline points required for draw_ladders")
        return None

    if shift_distances is None or len(shift_distances) < 2:
        session.logger.error("shift_distances must contain at least two values")
        return None
    
    if periodicity <= 0:
        session.logger.error("Periodicity must be positive")
        return None

    # Calculate shifted centerlines
    shift1_dist = shift_distances[0]
    shift2_dist = shift_distances[1]

    shift_line_1 = shift_path_perpendicular(centerline_points, shift1_dist, angle)
    shift_line_2 = shift_path_perpendicular(centerline_points, shift2_dist, angle)
    
    # Find rung anchor points based on periodicity
    path_segments = np.linalg.norm(centerline_points[1:] - centerline_points[:-1], axis=1)
    cumulative_distance = np.insert(np.cumsum(path_segments), 0, 0.0)
    
    # Start first rung one period from base
    rung_distances = np.arange(radius*2, cumulative_distance[-1]-radius*2, periodicity)
    
    if len(rung_distances) == 0:
        session.logger.info("No rungs drawn (centerline too short for periodicity)")
        return None

    rung_indices = []
    for dist in rung_distances:
        idx = np.argmin(np.abs(cumulative_distance - dist))
        if not rung_indices or idx != rung_indices[-1]:
            rung_indices.append(idx)
           
    # Generate rung geometries
    rung_geometries = []

    for idx in rung_indices:
        start_pt = shift_line_1[idx]
        end_pt = shift_line_2[idx]
        rung_centerline = np.array([start_pt, end_pt], dtype=np.float32)

        rung_verts, rung_tris = _create_tube_geometry(
            rung_centerline, radius, segments, capped=True
        )

        if rung_verts is not None and rung_verts.shape[0] > 0:
            rung_geometries.append((rung_verts, rung_tris))
    
    if not rung_geometries:
        session.logger.info("No rungs successfully generated")
        return None

    # Combine geometries
    final_vertices, final_triangles = _combine_geometries(rung_geometries)
    
    # Create surface
    surface = Surface(name, session)
    final_normals = calculate_vertex_normals(final_vertices, final_triangles)
    surface.set_geometry(final_vertices, final_normals, final_triangles)

    if color:
        surface.color = np.array(color, dtype=np.uint8)
        
    if add_to_session:
        session.models.add([surface])

    session.logger.info(f"Generated ladder '{name}' with {len(rung_indices)} rungs")
    
    return surface