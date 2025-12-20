"""
General function to draw multiple tubules (doublet or triplet)
"""
import numpy as np
from chimerax.core.models import Surface
from chimerax.surface import calculate_vertex_normals

# --- Helper Function to Replicate Tube Frame Calculation (For Seamlessness) ---
def _calculate_local_frame(tangent):
    """Replicates the frame calculation inside create_tube_geometry."""
    tangent = tangent / np.linalg.norm(tangent)
    
    up = np.array([0.0, 1.0, 0.0]) 
    
    # If the tangent is near parallel to 'up', choose X-axis instead
    if np.linalg.norm(np.cross(tangent, up)) < 1e-6:
         up = np.array([1.0, 0.0, 0.0])
    
    normal = np.cross(tangent, up)
    normal = normal / np.linalg.norm(normal)
    binormal = np.cross(tangent, normal)
    
    return normal, binormal

# --- Helper Function to Generate Sphere Geometry ---
def create_sphere_geometry(center, radius=10.0, segments_u=32, segments_v=16):
    """
    Generate vertices and triangles for a sphere using spherical coordinates.
    ... (Existing code for create_sphere_geometry) ...
    """
    center = np.array(center)
    
    # Angles for latitude (phi, 0 to pi) and longitude (theta, 0 to 2*pi)
    phi = np.linspace(0, np.pi, segments_v + 1)
    theta = np.linspace(0, 2 * np.pi, segments_u, endpoint=False)

    vertices = []
    
    # Generate vertices
    for p in phi:
        for t in theta:
            x = center[0] + radius * np.sin(p) * np.cos(t)
            y = center[1] + radius * np.sin(p) * np.sin(t)
            z = center[2] + radius * np.cos(p)
            vertices.append([x, y, z])

    vertices = np.array(vertices, dtype=np.float32)
    
    # Generate triangles (quads)
    triangles = []
    
    # Number of vertices per 'ring' (segments_u)
    v_ring = segments_u
    
    # Iterate through each quad
    for i in range(segments_v):
        for j in range(segments_u):
            
            p1 = i * v_ring + j
            p2 = i * v_ring + (j + 1) % v_ring
            p3 = (i + 1) * v_ring + (j + 1) % v_ring
            p4 = (i + 1) * v_ring + j
            
            # Triangle 1
            triangles.append([p1, p4, p3])
            # Triangle 2
            triangles.append([p1, p3, p2])

    return vertices, np.array(triangles, dtype=np.int32)

def create_tube_geometry(path_points, radius=10.0, segments=16, capped=True):
    """
    Create tube geometry from a path of points.
    ... (Existing code for create_tube_geometry) ...
    """
    n_points = len(path_points)
    
    # --- ERROR PREVENTION FIX ---
    if n_points < 2:
        return np.empty((0, 3), dtype=np.float32), np.empty((0, 3), dtype=np.int32)
    # ----------------------------

    # Create circle points around the path
    theta = np.linspace(0, 2*np.pi, segments, endpoint=False)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)
    
    vertices = []
    
    # For each point along the path
    for i in range(n_points):
        # Calculate tangent vector
        if i == 0:
            tangent = path_points[i+1] - path_points[i]
        elif i == n_points - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        
        # Calculate local frame (same logic as in _calculate_local_frame)
        normal, binormal = _calculate_local_frame(tangent)
        
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

# --- NEW HELPER FUNCTION for Combining Geometry ---
def combine_geometries(geom_list):
    """Utility to combine multiple (vertices, triangles) tuples."""
    all_vertices = []
    all_triangles = []
    vertex_count = 0
    
    for vertices, triangles in geom_list:
        if len(vertices) > 0:
            all_vertices.append(vertices)
            # Offset triangle indices by the current total vertex count
            all_triangles.append(triangles + vertex_count)
            vertex_count += len(vertices)
            
    if len(all_vertices) == 0:
        return np.empty((0, 3), dtype=np.float32), np.empty((0, 3), dtype=np.int32)
        
    return np.vstack(all_vertices), np.vstack(all_triangles)


# --- NEW HELPER FUNCTION for Flat Cap Geometry (Modified to use frame) ---
def create_disk_geometry(center, radius=10.0, segments=32, frame=None):
    """
    Generate vertices and triangles for a flat circular disk, using a frame 
    to ensure alignment with the cylinder.
    
    Parameters:
    ...
    frame (tuple): The (normal, binormal) vectors defining the disk's plane and orientation.
    """
    center = np.array(center)
    
    if frame is None:
        # Default frame if called standalone (should not happen in capsule function)
        normal = np.array([1.0, 0.0, 0.0])
        binormal = np.array([0.0, 1.0, 0.0])
    else:
        normal, binormal = frame
        
    # Generate ring vertices
    theta = np.linspace(0, 2*np.pi, segments, endpoint=False)
    
    vertices = [center] # The center point is the first vertex (index 0)
    
    # Outer circle vertices (must use the same cosine/sine sequence as the cylinder)
    for angle in theta:
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        vertex = center + x * normal + y * binormal
        vertices.append(vertex)
        
    vertices = np.array(vertices, dtype=np.float32)
    
    # Create triangles (fan pattern from the center vertex)
    triangles = []
    center_idx = 0
    
    for i in range(segments):
        v1 = i + 1  # Index of current point on the circumference
        v2 = (i + 1) % segments + 1 # Index of next point (loops back)
        
        # Triangle (Center, V1, V2)
        # Note: Winding should be consistent (e.g., counter-clockwise from the center's perspective)
        triangles.append([center_idx, v1, v2])
        
    return vertices, np.array(triangles, dtype=np.int32)


# --- MODIFIED FUNCTION: CAPSULE SURFACE GENERATOR (Fixed for seamless flat end) ---
def generate_capsule_surface(session, center=(0.0, 0.0, 0.0), capsule_length=25.0, radius=10.0, 
                             segments=32, color=(255, 255, 0, 255), name="capsule", 
                             flat_end_z='low', add_to_session=True):
    """
    Generate a capsule surface model (cylinder with two hemispherical caps) 
    in ChimeraX, centered around the Z-axis, with an option for one end to be flat.
    """
    center = np.array(center)
    
    # --- 1. Calculate dimensions and centers ---
    
    if capsule_length < 2 * radius:
        session.logger.warning(f"Capsule length ({capsule_length}) is less than 2*radius ({2*radius}). Generating a sphere.")
        # Assuming generate_sphere_surface takes segments_v=segments//2
        return generate_sphere_surface(session, center, radius, segments, segments // 2, color, name, add_to_session)

    cylinder_centerline_length = capsule_length - 2 * radius
    half_cyl_centerline = cylinder_centerline_length / 2.0
    
    # Sphere centers
    sphere_center_a = center + np.array([0, 0, -capsule_length / 2.0 + radius]) 
    sphere_center_b = center + np.array([0, 0, capsule_length / 2.0 - radius])  
    
    # Cylinder Path Points (Cap attachment points)
    cyl_start_z = center[2] - half_cyl_centerline
    cyl_end_z = center[2] + half_cyl_centerline
    cap_center_low_z = np.array([center[0], center[1], cyl_start_z])
    cap_center_high_z = np.array([center[0], center[1], cyl_end_z])
    
    # Generate the centerline path for the cylinder
    path_points = generate_centerline_points(
        length=cylinder_centerline_length, 
        interval=cylinder_centerline_length, 
        start_point=cap_center_low_z
    )
    
    # --- Calculate cylinder frame (TNB) for seamless flat caps ---
    path_tangent = path_points[1] - path_points[0]
    cylinder_frame = _calculate_local_frame(path_tangent)
    
    # --- 2. Generate Geometries ---
    
    geom_list = []
    segments_v = max(4, segments // 2)
    
    # Low Z End Cap (A)
    if flat_end_z == 'low':
        vertices_a, triangles_a = create_disk_geometry(
            cap_center_low_z, 
            radius, 
            segments, 
            frame=cylinder_frame # Use cylinder's frame for seamless connection
        )
    else:
        # Rounded sphere cap
        vertices_a, triangles_a = create_sphere_geometry(
            sphere_center_a, radius, segments_u=segments, segments_v=segments_v
        )
        
    geom_list.append((vertices_a, triangles_a))
    
    # High Z End Cap (B)
    if flat_end_z == 'high':
        vertices_b, triangles_b = create_disk_geometry(
            cap_center_high_z, 
            radius, 
            segments, 
            frame=cylinder_frame # Use cylinder's frame for seamless connection
        )
    else:
        # Rounded sphere cap
        vertices_b, triangles_b = create_sphere_geometry(
            sphere_center_b, radius, segments_u=segments, segments_v=segments_v
        )

    geom_list.append((vertices_b, triangles_b))
    
    # Generate Cylinder Body
    if cylinder_centerline_length > 0.0:
        vertices_cyl, triangles_cyl = create_tube_geometry(path_points, radius, segments, capped=False)
        geom_list.append((vertices_cyl, triangles_cyl))

    # --- 3. Combine Geometries and Create Surface ---
    
    vertices, triangles = combine_geometries(geom_list)
    
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
        
    session.logger.info(f"Created capsule \"{name}\" (Total Length: {capsule_length}, Radius: {radius}, Flat End: {flat_end_z})")
    
    return surf

# --- Main Sphere Surface Generator ---
def generate_sphere_surface(session, center=(0.0, 0.0, 0.0), radius=10.0, segments_u=32, segments_v=16,
                            color=(255, 255, 0, 255), name="sphere", add_to_session=True):
    """
    Generate a sphere surface model for a visualization session (replicates tube structure).
    
    Note: Sphere vertex normals are calculated directly from the center, 
    simplifying the normal calculation step.
    """
    center = np.array(center)

    # Create geometry
    vertices, triangles = create_sphere_geometry(center, radius, segments_u, segments_v)
    
    # Check for valid geometry (always valid for sphere, but follows tube structure)
    if len(vertices) == 0:
        session.logger.warning(f"Skipping surface creation for '{name}': No vertices generated.")
        return None
    
    # Calculate normals (optimized for sphere: normal is the normalized vector from center to vertex)
    # The 'calculate_vertex_normals' helper is not needed
    vectors_from_center = vertices - center
    normals = vectors_from_center / np.linalg.norm(vectors_from_center, axis=1)[:, np.newaxis]
    
    # Create surface model
    surf = Surface(name, session)
    surf.set_geometry(vertices, normals, triangles)
    
    # Set color (convert to 0-255 uint8 array)
    color_array = np.array(color, dtype=np.uint8)
    surf.color = color_array
    
    # Add to session if requested
    if add_to_session:
        # Assumes session.models is a list-like object that accepts new models
        session.models.add([surf])
    
    return surf
    
def generate_tube_surface(session, path_points, radius=10.0, segments=16, 
                         color=(255, 255, 0, 255), name="tube", capped=True, add_to_session=True):
    """
    Generate a tube surface model in ChimeraX.
    """
    # Create geometry
    vertices, triangles = create_tube_geometry(path_points, radius, segments, capped)
    
    # If geometry creation failed due to too few points, skip surface creation
    if len(vertices) == 0:
        session.logger.warning(f"Skipping surface creation for '{name}': Path must have at least 2 points.")
        return None
    
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
    
    This function uses a stable reference vector (Y-axis) to calculate the Normal/Binormal 
    frame, which prevents sudden kinks/twists.
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
        
        # --- Kink Fix: Use a stable reference vector (Y-axis) ---
        up = np.array([0.0, 1.0, 0.0]) 
        
        # If the tangent is near parallel to 'up', choose X-axis instead
        if np.linalg.norm(np.cross(tangent, up)) < 1e-6:
             up = np.array([1.0, 0.0, 0.0])
        # ---------------------------------------------------------

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
    """
    
    # Validate input: need either length or centerline_points
    if length is None and centerline_points is None:
        raise ValueError("Must provide either 'length' or 'centerline_points'")
    
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
        # Determine centerline for this tubule
        center_path_to_shift = None
        
        if centerline_points is not None:
            # --- FIX: Truncate the provided centerline based on length_diffs ---
            center_path = centerline_points
            
            # This logic enables length_diffs when centerline_points is passed
            if length is not None and length_diffs is not None:
                tubule_length = length + length_diffs[i]
    
                # Calculate cumulative arc length along the centerline
                if len(center_path) >= 2:
                    segment_lengths = np.linalg.norm(np.diff(center_path, axis=0), axis=1)
                    cumulative_length = np.concatenate(([0], np.cumsum(segment_lengths)))
        
                    # Find the index where cumulative length exceeds tubule_length
                    n_keep = np.searchsorted(cumulative_length, tubule_length, side='right')
                    n_keep = max(2, min(n_keep, len(center_path)))  # Ensure valid range
                else:
                    n_keep = len(center_path)

                # Truncate the path to the correct length
                center_path_to_shift = center_path[:n_keep]
            else:
                center_path_to_shift = center_path
            # -----------------------------------------------------------------
        else:
            # Original logic: Calculate centerline from scratch (only used if centerline_points=None)
            tubule_length = length + length_diffs[i]
            center_path_to_shift = generate_centerline_points(tubule_length, interval, start_point=(0, 0, 0))
        
        # Shift the path perpendicular to axis using the single angle
        path = shift_path_perpendicular(center_path_to_shift, shift_distances[i], angle)
        
        # Create tubule surface
        tube = generate_tube_surface(session, path,
                                     radius=radii[i],
                                     segments=32,
                                     color=colors[i],
                                     name=tubule_names[i],
                                     capped=True,
                                     add_to_session=False)
        if tube is not None:
            surfs.append(tube)
    
    # Add all as a group if requested
    if surfs and add_to_session:
        session.models.add_group(surfs, name=group_name)
    
    # Log information
    if centerline_points is not None:
        info_lines = [f"Created {len(surfs)}-tubule group \"{group_name}\" from centerline at angle {angle}°:"]
    else:
        info_lines = [f"Created {len(surfs)}-tubule group \"{group_name}\" at angle {angle}°:"]
    
    if surfs:
        for i in range(len(surfs)):
            if centerline_points is not None:
                info_lines.append(f"  {tubule_names[i]}: radius={radii[i]}, shift={shift_distances[i]:+.1f}")
            else:
                info_lines.append(f"  {tubule_names[i]}: length={length + length_diffs[i]}, "
                                 f"radius={radii[i]}, shift={shift_distances[i]:+.1f}")
    session.logger.info("\n".join(info_lines))
    
    return surfs



def draw_membrane(session, centerline_points, radius=1100.0, segments=32, 
                  color=(105, 105, 105, 255), name="membrane", capped=True,
                  membrane_thickness=40.0):
    """
    Draw a double membrane (inner and outer) with optional cylindrical ring caps.
    
    Parameters:
    -----------
    session : ChimeraX session
        The ChimeraX session object
    centerline_points : numpy.ndarray
        Array of shape (N, 3) containing the path points for the membrane centerline
    radius : float
        Radius of the inner membrane tube (default: 1100.0)
    segments : int
        Number of segments around the tube circumference (default: 32)
    color : tuple
        RGBA color tuple (0-255 values) (default: semi-transparent gray)
    name : str
        Name for the membrane surface group (default: "membrane")
    capped : bool
        Whether to add cylindrical ring caps at the ends (default: False)
    membrane_thickness : float
        Thickness between inner and outer membranes in Angstroms (default: 40.0)
    
    Returns:
    --------
    surfaces : list
        List of created membrane surface models [inner, outer, start_cap, end_cap]
    """
    
    session.logger.info(f"draw_membrane received {len(centerline_points)} points, Z range: {centerline_points[0][2]:.1f} to {centerline_points[-1][2]:.1f}")

    surfaces = []
    
    # Create inner membrane (uncapped tube)
    inner_membrane = generate_tube_surface(
        session, 
        centerline_points,
        radius=radius,
        segments=segments,
        color=color,
        name=f"{name}_inner",
        capped=False,
        add_to_session=False
    )
    
    if inner_membrane is not None:
        surfaces.append(inner_membrane)
        session.logger.info(f"Created inner membrane \"{name}_inner\" (radius={radius})")
    
    # Create outer membrane (uncapped tube)
    outer_radius = radius + membrane_thickness
    outer_membrane = generate_tube_surface(
        session, 
        centerline_points,
        radius=outer_radius,
        segments=segments,
        color=color,
        name=f"{name}_outer",
        capped=False,
        add_to_session=False
    )
    
    if outer_membrane is not None:
        surfaces.append(outer_membrane)
        session.logger.info(f"Created outer membrane \"{name}_outer\" (radius={outer_radius})")
    
    # Add cylindrical ring caps if requested
    if capped and len(centerline_points) >= 2:
        # Start cap (ring at the beginning)
        start_cap_vertices, start_cap_triangles = create_ring_cap_geometry(
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
            session.logger.info(f"Created start cap \"{name}_start_cap\"")
        
        # End cap (ring at the end)
        end_cap_vertices, end_cap_triangles = create_ring_cap_geometry(
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
            session.logger.info(f"Created end cap \"{name}_end_cap\"")
    
    return surfaces


def create_ring_cap_geometry(center, direction, inner_radius, outer_radius, segments=32):
    """
    Create a cylindrical ring (annulus) cap geometry.
    
    Parameters:
    -----------
    center : numpy.ndarray
        Center point of the ring (shape: (3,))
    direction : numpy.ndarray
        Direction vector for the ring normal (shape: (3,))
    inner_radius : float
        Inner radius of the ring
    outer_radius : float
        Outer radius of the ring
    segments : int
        Number of segments around the ring circumference
    
    Returns:
    --------
    vertices : numpy.ndarray
        Array of vertex positions
    triangles : numpy.ndarray
        Array of triangle indices
    """
    # Normalize direction vector
    direction = direction / np.linalg.norm(direction)
    
    # Create orthonormal basis for the ring plane
    up = np.array([0.0, 1.0, 0.0])
    if np.linalg.norm(np.cross(direction, up)) < 1e-6:
        up = np.array([1.0, 0.0, 0.0])
    
    normal = np.cross(direction, up)
    normal = normal / np.linalg.norm(normal)
    binormal = np.cross(direction, normal)
    
    # Generate ring vertices
    theta = np.linspace(0, 2*np.pi, segments, endpoint=False)
    vertices = []
    
    # Inner ring vertices
    for angle in theta:
        x = inner_radius * np.cos(angle)
        y = inner_radius * np.sin(angle)
        vertex = center + x * normal + y * binormal
        vertices.append(vertex)
    
    # Outer ring vertices
    for angle in theta:
        x = outer_radius * np.cos(angle)
        y = outer_radius * np.sin(angle)
        vertex = center + x * normal + y * binormal
        vertices.append(vertex)
    
    vertices = np.array(vertices, dtype=np.float32)
    
    # Create triangles connecting inner and outer rings
    triangles = []
    for i in range(segments):
        next_i = (i + 1) % segments
        
        # Inner ring index
        v0 = i
        v1 = next_i
        
        # Outer ring index
        v2 = segments + i
        v3 = segments + next_i
        
        # Two triangles per quad (facing outward along direction)
        triangles.append([v0, v2, v1])
        triangles.append([v1, v2, v3])
    
    triangles = np.array(triangles, dtype=np.int32)
    
    return vertices, triangles



# Example usage:
"""
# Singlet (single microtubule along centerline)
draw_mt(session, length=1500, radius=125, name="center_mt")


# Custom configuration - 3 tubules all same length
draw_tubules(session, 
            length=2000,
            angle=0,
            radii=[125, 130, 135],
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