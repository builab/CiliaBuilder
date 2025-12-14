"""
Draw tube for simulation
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


COLOR_CENTER = (255, 255, 0, 255)
COLOR_A_TUBULE = (255, 100, 100, 255)  # Reddish
COLOR_B_TUBULE = (100, 100, 255, 255)  # Blueish
COLOR_CP_TUBULE = (100, 255, 255, 255)  # Blueish

NAME_CENTER = 'center'


def draw_mt(session, length=1500, interval=80, radius=125, name=NAME_CENTER, color=COLOR_CENTER):
    """
    Draw a singlet microtubule.
    """
    # Calculate path points
    path_points = generate_centerline_points(length, interval, start_point=(0, 0, 0))
    
    # Create tube
    tube = generate_tube_surface(session, path_points,
                                 radius=radius,
                                 segments=32,
                                 color=color,
                                 name=name,
                                 capped=True)
    session.logger.info(f"Created tube \"{name}\"")
    return tube


def draw_doublet(session, length=1500, interval=80, angle=0, 
                radius_a=125, radius_b=130, shift_distance=70,
                length_diff=5,
                name="doublet", color_a=COLOR_A_TUBULE, color_b=COLOR_B_TUBULE):
    """
    Draw a microtubule doublet (A-tubule and B-tubule).
    
    Parameters:
    -----------
    session : chimerax.core.session.Session
        ChimeraX session
    length : float
        Length of the A-tubule (default: 1500)
    interval : float
        Distance between path points (default: 80)
    angle : float
        Angle in degrees for doublet orientation (default: 0)
    radius_a : float
        Radius of A-tubule (default: 125)
    radius_b : float
        Radius of B-tubule (default: 130)
    shift_distance : float
        Distance between center and tubules (default: 70)
    length_diff : float
        Length difference: B-tubule will be shorter by this amount (default: 5)
    name : str
        Base name for the doublet (default: "doublet")
    color_a : tuple
        RGBA color for A-tubule (default: reddish)
    color_b : tuple
        RGBA color for B-tubule (default: blueish)
    
    Returns:
    --------
    tuple : (tube_a, tube_b)
        The two tube models
    """
    # Generate center line path for A-tubule
    center_path_a = generate_centerline_points(length, interval, start_point=(0, 0, 0))
    
    # Generate center line path for B-tubule (shorter)
    center_path_b = generate_centerline_points(length - length_diff, interval, start_point=(0, 0, 0))
    
    # Generate A-tubule path (shifted by +shift_distance)
    path_a = shift_path_perpendicular(center_path_a, shift_distance, angle)
    
    # Generate B-tubule path (shifted by -shift_distance, opposite direction)
    path_b = shift_path_perpendicular(center_path_b, -shift_distance, angle)
    
    # Create list to hold surfaces
    surfs = []
    
    # Create A-tubule (don't add to session yet)
    tube_a = generate_tube_surface(session, path_a,
                                   radius=radius_a,
                                   segments=32,
                                   color=color_a,
                                   name=f"{name}_A",
                                   capped=True,
                                   add_to_session=False)
    surfs.append(tube_a)
    
    # Create B-tubule (don't add to session yet)
    tube_b = generate_tube_surface(session, path_b,
                                   radius=radius_b,
                                   segments=32,
                                   color=color_b,
                                   name=f"{name}_B",
                                   capped=True,
                                   add_to_session=False)
    surfs.append(tube_b)
    
    # Add both as a group
    session.models.add_group(surfs, name=name)
    
    session.logger.info(f"Created doublet \"{name}\" with A-tubule (length={length}, radius={radius_a}) and B-tubule (length={length - length_diff}, radius={radius_b})")
    
    return tube_a, tube_b

def draw_cp(session, length=1500, interval=80, angle=0, 
                radius_a=125, radius_b=125, shift_distance=200,
                length_diff=0,
                name="central_pair", color_a=COLOR_CP_TUBULE, color_b=COLOR_CP_TUBULE):
    # Draw central pair but calling draw_doublet
    cp_tubes = draw_doublet(session, 
                           length=length, 
                           interval=interval, 
                           angle=angle, 
                           radius_a=radius_a, 
                           radius_b=radius_b, 
                           shift_distance=shift_distance,
                           length_diff=length_diff,
                           name=name, 
                           color_a=color_a, 
                           color_b=color_b)
                
    return cp_tubes

def register_command(logger):
    from chimerax.core.commands import CmdDesc, register, FloatArg, StringArg, Color8Arg
    
    # Register draw_mt command
    desc_mt = CmdDesc(
        keyword = [('length', FloatArg),
                   ('interval', FloatArg),
                   ('radius', FloatArg),
                   ('name', StringArg),
                   ('color', Color8Arg)],
        synopsis = 'Draw single microtubule tube'
    )
    register('draw_mt', desc_mt, draw_mt, logger=logger)
    
    # Register draw_doublet command
    desc_doublet = CmdDesc(
        keyword = [('length', FloatArg),
                   ('interval', FloatArg),
                   ('angle', FloatArg),
                   ('radius_a', FloatArg),
                   ('radius_b', FloatArg),
                   ('shift_distance', FloatArg),
                   ('length_diff', FloatArg),
                   ('name', StringArg),
                   ('color_a', Color8Arg),
                   ('color_b', Color8Arg)],
        synopsis = 'Draw doublet (A and B tubules)'
    )
    register('draw_doublet', desc_doublet, draw_doublet, logger=logger)
    
    # Register draw_cp command
    desc_cp = CmdDesc(
        keyword = [('length', FloatArg),
                   ('interval', FloatArg),
                   ('angle', FloatArg),
                   ('radius_a', FloatArg),
                   ('radius_b', FloatArg),
                   ('shift_distance', FloatArg),
                   ('length_diff', FloatArg),
                   ('name', StringArg),
                   ('color_a', Color8Arg),
                   ('color_b', Color8Arg)],
        synopsis = 'Draw cp (C1 and C2 tubules)'
    )
    register('draw_cp', desc_cp, draw_cp, logger=logger)


# Register commands
#register_command(session.logger)

# Example usage:
# draw_mt length 1500 interval 80 radius 125 name centerline
# draw_doublet length 1500 interval 80 angle 45 radius_a 125 radius_b 130 length_diff 5 name mt_doublet