# Starting point
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
                         color=(255, 255, 0, 255), name="tube", capped=True):
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
    """
    # Create geometry
    vertices, triangles = create_tube_geometry(path_points, radius, segments, capped)
    
    # Calculate normals
    normals = calculate_vertex_normals(vertices, triangles)
    
    # Create surface model
    from chimerax.core.models import Surface
    surf = Surface(name, session)
    surf.set_geometry(vertices, normals, triangles)
    
    # Set color (convert to 0-255 uint8 array)
    color_array = np.array(color, dtype=np.uint8)
    surf.color = color_array
    
    # Add to session
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
    
COLOR_CENTER = (255,255,0,255)
NAME_CENTER = 'center'
str_COLOR_CENTER = ",".join(tuple(str(x) for x in COLOR_CENTER))

# To use in ChimeraX, run:
def drawmt(session, length=1500, interval=80, radius=125, name=NAME_CENTER, color=COLOR_CENTER):
    # Calculate number of points
    path_points = generate_centerline_points(length, interval, start_point=(0, 0, 0))

    # Create tube
    tube = generate_tube_surface(session, path_points,
                                 radius=radius,
                                 segments=num_points,
                                 color=COLOR_CENTER,
                                 name=name,
                                 capped=True)
    session.logger.info(f"Created tube \"{name}\"")
    return tube

def register_command(logger):
    from chimerax.core.commands import CmdDesc, register, FloatArg, StringArg, Color8Arg
    desc = CmdDesc(
        keyword = [('length', FloatArg),
                   ('interval', FloatArg),
                   ('radius', FloatArg),
                   ('name', StringArg),
                   ('color', Color8Arg)],
        synopsis = 'Draw tube based on a set of points'
    )
    register('drawmt', desc, drawmt, logger=logger)

register_command(session.logger)
