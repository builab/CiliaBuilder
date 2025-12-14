# This script generate line at high-level by generating marker, then tube
# Easy to understand but not so practical

import numpy as np

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
    
    
def idstr2tuple(model_id_string):
    """
    Convert model ID string to tuple
    """
    model_id = tuple(int(i) for i in model_id_string.lstrip("#").split('.'))
    return model_id
    
def model_id_to_string(model_id):
    """
    Convert model ID tuple to string format.
    Examples:
    (1,) -> "#1"
    (1, 2) -> "#1.2"
    (1, 2, 3) -> "#1.2.3"
    """
    return ".".join(str(i) for i in model_id)
    
def pts2markerset(session, set_name, points, radius = 20, color = (255,255,0,255)):
    """
    Making a marker set from arrays of points
    Return model_index
    """
    from chimerax.markers import MarkerSet
    marker_set = MarkerSet(session, name = set_name)
    markers = [marker_set.create_marker(pts, color, radius) for pts in points]
    session.models.add([marker_set])
    # Get the last model
    model = session.models[-1]
    #print(model.color)
    session.logger.status('Create %s markerset with %d points' % (set_name, len(markers)), log = True)
    return model.id

def get_model_id(session, model_name):
# Script to get model name from model ID string like #1.2.1
    #print(model_id)
    for model in session.models:
        print(model.name)
        if model.name == model_name:
            print(f'Found model name {model.name}')
            return model.id
    print("ERROR: No model id found")
    return None

def get_model_from_id(session, model_id):
# Script to get model name from model.id tuple
    #print(model_id)
    for model in session.models:
        if model.id == model_id:
            return model
    print("ERROR: No model found with #%s" % model_id_to_string(model_id))
    return None
    
COLOR_CENTER = (255,255,0,255)
NAME_CENTER = 'center'
str_COLOR_CENTER = ",".join(tuple(str(x) for x in COLOR_CENTER))


def generate_tube(session, set_name, model_id, tube_radius=125, color = (255,255,0,255)):
    """
    Making the tube base on the markerset
    model_id is the tuple
    """
    from chimerax.core.commands import run
    str_color = ",".join(tuple(str(x) for x in color))
    tube_name = f"{set_name}_tube"
    run(session, f'shape tube #%s radius %f color %s name %s' % (model_id_to_string(model_id), tube_radius, str_color, tube_name), log=False)
    #run(session, f'hide #%s target a' % model_id_to_string(model_id), log=False)
    model = get_model_from_id(session, model_id)
    model.display = False
  

def drawmt(session, length=1500, interval=80, radius=125, name=NAME_CENTER, color=COLOR_CENTER):
    points = generate_centerline_points(length, interval, start_point=(0, 0, 0))
    model_id = pts2markerset(session, name, points, color=color)
    generate_tube(session, name, model_id, radius, color=color)


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

#drawmt length 1500 interval 80 radius 125 name centerline 