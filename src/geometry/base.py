# geometry/base.py
def calculate_doublet_angle(doublet_index, total_doublets=9, offset_angle=0):
    """Calculate angle for a doublet (evenly distributed around 360°)"""
    return (360.0 / total_doublets) * doublet_index + offset_angle

def calculate_radial_position(centerline_points, doublet_index, 
                                total_doublets=9, cilia_radius=875.0):
    """Calculate angle for this doublet (evenly distributed around 360°)"""
    angle = (360.0 / total_doublets) * doublet_index
    
    return angle, cilia_radius