# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc, IntArg, FloatArg, StringArg, BoolArg

from .curve import generate_cilia_structure, get_doublet_centerline
from .draw import draw_tubules

# Define the optimal sampling interval for smoothness (must match curve.py)
MAX_INTERVAL = 10.0 # 10 Angstroms interval for centerline points
CILIA_OFFSET_ANGLE = 90.0 # Offset angle for cilia doublets


def ciliasim(session, 
            length=15000, 
            line='straight', # Simplified from centerline_type
            curve_radius=10000.0, 
            sine_frequency=2.0, sine_amplitude=500.0,
            # Cilia Structure Defaults
            num_doublets=9, 
            cilia_radius=875.0,
            draw_central_pair=True,
            # Doublet Geometry Defaults
            doublet_a_radius=125.0, # A-tubule radius
            doublet_b_radius=145.0, # B-tubule radius
            doublet_shift=70.0,     # A/B tubule shift from doublet centerline
            doublet_length_diff=250.0, # Length difference between A and B tubules (A - B)
            # Central Pair Geometry Defaults
            cp_radius=125.0,        # C1/C2 tubule radius
            cp_shift=160.0          # C1/C2 shift distance from cilia center
            ):
    """
    Generate and draw a complete cilia structure with doublet microtubules.
    
    Parameters:
    -----------
    session : chimerax session
        ChimeraX session
    length : float
        Length of the cilia in Angstroms (default: 15000). 
    line : str 
        Type of centerline: 'straight', 'curve', or 'sinusoidal' (default: 'straight')
    curve_radius : float
        Radius of curvature for 'curve' type (default: 10000.0)
    sine_frequency : float
        Frequency of sinusoidal oscillation (default: 2.0)
    sine_amplitude : float
        Amplitude of sinusoidal oscillation (default: 500.0)
        
    num_doublets : int
        Number of doublet microtubules (default: 9)
    cilia_radius : float
        Radius from center to doublets in Angstroms (default: 875.0)
    draw_central_pair : bool
        Whether to draw the central pair (default: True)
        
    doublet_a_radius : float
        Radius of the A-tubule (default: 125.0)
    doublet_b_radius : float
        Radius of the B-tubule (default: 145.0)
    doublet_shift : float
        Radial distance of A and B tubules from the doublet centerline (default: 70.0)
    doublet_length_diff : float
        Length difference between A and B tubules (A - B length, default: 250.0)
        
    cp_radius : float
        Radius of central pair (C1/C2) tubules (default: 125.0)
    cp_shift : float
        Distance of C1/C2 tubules from the cilia center line (default: 160.0)
    """
    
    centerline_type = line 
    
    session.logger.info(f"Generating cilia structure with {centerline_type} centerline...")
    
    # Calculate cilia structure points (num_points calculated internally in curve.py)
    structure = generate_cilia_structure(
        length=length,
        centerline_type=centerline_type,
        curve_radius=curve_radius,
        sine_frequency=sine_frequency,
        sine_amplitude=sine_amplitude,
        num_doublets=num_doublets,
        cilia_radius=cilia_radius
    )
    
    # Draw central pair if requested
    if draw_central_pair:
        session.logger.info("Drawing central pair (C1, C2)...")
        draw_tubules(
            session=session,
            length=length, 
            interval=MAX_INTERVAL,
            centerline_points=structure['centerline'],
            angle=0,
            radii=[cp_radius, cp_radius],
            shift_distances=[cp_shift, -cp_shift],
            length_diffs=[0.0, 0.0], 
            tubule_names=["C1", "C2"],
            colors=[(100, 255, 255, 255), (100, 255, 255, 255)],
            group_name="central_pair"
        )
    
    # Draw each doublet microtubule
    session.logger.info(f"Drawing {num_doublets} doublet microtubules...")
    for doublet_info in structure['doublets']:
        # Get the shifted centerline for this doublet
        doublet_centerline = get_doublet_centerline(
            structure['centerline'],
            doublet_info['angle'],
            doublet_info['shift_distance']
        )
        
        # Draw the doublet (A and B tubules)
        draw_tubules(
            session=session,
            length=length,
            interval=MAX_INTERVAL,
            centerline_points=doublet_centerline,
            angle=doublet_info['angle'] + CILIA_OFFSET_ANGLE,  # Orientation of A-B pair
            radii=[doublet_a_radius, doublet_b_radius],
            shift_distances=[-doublet_shift, doublet_shift],
            length_diffs=[0.0, -doublet_length_diff],
            tubule_names=[f"MT{doublet_info['index']+1}_A", 
                         f"MT{doublet_info['index']+1}_B"],
            colors=[(100, 100, 255, 255), (100, 100, 255, 255)],
            group_name=doublet_info['name']
        )
    
    session.logger.info(f"Cilia model generated successfully!")
    session.logger.info(f"  Type: {centerline_type}")
    session.logger.info(f"  Length: {length} Å")
    session.logger.info(f"  Doublets: {num_doublets}")
    session.logger.info(f"  Cilia radius: {cilia_radius} Å")


def centriolesim(session,
                length=5000,
                line='straight',
                curve_radius=10000.0,
                sine_frequency=2.0, sine_amplitude=500.0,
                # Centriole Structure Defaults
                num_triplets=9,
                centriole_radius=1100.0,
                centriole_angle_offset=60.0,  # Tunable angle offset for triplets
                # Triplet Geometry Defaults
                triplet_a_radius=125.0,
                triplet_b_radius=135.0,
                triplet_c_radius=140.0,
                triplet_ab_shift=70.0,  # A and B shift from triplet centerline
                triplet_c_shift=200.0,   # C shift from triplet centerline
                triplet_b_length_diff=5.0,  # B shorter than A
                triplet_c_length_diff=300.0   # C shorter than A
                ):
    """
    Generate and draw a complete centriole structure with triplet microtubules.
    
    Parameters:
    -----------
    session : chimerax session
        ChimeraX session
    length : float
        Length of the centriole in Angstroms (default: 5000)
    line : str
        Type of centerline: 'straight', 'curve', or 'sinusoidal' (default: 'straight')
    curve_radius : float
        Radius of curvature for 'curve' type (default: 10000.0)
    sine_frequency : float
        Frequency of sinusoidal oscillation (default: 2.0)
    sine_amplitude : float
        Amplitude of sinusoidal oscillation (default: 500.0)
        
    num_triplets : int
        Number of triplet microtubules (default: 9)
    centriole_radius : float
        Radius from center to triplets in Angstroms (default: 1100.0)
    centriole_angle_offset : float
        Angle offset for triplet orientation in degrees (default: 0.0)
        This controls the A-B-C orientation relative to the radial direction
        
    triplet_a_radius : float
        Radius of the A-tubule (default: 125.0)
    triplet_b_radius : float
        Radius of the B-tubule (default: 145.0)
    triplet_c_radius : float
        Radius of the C-tubule (default: 145.0)
    triplet_ab_shift : float
        Radial distance of A and B tubules from triplet centerline (default: 70.0)
    triplet_c_shift : float
        Radial distance of C tubule from triplet centerline (default: 140.0)
    triplet_b_length_diff : float
        Length difference: B shorter than A (default: 250.0)
    triplet_c_length_diff : float
        Length difference: C shorter than A (default: 500.0)
    """
    
    centerline_type = line
    
    session.logger.info(f"Generating centriole structure with {centerline_type} centerline...")
    
    # Calculate centriole structure points using the same function
    structure = generate_cilia_structure(
        length=length,
        centerline_type=centerline_type,
        curve_radius=curve_radius,
        sine_frequency=sine_frequency,
        sine_amplitude=sine_amplitude,
        num_doublets=num_triplets,  # Reuse parameter name
        cilia_radius=centriole_radius
    )
    
    # Draw each triplet microtubule
    session.logger.info(f"Drawing {num_triplets} triplet microtubules...")
    for triplet_info in structure['doublets']:  # Reuse 'doublets' key
        # Get the shifted centerline for this triplet
        triplet_centerline = get_doublet_centerline(
            structure['centerline'],
            triplet_info['angle'],
            triplet_info['shift_distance']
        )
        
        # Draw the triplet (A, B, and C tubules)
        draw_tubules(
            session=session,
            length=length,
            interval=MAX_INTERVAL,
            centerline_points=triplet_centerline,
            angle=triplet_info['angle'] + centriole_angle_offset,  # Tunable offset
            radii=[triplet_a_radius, triplet_b_radius, triplet_c_radius],
            shift_distances=[-triplet_ab_shift, triplet_ab_shift, triplet_c_shift],
            length_diffs=[0.0, -triplet_b_length_diff, -triplet_c_length_diff],
            tubule_names=[f"MT{triplet_info['index']+1}_A",
                         f"MT{triplet_info['index']+1}_B",
                         f"MT{triplet_info['index']+1}_C"],
            colors=[(100, 255, 255, 255), (100, 255, 255, 255), (100, 100, 255, 255)],
            group_name=f"triplet_{triplet_info['index']+1}"
        )
    
    session.logger.info(f"Centriole model generated successfully!")
    session.logger.info(f"  Type: {centerline_type}")
    session.logger.info(f"  Length: {length} Å")
    session.logger.info(f"  Triplets: {num_triplets}")
    session.logger.info(f"  Centriole radius: {centriole_radius} Å")
    session.logger.info(f"  Angle offset: {centriole_angle_offset}°")


# Command description for ciliasim
ciliasim_desc = CmdDesc(
    keyword=[
        ('length', FloatArg),
        ('line', StringArg),
        ('curve_radius', FloatArg),
        ('sine_frequency', FloatArg),
        ('sine_amplitude', FloatArg),
        ('num_doublets', IntArg),
        ('cilia_radius', FloatArg),
        ('draw_central_pair', BoolArg),
        ('doublet_a_radius', FloatArg),
        ('doublet_b_radius', FloatArg),
        ('doublet_shift', FloatArg),
        ('doublet_length_diff', FloatArg),
        ('cp_radius', FloatArg),
        ('cp_shift', FloatArg)
    ],
    synopsis='Generate complete cilia structure with customizable geometry'
)

# Command description for centriolesim
centriolesim_desc = CmdDesc(
    keyword=[
        ('length', FloatArg),
        ('line', StringArg),
        ('curve_radius', FloatArg),
        ('sine_frequency', FloatArg),
        ('sine_amplitude', FloatArg),
        ('num_triplets', IntArg),
        ('centriole_radius', FloatArg),
        ('centriole_angle_offset', FloatArg),
        ('triplet_a_radius', FloatArg),
        ('triplet_b_radius', FloatArg),
        ('triplet_c_radius', FloatArg),
        ('triplet_ab_shift', FloatArg),
        ('triplet_c_shift', FloatArg),
        ('triplet_b_length_diff', FloatArg),
        ('triplet_c_length_diff', FloatArg)
    ],
    synopsis='Generate complete centriole structure with triplet microtubules'
)