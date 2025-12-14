# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc, IntArg, FloatArg, StringArg, BoolArg

from .curve import generate_cilia_structure, get_doublet_centerline
from .draw import draw_tubules

# Define the optimal sampling interval for smoothness (must match curve.py)
MAX_INTERVAL = 10.0 # 10 Angstroms interval for centerline points

def ciliasim(session, 
            length=5000, 
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
    Generate and draw a complete cilia structure.
    
    Parameters:
    -----------
    session : chimerax session
        ChimeraX session
    length : float
        Length of the cilia in Angstroms (default: 5000). 
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
        Radius from center to doublets in Angstroms (default: 190.0)
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
            interval=MAX_INTERVAL, # Use defined interval for truncation
            centerline_points=structure['centerline'],
            angle=0,
            radii=[cp_radius, cp_radius], # Use argument
            shift_distances=[cp_shift, -cp_shift], # Use argument
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
            length=length, # Pass the original nominal length
            interval=MAX_INTERVAL, # Pass the sampling interval
            centerline_points=doublet_centerline,
            angle=doublet_info['angle'] + 90,  # Orientation of A-B pair
            radii=[doublet_a_radius, doublet_b_radius], # Use arguments
            shift_distances=[-doublet_shift, doublet_shift], # Use argument, B is usually shifted "outward" 
            length_diffs=[0.0, -doublet_length_diff], # Use argument
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


# Command description - Register ALL keywords with their default types
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
        # Doublet Geometry Arguments
        ('doublet_a_radius', FloatArg),
        ('doublet_b_radius', FloatArg),
        ('doublet_shift', FloatArg),
        ('doublet_length_diff', FloatArg),
        # Central Pair Geometry Arguments
        ('cp_radius', FloatArg),
        ('cp_shift', FloatArg)
    ],
    synopsis='Generate complete cilia structure with customizable geometry'
)