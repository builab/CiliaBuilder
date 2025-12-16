# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc, IntArg, FloatArg, StringArg, BoolArg, Color8Arg

import numpy as np
from chimerax.core.models import Surface
from .curve import generate_cilia_structure, get_doublet_centerline
from .draw import draw_tubules, draw_membrane

# Import default value from default_config.py
from . import default_config




def ciliabuild(session, 
            length=default_config.CILIA_LENGTH, 
            line=default_config.CILIA_LINE,
            curve_radius=default_config.CILIA_CURVE_RADIUS, 
            sine_frequency=default_config.CILIA_SINE_FREQUENCY, 
            sine_amplitude=default_config.CILIA_SINE_AMPLITUDE,
            # Cilia Structure Defaults
            num_doublets=default_config.CILIA_NUM_DOUBLETS, 
            cilia_radius=default_config.CILIA_RADIUS,
            draw_central_pair=default_config.CILIA_DRAW_CENTRAL_PAIR,
            membrane=default_config.CILIA_MEMBRANE,
            membrane_fraction=default_config.CILIA_MEMBRANE_FRACTION,
            membrane_radius=default_config.CILIA_MEMBRANE_RADIUS,
            # Doublet Geometry Defaults
            doublet_a_radius=default_config.CILIA_DOUBLET_A_RADIUS,
            doublet_b_radius=default_config.CILIA_DOUBLET_B_RADIUS,
            doublet_shift=default_config.CILIA_DOUBLET_SHIFT,
            doublet_length_diff=default_config.CILIA_DOUBLET_LENGTH_DIFF,
            cp_doublet_length_diff=default_config.CILIA_CP_DOUBLET_LENGTH_DIFF,
            # Central Pair Geometry Defaults
            cp_radius=default_config.CILIA_CP_RADIUS,
            cp_shift=default_config.CILIA_CP_SHIFT,
            # Color parameters
            doublet_a_color=default_config.CILIA_DOUBLET_A_COLOR,
            doublet_b_color=default_config.CILIA_DOUBLET_B_COLOR,
            cp_color=default_config.CILIA_CP_COLOR
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
    membrane : bool
        Whether to draw the membrane (default: True)
    membrane_fraction : float
        Fraction of cilia length covered by membrane, from base (0.0-1.0) (default: 0.7)
    membrane_radius : float
        Radius of the membrane in Angstroms (default: 1000)
        
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
    
    doublet_a_color : tuple
        RGBA color for A-tubules (default: (100, 100, 255, 255))
    doublet_b_color : tuple
        RGBA color for B-tubules (default: (100, 100, 255, 255))
    cp_color : tuple
        RGBA color for central pair tubules (default: (255, 255, 100, 255))
        
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
    
    cilia_root = Surface("Cilia", session)
    session.models.add([cilia_root])
            
    # Draw each doublet microtubule
    session.logger.info(f"Drawing {num_doublets} doublet microtubules...")
    
    # For B-tubule: 
    # - Find START index (offset from beginning for visual separation)
    idx_b_start = 1
    
    for doublet_info in structure['doublets']:
        # Get the shifted centerline for this doublet
        doublet_centerline = get_doublet_centerline(
            structure['centerline'],
            doublet_info['angle'],
            doublet_info['shift_distance']
        )
        
        # Ensure doublet centerline has same number of points as main centerline
        if len(doublet_centerline) != len(structure['centerline']):
            doublet_centerline = doublet_centerline[:len(structure['centerline'])]
        
        # Calculate arc length along this doublet's centerline
        doublet_segment_lengths = np.linalg.norm(np.diff(doublet_centerline, axis=0), axis=1)
        doublet_cumulative_length = np.concatenate(([0], np.cumsum(doublet_segment_lengths)))
        
        # For A-tubule: shorten if CP should be longer
        total_doublet_length = doublet_cumulative_length[-1]
        if cp_doublet_length_diff > 0:
            target_a_length = total_doublet_length - cp_doublet_length_diff
            idx_a = np.searchsorted(doublet_cumulative_length, target_a_length, side='right')
            idx_a = max(2, min(idx_a, len(doublet_centerline)))
        else:
            idx_a = len(doublet_centerline)

        # For B-tubule: find index relative to A-tubule length
        doublet_centerline_a = doublet_centerline[:idx_a]
        a_segment_lengths = np.linalg.norm(np.diff(doublet_centerline_a, axis=0), axis=1)
        a_cumulative_length = np.concatenate(([0], np.cumsum(a_segment_lengths)))
        total_a_length = a_cumulative_length[-1]
        target_b_length = total_a_length - doublet_length_diff
        idx_b = np.searchsorted(a_cumulative_length, target_b_length, side='right')
        idx_b = max(2, min(idx_b, len(doublet_centerline_a)))

        # Create separate centerlines for A and B
        doublet_centerline_b = doublet_centerline_a[idx_b_start:idx_b]  # Shortened relative to A
        
        # Create empty doublet surfs
        doublet_surfs = []
        
        # Draw A-tubule (full length)
        a_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=doublet_centerline_a,
            angle=doublet_info['angle'] + default_config.CILIA_OFFSET_ANGLE,
            radii=[doublet_a_radius],
            shift_distances=[-doublet_shift],
            length_diffs=None,
            tubule_names=[f"A_tubule"],
            colors=[doublet_a_color],
            group_name=f"A_tubule",
            add_to_session=False  # Don't add yet
        )
        if a_surfs:
            doublet_surfs.extend(a_surfs)
        
        # Draw B-tubule (shortened)
        b_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=doublet_centerline_b,
            angle=doublet_info['angle'] + default_config.CILIA_OFFSET_ANGLE,
            radii=[doublet_b_radius],
            shift_distances=[doublet_shift],
            length_diffs=None,
            tubule_names=[f"B_tubule"],
            colors=[doublet_b_color],
            group_name=f"B_tubule",
            add_to_session=False  # Don't add yet
        )
        if b_surfs:
            doublet_surfs.extend(b_surfs)
        
        session.models.add_group(doublet_surfs, parent=cilia_root, name=f"DMT{doublet_info['index']+1}")
        session.logger.info(f"Added doublet {doublet_info['index']+1}")
       
    # Draw central pair if requested
    if draw_central_pair:
        session.logger.info("Drawing central pair (C1, C2)...")
        cp_surfs = draw_tubules(
            session=session,
            length=None,  # Use full centerline
            interval=default_config.MAX_INTERVAL,
            centerline_points=structure['centerline'],
            angle=0,
            radii=[cp_radius, cp_radius],
            shift_distances=[cp_shift, -cp_shift],
            length_diffs=None,
            tubule_names=["C1", "C2"],
            colors=[cp_color, cp_color],
            group_name="central_pair"
        )
        session.models.add_group(cp_surfs, parent=cilia_root, name="Central Pair")
        session.logger.info(f"Added central pair")

    # Draw membrane if requested
    if membrane:
        # Calculate membrane path - use only the fraction from the base
        
        # Clamp membrane_fraction to valid range
        membrane_fraction = max(0.0, min(1.0, membrane_fraction))
        
        total_points = len(structure['centerline'])
        # Calculate number of points for the membrane (at least 2 for valid geometry)
        # Use floor division to ensure we don't exceed the fraction
        membrane_points = max(2, int(total_points * membrane_fraction))
        
        # Extract the membrane path (from base, index 0)
        membrane_path = structure['centerline'][:membrane_points]
        
        # Calculate actual membrane length from the actual path
        if len(membrane_path) >= 2:
            # Sum up distances between consecutive points
            dists = np.linalg.norm(np.diff(membrane_path, axis=0), axis=1)
            actual_membrane_length = np.sum(dists)
        else:
            actual_membrane_length = 0.0
        
        session.logger.info(f"Drawing membrane:")
        session.logger.info(f"  Cilia length: {length} Å")
        session.logger.info(f"  Membrane fraction: {membrane_fraction*100:.1f}%")
        session.logger.info(f"  Membrane length: {actual_membrane_length:.0f} Å ({actual_membrane_length/length*100:.1f}% of cilia)")
        session.logger.info(f"  Points used: {membrane_points}/{total_points}")
        session.logger.info(f"  First point Z: {membrane_path[0][2]:.1f}, Last point Z: {membrane_path[-1][2]:.1f}")
        
        membrane_surfs = draw_membrane(
            session=session,
            centerline_points=membrane_path,
            radius=membrane_radius,
            segments=32, 
            color=(105, 105, 105, 255),
            name="Membrane"
        )
        if membrane_surfs:
            session.models.add_group(membrane_surfs, parent=cilia_root, name="Membrane")
            session.logger.info(f"Added membrane")
    
    # Get the model ID and update the name
    model_id = cilia_root.id_string
    cilia_root.name = f"Cilia {model_id}"
    
    session.logger.info(f"Cilia model generated successfully!")
    session.logger.info(f"  Type: {centerline_type}")
    session.logger.info(f"  Length: {length} Å")
    session.logger.info(f"  Doublets: {num_doublets}")
    session.logger.info(f"  Cilia radius: {cilia_radius} Å")
    if membrane:
        session.logger.info(f"  Membrane: {membrane_fraction*100:.1f}% coverage, radius {membrane_radius} Å")
    
    return cilia_root


def centriolebuild(session,
                length=default_config.CENTRIOLE_LENGTH,
                line=default_config.CENTRIOLE_LINE,
                curve_radius=default_config.CENTRIOLE_CURVE_RADIUS,
                sine_frequency=default_config.CENTRIOLE_SINE_FREQUENCY,
                sine_amplitude=default_config.CENTRIOLE_SINE_AMPLITUDE,
                # Centriole Structure Defaults
                num_triplets=default_config.CENTRIOLE_NUM_TRIPLETS,
                centriole_radius=default_config.CENTRIOLE_RADIUS,
                centriole_angle_offset=default_config.CILIA_OFFSET_ANGLE - default_config.CENTRIOLE_OFFSET_ANGLE,
                # Triplet Geometry Defaults
                triplet_a_radius=default_config.CENTRIOLE_TRIPLET_A_RADIUS,
                triplet_b_radius=default_config.CENTRIOLE_TRIPLET_B_RADIUS,
                triplet_c_radius=default_config.CENTRIOLE_TRIPLET_C_RADIUS,
                triplet_ab_shift=default_config.CENTRIOLE_TRIPLET_AB_SHIFT,
                triplet_c_shift=default_config.CENTRIOLE_TRIPLET_C_SHIFT,
                triplet_b_length_diff=default_config.CENTRIOLE_TRIPLET_B_LENGTH_DIFF,
                triplet_c_length_diff=default_config.CENTRIOLE_TRIPLET_C_LENGTH_DIFF,
                # NEW PARAMETER FOR ALIGNMENT
                z_offset_end=default_config.CENTRIOLE_Z_OFFSET_END,
                # Color parameters
                triplet_a_color=default_config.CENTRIOLE_TRIPLET_A_COLOR,
                triplet_b_color=default_config.CENTRIOLE_TRIPLET_B_COLOR,
                triplet_c_color=default_config.CENTRIOLE_TRIPLET_C_COLOR
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
        Angle offset for triplet orientation in degrees (default: 60.0)
        This controls the A-B-C orientation relative to the radial direction
        
    triplet_a_radius : float
        Radius of the A-tubule (default: 125.0)
    triplet_b_radius : float
        Radius of the B-tubule (default: 135.0)
    triplet_c_radius : float
        Radius of the C-tubule (default: 135.0)
    triplet_ab_shift : float
        Radial distance of A and B tubules from triplet centerline (default: 140.0)
    triplet_c_shift : float
        Radial distance of C tubule from triplet centerline (default: -160.0)
    triplet_b_length_diff : float
        Length difference: B shorter than A at the END (default: 0.1)
    triplet_c_length_diff : float
        Length difference: C shorter than A at the END (default: 0.2)
    z_offset_end : float
        The target Z-coordinate for the end (tip) of the centriole structure (default: 20.0)
    triplet_a_color : tuple
        RGBA color for A-tubules (default: (100, 100, 255, 255))
    triplet_b_color : tuple
        RGBA color for B-tubules (default: (100, 100, 255, 255))
    triplet_c_color : tuple
        RGBA color for C-tubules (default: (179, 179, 255, 255))
    """
    
    centerline_type = line
    centriole_root = Surface("Centriole", session)
    session.models.add([centriole_root])
    
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
    
    # --- Centriole Alignment Shift to match Z_offset_end ---
    # Centriole is generated from Z=0 to Z=length. Calculate the Z-shift needed to place 
    # the end of the Centriole (originally at Z=length) at the desired z_offset_end.
    z_shift = z_offset_end - length
    
    # Apply the shift to all Z-coordinates in the centerline.
    # We assume the centerline is generated roughly along the Z-axis (index 2).
    structure['centerline'][:, 2] += z_shift
    
    session.logger.info(f"Centriole length {length} Å shifted by {z_shift:.1f} Å to end at Z={z_offset_end:.1f} Å.")
    session.logger.info(f"Centriole now spans approximately Z={0.0 + z_shift:.1f} to Z={length + z_shift:.1f}.")
    # --- END Shift ---
    
    # Draw each triplet microtubule
    session.logger.info(f"Drawing {num_triplets} triplet microtubules...")
    for triplet_info in structure['doublets']:  # Reuse 'doublets' key
        # Get the shifted centerline for this triplet
        triplet_centerline = get_doublet_centerline(
            structure['centerline'],
            triplet_info['angle'],
            triplet_info['shift_distance']
        )
        
        # Ensure triplet centerline has same number of points as main centerline
        if len(triplet_centerline) != len(structure['centerline']):
            triplet_centerline = triplet_centerline[:len(structure['centerline'])]
        
        # Calculate arc length along this triplet's centerline
        triplet_segment_lengths = np.linalg.norm(np.diff(triplet_centerline, axis=0), axis=1)
        triplet_cumulative_length = np.concatenate(([0], np.cumsum(triplet_segment_lengths)))
        
        total_triplet_length = triplet_cumulative_length[-1]
        
        # For B-tubule: 
        # - Find START index (offset from beginning for visual separation)
        idx_b_start = 1
        
        # - Find END index (shortened from end)
        target_b_length = total_triplet_length - triplet_b_length_diff
        idx_b_end = np.searchsorted(triplet_cumulative_length, target_b_length, side='right')
        idx_b_end = max(2, min(idx_b_end, len(triplet_centerline)))
        
        # For C-tubule:
        # - Find START index (offset from beginning for visual separation)
        idx_c_start = 2
        
        # - Find END index (shortened from end)
        target_c_length = total_triplet_length - triplet_c_length_diff
        idx_c_end = np.searchsorted(triplet_cumulative_length, target_c_length, side='right')
        idx_c_end = max(2, min(idx_c_end, len(triplet_centerline)))
        
        # Create separate centerlines for A, B, and C
        triplet_centerline_a = triplet_centerline  # Full length, starts at beginning
        triplet_centerline_b = triplet_centerline[idx_b_start:idx_b_end]  # Staggered start AND shortened end
        triplet_centerline_c = triplet_centerline[idx_c_start:idx_c_end]  # Staggered start AND shortened end
        
        # Create triplet surface
        triplet_surfs = []
        
        # Draw A-tubule (full length)
        a_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=triplet_centerline_a,
            angle=triplet_info['angle'] + centriole_angle_offset,
            radii=[triplet_a_radius],
            shift_distances=[-triplet_ab_shift],
            length_diffs=None,
            tubule_names=[f"A_tubule"],
            colors=[triplet_a_color],
            group_name=f"A_tubule",
            add_to_session=False
        )
        if a_surfs:
            triplet_surfs.extend(a_surfs)
        
        # Draw B-tubule (staggered start + shortened end)
        b_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=triplet_centerline_b,
            angle=triplet_info['angle'] + centriole_angle_offset,
            radii=[triplet_b_radius],
            shift_distances=[triplet_ab_shift],
            length_diffs=None,
            tubule_names=[f"B_tubule"],
            colors=[triplet_b_color],
            group_name=f"B_tubule",
            add_to_session=False
        )
        if b_surfs:
            triplet_surfs.extend(b_surfs)
        
        # Draw C-tubule (staggered start + shortened end)
        c_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=triplet_centerline_c,
            angle=triplet_info['angle'] + centriole_angle_offset,
            radii=[triplet_c_radius],
            shift_distances=[triplet_c_shift],
            length_diffs=None,
            tubule_names=[f"C_tubule"],
            colors=[triplet_c_color],
            group_name=f"C_tubule",
            add_to_session=False
        )
        if c_surfs:
            triplet_surfs.extend(c_surfs)
    
        # Add all surfaces as a single centriole group
        if triplet_surfs:
            session.models.add_group(triplet_surfs, parent=centriole_root, name=f"TMT{triplet_info['index']+1}")
            session.logger.info(f"Added triplet {triplet_info['index']+1} to 'Centriole' group")
    
    # Get the model ID and update the name
    model_id = centriole_root.id_string
    centriole_root.name = f"Centriole {model_id}"
    
    session.logger.info(f"Centriole model generated successfully!")
    session.logger.info(f"  Type: {centerline_type}")
    session.logger.info(f"  Length: {length} Å")
    session.logger.info(f"  Triplets: {num_triplets}")
    session.logger.info(f"  Centriole radius: {centriole_radius} Å")
    session.logger.info(f"  Angle offset: {centriole_angle_offset}°")
    
    return centriole_root


# Command description for ciliabuild
ciliabuild_desc = CmdDesc(
    keyword=[
        ('length', FloatArg),
        ('line', StringArg),
        ('curve_radius', FloatArg),
        ('sine_frequency', FloatArg),
        ('sine_amplitude', FloatArg),
        ('num_doublets', IntArg),
        ('cilia_radius', FloatArg),
        ('draw_central_pair', BoolArg),
        ('membrane', BoolArg),
        ('membrane_fraction', FloatArg),
        ('membrane_radius', FloatArg),
        ('doublet_a_radius', FloatArg),
        ('doublet_b_radius', FloatArg),
        ('doublet_shift', FloatArg),
        ('doublet_length_diff', FloatArg),
        ('cp_doublet_length_diff', FloatArg),
        ('cp_radius', FloatArg),
        ('cp_shift', FloatArg),
        ('doublet_a_color', Color8Arg),
        ('doublet_b_color', Color8Arg),
        ('cp_color', Color8Arg)
    ],
    synopsis='Generate complete cilia structure with customizable geometry'
)

centriolebuild_desc = CmdDesc(
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
        ('triplet_c_length_diff', FloatArg),
        ('z_offset_end', FloatArg), # Added new keyword argument
        ('triplet_a_color', Color8Arg),
        ('triplet_b_color', Color8Arg),
        ('triplet_c_color', Color8Arg)
    ],
    synopsis='Generate complete centriole structure with triplet microtubules'
)