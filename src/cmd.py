# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc, IntArg, FloatArg, StringArg, BoolArg, Color8Arg

import pandas as pd
import numpy as np
from chimerax.core.models import Surface
from .draw import draw_tubules, draw_membrane
from .geometry.centerline import generate_cilia_structure, get_doublet_centerline
from .io import generate_cilia_with_tip_csv

# Experimental
import csv

# Import default value from default_config.py
from . import default_config

def ciliabuild(session, 
            length=default_config.CILIA_LENGTH, 
            line=default_config.CILIA_LINE,
            curve_radius=default_config.CILIA_CURVE_RADIUS, 
            sine_frequency=default_config.CILIA_SINE_FREQUENCY, 
            sine_amplitude=default_config.CILIA_SINE_AMPLITUDE,
            template_file=default_config.TEMPLATE_FILE,
            tip_length=default_config.TIP_LENGTH,
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
            cp_color=default_config.CILIA_CP_COLOR,
            membrane_color=default_config.CILIA_MEMBRANE_COLOR,
            write_csv=False
            ):
    """
    Generate and draw a complete cilia structure with doublet microtubules.
    """
    
    centerline_type = line 
    
    if centerline_type == 'tip':
        session.logger.info(f"Generating cilia structure with {centerline_type} geometry...")
        
        # Generate tip CSV data
        cilia_data_df = generate_cilia_with_tip_csv(
            cilia_length=length,
            cilia_radius=cilia_radius,
            tip_length=tip_length,
            num_doublets=num_doublets,
            doublet_shift=doublet_shift,
            cp_doublet_length_diff=cp_doublet_length_diff,
            cp_shift=cp_shift,
            membrane_radius=membrane_radius,
            membrane_fraction=membrane_fraction,
            max_interval=default_config.MAX_INTERVAL
        )
        # Only turn on writing by default for this
        write_csv = True
        
    else:
        # Generate structure for 'straight', 'curve', 'sinusoidal', 'template'
        session.logger.info(f"Generating cilia structure with {centerline_type} centerline...")
        
        structure = generate_cilia_structure(
            length=length,
            centerline_type=centerline_type,
            curve_radius=curve_radius,
            sine_frequency=sine_frequency,
            sine_amplitude=sine_amplitude,
            template_file=template_file,
            num_doublets=num_doublets,
            cilia_radius=cilia_radius
        )
        
        # Build DataFrame from structure
        session.logger.info(f"Building DataFrame for {num_doublets} doublets...")
        
        all_data = []
        idx_b_start = 1
        
        # Process each doublet
        for doublet_info in structure['doublets']:
            # Get the shifted centerline for this doublet
            doublet_centerline = get_doublet_centerline(
                structure['centerline'],
                doublet_info['angle'],
                doublet_info['shift_distance']
            )
            
            if len(doublet_centerline) != len(structure['centerline']):
                doublet_centerline = doublet_centerline[:len(structure['centerline'])]
            
            # Calculate arc length
            doublet_segment_lengths = np.linalg.norm(np.diff(doublet_centerline, axis=0), axis=1)
            doublet_cumulative_length = np.concatenate(([0], np.cumsum(doublet_segment_lengths)))
            
            # Calculate A-tubule end index
            total_doublet_length = doublet_cumulative_length[-1]
            if cp_doublet_length_diff > 0:
                target_a_length = total_doublet_length - cp_doublet_length_diff
                idx_a = np.searchsorted(doublet_cumulative_length, target_a_length, side='right')
                idx_a = max(2, min(idx_a, len(doublet_centerline)))
            else:
                idx_a = len(doublet_centerline)
            
            # Calculate B-tubule end index
            doublet_centerline_a = doublet_centerline[:idx_a]
            a_segment_lengths = np.linalg.norm(np.diff(doublet_centerline_a, axis=0), axis=1)
            a_cumulative_length = np.concatenate(([0], np.cumsum(a_segment_lengths)))
            total_a_length = a_cumulative_length[-1]
            target_b_length = total_a_length - doublet_length_diff
            idx_b = np.searchsorted(a_cumulative_length, target_b_length, side='right')
            idx_b = max(2, min(idx_b, len(doublet_centerline_a)))
            
            # Add data for this doublet
            doublet_number = doublet_info['index'] + 1
            angle = doublet_info['angle'] + default_config.CILIA_OFFSET_ANGLE
            
            for i, point in enumerate(doublet_centerline_a):
                x, y, z = point
                idx_a_val = 1
                idx_b_val = 1 if (idx_b_start <= i < idx_b) else 0
                
                all_data.append([
                    doublet_number, x, y, z, 
                    idx_a_val, idx_b_val, angle, 
                    -doublet_shift, doublet_shift
                ])
        
        # Add central pair data (DoubletNumber = -1)
        if draw_central_pair:
            for point in structure['centerline']:
                x, y, z = point
                all_data.append([
                    -1, x, y, z, 
                    1, 1, 0, 
                    cp_shift, -cp_shift
                ])
        
        # Add membrane data (DoubletNumber = 0)
        if membrane:
            membrane_fraction = max(0.0, min(1.0, membrane_fraction))
            total_points = len(structure['centerline'])
            membrane_points = max(2, int(total_points * membrane_fraction))
            membrane_path = structure['centerline'][:membrane_points]
            
            for point in membrane_path:
                x, y, z = point
                all_data.append([
                    0, x, y, z, 
                    1, 1, 0, 
                    0, 0
                ])
        
        # Create DataFrame
        columns = ['DoubletNumber', 'X', 'Y', 'Z', 'Idx_A', 'Idx_B', 'Angle', 'A_Shift', 'B_Shift']
        cilia_data_df = pd.DataFrame(all_data, columns=columns)
        
    
    # Draw the structure using _ciliabuild_from_df
    cilia_root = _ciliabuild_from_df(
        session=session,
        df=cilia_data_df,
        draw_central_pair=draw_central_pair,
        membrane=membrane,
        membrane_fraction=membrane_fraction,
        membrane_radius=membrane_radius,
        doublet_a_radius=doublet_a_radius,
        doublet_b_radius=doublet_b_radius,
        doublet_shift=doublet_shift,
        cp_radius=cp_radius,
        cp_shift=cp_shift,
        doublet_a_color=doublet_a_color,
        doublet_b_color=doublet_b_color,
        cp_color=cp_color,
        membrane_color=membrane_color
    )
    
    # Write CSV if requested
    if write_csv:
        csv_filename = 'cilia.csv'
        cilia_data_df['X'] = cilia_data_df['X'].map('{:.2f}'.format)
        cilia_data_df['Y'] = cilia_data_df['Y'].map('{:.2f}'.format)
        cilia_data_df['Z'] = cilia_data_df['Z'].map('{:.2f}'.format)
        cilia_data_df.to_csv(csv_filename, index=False)
        session.logger.info(f"Saved structure to {csv_filename}")
        
    # Update root name based on type
    if cilia_root:
        model_id = cilia_root.id_string
        if centerline_type == 'tip':
            cilia_root.name = f"Cilia_Tip {model_id}"
        else:
            cilia_root.name = f"Cilia_{centerline_type} {model_id}"
        
        session.logger.info(f"Cilia model generated successfully!")
        session.logger.info(f"  Type: {centerline_type}")
        session.logger.info(f"  Length: {length if centerline_type != 'tip' else tip_length} Å")
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
                template_file=default_config.TEMPLATE_FILE,
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
        template_file=default_config.TEMPLATE_FILE,
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

def _ciliabuild_from_df(session, df,
            # Cilia Structure Defaults
            draw_central_pair=default_config.CILIA_DRAW_CENTRAL_PAIR,
            membrane=default_config.CILIA_MEMBRANE,
            membrane_fraction=default_config.CILIA_MEMBRANE_FRACTION,
            membrane_radius=default_config.CILIA_MEMBRANE_RADIUS,
            # Doublet Geometry Defaults
            doublet_a_radius=default_config.CILIA_DOUBLET_A_RADIUS,
            doublet_b_radius=default_config.CILIA_DOUBLET_B_RADIUS,
            doublet_shift=default_config.CILIA_DOUBLET_SHIFT,
            # Central Pair Geometry Defaults
            cp_radius=default_config.CILIA_CP_RADIUS,
            cp_shift=default_config.CILIA_CP_SHIFT,
            # Color parameters
            doublet_a_color=default_config.CILIA_DOUBLET_A_COLOR,
            doublet_b_color=default_config.CILIA_DOUBLET_B_COLOR,
            cp_color=default_config.CILIA_CP_COLOR,
            membrane_color=default_config.CILIA_MEMBRANE_COLOR
            ):
    """
    Internal function to generate and draw cilia structure from a DataFrame.
    
    Parameters:
    -----------
    session : chimerax session
        ChimeraX session
    df : pandas.DataFrame
        DataFrame with columns: DoubletNumber,X,Y,Z,Idx_A,Idx_B,Angle,A_shift,B_shift
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
    membrane_color : tuple
        RGBA color for membrane (default: (105, 105, 105, 255))
        
    Returns:
    --------
    cilia_root : Surface
        The root surface model containing all cilia components
    """
    
    # Validate DataFrame columns
    required_columns = ['DoubletNumber', 'X', 'Y', 'Z', 'Idx_A', 'Idx_B', 'Angle', 'A_Shift', 'B_Shift']
    if not all(col in df.columns for col in required_columns):
        session.logger.error(f"DataFrame must contain columns: {required_columns}")
        return None
    
    cilia_root = Surface("Cilia_from_template", session)
    session.models.add([cilia_root])
    
    # Get unique doublet numbers
    doublet_numbers = sorted(df.loc[df['DoubletNumber'] > 0, 'DoubletNumber'].unique())
    num_doublets = len(doublet_numbers)
    
    session.logger.info(f"Found {num_doublets} doublets in data")
    
    # Process each doublet
    for doublet_num in doublet_numbers:
        doublet_data = df[df['DoubletNumber'] == doublet_num]
        
        # Extract points where Idx_A = 1 (A-tubule)
        a_data = doublet_data[doublet_data['Idx_A'] == 1]
        doublet_centerline_a = a_data[['X', 'Y', 'Z']].values
        
        # Extract points where Idx_B = 1 (B-tubule)
        b_data = doublet_data[doublet_data['Idx_B'] == 1]
        doublet_centerline_b = b_data[['X', 'Y', 'Z']].values
        
        angle = doublet_data['Angle'].iloc[0]
        a_shift = doublet_data['A_Shift'].iloc[0]
        b_shift = doublet_data['B_Shift'].iloc[0]
        
        if len(doublet_centerline_a) < 2:
            session.logger.warning(f"Doublet {doublet_num} has insufficient A-tubule points, skipping")
            continue
        
        # Create empty doublet surfs
        doublet_surfs = []
        
        # Draw A-tubule (full length)
        a_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=doublet_centerline_a,
            angle=angle,
            radii=[doublet_a_radius],
            shift_distances=[a_shift],
            length_diffs=None,
            tubule_names=[f"A_tubule"],
            colors=[doublet_a_color],
            group_name=f"A_tubule",
            add_to_session=False
        )
        if a_surfs:
            doublet_surfs.extend(a_surfs)
        
        # Draw B-tubule if it has points
        if len(doublet_centerline_b) >= 2:
            b_surfs = draw_tubules(
                session=session,
                length=None,
                interval=default_config.MAX_INTERVAL,
                centerline_points=doublet_centerline_b,
                angle=angle,
                radii=[doublet_b_radius],
                shift_distances=[b_shift],
                length_diffs=None,
                tubule_names=[f"B_tubule"],
                colors=[doublet_b_color],
                group_name=f"B_tubule",
                add_to_session=False
            )
            if b_surfs:
                doublet_surfs.extend(b_surfs)
        
        session.models.add_group(doublet_surfs, parent=cilia_root, name=f"DMT{int(doublet_num)}")
        session.logger.info(f"Added doublet {int(doublet_num)}")
        
        # Draw central pair if requested (using the DoubletNumber = -1)
    cp_data = df[df['DoubletNumber'] == -1]
    if draw_central_pair and len(cp_data) > 0:
        session.logger.info("Drawing central pair (C1, C2)...")
        cp_centerline = cp_data[['X', 'Y', 'Z']].values
        cp_shift1 = cp_data['A_Shift'].iloc[0]
        cp_shift2 = cp_data['B_Shift'].iloc[0]

        cp_surfs = draw_tubules(
            session=session,
            length=None,  # Use full centerline
            interval=default_config.MAX_INTERVAL,
            centerline_points=cp_centerline,
            angle=0,
            radii=[cp_radius, cp_radius],
            shift_distances=[cp_shift1, cp_shift2],
            length_diffs=None,
            tubule_names=["C1", "C2"],
            colors=[cp_color, cp_color],
            group_name="central_pair"
        )
        session.models.add_group(cp_surfs, parent=cilia_root, name="Central Pair")
        session.logger.info(f"Added central pair")
            
    # Draw membrane if requested (using DoubletNumber = 0) 
    membrane_data = df[df['DoubletNumber'] == 0]
    if draw_central_pair and len(membrane_data) > 0:
        session.logger.info("Drawing membrane ...")
        membrane_centerline = membrane_data[['X', 'Y', 'Z']].values
            
        membrane_surfs = draw_membrane(
            session=session,
            centerline_points=membrane_centerline,
            radius=membrane_radius,
            segments=32, 
            color=membrane_color,
            name="Membrane"
        )
        if membrane_surfs:
            session.models.add_group(membrane_surfs, parent=cilia_root, name="Membrane")
            session.logger.info(f"Added membrane")

    
    # Get the model ID and update the name
    model_id = cilia_root.id_string
    cilia_root.name = f"Cilia_from_template {model_id}"
    
    session.logger.info(f"Cilia model from template generated successfully!")
    session.logger.info(f"  Doublets: {num_doublets}")
    if membrane:
        session.logger.info(f"  Membrane: {membrane_fraction*100:.1f}% coverage, radius {membrane_radius} Å")
    
    return cilia_root


def ciliabuild_from_csv(session,
            template_csv='cilia.csv',
            # Cilia Structure Defaults
            draw_central_pair=default_config.CILIA_DRAW_CENTRAL_PAIR,
            membrane=default_config.CILIA_MEMBRANE,
            membrane_fraction=default_config.CILIA_MEMBRANE_FRACTION,
            membrane_radius=default_config.CILIA_MEMBRANE_RADIUS,
            # Doublet Geometry Defaults
            doublet_a_radius=default_config.CILIA_DOUBLET_A_RADIUS,
            doublet_b_radius=default_config.CILIA_DOUBLET_B_RADIUS,
            doublet_shift=default_config.CILIA_DOUBLET_SHIFT,
            # Central Pair Geometry Defaults
            cp_radius=default_config.CILIA_CP_RADIUS,
            cp_shift=default_config.CILIA_CP_SHIFT,
            # Color parameters
            doublet_a_color=default_config.CILIA_DOUBLET_A_COLOR,
            doublet_b_color=default_config.CILIA_DOUBLET_B_COLOR,
            cp_color=default_config.CILIA_CP_COLOR,
            membrane_color=default_config.CILIA_MEMBRANE_COLOR
            ):
    """
    Generate and draw a complete cilia structure from a template CSV file.
    
    Parameters:
    -----------
    session : chimerax session
        ChimeraX session
    template_csv : str
        Path to CSV file with columns: DoubletNumber,X,Y,Z,Idx_A,Idx_B,Angle,A_shift,B_shift
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
        
    Returns:
    --------
    cilia_root : Surface or None
        The root surface model containing all cilia components, or None if loading failed
    """
    
    session.logger.info(f"Loading cilia structure from template: {template_csv}")
    
    # Read the CSV file
    try:
        df = pd.read_csv(template_csv)
    except FileNotFoundError:
        session.logger.error(f"Template file not found: {template_csv}")
        return None
    except Exception as e:
        session.logger.error(f"Error reading template file: {e}")
        return None
    
    # Call the internal function with the loaded DataFrame
    return _ciliabuild_from_df(
        session=session,
        df=df,
        draw_central_pair=draw_central_pair,
        membrane=membrane,
        membrane_fraction=membrane_fraction,
        membrane_radius=membrane_radius,
        doublet_a_radius=doublet_a_radius,
        doublet_b_radius=doublet_b_radius,
        doublet_shift=doublet_shift,
        cp_radius=cp_radius,
        cp_shift=cp_shift,
        doublet_a_color=doublet_a_color,
        doublet_b_color=doublet_b_color,
        cp_color=cp_color,
        membrane_color=membrane_color
    )

# Command description for ciliabuild
ciliabuild_desc = CmdDesc(
    keyword=[
        ('length', FloatArg),
        ('line', StringArg),
        ('curve_radius', FloatArg),
        ('sine_frequency', FloatArg),
        ('sine_amplitude', FloatArg),
        ('template_file', StringArg),
        ('tip_length', FloatArg),
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
        ('cp_color', Color8Arg),
        ('membrane_color', Color8Arg),
        ('write_csv', BoolArg)
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
        ('template_file', StringArg),
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

# Command description for ciliabuild_from_template
ciliabuild_from_csv_desc = CmdDesc(
    keyword=[
        ('data_source', StringArg),
        ('draw_central_pair', BoolArg),
        ('membrane', BoolArg),
        ('membrane_fraction', FloatArg),
        ('membrane_radius', FloatArg),
        ('doublet_a_radius', FloatArg),
        ('doublet_b_radius', FloatArg),
        ('doublet_shift', FloatArg),
        ('cp_radius', FloatArg),
        ('cp_shift', FloatArg),
        ('doublet_a_color', Color8Arg),
        ('doublet_b_color', Color8Arg),
        ('cp_color', Color8Arg),
        ('membrane_color', Color8Arg)
    ],
    synopsis='Generate a cilia structure from a template CSV file or in-memory data'
)

