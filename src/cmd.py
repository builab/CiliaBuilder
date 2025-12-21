# vim: set expandtab shiftwidth=4 softtabstop=4:

"""
ChimeraX command interface for CiliaBuilder.

This module provides high-level commands for generating and visualizing
cilia and centriole structures in ChimeraX. It handles parameter parsing,
structure generation, and 3D rendering.

Commands:
    ciliabuild - Generate cilia structures (9+2 configuration)
    centriolebuild - Generate centriole structures (9x3 configuration)
    ciliabuild_from_csv - Load and render structures from CSV templates
"""

from chimerax.core.commands import CmdDesc, IntArg, FloatArg, StringArg, BoolArg, Color8Arg

import pandas as pd
import numpy as np
from chimerax.core.models import Surface
from .draw import draw_tubules, draw_membrane, generate_capsule_surface
from .geometry.centerline import generate_cilia_structure, get_doublet_centerline
from .io import read_3d_csv, write_3d_csv, load_template_data
from .geometry.primarycilia import generate_primary_cilia
from .geometry.tip import generate_cilia_with_tip

# Import default configuration values
from . import default_config

# Core column names for doublet/triplet data (internal use)
REQUIRED_COLUMNS = [
    "DoubletNumber", "X", "Y", "Z", 
    "Idx_A", "Idx_B", "Angle", 
    "A_Shift", "B_Shift"
]

# Extended column format including C-tubule data for triplet structures
EXTENDED_COLUMNS = [
    'DoubletNumber', 'X', 'Y', 'Z', 
    'Idx_A', 'Idx_B', 'Idx_C', 
    'Angle', 
    'A_Shift', 'B_Shift', 'C_Shift'
]


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

    centerline_type = line 
    csv_filename = 'cilia.csv'

    
    if centerline_type == 'tip':
        session.logger.info(f"Generating cilia with tip geometry...")
        
        # Generate tip CSV data
        cilia_data_df = generate_cilia_with_tip(
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
        csv_filename = 'ciliatip.csv'
        
    elif centerline_type == 'primarycilia':
        session.logger.info(f"Generating primary cilia ...")
        
        df_template = load_template_data(default_config.PRIMARYCILIA_TEMPLATE)
        # Generate tip CSV data
        cilia_data_df = generate_primary_cilia(
            df_template=df_template,
            cilia_length=length,
            membrane_radius=membrane_radius,
            membrane_fraction=membrane_fraction,
            interval=default_config.MAX_INTERVAL
        )
        # Only turn on writing by default for this
        write_csv = True
        csv_filename = 'primarycilia.csv'
        
    else:
        # Generate structure for 'straight', 'curve', 'sinusoidal', '2Dtemplate'
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
        write_3d_csv(cilia_data_df, csv_filename)
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
    Generate and draw a complete centriole structure with triplet microtubules 
    (typically 9 triplets) in a ChimeraX session.

    This function calculates a 3D centerline based on the `line` type, constructs 
    the triplet geometry DataFrame, applies specific length differences for the 
    B and C tubules (modeling the distal centriole), and renders the surfaces.
    The entire structure is shifted so its end aligns with `z_offset_end`.

    Args:
        session: The current ChimeraX session object.
        
        # --- Centerline Parameters ---
        length (float, optional): Total length of the centriole structure in Å.
                                  Defaults to `default_config.CENTRIOLE_LENGTH`.
        line (str, optional): The type of centerline geometry to use ('straight', 'curve', 'sine', 'template').
                              Defaults to `default_config.CENTRIOLE_LINE`.
        curve_radius (float, optional): Radius of curvature for the 'curve' line type in Å.
                                        Defaults to `default_config.CENTRIOLE_CURVE_RADIUS`.
        sine_frequency (float, optional): Frequency for the 'sine' line type.
                                          Defaults to `default_config.CENTRIOLE_SINE_FREQUENCY`.
        sine_amplitude (float, optional): Amplitude for the 'sine' line type in Å.
                                          Defaults to `default_config.CENTRIOLE_SINE_AMPLITUDE`.
        template_file (str, optional): Filename for the template CSV if `line` is 'template'.
                                       Defaults to `default_config.TEMPLATE_FILE`.

        # --- Structure & Alignment Parameters ---
        num_triplets (int, optional): The number of microtubule triplets (e.g., 9 for a standard centriole).
                                      Defaults to `default_config.CENTRIOLE_NUM_TRIPLETS`.
        centriole_radius (float, optional): The radius of the centriole barrel (distance from center to A-tubule centerline).
                                            Defaults to `default_config.CENTRIOLE_RADIUS`.
        centriole_angle_offset (float, optional): Angular offset in radians for the entire 9-fold structure.
                                                  Defaults to `default_config.CILIA_OFFSET_ANGLE - default_config.CENTRIOLE_OFFSET_ANGLE`.
        z_offset_end (float, optional): The target Z-coordinate where the *end* (most distal point) of the centriole should be placed.
                                        The entire structure is shifted to meet this value.
                                        Defaults to `default_config.CENTRIOLE_Z_OFFSET_END`.

        # --- Triplet Geometry Parameters ---
        triplet_a_radius (float, optional): Radius of the A-tubule in Å.
                                            Defaults to `default_config.CENTRIOLE_TRIPLET_A_RADIUS`.
        triplet_b_radius (float, optional): Radius of the B-tubule in Å.
                                            Defaults to `default_config.CENTRIOLE_TRIPLET_B_RADIUS`.
        triplet_c_radius (float, optional): Radius of the C-tubule in Å.
                                            Defaults to `default_config.CENTRIOLE_TRIPLET_C_RADIUS`.
        triplet_ab_shift (float, optional): The centerline shift magnitude for A and B tubules, defining their overlap.
                                            Defaults to `default_config.CENTRIOLE_TRIPLET_AB_SHIFT`.
        triplet_c_shift (float, optional): The centerline shift magnitude for the C tubule.
                                           Defaults to `default_config.CENTRIOLE_TRIPLET_C_SHIFT`.
        triplet_b_length_diff (float, optional): The length, in Å, by which the B-tubule is shorter than the A-tubule.
                                                 Defaults to `default_config.CENTRIOLE_TRIPLET_B_LENGTH_DIFF`.
        triplet_c_length_diff (float, optional): The length, in Å, by which the C-tubule is shorter than the A-tubule.
                                                 Defaults to `default_config.CENTRIOLE_TRIPLET_C_LENGTH_DIFF`.

        # --- Color Parameters ---
        triplet_a_color (tuple, optional): RGBA color tuple for A-tubules.
                                           Defaults to `default_config.CENTRIOLE_TRIPLET_A_COLOR`.
        triplet_b_color (tuple, optional): RGBA color tuple for B-tubules.
                                           Defaults to `default_config.CENTRIOLE_TRIPLET_B_COLOR`.
        triplet_c_color (tuple, optional): RGBA color tuple for C-tubules.
                                           Defaults to `default_config.CENTRIOLE_TRIPLET_C_COLOR`.

    Returns:
        ChimeraX.core.models.Model | None: The root ChimeraX Model object containing the generated surfaces, or None if generation failed.
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
        template_file=default_config.TEMPLATE_FILE,
        num_doublets=num_triplets,  # Reuse parameter name
        cilia_radius=centriole_radius
    )
    
    # --- Centriole Alignment Shift to match Z_offset_end ---
    z_shift = z_offset_end - length
    structure['centerline'][:, 2] += z_shift
    session.logger.info(f"Centriole length {length} Å shifted by {z_shift:.1f} Å to end at Z={z_offset_end:.1f} Å.")
    # --- END Shift ---
    
    # --- Triplet DataFrame Construction (New Logic) ---
    session.logger.info(f"Building DataFrame for {num_triplets} triplets...")
    all_data = []
    idx_b_start = 1
    idx_c_start = 2

    for triplet_info in structure['doublets']:  # Reuse 'doublets' key
        
        # Get the shifted centerline for this triplet
        triplet_centerline = get_doublet_centerline(
            structure['centerline'],
            triplet_info['angle'],
            triplet_info['shift_distance']
        )
        
        if len(triplet_centerline) != len(structure['centerline']):
            triplet_centerline = triplet_centerline[:len(structure['centerline'])]
        
        # Calculate arc length along this triplet's centerline
        triplet_segment_lengths = np.linalg.norm(np.diff(triplet_centerline, axis=0), axis=1)
        triplet_cumulative_length = np.concatenate(([0], np.cumsum(triplet_segment_lengths)))
        total_triplet_length = triplet_cumulative_length[-1]
        
        # Calculate B-tubule END index
        target_b_length = total_triplet_length - triplet_b_length_diff
        idx_b_end = np.searchsorted(triplet_cumulative_length, target_b_length, side='right')
        idx_b_end = max(2, min(idx_b_end, len(triplet_centerline)))
        
        # Calculate C-tubule END index
        target_c_length = total_triplet_length - triplet_c_length_diff
        idx_c_end = np.searchsorted(triplet_cumulative_length, target_c_length, side='right')
        idx_c_end = max(2, min(idx_c_end, len(triplet_centerline)))
        
        # Add data for this triplet
        triplet_number = triplet_info['index'] + 1
        angle = triplet_info['angle'] + centriole_angle_offset
        
        # Populate data point by point
        for i, point in enumerate(triplet_centerline):
            x, y, z = point
            
            idx_a_val = 1 # A-tubule is full length
            idx_b_val = 1 if (idx_b_start <= i < idx_b_end) else 0
            idx_c_val = 1 if (idx_c_start <= i < idx_c_end) else 0
            
            # Use triplet_ab_shift for A/B magnitude, and triplet_c_shift for C magnitude
            all_data.append([
                triplet_number, x, y, z, 
                idx_a_val, idx_b_val, idx_c_val, # Idx_A, Idx_B, Idx_C
                angle, 
                -triplet_ab_shift, triplet_ab_shift, triplet_c_shift # A_Shift, B_Shift, C_Shift
            ])
            
    # Create DataFrame
    centriole_data_df = pd.DataFrame(all_data, columns=EXTENDED_COLUMNS)
    
    # --- Centriole Drawing (Now calls _ciliabuild_from_df) ---
    centriole_root = _ciliabuild_from_df(
        session=session,
        df=centriole_data_df,
        draw_central_pair=False,  # Centrioles typically do not have a central pair
        membrane=False,           # Centrioles typically do not have a membrane
        # Doublet parameters are used for A/B, but triplet-specific values are passed
        doublet_a_radius=triplet_a_radius,
        doublet_b_radius=triplet_b_radius,
        doublet_shift=triplet_ab_shift, 
        doublet_a_color=triplet_a_color,
        doublet_b_color=triplet_b_color,
        # Triplet-specific C-tubule parameters (used as the 'flag' for Triplet mode)
        triplet_c_radius=triplet_c_radius, 
        triplet_c_shift=triplet_c_shift,
        triplet_c_color=triplet_c_color
    )
    
    if centriole_root is None:
        return None
    
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
            membrane_color=default_config.CILIA_MEMBRANE_COLOR,
            # Triplet Parameters (for centriole CSV files)
            triplet_c_radius=default_config.CENTRIOLE_TRIPLET_C_RADIUS, 
            triplet_c_shift=default_config.CENTRIOLE_TRIPLET_C_SHIFT,
            triplet_c_color=default_config.CENTRIOLE_TRIPLET_C_COLOR
            ):
    """
    Generate and draw a complete cilia/centriole structure from a template CSV file.
    
    Loads a pre-computed structure from a CSV file and renders it in 3D. 
    Automatically detects whether the file contains doublet (cilia) or triplet 
    (centriole) data based on the presence of Idx_C column values.
    
    Parameters
    ----------
    session : chimerax.core.session.Session
        ChimeraX session object for model management and logging
    template_csv : str, optional
        Path to CSV file containing structure data (default: 'cilia.csv')
    draw_central_pair : bool, optional
        Whether to render central pair microtubules
    membrane : bool, optional
        Whether to render ciliary membrane
    membrane_fraction : float, optional
        Unused parameter (kept for API compatibility)
    membrane_radius : float, optional
        Outer radius of membrane surface (Å)
    doublet_a_radius : float, optional
        Outer radius of A-tubules (Å)
    doublet_b_radius : float, optional
        Outer radius of B-tubules (Å)
    doublet_shift : float, optional
        Lateral displacement of A/B tubules from unit center (Å)
    cp_radius : float, optional
        Outer radius of central pair tubules (Å)
    cp_shift : float, optional
        Distance of each central tubule from central axis (Å)
    doublet_a_color : tuple, optional
        RGBA color for A-tubules (0-255)
    doublet_b_color : tuple, optional
        RGBA color for B-tubules (0-255)
    cp_color : tuple, optional
        RGBA color for central pair (0-255)
    membrane_color : tuple, optional
        RGBA color for membrane (0-255)
    triplet_c_radius : float, optional
        Outer radius of C-tubules for triplet structures (Å)
    triplet_c_shift : float, optional
        Lateral displacement of C-tubule from unit center (Å)
    triplet_c_color : tuple, optional
        RGBA color for C-tubules (0-255)
    
    Returns
    -------
    chimerax.core.models.Surface
        Root surface model containing all structure components, or None if failed
    
    Notes
    -----
    CSV File Format:
        Required columns: DoubletNumber, X, Y, Z, Idx_A, Idx_B, Angle, 
                         A_Shift, B_Shift
        Optional columns: Idx_C, C_Shift (for triplet structures)
        
    CSV files generated by ciliabuild() with write_csv=True are compatible
    with this function.
    
    Examples
    --------
    Load a cilia structure:
        ciliabuild_from_csv template_csv 'cilia.csv'
    
    Load a centriole structure with custom colors:
        ciliabuild_from_csv template_csv 'centriole.csv' \\
                            triplet_a_color 255,100,100,255 \\
                            triplet_c_color 100,255,100,255
    """
    
    session.logger.info(f"Loading cilia structure from template: {template_csv}")
    
    # Load CSV file
    try:
        df = pd.read_csv(template_csv)
    except FileNotFoundError:
        session.logger.error(f"Template file not found: {template_csv}")
        return None
    except Exception as e:
        session.logger.error(f"Error reading template file: {e}")
        return None
    
    # Render structure using internal function
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
        membrane_color=membrane_color,
        triplet_c_radius=triplet_c_radius, 
        triplet_c_shift=triplet_c_shift,
        triplet_c_color=triplet_c_color
    )


# ============================================================================
# CHIMERAX COMMAND DESCRIPTORS
# ============================================================================
# Define command-line argument specifications for ChimeraX integration

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
        ('z_offset_end', FloatArg),
        ('triplet_a_color', Color8Arg),
        ('triplet_b_color', Color8Arg),
        ('triplet_c_color', Color8Arg)
    ],
    synopsis='Generate complete centriole structure with triplet microtubules'
)

ciliabuild_from_csv_desc = CmdDesc(
    keyword=[
        ('template_csv', StringArg),
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
        ('membrane_color', Color8Arg),
        ('triplet_c_radius', FloatArg), 
        ('triplet_c_shift', FloatArg),
        ('triplet_c_color', Color8Arg)
    ],
    synopsis='Generate a cilia/centriole structure from a template CSV file'
)

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
            membrane_color=default_config.CILIA_MEMBRANE_COLOR,
            # Triplet Parameters (optional - presence determines triplet mode)
            triplet_c_radius=None, 
            triplet_c_shift=None,
            triplet_c_color=None
            ):
    """
    Internal function to render cilia/centriole structures from DataFrame geometry data.
    
    This is the core rendering engine that converts tabular geometric data into 3D 
    ChimeraX models. It automatically detects structure type (doublet vs triplet) and 
    creates a hierarchical model with proper grouping. Used by ciliabuild(), 
    centriolebuild(), and ciliabuild_from_csv().
    
    Parameters
    ----------
    session : chimerax.core.session.Session
        ChimeraX session for model management and logging
    df : pandas.DataFrame
        Structure geometry data containing centerline points and tubule indices
        Required columns: DoubletNumber, X, Y, Z, Idx_A, Idx_B, Angle, A_Shift, B_Shift
        Optional columns: Idx_C, C_Shift (for triplet structures)
    draw_central_pair : bool, optional
        Whether to render central pair tubules (DoubletNumber = -1)
    membrane : bool, optional
        Whether to render ciliary membrane (DoubletNumber = 0)
    membrane_fraction : float, optional
        Legacy parameter, unused (kept for API compatibility)
    membrane_radius : float, optional
        Outer radius of membrane surface (Å)
    doublet_a_radius : float, optional
        Outer radius of A-tubules (Å)
    doublet_b_radius : float, optional
        Outer radius of B-tubules (Å)
    doublet_shift : float, optional
        Lateral displacement of A/B tubules from unit centerline (Å)
    cp_radius : float, optional
        Outer radius of central pair tubules (Å)
    cp_shift : float, optional
        Distance of each central tubule from central axis (Å)
    doublet_a_color : tuple of int, optional
        RGBA color for A-tubules (values 0-255)
    doublet_b_color : tuple of int, optional
        RGBA color for B-tubules (values 0-255)
    cp_color : tuple of int, optional
        RGBA color for central pair tubules (values 0-255)
    membrane_color : tuple of int, optional
        RGBA color for membrane surface (values 0-255)
    triplet_c_radius : float or None, optional
        Outer radius of C-tubules (Å). If not None, enables triplet rendering mode
    triplet_c_shift : float or None, optional
        Lateral displacement of C-tubule from unit centerline (Å)
    triplet_c_color : tuple of int or None, optional
        RGBA color for C-tubules (values 0-255)
    
    Returns
    -------
    chimerax.core.models.Surface or None
        Root surface model containing all rendered components in a hierarchy, 
        or None if rendering fails
    
    Structure Type Detection
    ------------------------
    The function automatically determines rendering mode:
    - **Doublet mode**: When Idx_C column is absent or sum(Idx_C) == 0
    - **Triplet mode**: When sum(Idx_C) > 0 in the DataFrame
    
    DoubletNumber Encoding
    ----------------------
    Special values in the DoubletNumber column:
        -2 : Cap complex (terminal structure for tip geometry)
        -1 : Central pair (C1 and C2 singlet tubules)
         0 : Ciliary membrane
        >0 : Doublet/Triplet unit number (1-based indexing)
    
    Model Hierarchy
    ---------------
    Creates a hierarchical structure in ChimeraX::
    
        Root Model (Cilia_from_template or Centriole_from_template)
        ├── DMT1 or TMT1 (Doublet/Triplet Microtubule 1)
        │   ├── A_tubule (always present)
        │   ├── B_tubule (if Idx_B > 0)
        │   └── C_tubule (triplets only, if Idx_C > 0)
        ├── DMT2 or TMT2
        │   └── ...
        ├── ...
        ├── Central Pair (if draw_central_pair=True and DoubletNumber=-1 exists)
        │   ├── C1 (first central tubule)
        │   ├── C2 (second central tubule)
        │   └── CapComplex (if DoubletNumber=-2 exists)
        └── Membrane (if membrane=True and DoubletNumber=0 exists)
    
    Naming Convention
    -----------------
    - DMT = Doublet Microtubule (for cilia)
    - TMT = Triplet Microtubule (for centrioles)
    
    Notes
    -----
    - This function does NOT modify the input DataFrame
    - Tubules are only rendered when sufficient points exist (≥2)
    - Units with insufficient A-tubule points are skipped with a warning
    - Color tuples must contain 4 integers (RGBA) with values 0-255
    - All distance measurements are in Ångströms
    - The function relies on draw_tubules() and draw_membrane() for actual geometry generation
    
    Examples
    --------
    Render doublet structure from DataFrame:
        >>> df = pd.read_csv('cilia_geometry.csv')
        >>> model = _ciliabuild_from_df(session, df, 
        ...                              draw_central_pair=True,
        ...                              membrane=True)
    
    Render triplet structure with custom C-tubule parameters:
        >>> df = pd.read_csv('centriole_geometry.csv')
        >>> model = _ciliabuild_from_df(session, df,
        ...                              triplet_c_radius=135.0,
        ...                              triplet_c_shift=200.0,
        ...                              triplet_c_color=(179, 179, 255, 255))
    
    See Also
    --------
    ciliabuild : High-level function for generating cilia structures
    centriolebuild : High-level function for generating centriole structures
    ciliabuild_from_csv : Load and render structures from CSV files
    draw_tubules : Low-level function for rendering individual tubules
    draw_membrane : Low-level function for rendering membrane surfaces
    """
    
    # Validate DataFrame columns
    if not all(col in df.columns for col in REQUIRED_COLUMNS):
        session.logger.error(f"DataFrame must contain columns: {REQUIRED_COLUMNS}")
        return None
        
    # --- GLOBAL STRUCTURE TYPE DETERMINATION ---
    has_c_tubule_column = 'Idx_C' in df.columns
    is_triplet_structure = False
    
    if has_c_tubule_column:
        # Check if the sum of all Idx_C values across the whole DataFrame is > 0
        if df['Idx_C'].sum() > 0:
            is_triplet_structure = True
            session.logger.info("Global Triplet structure detected (sum(Idx_C) > 0).")

    # Set root name and prefix based on global determination
    if is_triplet_structure:
        root_name = "Centriole_from_template"
        name_prefix = "TMT"
    else:
        root_name = "Cilia_from_template"
        name_prefix = "DMT"
        
    cilia_root = Surface(root_name, session)
    session.models.add([cilia_root])
    
    # Get unique doublet/triplet numbers
    unit_numbers = sorted(df.loc[df['DoubletNumber'] > 0, 'DoubletNumber'].unique())
    num_units = len(unit_numbers)
    
    session.logger.info(f"Found {num_units} {name_prefix} units in data")
    
    # Process each unit
    for unit_num in unit_numbers:
        unit_data = df[df['DoubletNumber'] == unit_num]
        
        # Extract points where Idx_A = 1 (A-tubule)
        a_data = unit_data[unit_data['Idx_A'] == 1]
        unit_centerline_a = a_data[['X', 'Y', 'Z']].values
        
        # Extract points where Idx_B = 1 (B-tubule)
        b_data = unit_data[unit_data['Idx_B'] == 1]
        unit_centerline_b = b_data[['X', 'Y', 'Z']].values
        
        angle = unit_data['Angle'].iloc[0]
        a_shift = unit_data['A_Shift'].iloc[0]
        b_shift = unit_data['B_Shift'].iloc[0]
        
        if len(unit_centerline_a) < 2:
            session.logger.warning(f"Unit {unit_num} has insufficient A-tubule points, skipping")
            continue
        
        unit_surfs = []
        
        # Draw A-tubule 
        a_surfs = draw_tubules(
            session=session,
            length=None,
            interval=default_config.MAX_INTERVAL,
            centerline_points=unit_centerline_a,
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
            unit_surfs.extend(a_surfs)
        
        # Draw B-tubule if it has points
        if len(unit_centerline_b) >= 2:
            b_surfs = draw_tubules(
                session=session,
                length=None,
                interval=default_config.MAX_INTERVAL,
                centerline_points=unit_centerline_b,
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
                unit_surfs.extend(b_surfs)

        # Draw C-tubule if it's a global triplet structure AND this specific unit has C-tubule points
        if is_triplet_structure:
            c_data = unit_data[unit_data['Idx_C'] == 1] 
            unit_centerline_c = c_data[['X', 'Y', 'Z']].values
            
            # Draw C-tubule only if we have enough points for a centerline
            if len(unit_centerline_c) >= 2:
                session.logger.info(f"Unit {unit_num} drawing C-tubule")
                
                c_surfs = draw_tubules(
                    session=session,
                    length=None,
                    interval=default_config.MAX_INTERVAL,
                    centerline_points=unit_centerline_c,
                    angle=angle,
                    radii=[triplet_c_radius],
                    shift_distances=[triplet_c_shift],
                    length_diffs=None,
                    tubule_names=[f"C_tubule"],
                    colors=[triplet_c_color],
                    group_name=f"C_tubule",
                    add_to_session=False
                )
                if c_surfs:
                    unit_surfs.extend(c_surfs)
        
        # Add all surfaces as a single group, using the globally determined prefix
        if unit_surfs:
            session.models.add_group(unit_surfs, parent=cilia_root, name=f"{name_prefix}{int(unit_num)}")
            session.logger.info(f"Added {name_prefix}{int(unit_num)}")
        
    # Draw central pair if requested (using the DoubletNumber = -1)
    cp_data = df[df['DoubletNumber'] == -1]
    cap_data = df[df['DoubletNumber'] == -2]
    cap_surf = None

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
        
        if len(cap_data) > 0:
            cap_surf = generate_capsule_surface(session, tuple(cap_data[['X', 'Y', 'Z']].iloc[0]), capsule_length=default_config.TIP_CAP_LENGTH,
                            radius=default_config.TIP_CAP_RADIUS, color=default_config.CILIA_CAP_COLOR, name="CapComplex", add_to_session=False)
        
        if cap_surf:
            cp_surfs.append(cap_surf)

        session.models.add_group(cp_surfs, parent=cilia_root, name="Central Pair")
        session.logger.info(f"Added central pair")
            
    # Draw membrane if requested (using DoubletNumber = 0) 
    membrane_data = df[df['DoubletNumber'] == 0]
    if membrane and len(membrane_data) > 0:
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
    cilia_root.name = f"{root_name} {model_id}"
    
    session.logger.info(f"Model from template generated successfully!")
    session.logger.info(f"  Units: {num_units}")
    
    return cilia_root