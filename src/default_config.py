"""
Default configuration parameters for CiliaBuilder.

This module defines all default values for cilia and centriole structure
generation, including geometry, colors, and rendering parameters.

Units:
    - All linear measurements are in Ångströms (Å)
    - All angles are in degrees (°)
    - All colors are RGBA tuples with integer values 0-255

Structure:
    - GENERAL PARAMETERS: Core settings used across all structures
    - CILIA DEFAULTS: Parameters for 9+2 cilia (9 doublets + 2 central pair)
    - CENTRIOLE DEFAULTS: Parameters for 9x3 centriole (9 triplets)
"""

# ============================================================================
# GENERAL PARAMETERS
# ============================================================================

# Sampling interval for centerline point generation
# This determines the spacing between consecutive points along the centerline
# Smaller values = smoother curves but more computation
MAX_INTERVAL = 20.0  # Ångströms

# Angular offset applied to cilia doublet positioning
# Rotates the entire doublet array around the central axis
CILIA_OFFSET_ANGLE = 90.0  # Degrees

# Z-coordinate offset for structure alignment
# Used to position centriole and cilia relative to coordinate origin
Z_OFFSET_END = 2 * MAX_INTERVAL  # Ångströms (40 Å)

# Default template file for custom centerline geometry
TEMPLATE_FILE = 'template.csv'

# Primary cilia template file (contains singlet configuration data)
PRIMARYCILIA_TEMPLATE = 'primarycilia_template.csv'


# ============================================================================
# CILIA DEFAULTS (9+2 Configuration)
# ============================================================================
# Standard motile cilia with 9 peripheral doublets and 2 central singlets

# --- Structure Parameters ---
CILIA_LENGTH = 10000  # Total length of cilia structure
CILIA_NUM_DOUBLETS = 9  # Number of peripheral doublets
CILIA_RADIUS = 875.0  # Radial distance from center to doublet centers
CILIA_DRAW_CENTRAL_PAIR = True  # Whether to include C1/C2 central pair
CILIA_MEMBRANE = True  # Whether to draw ciliary membrane
CILIA_MEMBRANE_FRACTION = 0.7  # Fraction of cilia length covered by membrane (0.0-1.0)
CILIA_MEMBRANE_RADIUS = 1100  # Outer radius of ciliary membrane

# --- Centerline Parameters ---
# Define the path that the cilia structure follows
CILIA_LINE = 'straight'  # Options: 'straight', 'curve', 'sinusoidal', 'tip', 'primarycilia'
CILIA_CURVE_RADIUS = 10000.0  # Radius of curvature for 'curve' centerline type
CILIA_SINE_FREQUENCY = 1.5  # Number of oscillations per length unit (for 'sinusoidal')
CILIA_SINE_AMPLITUDE = 1000.0  # Maximum displacement from center (for 'sinusoidal')

# --- Doublet Microtubule Geometry ---
CILIA_DOUBLET_A_RADIUS = 125.0  # Outer radius of A-tubule
CILIA_DOUBLET_B_RADIUS = 135.0  # Outer radius of B-tubule (slightly larger due to incomplete structure)
CILIA_DOUBLET_SHIFT = 70.0  # Lateral shift distance for A and B tubules from doublet center
CILIA_DOUBLET_LENGTH_DIFF = 5.0  # A-tubule extends this much beyond B-tubule
CILIA_CP_DOUBLET_LENGTH_DIFF = 300.0  # Central pair extends this much beyond doublet tips

# --- Central Pair Geometry ---
CILIA_CP_RADIUS = 125.0  # Outer radius of C1 and C2 singlet tubules
CILIA_CP_SHIFT = 160.0  # Distance of each central tubule from the central axis
CILIA_CP_BRIDGE_RADIUS = 30 # Radius of rung linking C1 & C2
CILIA_CP_RUNG_PERIODICITY = 320 # Run periodicity 160 or 320

# --- Tip Geometry ---
# Parameters for generating tapered cilia tips
TIP_FINAL_RADIUS = CILIA_CP_RADIUS + CILIA_CP_SHIFT + 120  # Final tip radius (405 Å)
TIP_TRANSITION_LENGTH = 1500  # Length over which taper occurs
TIP_INITIAL_LENGTH = 300  # Straight section before taper begins
TIP_LENGTH = 5000  # Total length of tip structure
TIP_CAP_RADIUS = CILIA_CP_RADIUS + CILIA_CP_SHIFT + 100  # Radius of terminal cap (385 Å)
TIP_CAP_LENGTH = TIP_CAP_RADIUS * 2.5  # Length of terminal cap (962.5 Å)

# --- Colors (RGBA: 0-255) ---
CILIA_DOUBLET_A_COLOR = (100, 100, 255, 255)  # Blue for A-tubules
CILIA_DOUBLET_B_COLOR = (100, 100, 255, 255)  # Blue for B-tubules
CILIA_CP_COLOR = (255, 255, 100, 255)  # Yellow for central pair
CILIA_MEMBRANE_COLOR = (105, 105, 105, 255)  # Gray for membrane
CILIA_CAP_COLOR = (180, 180, 30, 255)  # Dark yellow for tip cap


# ============================================================================
# CENTRIOLE DEFAULTS (9x3 Configuration)
# ============================================================================
# Centriole structure with 9 peripheral triplet microtubules

# --- Structure Parameters ---
CENTRIOLE_LENGTH = 3000  # Total length of centriole structure
CENTRIOLE_NUM_TRIPLETS = 9  # Number of triplet microtubules
CENTRIOLE_RADIUS = 875.0  # Radial distance from center to triplet centers
CENTRIOLE_OFFSET_ANGLE = 30.0  # Angular offset relative to cilia orientation

# --- Centerline Parameters ---
CENTRIOLE_LINE = 'straight'  # Options: 'straight', 'curve', 'sinusoidal'
CENTRIOLE_CURVE_RADIUS = 10000.0  # Radius of curvature for 'curve' centerline type
CENTRIOLE_SINE_FREQUENCY = 1.0  # Number of oscillations (for 'sinusoidal')
CENTRIOLE_SINE_AMPLITUDE = 500.0  # Maximum displacement (for 'sinusoidal')

# --- Triplet Microtubule Geometry ---
CENTRIOLE_TRIPLET_A_RADIUS = 125.0  # Outer radius of A-tubule (complete tubule)
CENTRIOLE_TRIPLET_B_RADIUS = 135.0  # Outer radius of B-tubule (incomplete tubule)
CENTRIOLE_TRIPLET_C_RADIUS = 135.0  # Outer radius of C-tubule (incomplete tubule)
CENTRIOLE_TRIPLET_AB_SHIFT = 70.0  # Lateral shift for A and B tubules from triplet center
CENTRIOLE_TRIPLET_C_SHIFT = 200.0  # Lateral shift for C-tubule from triplet center
CENTRIOLE_TRIPLET_B_LENGTH_DIFF = MAX_INTERVAL  # B-tubule is shorter than A by this amount (20 Å)
CENTRIOLE_TRIPLET_C_LENGTH_DIFF = 2 * MAX_INTERVAL  # C-tubule is shorter than A by this amount (40 Å)
CENTRIOLE_Z_OFFSET_END = 2 * MAX_INTERVAL  # Z-coordinate where centriole end is positioned (40 Å)

# --- Colors (RGBA: 0-255) ---
CENTRIOLE_TRIPLET_A_COLOR = (100, 100, 255, 255)  # Blue for A-tubules
CENTRIOLE_TRIPLET_B_COLOR = (100, 100, 255, 255)  # Blue for B-tubules
CENTRIOLE_TRIPLET_C_COLOR = (179, 179, 255, 255)  # Light blue for C-tubules