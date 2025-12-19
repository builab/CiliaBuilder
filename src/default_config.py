#default_config.py
# ============================================================================
# GENERAL PARAMETERS
# ============================================================================
# Define the optimal sampling interval for smoothness (must match curve.py)
MAX_INTERVAL = 20.0  # Angstroms interval for centerline points
CILIA_OFFSET_ANGLE = 90.0  # Offset angle for cilia doublets
Z_OFFSET_END = 2 * MAX_INTERVAL  # Offset to draw centriole & cilia

# ============================================================================
# CILIA DEFAULTS (9x2 + 2)
# ============================================================================
# Structure Parameters
CILIA_LENGTH = 15000  # Angstroms
CILIA_NUM_DOUBLETS = 9
CILIA_RADIUS = 875.0  # Angstroms
CILIA_DRAW_CENTRAL_PAIR = True
CILIA_MEMBRANE = True
CILIA_MEMBRANE_FRACTION = 0.7
CILIA_MEMBRANE_RADIUS = 1100  # Angstroms

# Centerline Parameters
CILIA_LINE = 'straight'
CILIA_CURVE_RADIUS = 10000.0  # Angstroms
CILIA_SINE_FREQUENCY = 1.5
CILIA_SINE_AMPLITUDE = 1000.0  # Angstroms
TEMPLATE_FILE = 'template.csv'

# Doublet Geometry
CILIA_DOUBLET_A_RADIUS = 125.0  # Angstroms
CILIA_DOUBLET_B_RADIUS = 135.0  # Angstroms
CILIA_DOUBLET_SHIFT = 70.0  # Angstroms
CILIA_DOUBLET_LENGTH_DIFF = 5.0  # A - B length difference (Angstroms)
CILIA_CP_DOUBLET_LENGTH_DIFF = 0.0  # CP - doublet length difference (Angstroms)


# Central Pair Geometry
CILIA_CP_RADIUS = 125.0  # Angstroms
CILIA_CP_SHIFT = 160.0  # Angstroms

# Tip Geometry
TIP_TRANSITION_RADIUS = 0.75*CILIA_RADIUS
TIP_FINAL_RADIUS = CILIA_CP_RADIUS*2 + CILIA_CP_SHIFT + 5
TIP_LENGTH = 10000

# Colors (RGBA, 0-255)
CILIA_DOUBLET_A_COLOR = (100, 100, 255, 255)
CILIA_DOUBLET_B_COLOR = (100, 100, 255, 255)
CILIA_CP_COLOR = (255, 255, 100, 255)
CILIA_MEMBRANE_COLOR = (105, 105, 105, 255)

# ============================================================================
# CENTRIOLE DEFAULTS (9x3)
# ============================================================================
# Structure Parameters
CENTRIOLE_LENGTH = 5000  # Angstroms
CENTRIOLE_NUM_TRIPLETS = 9
CENTRIOLE_RADIUS = 875.0  # Angstroms
CENTRIOLE_OFFSET_ANGLE = 30.0  # Degrees

# Centerline Parameters
CENTRIOLE_LINE = 'straight'
CENTRIOLE_CURVE_RADIUS = 10000.0  # Angstroms
CENTRIOLE_SINE_FREQUENCY = 1.0
CENTRIOLE_SINE_AMPLITUDE = 500.0  # Angstroms

# Triplet Geometry
CENTRIOLE_TRIPLET_A_RADIUS = 125.0  # Angstroms
CENTRIOLE_TRIPLET_B_RADIUS = 135.0  # Angstroms
CENTRIOLE_TRIPLET_C_RADIUS = 135.0  # Angstroms
CENTRIOLE_TRIPLET_AB_SHIFT = 70.0  # Angstroms
CENTRIOLE_TRIPLET_C_SHIFT =  200.0  # Angstroms
CENTRIOLE_TRIPLET_B_LENGTH_DIFF = MAX_INTERVAL  # Angstroms
CENTRIOLE_TRIPLET_C_LENGTH_DIFF = 2*MAX_INTERVAL  # Angstroms
CENTRIOLE_Z_OFFSET_END = 2*MAX_INTERVAL  # Target Z-coordinate for centriole end (Angstroms)

# Colors (RGBA, 0-255)
CENTRIOLE_TRIPLET_A_COLOR = (100, 100, 255, 255)
CENTRIOLE_TRIPLET_B_COLOR = (100, 100, 255, 255)
CENTRIOLE_TRIPLET_C_COLOR = (179, 179, 255, 255)