# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QLabel, QLineEdit, QPushButton,
                          QComboBox, QCheckBox, QGroupBox)
from Qt.QtCore import Qt

# Define the optimal sampling interval for smoothness (10 Angstroms)
# This is used for calculating point counts for length truncation in draw.py
MAX_INTERVAL = 10.0 

class CiliaSim(ToolInstance):
    """
    CiliaSim tool for generating cilia microtubule models
    """

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        # Initialize base class
        super().__init__(session, tool_name)

        # Set name displayed on title bar
        self.display_name = "CiliaSim"

        # Create the main window
        self.tool_window = MainToolWindow(self)

        # Build the user interface
        self._build_ui()

    def _build_ui(self):
        """Build the user interface"""
        
        # Main layout
        main_layout = QVBoxLayout()
        main_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # Title
        title_label = QLabel("Cilia Microtubule Model Generator")
        title_label.setStyleSheet("font-size: 14pt; font-weight: bold; margin-bottom: 10px;")
        main_layout.addWidget(title_label)

        # --- 1. General Settings Group ---
        general_group = QGroupBox("General & Doublet Structure")
        general_layout = QGridLayout()

        # Cilia Length
        general_layout.addWidget(QLabel("Total Length (Å):"), 0, 0)
        self.length_input = QLineEdit("15000")
        general_layout.addWidget(self.length_input, 0, 1)

        # Number of Doublets
        general_layout.addWidget(QLabel("Number of Doublets (9):"), 1, 0)
        self.num_doublets_input = QLineEdit("9")
        general_layout.addWidget(self.num_doublets_input, 1, 1)

        # Cilia Radius (Distance from center to doublet centerline)
        general_layout.addWidget(QLabel("Doublet Ring Radius (Å):"), 2, 0)
        self.cilia_radius_input = QLineEdit("875.0")
        general_layout.addWidget(self.cilia_radius_input, 2, 1)
        
        # Draw Central Pair Checkbox
        self.draw_cp_check = QCheckBox("Draw Central Pair (C1/C2)")
        self.draw_cp_check.setChecked(True)
        general_layout.addWidget(self.draw_cp_check, 3, 0, 1, 2)
        
        # Doublet Colors (existing)
        general_layout.addWidget(QLabel("Doublet Color (R,G,B,A):"), 4, 0)
        self.color_a_input = QLineEdit("100,100,255,255")
        self.color_a_input.setToolTip("Color for both A and B tubules")
        general_layout.addWidget(self.color_a_input, 4, 1)
        
        general_group.setLayout(general_layout)
        main_layout.addWidget(general_group)

        # --- 2. Centerline Type Group ---
        centerline_group = QGroupBox("Centerline Type & Geometry")
        centerline_layout = QGridLayout()
        
        # Centerline Type (Dropdown)
        centerline_layout.addWidget(QLabel("Centerline Type:"), 0, 0)
        self.line_type_combo = QComboBox()
        self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal'])
        self.line_type_combo.setCurrentText('straight')
        self.line_type_combo.currentIndexChanged.connect(self._toggle_centerline_inputs)
        centerline_layout.addWidget(self.line_type_combo, 0, 1)
        
        # Curve Input (Angle is now calculated from length/radius)
        centerline_layout.addWidget(QLabel("Curve Radius (Å):"), 1, 0)
        self.curve_radius_input = QLineEdit("10000.0")
        self.curve_radius_input.setToolTip("R in s = R*theta arc length formula")
        centerline_layout.addWidget(self.curve_radius_input, 1, 1)

        # Sinusoidal Inputs
        centerline_layout.addWidget(QLabel("Sine Frequency:"), 2, 0)
        self.sine_frequency_input = QLineEdit("2.0")
        centerline_layout.addWidget(self.sine_frequency_input, 2, 1)
        centerline_layout.addWidget(QLabel("Sine Amplitude (Å):"), 3, 0)
        self.sine_amplitude_input = QLineEdit("500.0")
        centerline_layout.addWidget(self.sine_amplitude_input, 3, 1)
        
        centerline_group.setLayout(centerline_layout)
        main_layout.addWidget(centerline_group)
        
        # Initial state of centerline inputs
        self._toggle_centerline_inputs()


        # --- 3. Microtubule Geometry Group ---
        mt_group = QGroupBox("Microtubule Dimensions")
        mt_layout = QGridLayout()

        # Doublet A/B Radii
        mt_layout.addWidget(QLabel("A-Tubule Radius (Å):"), 0, 0)
        self.d_a_radius_input = QLineEdit("125.0")
        mt_layout.addWidget(self.d_a_radius_input, 0, 1)
        
        mt_layout.addWidget(QLabel("B-Tubule Radius (Å):"), 1, 0)
        self.d_b_radius_input = QLineEdit("145.0")
        mt_layout.addWidget(self.d_b_radius_input, 1, 1)
        
        # Doublet Shift
        mt_layout.addWidget(QLabel("Doublet A/B Shift (Å):"), 2, 0)
        self.doublet_shift_input = QLineEdit("70.0")
        self.doublet_shift_input.setToolTip("Distance from the doublet centerline for both A and B tubules")
        mt_layout.addWidget(self.doublet_shift_input, 2, 1)
        
        # Doublet Length Difference
        mt_layout.addWidget(QLabel("A-B Length Diff (Å):"), 3, 0)
        self.doublet_length_diff_input = QLineEdit("250.0")
        mt_layout.addWidget(self.doublet_length_diff_input, 3, 1)
        
        # Central Pair Radius
        mt_layout.addWidget(QLabel("CP Tubule Radius (Å):"), 4, 0)
        self.cp_radius_input = QLineEdit("125.0")
        mt_layout.addWidget(self.cp_radius_input, 4, 1)
        
        # Central Pair Shift
        mt_layout.addWidget(QLabel("CP Tubule Shift (Å):"), 5, 0)
        self.cp_shift_input = QLineEdit("160.0")
        mt_layout.addWidget(self.cp_shift_input, 5, 1)

        mt_group.setLayout(mt_layout)
        main_layout.addWidget(mt_group)


        # Generate button
        generate_button = QPushButton("Generate Cilia Model")
        generate_button.setStyleSheet("font-size: 12pt; padding: 10px; margin-top: 10px;")
        generate_button.clicked.connect(self._generate_model)
        main_layout.addWidget(generate_button)

        # Status label
        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color: blue; font-style: italic;")
        main_layout.addWidget(self.status_label)

        # Add stretch to push everything to the top
        main_layout.addStretch()

        # Set the layout
        self.tool_window.ui_area.setLayout(main_layout)

        # Show the window
        self.tool_window.manage('side')
        
    def _toggle_centerline_inputs(self):
        """Enable/disable curve or sine inputs based on line type dropdown"""
        line_type = self.line_type_combo.currentText()
        
        is_curve = line_type == 'curve'
        self.curve_radius_input.setEnabled(is_curve)
        
        is_sine = line_type == 'sinusoidal'
        self.sine_frequency_input.setEnabled(is_sine)
        self.sine_amplitude_input.setEnabled(is_sine)

    def _parse_color(self, color_string):
        """Parse color string 'R,G,B,A' to tuple of integers"""
        try:
            parts = [int(x.strip()) for x in color_string.split(',')]
            if len(parts) != 4:
                raise ValueError("Color must have 4 components (R,G,B,A)")
            if not all(0 <= x <= 255 for x in parts):
                raise ValueError("Color values must be between 0 and 255")
            return tuple(parts)
        except Exception as e:
            raise ValueError(f"Invalid color format: {e}")

    def _generate_model(self):
        """Generate the cilia model with current parameters"""
        try:
            # --- 1. Get and Validate Parameters ---
            
            # General
            length = float(self.length_input.text())
            num_doublets = int(self.num_doublets_input.text())
            cilia_radius = float(self.cilia_radius_input.text())
            draw_central_pair = self.draw_cp_check.isChecked()
            
            # Use color_a for both A and B, or you can add a separate input later
            doublet_color = self._parse_color(self.color_a_input.text())
            
            # Centerline
            centerline_type = self.line_type_combo.currentText()
            curve_radius = float(self.curve_radius_input.text())
            sine_frequency = float(self.sine_frequency_input.text())
            sine_amplitude = float(self.sine_amplitude_input.text())
            
            # MT Geometry
            doublet_a_radius = float(self.d_a_radius_input.text())
            doublet_b_radius = float(self.d_b_radius_input.text())
            doublet_shift = float(self.doublet_shift_input.text())
            doublet_length_diff = float(self.doublet_length_diff_input.text())
            cp_radius = float(self.cp_radius_input.text())
            cp_shift = float(self.cp_shift_input.text())
            
            if length <= 0 or num_doublets < 0 or cilia_radius <= 0:
                 raise ValueError("Length, num doublets, and cilia radius must be positive.")

            # Update status
            self.status_label.setText("Generating model...")
            self.status_label.setStyleSheet("color: blue; font-style: italic;")

            # --- 2. Core Generation Logic ---
            
            from .curve import generate_cilia_structure, get_doublet_centerline
            from .draw import draw_tubules

            # 2a. Calculate cilia structure points
            # NOTE: num_points is NOT passed here. It's calculated inside generate_cilia_structure 
            # to ensure the correct MAX_INTERVAL (10.0 A) is used.
            structure = generate_cilia_structure(
                length=length,
                centerline_type=centerline_type,
                curve_radius=curve_radius,
                # curve_angle removed, calculated from length/radius internally
                sine_frequency=sine_frequency,
                sine_amplitude=sine_amplitude,
                num_doublets=num_doublets,
                cilia_radius=cilia_radius
            )
            
            # 2b. Draw central pair if requested
            if draw_central_pair:
                draw_tubules(
                    session=self.session,
                    length=length, # Pass length for truncation logic in draw_tubules
                    interval=MAX_INTERVAL, # Pass interval for truncation logic
                    centerline_points=structure['centerline'],
                    angle=0,
                    radii=[cp_radius, cp_radius], 
                    shift_distances=[cp_shift, -cp_shift], 
                    length_diffs=[0.0, 0.0], 
                    tubule_names=["C1", "C2"],
                    colors=[(100, 255, 255, 255), (100, 255, 255, 255)], # Default CP color
                    group_name="central_pair"
                )
            
            # 2c. Draw each doublet microtubule
            for doublet_info in structure['doublets']:
                # Get the shifted centerline for this doublet
                doublet_centerline = get_doublet_centerline(
                    structure['centerline'],
                    doublet_info['angle'],
                    doublet_info['shift_distance']
                )
                
                # Draw the doublet (A and B tubules)
                draw_tubules(
                    session=self.session,
                    length=length, # Pass length for truncation logic in draw_tubules
                    interval=MAX_INTERVAL, # Pass interval for truncation logic
                    centerline_points=doublet_centerline,
                    angle=doublet_info['angle'] + 90,  # Orientation of A-B pair
                    radii=[doublet_a_radius, doublet_b_radius],
                    shift_distances=[-doublet_shift, doublet_shift],
                    length_diffs=[0.0, -doublet_length_diff], # A-tubule is nominal length, B-tubule is shorter
                    tubule_names=[f"MT{doublet_info['index']+1}_A", 
                                 f"MT{doublet_info['index']+1}_B"],
                    colors=[doublet_color, doublet_color],
                    group_name=doublet_info['name']
                )


            # Update status
            self.status_label.setText(f"✓ Model generated successfully! {num_doublets} Doublets, Type: {centerline_type}")
            self.status_label.setStyleSheet("color: green; font-weight: bold;")
            
            # Log to session
            self.session.logger.info(f"Generated cilia model (Type: {centerline_type}, Length: {length} Å)")

        except ValueError as e:
            error_msg = f"Invalid input: {str(e)}"
            self.status_label.setText(f"✗ {error_msg}")
            self.status_label.setStyleSheet("color: red; font-weight: bold;")
            self.session.logger.error(error_msg)
        except Exception as e:
            error_msg = f"Error generating model: {type(e).__name__}: {str(e)}"
            self.status_label.setText(f"✗ {error_msg}")
            self.status_label.setStyleSheet("color: red; font-weight: bold;")
            self.session.logger.error(error_msg)