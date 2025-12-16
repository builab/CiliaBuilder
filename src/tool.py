# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.core.models import Model


from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QLabel, QLineEdit, QPushButton,
                          QComboBox, QCheckBox, QGroupBox, QSpacerItem, QSizePolicy)
from Qt.QtCore import Qt

# Import the command functions from cmd.py
from .cmd import ciliabuild, centriolebuild 

# Define the optimal sampling interval for smoothness (10 Angstroms)
MAX_INTERVAL = 20.0 

# Fixed Microtubule Dimensions (kept for reference, but cmd.py handles the actual values)
# Cilia Doublet
D_A_RADIUS = 125.0
D_B_RADIUS = 140.0
D_SHIFT = 70.0
CP_RADIUS = 125.0
CP_SHIFT = 160.0
# Centriole Triplet
T_A_RADIUS = 125.0
T_B_RADIUS = 135.0
T_C_RADIUS = 135.0
T_SHIFTS = [140.0, 0.0, -160.0] 

class CiliaBuilder(ToolInstance):
    """
    CiliaBuilder tool for generating cilia and centriole microtubule models
    """

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        
        # --- NEW: Track the last generated model ---
        self.last_generated_model = None

        # Set name displayed on title bar
        self.display_name = "Cilia/Centriole Builder"
        
        self.tool_window = MainToolWindow(self)
        
        # Build the user interface
        self._build_ui()
        
        # --- PLACE THE COMMAND HERE ---
        # This runs the command to set the lighting as soon as the tool is instantiated.
        session.logger.info(f"Starting Cilia Builder - Builab 2025")
        session.logger.info("Default setting:\n\tset bgColor white\n\tlighting depthCue false\n\tlighting soft\n\tgraphics silhouettes true\n\tlighting shadows true")
        run(session, "set bgColor white", log=False)
        run(session, "lighting depthCue false", log=False)
        run(session, "lighting soft", log=False) 
        run(session, "graphics silhouettes true", log=False)
        run(session, "lighting shadows true", log=False)
        # ------------------------------


    def _build_ui(self):
        """Build the user interface with mode selector and dynamic visibility"""
        
        main_layout = QVBoxLayout()
        main_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # --- Mode Selector ---
        mode_h_layout = QHBoxLayout()
        mode_h_layout.addWidget(QLabel("Structure Type:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(['Cilia (9x2 + 2)', 'Centriole (9x3)'])
        self.mode_combo.setCurrentText('Cilia (9x2 + 2)')
        self.mode_combo.currentIndexChanged.connect(self._update_ui_visibility)
        mode_h_layout.addWidget(self.mode_combo)
        main_layout.addLayout(mode_h_layout)

        # --- 1. General Settings Group (Common to both) ---
        general_group = QGroupBox("General & Centerline Settings")
        general_layout = QGridLayout()

        # Row 0: Total Length
        general_layout.addWidget(QLabel("Total Length (Å):"), 0, 0)
        self.length_input = QLineEdit("15000")
        general_layout.addWidget(self.length_input, 0, 1)

        # Row 1: Structure Radius (Dynamically labeled)
        self.radius_label = QLabel("Doublet Ring Radius (Å):")
        general_layout.addWidget(self.radius_label, 1, 0)
        self.cilia_radius_input = QLineEdit("875.0")
        general_layout.addWidget(self.cilia_radius_input, 1, 1)
        
        # Row 2: Doublets/Triplets Count (Dynamically labeled)
        self.count_label = QLabel("Number of Doublets:")
        general_layout.addWidget(self.count_label, 2, 0)
        self.num_units_input = QLineEdit("9")
        general_layout.addWidget(self.num_units_input, 2, 1)

        # Row 3: Doublet/Triplet Color (NOTE: This is not used by cmd.py's logic)
        general_layout.addWidget(QLabel("MT Color (R,G,B,A):"), 3, 0)
        self.color_input = QLineEdit("100,100,255,255")
        general_layout.addWidget(self.color_input, 3, 1)
        
        # Row 4: Centerline Type (Dropdown)
        general_layout.addWidget(QLabel("Centerline Type:"), 4, 0)
        self.line_type_combo = QComboBox()
        self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal'])
        self.line_type_combo.setCurrentText('straight')
        self.line_type_combo.currentIndexChanged.connect(self._toggle_centerline_inputs)
        general_layout.addWidget(self.line_type_combo, 4, 1)
        
        # Row 5: Curve Input 
        self.curve_radius_label = QLabel("Curve Radius (Å):")
        general_layout.addWidget(self.curve_radius_label, 5, 0)
        self.curve_radius_input = QLineEdit("10000.0")
        general_layout.addWidget(self.curve_radius_input, 5, 1)

        # Row 6: Sinusoidal Inputs (Combined)
        sine_h_layout = QHBoxLayout()
        
        self.sine_freq_label = QLabel("Sine Frequency:")
        sine_h_layout.addWidget(self.sine_freq_label)
        self.sine_frequency_input = QLineEdit("2.0")
        sine_h_layout.addWidget(self.sine_frequency_input)
        
        self.sine_amp_label = QLabel("Amplitude (Å):")
        sine_h_layout.addWidget(self.sine_amp_label)
        self.sine_amplitude_input = QLineEdit("2000.0")
        sine_h_layout.addWidget(self.sine_amplitude_input)
        
        general_layout.addLayout(sine_h_layout, 6, 0, 1, 2) # Span 2 columns

        general_group.setLayout(general_layout)
        main_layout.addWidget(general_group)
        self._toggle_centerline_inputs() # Initialize centerline state

        # --- 2. Cilia-Specific Group ---
        self.cilia_group = QGroupBox("Cilia-Specific Parameters (9x2 + 2)")
        cilia_layout = QGridLayout()
        
        # Row 1: Draw Central Pair and Draw Membrane (Combined)
        draw_h_layout = QHBoxLayout()
        self.draw_cp_check = QCheckBox("Draw Central Pair (C1/C2)")
        self.draw_cp_check.setChecked(True)
        draw_h_layout.addWidget(self.draw_cp_check)
        
        self.draw_membrane_check = QCheckBox("Draw Membrane")
        self.draw_membrane_check.setChecked(True)
        draw_h_layout.addWidget(self.draw_membrane_check)
        draw_h_layout.addStretch()
        cilia_layout.addLayout(draw_h_layout, 1, 0, 1, 2)
        
        # Row 2: Doublet Length Difference
        cilia_layout.addWidget(QLabel("A-B Length Diff (Å):"), 2, 0)
        self.cilia_doublet_length_diff_input = QLineEdit("1.0")
        self.cilia_doublet_length_diff_input.setToolTip("Length of A-tubule - Length of B-tubule")
        cilia_layout.addWidget(self.cilia_doublet_length_diff_input, 2, 1)
        
        # Row 3: Membrane Radius and Fraction (Combined)
        membrane_param_h_layout = QHBoxLayout()
        membrane_param_h_layout.addWidget(QLabel("Membrane Radius (Å):"))
        self.membrane_radius_input = QLineEdit("1100.0")
        membrane_param_h_layout.addWidget(self.membrane_radius_input)

        membrane_param_h_layout.addWidget(QLabel("Fraction (0-1):"))
        self.membrane_fraction_input = QLineEdit("0.5")
        membrane_param_h_layout.addWidget(self.membrane_fraction_input)
        
        cilia_layout.addLayout(membrane_param_h_layout, 3, 0, 1, 2)
        
        self.cilia_group.setLayout(cilia_layout)
        main_layout.addWidget(self.cilia_group)


        # --- 3. Centriole-Specific Group ---
        self.centriole_group = QGroupBox("Centriole-Specific Parameters (9x3)")
        centriole_layout = QGridLayout()
        
        # Row 1: Centriole Angle Offset (Default 60.0)
        centriole_layout.addWidget(QLabel("Triplet Angle Offset (°):"), 1, 0)
        self.centriole_angle_offset_input = QLineEdit("60.0") 
        centriole_layout.addWidget(self.centriole_angle_offset_input, 1, 1)

        # Row 2: Triplet A-B Length Difference (Default 1.0)
        centriole_layout.addWidget(QLabel("A-B Length Diff (Å):"), 2, 0)
        self.centriole_ab_length_diff_input = QLineEdit("1.0") 
        centriole_layout.addWidget(self.centriole_ab_length_diff_input, 2, 1)
        
        # Row 3: Triplet B-C Length Difference (Default 300.0)
        centriole_layout.addWidget(QLabel("B-C Length Diff (Å):"), 3, 0)
        self.centriole_bc_length_diff_input = QLineEdit("300.0") 
        centriole_layout.addWidget(self.centriole_bc_length_diff_input, 3, 1)
        
        self.centriole_group.setLayout(centriole_layout)
        main_layout.addWidget(self.centriole_group)
        
        # Initial visibility update
        self._update_ui_visibility()

        # --- Control Button Group ---
        control_h_layout = QHBoxLayout()
        
        # NEW: Close Old Model Checkbox
        self.close_old_model_check = QCheckBox("Close old model")
        self.close_old_model_check.setChecked(True)
        control_h_layout.addWidget(self.close_old_model_check)
        
        # Generate button
        generate_button = QPushButton("Generate Model")
        generate_button.setStyleSheet("font-size: 12pt; padding: 10px; margin-top: 10px;")
        generate_button.clicked.connect(self._generate_model)
        control_h_layout.addWidget(generate_button)
        
        main_layout.addLayout(control_h_layout)

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

    def _update_ui_visibility(self):
        """Update labels and hide/show groups based on Cilia/Centriole mode"""
        mode = self.mode_combo.currentText()
        is_cilia = mode.startswith('Cilia')
        
        # Update General Labels and Defaults
        if is_cilia:
            self.radius_label.setText("Doublet Ring Radius (Å):")
            self.count_label.setText("Number of Doublets:")
            self.length_input.setText("15000")
            self.num_units_input.setText("9")
        else:
            self.radius_label.setText("Triplet Ring Radius (Å):")
            self.count_label.setText("Number of Triplets:")
            self.length_input.setText("5000")
            self.num_units_input.setText("9")
        
        # Toggle Specific Groups
        self.cilia_group.setVisible(is_cilia)
        self.centriole_group.setVisible(not is_cilia)
        
        # Ensure centerline inputs reflect current centerline type
        self._toggle_centerline_inputs()


    def _toggle_centerline_inputs(self):
        """Enable/disable curve or sine inputs based on line type dropdown"""
        line_type = self.line_type_combo.currentText()
        
        is_curve = line_type == 'curve'
        self.curve_radius_input.setEnabled(is_curve)
        self.curve_radius_label.setEnabled(is_curve)
        
        is_sine = line_type == 'sinusoidal'
        self.sine_frequency_input.setEnabled(is_sine)
        self.sine_amplitude_input.setEnabled(is_sine)
        self.sine_freq_label.setEnabled(is_sine)
        self.sine_amp_label.setEnabled(is_sine)


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
        """Generate the cilia or centriole model by calling cmd.py functions"""
        try:
            # --- 0. Pre-Generation Cleanup ---
            if self.close_old_model_check.isChecked() and self.last_generated_model:
                try:
                    self.session.models.remove([self.last_generated_model])
                    self.session.logger.info(f"Closed previously generated model: {self.last_generated_model.name}")
                    self.last_generated_model = None
                except Exception as close_e:
                    self.session.logger.warning(f"Could not close old model: {close_e}")
                    
            # --- 1. Get Common Parameters ---            
            length = float(self.length_input.text())
            num_units = int(self.num_units_input.text())
            ring_radius = float(self.cilia_radius_input.text())
            
            centerline_type = self.line_type_combo.currentText()
            curve_radius = float(self.curve_radius_input.text())
            sine_frequency = float(self.sine_frequency_input.text())
            sine_amplitude = float(self.sine_amplitude_input.text())

            if length <= 0 or num_units < 0 or ring_radius <= 0:
                 raise ValueError("Length, number of units, and ring radius must be positive.")

            # Update status
            self.status_label.setText("Generating model...")
            self.status_label.setStyleSheet("color: blue; font-style: italic;")

            is_cilia = self.mode_combo.currentText().startswith('Cilia')
            new_model = None
            
            if is_cilia:
                # --- CILIA (9x2 + 2) Logic - Call ciliabuild ---
                
                # Cilia-specific inputs
                draw_central_pair = self.draw_cp_check.isChecked()
                should_draw_membrane = self.draw_membrane_check.isChecked() 
                doublet_length_diff = float(self.cilia_doublet_length_diff_input.text())
                membrane_radius = float(self.membrane_radius_input.text())
                membrane_fraction = float(self.membrane_fraction_input.text())
                
                # Call the command function and get the returned model
                new_model = ciliabuild(
                    session=self.session,
                    length=length, 
                    line=centerline_type,
                    curve_radius=curve_radius,
                    sine_frequency=sine_frequency,
                    sine_amplitude=sine_amplitude,
                    num_doublets=num_units,
                    cilia_radius=ring_radius,
                    draw_central_pair=draw_central_pair,
                    membrane=should_draw_membrane,
                    membrane_radius=membrane_radius,
                    membrane_fraction=membrane_fraction,
                    doublet_length_diff=doublet_length_diff,
                )
                
            else:
                # --- CENTRIOLE (9x3) Logic - Call centriolebuild ---
                
                # Centriole-specific inputs
                angle_offset = float(self.centriole_angle_offset_input.text())
                ab_length_diff = float(self.centriole_ab_length_diff_input.text())
                bc_length_diff = float(self.centriole_bc_length_diff_input.text())
                
                # Call the command function and get the returned model
                new_model = centriolebuild(
                    session=self.session,
                    length=length,
                    line=centerline_type,
                    curve_radius=curve_radius,
                    sine_frequency=sine_frequency,
                    sine_amplitude=sine_amplitude,
                    num_triplets=num_units,
                    centriole_radius=ring_radius,
                    centriole_angle_offset=angle_offset,
                    triplet_b_length_diff=ab_length_diff,
                    triplet_c_length_diff=bc_length_diff
                )
            
            # --- 2. Post-Generation Tracking ---
            if new_model is not None:
                self.last_generated_model = new_model
            
            # Update status
            model_type_name = self.mode_combo.currentText().split('(')[0].strip()
            self.status_label.setText(f"✓ {model_type_name} Model generated successfully! {num_units} Units, Type: {centerline_type}")
            self.status_label.setStyleSheet("color: green; font-weight: bold;")            
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