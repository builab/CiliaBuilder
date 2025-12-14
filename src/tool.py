# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QLabel, QLineEdit, QPushButton,
                          QComboBox, QCheckBox, QGroupBox, QSpacerItem, QSizePolicy)
from Qt.QtCore import Qt

# Define the optimal sampling interval for smoothness (10 Angstroms)
MAX_INTERVAL = 10.0 

# Fixed Microtubule Dimensions (based on previous request)
# Cilia Doublet
D_A_RADIUS = 125.0
D_B_RADIUS = 145.0
D_SHIFT = 70.0
CP_RADIUS = 125.0
CP_SHIFT = 160.0
# Centriole Triplet
T_A_RADIUS = 125.0
T_B_RADIUS = 135.0
T_C_RADIUS = 135.0
T_SHIFTS = [140.0, 0.0, -160.0] # A, B, C distance from triplet centerline

class CiliaSim(ToolInstance):
    """
    CiliaSim tool for generating cilia and centriole microtubule models
    """

    SESSION_ENDURING = False
    SESSION_SAVE = False

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)

        # Set name displayed on title bar
        self.display_name = "Cilia/Centriole Sim"
        
        self.tool_window = MainToolWindow(self)
        
        # Build the user interface
        self._build_ui()


    def _build_ui(self):
        """Build the user interface with mode selector and dynamic visibility"""
        
        # Need to import QSpacerItem and QSizePolicy if not already imported
        # The import list above already contains QSpacerItem, QSizePolicy 
        
        main_layout = QVBoxLayout()
        main_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # NOTE: Removed the main Title Label as requested

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
        self.length_input = QLineEdit("5000") # Default for Cilia
        general_layout.addWidget(self.length_input, 0, 1)

        # Row 1: Structure Radius (Dynamically labeled)
        self.radius_label = QLabel("Doublet Ring Radius (Å):") # Cilia Default
        general_layout.addWidget(self.radius_label, 1, 0)
        self.cilia_radius_input = QLineEdit("875.0")
        general_layout.addWidget(self.cilia_radius_input, 1, 1)
        
        # Row 2: Doublets/Triplets Count (Dynamically labeled)
        # NOTE: Removed (9-0) from label as requested
        self.count_label = QLabel("Number of Doublets:") # Cilia Default
        general_layout.addWidget(self.count_label, 2, 0)
        self.num_units_input = QLineEdit("9")
        general_layout.addWidget(self.num_units_input, 2, 1)

        # Row 3: Doublet/Triplet Color
        general_layout.addWidget(QLabel("MT Color (R,G,B,A):"), 3, 0)
        self.color_input = QLineEdit("100,100,255,255")
        general_layout.addWidget(self.color_input, 3, 1)
        
        # Centerline Type (Dropdown)
        general_layout.addWidget(QLabel("Centerline Type:"), 4, 0)
        self.line_type_combo = QComboBox()
        self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal'])
        self.line_type_combo.setCurrentText('straight')
        self.line_type_combo.currentIndexChanged.connect(self._toggle_centerline_inputs)
        general_layout.addWidget(self.line_type_combo, 4, 1)
        
        # Curve Input 
        self.curve_radius_label = QLabel("Curve Radius (Å):")
        general_layout.addWidget(self.curve_radius_label, 5, 0)
        self.curve_radius_input = QLineEdit("10000.0")
        general_layout.addWidget(self.curve_radius_input, 5, 1)

        # Sinusoidal Inputs
        self.sine_freq_label = QLabel("Sine Frequency:")
        general_layout.addWidget(self.sine_freq_label, 6, 0)
        self.sine_frequency_input = QLineEdit("2.0")
        general_layout.addWidget(self.sine_frequency_input, 6, 1)
        
        self.sine_amp_label = QLabel("Sine Amplitude (Å):")
        general_layout.addWidget(self.sine_amp_label, 7, 0)
        self.sine_amplitude_input = QLineEdit("500.0")
        general_layout.addWidget(self.sine_amplitude_input, 7, 1)

        general_group.setLayout(general_layout)
        main_layout.addWidget(general_group)
        self._toggle_centerline_inputs() # Initialize centerline state

        # --- 2. Cilia-Specific Group ---
        self.cilia_group = QGroupBox("Cilia-Specific Parameters (9x2 + 2)")
        cilia_layout = QGridLayout()
        
        # Fixed Dimensions Summary
        fixed_cilia_label = QLabel(f"Fixed Dims: A R={D_A_RADIUS}Å, B R={D_B_RADIUS}Å, D Shift={D_SHIFT}Å, CP Shift={CP_SHIFT}Å")
        fixed_cilia_label.setStyleSheet("font-size: 8pt; color: gray;")
        cilia_layout.addWidget(fixed_cilia_label, 0, 0, 1, 2)
        
        # Draw Central Pair Checkbox
        self.draw_cp_check = QCheckBox("Draw Central Pair (C1/C2)")
        self.draw_cp_check.setChecked(True)
        cilia_layout.addWidget(self.draw_cp_check, 1, 0, 1, 2)
        
        # Doublet Length Difference
        cilia_layout.addWidget(QLabel("A-B Length Diff (Å):"), 2, 0)
        self.cilia_doublet_length_diff_input = QLineEdit("250.0")
        self.cilia_doublet_length_diff_input.setToolTip("Length of A-tubule - Length of B-tubule")
        cilia_layout.addWidget(self.cilia_doublet_length_diff_input, 2, 1)
        
        self.cilia_group.setLayout(cilia_layout)
        main_layout.addWidget(self.cilia_group)


        # --- 3. Centriole-Specific Group ---
        self.centriole_group = QGroupBox("Centriole-Specific Parameters (9x3)")
        centriole_layout = QGridLayout()
        
        # Fixed Dimensions Summary
        fixed_centriole_label = QLabel(f"Fixed Triplet R={T_A_RADIUS}Å, Shifts={T_SHIFTS[0]}, {T_SHIFTS[1]}, {T_SHIFTS[2]}Å")
        fixed_centriole_label.setStyleSheet("font-size: 8pt; color: gray;")
        centriole_layout.addWidget(fixed_centriole_label, 0, 0, 1, 2)

        # Centriole Angle Offset (Default 60.0)
        centriole_layout.addWidget(QLabel("Triplet Angle Offset (°):"), 1, 0)
        self.centriole_angle_offset_input = QLineEdit("60.0") # NEW DEFAULT
        self.centriole_angle_offset_input.setToolTip("Angular offset for the triplet unit relative to the tangent frame.")
        centriole_layout.addWidget(self.centriole_angle_offset_input, 1, 1)

        # Triplet A-B Length Difference (Default 5.0)
        centriole_layout.addWidget(QLabel("A-B Length Diff (Å):"), 2, 0)
        self.centriole_ab_length_diff_input = QLineEdit("5.0") # NEW DEFAULT
        self.centriole_ab_length_diff_input.setToolTip("Length of A-tubule - Length of B-tubule")
        centriole_layout.addWidget(self.centriole_ab_length_diff_input, 2, 1)
        
        # Triplet B-C Length Difference (Default 300.0)
        centriole_layout.addWidget(QLabel("B-C Length Diff (Å):"), 3, 0)
        self.centriole_bc_length_diff_input = QLineEdit("300.0") # NEW DEFAULT
        self.centriole_bc_length_diff_input.setToolTip("Length of B-tubule - Length of C-tubule")
        centriole_layout.addWidget(self.centriole_bc_length_diff_input, 3, 1)
        
        self.centriole_group.setLayout(centriole_layout)
        main_layout.addWidget(self.centriole_group)
        
        # Initial visibility update
        self._update_ui_visibility()


        # Generate button
        generate_button = QPushButton("Generate Model")
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
        """Generate the cilia or centriole model with current parameters"""
        try:
            # --- 1. Get Common Parameters ---
            
            length = float(self.length_input.text())
            num_units = int(self.num_units_input.text())
            ring_radius = float(self.cilia_radius_input.text())
            mt_color = self._parse_color(self.color_input.text())
            
            centerline_type = self.line_type_combo.currentText()
            curve_radius = float(self.curve_radius_input.text())
            sine_frequency = float(self.sine_frequency_input.text())
            sine_amplitude = float(self.sine_amplitude_input.text())

            if length <= 0 or num_units < 0 or ring_radius <= 0:
                 raise ValueError("Length, number of units, and ring radius must be positive.")

            # Update status
            self.status_label.setText("Generating model...")
            self.status_label.setStyleSheet("color: blue; font-style: italic;")

            # --- 2. Core Generation Logic ---
            
            from .curve import generate_cilia_structure, get_doublet_centerline
            from .draw import draw_tubules

            is_cilia = self.mode_combo.currentText().startswith('Cilia')
            
            # 2a. Calculate structure points (Applies to both Cilia and Centriole)
            structure = generate_cilia_structure(
                length=length,
                centerline_type=centerline_type,
                curve_radius=curve_radius,
                sine_frequency=sine_frequency,
                sine_amplitude=sine_amplitude,
                num_doublets=num_units,
                cilia_radius=ring_radius
            )
            
            centerline = structure['centerline']
            model_type_name = self.mode_combo.currentText().split('(')[0].strip()

            if is_cilia:
                # --- CILIA (9x2 + 2) Logic ---
                
                # Cilia-specific inputs
                draw_central_pair = self.draw_cp_check.isChecked()
                doublet_length_diff = float(self.cilia_doublet_length_diff_input.text())

                # Central Pair (C1, C2)
                if draw_central_pair:
                    draw_tubules(
                        session=self.session,
                        length=length, 
                        interval=MAX_INTERVAL, 
                        centerline_points=centerline,
                        angle=0,
                        radii=[CP_RADIUS, CP_RADIUS], 
                        shift_distances=[CP_SHIFT, -CP_SHIFT], 
                        length_diffs=[0.0, 0.0], 
                        tubule_names=["C1", "C2"],
                        colors=[(100, 255, 255, 255), (100, 255, 255, 255)],
                        group_name="central_pair"
                    )
                
                # Doublets (A, B)
                for doublet_info in structure['doublets']:
                    doublet_centerline = get_doublet_centerline(
                        centerline,
                        doublet_info['angle'],
                        doublet_info['shift_distance']
                    )
                    draw_tubules(
                        session=self.session,
                        length=length, 
                        interval=MAX_INTERVAL, 
                        centerline_points=doublet_centerline,
                        angle=doublet_info['angle'] + 90,  
                        radii=[D_A_RADIUS, D_B_RADIUS],
                        shift_distances=[-D_SHIFT, D_SHIFT],
                        length_diffs=[0.0, -doublet_length_diff], 
                        tubule_names=[f"MT{doublet_info['index']+1}_A", f"MT{doublet_info['index']+1}_B"],
                        colors=[mt_color, mt_color],
                        group_name=doublet_info['name']
                    )
                
            else:
                # --- CENTRIOLE (9x3) Logic ---
                
                # Centriole-specific inputs
                angle_offset = float(self.centriole_angle_offset_input.text())
                ab_length_diff = float(self.centriole_ab_length_diff_input.text())
                bc_length_diff = float(self.centriole_bc_length_diff_input.text())
                
                # Triplet length diffs (A is nominal length)
                length_diffs = [
                    0.0, 
                    -ab_length_diff,
                    -(ab_length_diff + bc_length_diff) 
                ]

                # Triplets (A, B, C)
                for triplet_info in structure['doublets']:
                    triplet_centerline = get_doublet_centerline(
                        centerline,
                        triplet_info['angle'],
                        triplet_info['shift_distance']
                    )
                    
                    # Centriole uses Triplet Radii/Shifts/Colors
                    draw_tubules(
                        session=self.session,
                        length=length, 
                        interval=MAX_INTERVAL, 
                        centerline_points=triplet_centerline,
                        angle=triplet_info['angle'] + 90 + angle_offset, 
                        radii=[T_A_RADIUS, T_B_RADIUS, T_C_RADIUS],
                        shift_distances=T_SHIFTS,
                        length_diffs=length_diffs, 
                        tubule_names=[f"MT{triplet_info['index']+1}_A", f"MT{triplet_info['index']+1}_B", f"MT{triplet_info['index']+1}_C"],
                        colors=[mt_color, mt_color, mt_color],
                        group_name=triplet_info['name'].replace("doublet", "triplet")
                    )

            # Update status
            self.status_label.setText(f"✓ {model_type_name} Model generated successfully! {num_units} Units, Type: {centerline_type}")
            self.status_label.setStyleSheet("color: green; font-weight: bold;")
            
            # Log to session
            self.session.logger.info(f"Generated {model_type_name} model (Type: {centerline_type}, Length: {length} Å)")

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