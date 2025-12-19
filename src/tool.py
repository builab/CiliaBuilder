# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.core.models import Model

from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QLabel, QLineEdit, QPushButton,
                          QComboBox, QCheckBox, QGroupBox, QSpacerItem, QSizePolicy, QFileDialog) 
from Qt.QtCore import Qt

# Import the command functions from cmd.py
from .cmd import ciliabuild, centriolebuild 

# Import default_config
from . import default_config

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
        session.logger.info("Default setting:\n\tset bgColor white\n\tlighting depthCue false\n\tlighting soft\n\tgraphics silhouettes true\n\tlighting shadows true\nsurface cap true")
        run(session, "set bgColor white", log=False)
        run(session, "lighting depthCue false", log=False)
        run(session, "lighting soft", log=False) 
        run(session, "graphics silhouettes true", log=False)
        run(session, "lighting shadows true", log=False)
        run(session, "surface cap true", log=False)
        run(session, "sop cap on", log=False)
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

        # Row 0: Total Length (Starts with CILIA_LENGTH)
        general_layout.addWidget(QLabel("Total Length (Å):"), 0, 0)
        self.length_input = QLineEdit(str(default_config.CILIA_LENGTH)) # Use CILIA default
        general_layout.addWidget(self.length_input, 0, 1)

        # Row 1: Structure Radius (Dynamically labeled, starts with CILIA_RADIUS)
        self.radius_label = QLabel("Doublet Ring Radius (Å):")
        general_layout.addWidget(self.radius_label, 1, 0)
        self.cilia_radius_input = QLineEdit(str(default_config.CILIA_RADIUS)) # Use CILIA default
        general_layout.addWidget(self.cilia_radius_input, 1, 1)
        
        # Row 2: Doublets/Triplets Count (Dynamically labeled, starts with CILIA_NUM_DOUBLETS)
        self.count_label = QLabel("Number of Doublets:")
        general_layout.addWidget(self.count_label, 2, 0)
        self.num_units_input = QLineEdit(str(default_config.CILIA_NUM_DOUBLETS)) # Use CILIA default
        general_layout.addWidget(self.num_units_input, 2, 1)
        
        # Row 3: Centerline Type (Dropdown, starts with CILIA_LINE)
        general_layout.addWidget(QLabel("Line Type:"), 3, 0)
        self.line_type_combo = QComboBox()
        self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal', 'template']) 
        self.line_type_combo.setCurrentText(default_config.CILIA_LINE) # Use CILIA default
        self.line_type_combo.currentIndexChanged.connect(self._toggle_centerline_inputs)
        general_layout.addWidget(self.line_type_combo, 3, 1)
        
        # --- Row 4: Template File Input (New) ---
        self.template_file_h_layout = QHBoxLayout()
        
        self.template_file_input = QLineEdit("") # Empty default
        self.template_file_h_layout.addWidget(self.template_file_input)
        
        self.browse_button = QPushButton("Browse...")
        self.browse_button.clicked.connect(self._browse_template_file)
        self.template_file_h_layout.addWidget(self.browse_button)
        
        self.template_file_group = QWidget()
        self.template_file_group.setLayout(self.template_file_h_layout)
        # Position in the same row as curve input, will be hidden by default
        general_layout.addWidget(self.template_file_group, 4, 0, 1, 2)
        
        # Row 5: Curve Input (Starts with CILIA_CURVE_RADIUS)
        self.curve_radius_label = QLabel("Curve Radius (Å):")
        general_layout.addWidget(self.curve_radius_label, 5, 0)
        self.curve_radius_input = QLineEdit(str(default_config.CILIA_CURVE_RADIUS)) # Use CILIA default
        general_layout.addWidget(self.curve_radius_input, 5, 1)

        # Row 6: Sinusoidal Inputs (Combined, starts with CILIA SINE defaults)
        # --- FIX START: Wrap Sine controls in a QWidget container ---
        self.sine_controls_widget = QWidget()
        sine_h_layout = QHBoxLayout(self.sine_controls_widget) 
        sine_h_layout.setContentsMargins(0, 0, 0, 0) # Remove extra margin from the inner layout
        
        self.sine_freq_label = QLabel("Sine Frequency:")
        sine_h_layout.addWidget(self.sine_freq_label)
        self.sine_frequency_input = QLineEdit(str(default_config.CILIA_SINE_FREQUENCY)) 
        sine_h_layout.addWidget(self.sine_frequency_input)
        
        self.sine_amp_label = QLabel("Amplitude (Å):")
        sine_h_layout.addWidget(self.sine_amp_label)
        self.sine_amplitude_input = QLineEdit(str(default_config.CILIA_SINE_AMPLITUDE)) 
        sine_h_layout.addWidget(self.sine_amplitude_input)
        
        general_layout.addWidget(self.sine_controls_widget, 6, 0, 1, 2) # Add the QWidget container to the grid
        # --- FIX END ---

        general_group.setLayout(general_layout)
        main_layout.addWidget(general_group)
        self._toggle_centerline_inputs() # Initialize centerline state

        # --- 2. Cilia-Specific Group ---
        self.cilia_group = QGroupBox("Cilia-Specific Parameters (9x2 + 2)")
        cilia_layout = QGridLayout()
        
        # Helper to convert RGBA tuple to comma-separated string
        def color_to_string(color_tuple):
            return ",".join(map(str, color_tuple))

        # Row 0: Cilia Tubule Colors (A, B, CP on same line)
        cilia_color_h_layout = QHBoxLayout()
        cilia_color_h_layout.addWidget(QLabel("A-tubule:"))
        self.cilia_a_color_input = QLineEdit(color_to_string(default_config.CILIA_DOUBLET_A_COLOR))
        self.cilia_a_color_input.setMaximumWidth(120)
        cilia_color_h_layout.addWidget(self.cilia_a_color_input)
        
        cilia_color_h_layout.addWidget(QLabel("B-tubule:"))
        self.cilia_b_color_input = QLineEdit(color_to_string(default_config.CILIA_DOUBLET_B_COLOR))
        self.cilia_b_color_input.setMaximumWidth(120)
        cilia_color_h_layout.addWidget(self.cilia_b_color_input)
        
        cilia_color_h_layout.addWidget(QLabel("C1/C2:"))
        self.cilia_cp_color_input = QLineEdit(color_to_string(default_config.CILIA_CP_COLOR))
        self.cilia_cp_color_input.setMaximumWidth(120)
        cilia_color_h_layout.addWidget(self.cilia_cp_color_input)
        
        cilia_layout.addLayout(cilia_color_h_layout, 0, 0, 1, 2)
        
        # Row 1: Draw Central Pair and Draw Membrane (Combined)
        draw_h_layout = QHBoxLayout()
        self.draw_cp_check = QCheckBox("Draw Central Pair (C1/C2)")
        self.draw_cp_check.setChecked(default_config.CILIA_DRAW_CENTRAL_PAIR)
        draw_h_layout.addWidget(self.draw_cp_check)
        
        self.draw_membrane_check = QCheckBox("Draw Membrane")
        self.draw_membrane_check.setChecked(default_config.CILIA_MEMBRANE)
        draw_h_layout.addWidget(self.draw_membrane_check)
        draw_h_layout.addStretch()
        cilia_layout.addLayout(draw_h_layout, 1, 0, 1, 2)
        
        # Row 2: Length Differences (A-B and CP-Doublet on same line)
        length_diff_h_layout = QHBoxLayout()
        length_diff_h_layout.addWidget(QLabel("A-B Diff (Å):"))
        self.cilia_doublet_length_diff_input = QLineEdit(str(default_config.CILIA_DOUBLET_LENGTH_DIFF))
        self.cilia_doublet_length_diff_input.setToolTip("A-tubule minus B-tubule length")
        self.cilia_doublet_length_diff_input.setMaximumWidth(100)
        length_diff_h_layout.addWidget(self.cilia_doublet_length_diff_input)

        length_diff_h_layout.addWidget(QLabel("CP-Doublet Diff (Å):"))
        self.cilia_cp_doublet_length_diff_input = QLineEdit(str(default_config.CILIA_CP_DOUBLET_LENGTH_DIFF))
        self.cilia_cp_doublet_length_diff_input.setToolTip("CP minus doublet A-tubule length (positive = CP longer)")
        self.cilia_cp_doublet_length_diff_input.setMaximumWidth(100)
        length_diff_h_layout.addWidget(self.cilia_cp_doublet_length_diff_input)

        cilia_layout.addLayout(length_diff_h_layout, 2, 0, 1, 2)
        
        # Row 3: Membrane Radius and Fraction (Combined)
        membrane_param_h_layout = QHBoxLayout()
        membrane_param_h_layout.addWidget(QLabel("Membrane Radius (Å):"))
        self.membrane_radius_input = QLineEdit(str(default_config.CILIA_MEMBRANE_RADIUS))
        membrane_param_h_layout.addWidget(self.membrane_radius_input)

        membrane_param_h_layout.addWidget(QLabel("Fraction (0-1):"))
        self.membrane_fraction_input = QLineEdit(str(default_config.CILIA_MEMBRANE_FRACTION))
        membrane_param_h_layout.addWidget(self.membrane_fraction_input)
        
        cilia_layout.addLayout(membrane_param_h_layout, 3, 0, 1, 2)
        
        self.cilia_group.setLayout(cilia_layout)
        main_layout.addWidget(self.cilia_group)


        # --- 3. Centriole-Specific Group ---
        self.centriole_group = QGroupBox("Centriole-Specific Parameters (9x3)")
        centriole_layout = QGridLayout()
        
        # Row 0: Centriole Tubule Colors (A, B, C on same line)
        centriole_color_h_layout = QHBoxLayout()
        centriole_color_h_layout.addWidget(QLabel("A-tubule:"))
        self.centriole_a_color_input = QLineEdit(color_to_string(default_config.CENTRIOLE_TRIPLET_A_COLOR))
        self.centriole_a_color_input.setMaximumWidth(120)
        centriole_color_h_layout.addWidget(self.centriole_a_color_input)
        
        centriole_color_h_layout.addWidget(QLabel("B-tubule:"))
        self.centriole_b_color_input = QLineEdit(color_to_string(default_config.CENTRIOLE_TRIPLET_B_COLOR))
        self.centriole_b_color_input.setMaximumWidth(120)
        centriole_color_h_layout.addWidget(self.centriole_b_color_input)
        
        centriole_color_h_layout.addWidget(QLabel("C-tubule:"))
        self.centriole_c_color_input = QLineEdit(color_to_string(default_config.CENTRIOLE_TRIPLET_C_COLOR))
        self.centriole_c_color_input.setMaximumWidth(120)
        centriole_color_h_layout.addWidget(self.centriole_c_color_input)
        
        centriole_layout.addLayout(centriole_color_h_layout, 0, 0, 1, 2)
        
        # Row 1: Centriole Angle Offset (Default 60.0)
        centriole_layout.addWidget(QLabel("Triplet Angle Offset (°):"), 1, 0)
        self.centriole_angle_offset_input = QLineEdit(str(default_config.CENTRIOLE_OFFSET_ANGLE)) 
        centriole_layout.addWidget(self.centriole_angle_offset_input, 1, 1)

        # Row 2: Triplet B-Length Difference (Default 0.1)
        centriole_layout.addWidget(QLabel("B Length Diff (Å):"), 2, 0)
        self.centriole_b_length_diff_input = QLineEdit(str(default_config.CENTRIOLE_TRIPLET_B_LENGTH_DIFF)) # Changed var name
        centriole_layout.addWidget(self.centriole_b_length_diff_input, 2, 1)
        
        # Row 3: Triplet C-Length Difference (Default 0.2)
        centriole_layout.addWidget(QLabel("C Length Diff (Å):"), 3, 0)
        self.centriole_c_length_diff_input = QLineEdit(str(default_config.CENTRIOLE_TRIPLET_C_LENGTH_DIFF)) # Changed var name
        centriole_layout.addWidget(self.centriole_c_length_diff_input, 3, 1)
        
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
        
    def _browse_template_file(self):
        """Open a file dialog to select the template centerline file"""
        # FIX: Use self.tool_window.ui_area, which is the QWidget containing the tool's UI,
        # ensuring it satisfies the QFileDialog constructor's type requirement.
        file_dialog = QFileDialog(self.tool_window.ui_area)
        
        # Suggest opening common file types for coordinates (e.g., PDB, CSV, TXT)
        file_dialog.setNameFilter("Centerline Files (*.pdb *.csv *.txt);;All Files (*)")
        file_dialog.setWindowTitle("Select Centerline Template File")
        
        if file_dialog.exec():
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                self.template_file_input.setText(selected_files[0])


    def _update_ui_visibility(self):
        """Update labels and hide/show groups based on Cilia/Centriole mode"""
        mode = self.mode_combo.currentText()
        is_cilia = mode.startswith('Cilia')
        
        # Update General Labels and Defaults
        if is_cilia:
            self.radius_label.setText("Cilia Radius (Å):")
            self.count_label.setText("Number of Doublets:")
            
            # Update general settings to Cilia defaults
            self.length_input.setText(str(default_config.CILIA_LENGTH))
            self.cilia_radius_input.setText(str(default_config.CILIA_RADIUS))
            self.num_units_input.setText(str(default_config.CILIA_NUM_DOUBLETS))
            self.line_type_combo.setCurrentText(default_config.CILIA_LINE)
            self.curve_radius_input.setText(str(default_config.CILIA_CURVE_RADIUS))
            self.sine_frequency_input.setText(str(default_config.CILIA_SINE_FREQUENCY))
            self.sine_amplitude_input.setText(str(default_config.CILIA_SINE_AMPLITUDE))
            
        else:
            self.radius_label.setText("Centriole Radius (Å):")
            self.count_label.setText("Number of Triplets:")
            
            # Update general settings to Centriole defaults
            self.length_input.setText(str(default_config.CENTRIOLE_LENGTH))
            self.cilia_radius_input.setText(str(default_config.CENTRIOLE_RADIUS))
            self.num_units_input.setText(str(default_config.CENTRIOLE_NUM_TRIPLETS))
            self.line_type_combo.setCurrentText(default_config.CENTRIOLE_LINE)
            self.curve_radius_input.setText(str(default_config.CENTRIOLE_CURVE_RADIUS))
            self.sine_frequency_input.setText(str(default_config.CENTRIOLE_SINE_FREQUENCY))
            self.sine_amplitude_input.setText(str(default_config.CENTRIOLE_SINE_AMPLITUDE))
        
        # Toggle Specific Groups
        self.cilia_group.setVisible(is_cilia)
        self.centriole_group.setVisible(not is_cilia)
        
        # Ensure centerline inputs reflect current centerline type
        self._toggle_centerline_inputs()


    def _toggle_centerline_inputs(self):
        """Enable/disable curve, sine, or template inputs based on line type dropdown"""
        line_type = self.line_type_combo.currentText()
        
        is_curve = line_type == 'curve'
        is_sine = line_type == 'sinusoidal'
        is_template = line_type == 'template'
        
        # 1. Curve Inputs (Only visible for 'curve')
        self.curve_radius_label.setVisible(is_curve)
        self.curve_radius_input.setVisible(is_curve)
        self.curve_radius_label.setEnabled(is_curve)
        self.curve_radius_input.setEnabled(is_curve)
        
        # 2. Sine Inputs (Only visible for 'sinusoidal') - Target the QWidget container
        self.sine_controls_widget.setVisible(is_sine)
        # We set the overall widget visibility, so individual enable/disable is not strictly needed but included for thoroughness if needed later
        # for inputs inside the widget, e.g. self.sine_frequency_input.setEnabled(is_sine)
        
        # 3. Template Inputs (Only visible for 'template')
        self.template_file_group.setVisible(is_template)
        # Note: 'straight' line type will hide all three groups above (curve, sine, template).


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
            
            template_file = self.template_file_input.text() # Get the new template file path
            
            if centerline_type == 'template' and not template_file:
                raise ValueError("Template file path is required when line type is 'template'.")

            if length <= 0 or num_units < 0 or ring_radius <= 0:
                 raise ValueError("Length, number of units, and ring radius must be positive.")

            # Update status
            self.status_label.setText("Generating model...")
            self.status_label.setStyleSheet("color: blue; font-style: italic;")

            is_cilia = self.mode_combo.currentText().startswith('Cilia')
            new_model = None
            
            # NOTE: We assume the target cmd.py functions (ciliabuild and centriolebuild) 
            # have been updated to accept a 'template_file' keyword argument.
            
            if is_cilia:
                # --- CILIA (9x2 + 2) Logic - Call ciliabuild ---
                
                # Parse cilia colors
                cilia_a_color = self._parse_color(self.cilia_a_color_input.text())
                cilia_b_color = self._parse_color(self.cilia_b_color_input.text())
                cilia_cp_color = self._parse_color(self.cilia_cp_color_input.text())
                
                self.session.logger.info(f"Cilia colors - A: {cilia_a_color}, B: {cilia_b_color}, CP: {cilia_cp_color}")
                
                # Cilia-specific inputs
                draw_central_pair = self.draw_cp_check.isChecked()
                should_draw_membrane = self.draw_membrane_check.isChecked() 
                doublet_length_diff = float(self.cilia_doublet_length_diff_input.text())
                cp_doublet_length_diff = float(self.cilia_cp_doublet_length_diff_input.text())
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
                    template_file=template_file, # PASSED NEW ARGUMENT
                    num_doublets=num_units,
                    cilia_radius=ring_radius,
                    draw_central_pair=draw_central_pair,
                    membrane=should_draw_membrane,
                    membrane_radius=membrane_radius,
                    membrane_fraction=membrane_fraction,
                    doublet_length_diff=doublet_length_diff,
                    cp_doublet_length_diff=cp_doublet_length_diff,
                    doublet_a_color=cilia_a_color,
                    doublet_b_color=cilia_b_color,
                    cp_color=cilia_cp_color
                )
                
            else:
                # --- CENTRIOLE (9x3) Logic - Call centriolebuild ---
                
                # Parse centriole colors
                centriole_a_color = self._parse_color(self.centriole_a_color_input.text())
                centriole_b_color = self._parse_color(self.centriole_b_color_input.text())
                centriole_c_color = self._parse_color(self.centriole_c_color_input.text())
                
                self.session.logger.info(f"Centriole colors - A: {centriole_a_color}, B: {centriole_b_color}, C: {centriole_c_color}")
                
                # Centriole-specific inputs
                angle_offset = default_config.CILIA_OFFSET_ANGLE - float(self.centriole_angle_offset_input.text())
                b_length_diff = float(self.centriole_b_length_diff_input.text())
                c_length_diff = float(self.centriole_c_length_diff_input.text())
                
                # Call the command function and get the returned model
                new_model = centriolebuild(
                    session=self.session,
                    length=length,
                    line=centerline_type,
                    curve_radius=curve_radius,
                    sine_frequency=sine_frequency,
                    sine_amplitude=sine_amplitude,
                    template_file=template_file, # PASSED NEW ARGUMENT
                    num_triplets=num_units,
                    centriole_radius=ring_radius,
                    centriole_angle_offset=angle_offset,
                    triplet_b_length_diff=b_length_diff,
                    triplet_c_length_diff=c_length_diff,
                    triplet_a_color=centriole_a_color,
                    triplet_b_color=centriole_b_color,
                    triplet_c_color=centriole_c_color
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