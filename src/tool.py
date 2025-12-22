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
# ADDED ciliabuild_from_csv for the new '3Dtemplate' option
from .cmd import ciliabuild, centriolebuild, ciliabuild_from_csv 

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
        
        # --- Track the last generated model ---
        self.last_generated_model = None

        # Set name displayed on title bar
        self.display_name = "Cilia/Centriole Builder"
        
        self.tool_window = MainToolWindow(self)
        
        # Build the user interface
        self._build_ui()
        
        # Default lighting setup
        session.logger.info(f"Starting Cilia Builder - Builab 2025")
        session.logger.info("Default setting:\n\tset bgColor white\n\tlighting depthCue false\n\tlighting soft\n\tgraphics silhouettes true\n\tlighting shadows true\nsurface cap true")
        run(session, "set bgColor white", log=False)
        run(session, "lighting depthCue false", log=False)
        run(session, "lighting soft", log=False) 
        run(session, "graphics silhouettes true", log=False)
        run(session, "lighting shadows true", log=False)
        run(session, "surface cap true", log=False)
        run(session, "sop cap on", log=False)


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
        self.length_input = QLineEdit(str(default_config.CILIA_LENGTH))
        general_layout.addWidget(self.length_input, 0, 1)

        # Row 1: Structure Radius (Dynamically labeled, starts with CILIA_RADIUS)
        self.radius_label = QLabel("Doublet Ring Radius (Å):")
        general_layout.addWidget(self.radius_label, 1, 0)
        self.cilia_radius_input = QLineEdit(str(default_config.CILIA_RADIUS))
        general_layout.addWidget(self.cilia_radius_input, 1, 1)
        
        # Row 2: Doublets/Triplets Count (Dynamically labeled, starts with CILIA_NUM_DOUBLETS)
        self.count_label = QLabel("Number of Doublets:")
        general_layout.addWidget(self.count_label, 2, 0)
        self.num_units_input = QLineEdit(str(default_config.CILIA_NUM_DOUBLETS))
        general_layout.addWidget(self.num_units_input, 2, 1)
        
        # Row 3: Centerline Type (Dropdown, starts with CILIA_LINE)
        general_layout.addWidget(QLabel("Line Type:"), 3, 0)
        self.line_type_combo = QComboBox()
        # Updated Cilia line types: added '3Dtemplate'
        self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal', 'tip', 'primarycilia', '3Dtemplate', '2Dtemplate'])
        self.line_type_combo.setCurrentText(default_config.CILIA_LINE)
        self.line_type_combo.currentIndexChanged.connect(self._toggle_centerline_inputs)
        general_layout.addWidget(self.line_type_combo, 3, 1)
        
        # --- Row 4: Template File Input ---
        self.template_file_h_layout = QHBoxLayout()
        
        self.template_file_input = QLineEdit("")
        self.template_file_h_layout.addWidget(self.template_file_input)
        
        self.browse_button = QPushButton("Browse...")
        self.browse_button.clicked.connect(self._browse_template_file)
        self.template_file_h_layout.addWidget(self.browse_button)
        
        self.template_file_group = QWidget()
        self.template_file_group.setLayout(self.template_file_h_layout)
        general_layout.addWidget(self.template_file_group, 4, 0, 1, 2)
        
        # --- Row 5: Tip Length Input (NEW) ---
        self.tip_length_label = QLabel("Tip Length (Å):")
        general_layout.addWidget(self.tip_length_label, 5, 0)
        self.tip_length_input = QLineEdit(str(default_config.TIP_LENGTH))
        general_layout.addWidget(self.tip_length_input, 5, 1)
        
        # Row 6: Curve Input
        self.curve_radius_label = QLabel("Curve Radius (Å):")
        general_layout.addWidget(self.curve_radius_label, 6, 0)
        self.curve_radius_input = QLineEdit(str(default_config.CILIA_CURVE_RADIUS))
        general_layout.addWidget(self.curve_radius_input, 6, 1)

        # Row 7: Sinusoidal Inputs
        self.sine_controls_widget = QWidget()
        sine_h_layout = QHBoxLayout(self.sine_controls_widget) 
        sine_h_layout.setContentsMargins(0, 0, 0, 0)
        
        self.sine_freq_label = QLabel("Sine Frequency:")
        sine_h_layout.addWidget(self.sine_freq_label)
        self.sine_frequency_input = QLineEdit(str(default_config.CILIA_SINE_FREQUENCY)) 
        sine_h_layout.addWidget(self.sine_frequency_input)
        
        self.sine_amp_label = QLabel("Amplitude (Å):")
        sine_h_layout.addWidget(self.sine_amp_label)
        self.sine_amplitude_input = QLineEdit(str(default_config.CILIA_SINE_AMPLITUDE)) 
        sine_h_layout.addWidget(self.sine_amplitude_input)
        
        general_layout.addWidget(self.sine_controls_widget, 7, 0, 1, 2)

        general_group.setLayout(general_layout)
        main_layout.addWidget(general_group)
        self._toggle_centerline_inputs() # Initialize centerline state

        # --- 2. Cilia-Specific Group ---
        self.cilia_group = QGroupBox("Cilia-Specific Parameters (9x2 + 2)")
        cilia_layout = QGridLayout()
        
        # Helper to convert RGBA tuple to comma-separated string
        def color_to_string(color_tuple):
            return ",".join(map(str, color_tuple))

        # Row 0: Cilia Tubule Colors
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
        
        # Row 1: Draw Central Pair and Draw Membrane
        draw_h_layout = QHBoxLayout()
        self.draw_cp_check = QCheckBox("Draw Central Pair (C1/C2)")
        self.draw_cp_check.setChecked(default_config.CILIA_DRAW_CENTRAL_PAIR)
        draw_h_layout.addWidget(self.draw_cp_check)
        
        self.draw_membrane_check = QCheckBox("Draw Membrane")
        self.draw_membrane_check.setChecked(default_config.CILIA_MEMBRANE)
        draw_h_layout.addWidget(self.draw_membrane_check)
        draw_h_layout.addStretch()
        cilia_layout.addLayout(draw_h_layout, 1, 0, 1, 2)
        
        # Row 2: Length Differences
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
        
        # Row 3: Membrane Parameters
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
        
        # Row 0: Centriole Tubule Colors
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
        
        # Row 1: Centriole Angle Offset
        centriole_layout.addWidget(QLabel("Triplet Angle Offset (°):"), 1, 0)
        self.centriole_angle_offset_input = QLineEdit(str(default_config.CENTRIOLE_OFFSET_ANGLE)) 
        centriole_layout.addWidget(self.centriole_angle_offset_input, 1, 1)

        # Row 2: Triplet B-Length Difference
        centriole_layout.addWidget(QLabel("B Length Diff (Å):"), 2, 0)
        self.centriole_b_length_diff_input = QLineEdit(str(default_config.CENTRIOLE_TRIPLET_B_LENGTH_DIFF))
        centriole_layout.addWidget(self.centriole_b_length_diff_input, 2, 1)
        
        # Row 3: Triplet C-Length Difference
        centriole_layout.addWidget(QLabel("C Length Diff (Å):"), 3, 0)
        self.centriole_c_length_diff_input = QLineEdit(str(default_config.CENTRIOLE_TRIPLET_C_LENGTH_DIFF))
        centriole_layout.addWidget(self.centriole_c_length_diff_input, 3, 1)
        
        self.centriole_group.setLayout(centriole_layout)
        main_layout.addWidget(self.centriole_group)
        
        # Initial visibility update
        self._update_ui_visibility()

        # --- Control Button Group ---
        control_h_layout = QHBoxLayout()
        
        # Close Old Model Checkbox
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
        file_dialog = QFileDialog(self.tool_window.ui_area)
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
        
        # Update line type combo options based on mode
        current_line = self.line_type_combo.currentText()
        self.line_type_combo.clear()
        
        if is_cilia:
            # Cilia updated options: added '3Dtemplate'
            self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal', 'tip', 'primarycilia', '2Dtemplate', '3Dtemplate'])
            self.radius_label.setText("Cilia Radius (Å):")
            self.count_label.setText("Number of Doublets:")
            
            # Update general settings to Cilia defaults
            self.length_input.setText(str(default_config.CILIA_LENGTH))
            self.cilia_radius_input.setText(str(default_config.CILIA_RADIUS))
            self.num_units_input.setText(str(default_config.CILIA_NUM_DOUBLETS))
            self.curve_radius_input.setText(str(default_config.CILIA_CURVE_RADIUS))
            self.sine_frequency_input.setText(str(default_config.CILIA_SINE_FREQUENCY))
            self.sine_amplitude_input.setText(str(default_config.CILIA_SINE_AMPLITUDE))
            
            # Restore line type if valid for cilia (added '3Dtemplate')
            valid_cilia_lines = ['straight', 'curve', 'sinusoidal', 'template', 'tip', 'primarycilia', '2Dtemplate', '3Dtemplate']
            if current_line in valid_cilia_lines:
                self.line_type_combo.setCurrentText(current_line)
            else:
                self.line_type_combo.setCurrentText(default_config.CILIA_LINE)
            
        else:
            # Centriole updated options: 'template' removed
            self.line_type_combo.addItems(['straight', 'curve', 'sinusoidal'])
            self.radius_label.setText("Centriole Radius (Å):")
            self.count_label.setText("Number of Triplets:")
            
            # Update general settings to Centriole defaults
            self.length_input.setText(str(default_config.CENTRIOLE_LENGTH))
            self.cilia_radius_input.setText(str(default_config.CENTRIOLE_RADIUS))
            self.num_units_input.setText(str(default_config.CENTRIOLE_NUM_TRIPLETS))
            self.curve_radius_input.setText(str(default_config.CENTRIOLE_CURVE_RADIUS))
            self.sine_frequency_input.setText(str(default_config.CENTRIOLE_SINE_FREQUENCY))
            self.sine_amplitude_input.setText(str(default_config.CENTRIOLE_SINE_AMPLITUDE))
            
            # Restore line type if valid for centriole (excluding 'tip' and 'template' options)
            valid_centriole_lines = ['straight', 'curve', 'sinusoidal']
            if current_line in valid_centriole_lines:
                self.line_type_combo.setCurrentText(current_line)
            else:
                self.line_type_combo.setCurrentText(default_config.CENTRIOLE_LINE)
        
        # Toggle Specific Groups
        self.cilia_group.setVisible(is_cilia)
        self.centriole_group.setVisible(not is_cilia)
        
        # Ensure centerline inputs reflect current centerline type
        self._toggle_centerline_inputs()


    def _toggle_centerline_inputs(self):
        """Enable/disable curve, sine, template, or tip inputs based on line type dropdown"""
        line_type = self.line_type_combo.currentText()
        
        is_curve = line_type == 'curve'
        is_sine = line_type == 'sinusoidal'
        # Check for template types that require a file input box
        is_file_template = line_type in ('2Dtemplate', '3Dtemplate') 
        is_tip = line_type == 'tip'
        
        # 1. Curve Inputs
        self.curve_radius_label.setVisible(is_curve)
        self.curve_radius_input.setVisible(is_curve)
        self.curve_radius_label.setEnabled(is_curve)
        self.curve_radius_input.setEnabled(is_curve)
        
        # 2. Sine Inputs
        self.sine_controls_widget.setVisible(is_sine)
        
        # 3. Template File Inputs (Visible ONLY for 2Dtemplate and 3Dtemplate)
        self.template_file_group.setVisible(is_file_template)
        
        # 4. Tip Inputs 
        self.tip_length_label.setVisible(is_tip)
        self.tip_length_input.setVisible(is_tip)
        self.tip_length_label.setEnabled(is_tip)
        self.tip_length_input.setEnabled(is_tip)


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
            
            template_file = self.template_file_input.text()
            tip_length = float(self.tip_length_input.text())
            
            if length <= 0 or num_units < 0 or ring_radius <= 0:
                 raise ValueError("Length, number of units, and ring radius must be positive.")

            # Update status
            self.status_label.setText("Generating model...")
            self.status_label.setStyleSheet("color: blue; font-style: italic;")

            is_cilia = self.mode_combo.currentText().startswith('Cilia')
            new_model = None
            
            # --- CILIA (9x2 + 2) Logic ---
            if is_cilia:
                # Parse Cilia colors and other parameters needed for both generation methods
                cilia_a_color = self._parse_color(self.cilia_a_color_input.text())
                cilia_b_color = self._parse_color(self.cilia_b_color_input.text())
                cilia_cp_color = self._parse_color(self.cilia_cp_color_input.text())
                
                draw_central_pair = self.draw_cp_check.isChecked()
                should_draw_membrane = self.draw_membrane_check.isChecked() 
                membrane_radius = float(self.membrane_radius_input.text())
                membrane_fraction = float(self.membrane_fraction_input.text())
                
                if centerline_type == '3Dtemplate':
                    # --- EXECUTION BRANCH 1: 3Dtemplate (ciliabuild_from_csv) ---
                    if not template_file:
                        raise ValueError("Template file path is required when line type is '3Dtemplate'.")

                    self.session.logger.info(f"Generating cilia from 3D template file: {template_file}")
                    
                    # Call ciliabuild_from_csv, passing only the necessary rendering parameters
                    new_model = ciliabuild_from_csv(
                        session=self.session,
                        template_csv=template_file,
                        draw_central_pair=draw_central_pair,
                        membrane=should_draw_membrane,
                        membrane_radius=membrane_radius,
                        # Pass geometry defaults (not exposed in UI)
                        doublet_a_radius=default_config.CILIA_DOUBLET_A_RADIUS, 
                        doublet_b_radius=default_config.CILIA_DOUBLET_B_RADIUS,
                        cp_radius=default_config.CILIA_CP_RADIUS,
                        # Pass colors
                        doublet_a_color=cilia_a_color,
                        doublet_b_color=cilia_b_color,
                        cp_color=cilia_cp_color,
                        membrane_color=default_config.CILIA_MEMBRANE_COLOR # Use default as it's not exposed
                    )
                
                else:
                    # --- EXECUTION BRANCH 2: All other Cilia types (ciliabuild) ---
                    
                    # Enforce file path only for '2Dtemplate'
                    if centerline_type == '2Dtemplate' and not template_file:
                        raise ValueError("Template file path is required when line type is '2Dtemplate'.")

                    self.session.logger.info(f"Cilia colors - A: {cilia_a_color}, B: {cilia_b_color}, CP: {cilia_cp_color}")
                    
                    doublet_length_diff = float(self.cilia_doublet_length_diff_input.text())
                    cp_doublet_length_diff = float(self.cilia_cp_doublet_length_diff_input.text())
                    
                    # Call the command function and get the returned model
                    new_model = ciliabuild(
                        session=self.session,
                        length=length, 
                        type=centerline_type,
                        curve_radius=curve_radius,
                        sine_frequency=sine_frequency,
                        sine_amplitude=sine_amplitude,
                        # template_file is used by 2Dtemplate, but ignored by primarycilia
                        template_file=template_file, 
                        tip_length=tip_length,
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
                # --- CENTRIOLE Logic - Call centriolebuild ---
                
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
                    type=centerline_type,
                    curve_radius=curve_radius,
                    sine_frequency=sine_frequency,
                    sine_amplitude=sine_amplitude,
                    template_file=template_file,
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