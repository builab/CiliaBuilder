# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QLabel, QLineEdit, QPushButton)
from Qt.QtCore import Qt

class CiliaSim(ToolInstance):
    """
    CiliaSim tool for generating cilia microtubule doublet models
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
        title_label = QLabel("Cilia Microtubule Doublet Generator")
        title_label.setStyleSheet("font-size: 14pt; font-weight: bold;")
        main_layout.addWidget(title_label)

        # Parameters grid layout
        params_layout = QGridLayout()

        # Cilia Length
        params_layout.addWidget(QLabel("Cilia Length (Å):"), 0, 0)
        self.length_input = QLineEdit("5000")
        params_layout.addWidget(self.length_input, 0, 1)

        # Cilia Diameter (not used for now, but included for future)
        params_layout.addWidget(QLabel("Cilia Diameter (Å):"), 1, 0)
        self.diameter_input = QLineEdit("190")
        self.diameter_input.setEnabled(False)  # Disabled for now
        self.diameter_input.setToolTip("Not used in current version")
        params_layout.addWidget(self.diameter_input, 1, 1)

        # A-tubule length difference
        params_layout.addWidget(QLabel("A-tubule extra length (Å):"), 2, 0)
        self.length_diff_input = QLineEdit("5")
        params_layout.addWidget(self.length_diff_input, 2, 1)

        # Color of A-tubule
        params_layout.addWidget(QLabel("A-tubule color (R,G,B,A):"), 3, 0)
        self.color_a_input = QLineEdit("255,100,100,255")
        self.color_a_input.setToolTip("Enter color as: R,G,B,A (0-255)")
        params_layout.addWidget(self.color_a_input, 3, 1)

        # Color of B-tubule
        params_layout.addWidget(QLabel("B-tubule color (R,G,B,A):"), 4, 0)
        self.color_b_input = QLineEdit("100,100,255,255")
        self.color_b_input.setToolTip("Enter color as: R,G,B,A (0-255)")
        params_layout.addWidget(self.color_b_input, 4, 1)

        main_layout.addLayout(params_layout)

        # Generate button
        generate_button = QPushButton("Generate Cilia Model")
        generate_button.setStyleSheet("font-size: 12pt; padding: 10px;")
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
            # Get parameters
            length = float(self.length_input.text())
            length_diff = float(self.length_diff_input.text())
            
            # Parse colors
            color_a = self._parse_color(self.color_a_input.text())
            color_b = self._parse_color(self.color_b_input.text())

            # Update status
            self.status_label.setText("Generating model...")
            self.status_label.setStyleSheet("color: blue; font-style: italic;")

            # Import and call the drawdoublet function
            from .drawdoublet import drawdoublet
            
            # Generate the model
            tube_a, tube_b = drawdoublet(
                self.session,
                length=length,
                interval=80,
                angle=0,
                radius_a=125,
                radius_b=130,
                shift_distance=100,
                length_diff=length_diff,
                name="cilia_doublet",
                color_a=color_a,
                color_b=color_b
            )

            # Update status
            self.status_label.setText(f"✓ Model generated successfully! Length: {length} Å")
            self.status_label.setStyleSheet("color: green; font-weight: bold;")
            
            # Log to session
            self.session.logger.info(f"Generated cilia doublet model:")
            self.session.logger.info(f"  Length: {length} Å")
            self.session.logger.info(f"  A-tubule extra length: {length_diff} Å")
            self.session.logger.info(f"  A-tubule radius: 125 Å, B-tubule radius: 130 Å")

        except ValueError as e:
            error_msg = f"Invalid input: {str(e)}"
            self.status_label.setText(f"✗ {error_msg}")
            self.status_label.setStyleSheet("color: red; font-weight: bold;")
            self.session.logger.error(error_msg)
        except Exception as e:
            error_msg = f"Error generating model: {str(e)}"
            self.status_label.setText(f"✗ {error_msg}")
            self.status_label.setStyleSheet("color: red; font-weight: bold;")
            self.session.logger.error(error_msg)