CiliaBuilder/
├── bundle_info.xml			
├── license.txt			    
├── README.md			
├── structure.md		# Source script
└── src/
    ├── __init__.py
    ├── default_config.py           # ChimeraX commands (high-level interface)
    ├── cmd.py           # ChimeraX commands (high-level interface)
    ├── draw.py          # Rendering functions
    ├── geometry/
    │   ├── __init__.py
    │   ├── centerline.py    # Centerline generation (from curve.py)
    │   ├── tip.py           # Tip-specific geometry
    │   └── base.py          # Shared geometry utilities
    └── io.py        # CSV read/write functions
