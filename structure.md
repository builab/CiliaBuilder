CiliaBuilder/
├── bundle_info.xml			
├── license.txt			    
├── README.md			
├── structure.md		# Source script
├── example/
└── src/
    ├── __init__.py
    ├── default_config.py           # ChimeraX commands (high-level interface)
    ├── cmd.py           # ChimeraX commands (high-level interface)
    ├── io.py            # io but might killed in the future as we 
    ├── draw.py          # Rendering functions
    └── geometry/
        ├── __init__.py
        ├── centerline.py    # Centerline generation (from curve.py)
        └── tip.py          # Shared geometry utilities
