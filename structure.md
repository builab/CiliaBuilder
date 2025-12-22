CiliaBuilder/
├── bundle_info.xml			
├── license.txt			    
├── README.md			
├── structure.md		# Source script
├── example/
└── src/
    ├── __init__.py      # Export function
    ├── default_config.py           # Default values (to be modified before installing)
    ├── cmd.py           # ChimeraX commands
    ├── tool.py          # ChimeraX UI
    ├── io.py            # io function
    ├── draw.py          # Rendering functions
    └── geometry/
        ├── __init__.py
        ├── primarycilia_template.csv  # Primary cilia template
        ├── centerline.py    # Centerline generation
        ├── primarycilia.py  # Primary cilia generation
        └── tip.py           # Tip generation