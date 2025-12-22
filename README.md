# ChimeraX CiliaBuilder 

[![DOI](https://zenodo.org/badge/1115546661.svg)](https://zenodo.org/badge/latestdoi/1115546661)

ChimeraX Tool to build 3D model of cilia and centriole.

## Installation

Download the repository or git clone and install **CiliaBuilder**:

Download or git clone
```bash
git clone https://github.com/builab/CiliaBuilder.git
```

If you download, then unzip

Open ChimeraX
```bash
devel install [path_to_CiliaBuilder_folder]
```

## Interface usage
Tools > High-order structure > Cilia Builder

Then enjoy

## Commandline usage (NOT RECOMMENDED)

### Simple straight cilia
ciliabuild length 10000

### Curved cilia
ciliabuild length 10000 type curve curve_radius 10000

### Sinusoidal cilia
ciliabuild length 15000 type sinusoidal sine_frequency 1 sine_amplitude 1000

### Cilia with tip
ciliabuild length 10000 type tip tip_length 4000 

### Primary cilia
ciliabuild length 10000 type primarycilia

### Using 2D template file
ciliabuild length 10000 type 2Dtemplate template_csv

### Custom number of doublets
ciliabuild length 5000 num_doublets 12 cilia_radius 1000


### Basic centriole (offset 30 degree)
centriolebuild length 3000

### Centriole with custom angle offset to match with cilia
centriolebuild length 2000 centriole_angle_offset 0

### Custom triplet configuration
centriolebuild length 5000 num_triplets 12 centriole_radius 1200


TODO
- Movie part (low priority)
- For ArtiaX build, perhaps presupply a model of doublet & central pair volume at 8-nm repeat (expand curve.py).

DONE
- Curve and sinusoidal from same direction (Done)
- Make the UI version
- Make membrane.
- Group of group modification
- Name change to CiliaBuilder
- Membrane as two layers and cover (Much better) 0.4.4
- When doing sinuisoidal, the CP is actually much shorter. (Priority)
- Make 1 Angstrom different at the starting end as well for better visualization. (Low priority)
- Button to default close the old models & generate new models. 0.4.6
- Change default color
- Make coloring possible for doublet, triplet. 0.4.7
- Make a,b different for cilia as well to make central pair longer 0.4.8
- Generate centriole with offset so it is continuous with cilia
- Make default_config.py for easy controlling of the data. 0.5.2
- Add initial part for sinusoidal to align it. cos(90 - atan2(amplitude, length/(freq*4)))*(2*cilia_radius)
- Use a template curve for line type (csv: X,Y only) and interface update. 0.5.6
- Too many membrane, too many central pair
- generate_tip_curves not ideal yet
- The merged has a wrinkle 0.9.1
- Add to interface & cosine drop better 0.9.3
- Add the cap
- CP rung added
- Add primary to interface, change template to template2D
- Fix Idx_B = 0 for the 1st Z
- Gap in primary cilia
- Change the format of 3D csv to accomodate triplet
- Reorganize read_2D_csv, read_3D_csv to io
- Primary cilia (Using template first)
- Clean up & refactor code (Priority)




