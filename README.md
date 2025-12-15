# ChimeraX CiliaBuilder
Only testing internally now. Probably release in 6 months with a lot of more function planned.

## Installation

Clone the repository and install **CiliaBuilder**:

```bash
git clone https://github.com/builab/CiliaBuilder.git
cd [DownloadFolder]
devel install CiliaBuilder
```

### Simple straight cilia
ciliabuild length 15000 cilia_radius 875

### Curved cilia
ciliabuild length 15000 line curve curve_radius 10000

### Sinusoidal cilia
ciliabuild length 25000 line sinusoidal sine_frequency 1 sine_amplitude 2000

### Custom number of doublets
ciliabuild length 5000 num_doublets 12 cilia_radius 1000


### Basic centriole
centriolebuild length 5000

### Centriole with custom angle offset
centriolebuild length 5000 centriole_angle_offset 45

### Custom triplet configuration
centriolebuild length 5000 num_triplets 12 centriole_radius 1200


TODO
- When doing sinuisoidal, the CP is actually much shorter.
- Make a,b different for cilia as well to make central pair longer
- Make coloring possible for doublet, triplet
- Pre-calculated tip for Tetrahymena & Primary cilia.

DONE
- Curve and sinusoidal from same direction (Done)
- Make the UI version
- Make membrane.
- Group of group modification
- Name change to CiliaBuilder
- Membrane as two layers and cover (Much better) 0.4.4

