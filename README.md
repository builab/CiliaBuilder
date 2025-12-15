# CiliaSim or might change to CiliaBuilder later
Only testing internally now. Release in 6 months with a lot of more function planned.

## Installation

Clone the repository and install **CiliaSim**:

```bash
git clone https://github.com/builab/CiliaSim.git
cd [DownloadFolder]
devel install CiliaSim
```

### Simple straight cilia
ciliasim length 15000 cilia_radius 875

### Curved cilia
ciliasim length 15000 line curve curve_radius 10000

### Sinusoidal cilia
ciliasim length 25000 line sinusoidal sine_frequency 1 sine_amplitude 2000

### Custom number of doublets
ciliasim length 5000 num_doublets 12 cilia_radius 1000


### Basic centriole
centriolesim length 5000

### Centriole with custom angle offset
centriolesim length 5000 centriole_angle_offset 45

### Custom triplet configuration
centriolesim length 5000 num_triplets 12 centriole_radius 1200


TODO
- When doing sinuisoidal, the CP is actually much shorter.
- Make a,b different for cilia as well to make central pair longer
- Make coloring possible for doublet, triplet
- Membrane as two layers and cover.
- Pre-calculated tip for Tetrahymena & Primary cilia.

DONE
- Curve and sinusoidal from same direction (Done)
- Make the UI version
- Make membrane.
- Group of group modification
