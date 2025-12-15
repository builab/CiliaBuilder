INSTALLATION

devel install CiliaSim_extension_path

# Simple straight cilia
ciliasim length 5000 cilia_radius 875

# Curved cilia
ciliasim length 5000 line curve curve_radius 10000

# Sinusoidal cilia
ciliasim length 25000 line sinusoidal sine_frequency 1 sine_amplitude 2000

# Custom number of doublets
ciliasim length 5000 num_doublets 12 cilia_radius 1000


# Basic centriole
centriolesim length 5000

# Centriole with custom angle offset
centriolesim length 5000 centriole_angle_offset 45

# Curved centriole
centriolesim length 5000 line curve curve_radius 10000

# Custom triplet configuration
centriolesim length 5000 num_triplets 12 centriole_radius 1200


TODO
- Make coloring possible for doublet, triplet
- what about group of group modification?
- Make a,b different for cilia as well to make central pair longer
- Make membrane look nicer
- Can you subtract surface in chimeraX?

DONE
- Curve and sinusoidal from same direction (Done)
- Make the UI version
- Make membrane.