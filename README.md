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

PROGRESS
- First working version
- 0.2 now drawing cp, doublet, triplet


TODO
- Make the centriolesim

DONE
- Curve and sinusoidal from same direction (Done)
- Make the UI version
