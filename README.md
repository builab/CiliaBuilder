# ChimeraX CiliaBuilder 

[![DOI](https://zenodo.org/badge/1115546661.svg)](https://zenodo.org/badge/latestdoi/1115546661)

ChimeraX Tool to build 3D models of cilia and centriole.

## Installation

Download the repository or git clone and install **CiliaBuilder**:

Download or git clone
```bash
git clone https://github.com/builab/CiliaBuilder.git
```

If you download, then unzip

Open ChimeraX and type this in # ChimeraX commandline # (Tested to be compatible with ChimeraX 1.7 and up to 1.10 but not tested with earlier or later version)
```bash
devel install [path_to_CiliaBuilder_folder]
```

## Interface usage
Open Tools > High-order structure > Cilia Builder

### UI

[![Screenshot](example/CiliaBuilder_screenshot.png)](example/CiliaBuilder_screenshot.png)


### Video Tutorial
To Come



## Commandline usage

### Simple straight cilia
```bash
ciliabuild length 10000
```

### Curved cilia
```bash
ciliabuild length 10000 type curve curve_radius 10000
```

### Sinusoidal cilia
```bash
ciliabuild length 15000 type sinusoidal sine_frequency 1 sine_amplitude 1000
```

### Cilia with tip
```bash
ciliabuild length 10000 type tip tip_length 4000 
```

### Primary cilia
```bash
ciliabuild length 10000 type primarycilia
```

### Using 2D template file
```bash
ciliabuild length 10000 type 2Dtemplate template_csv
```

### Custom number of doublets
```bash
ciliabuild length 5000 num_doublets 12 cilia_radius 1000
```


### Basic centriole (offset 30 degree)
```bash
centriolebuild length 3000
```

### Centriole with custom angle offset to match with cilia
```bash
centriolebuild length 2000 centriole_angle_offset 0
```

### Custom triplet configuration
```bash
centriolebuild length 5000 num_triplets 12 centriole_radius 1200
```



