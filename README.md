# Lumen network
**Lumen network** is a Python-based code for simulations of coarsening in a 2D network of lumens. It was developed to describe the process of blastocoel formation in the mouse embryo [1]. The network is fixed and the lumen are represented by nodes, each described by their area, tension and contact tension, and connected by edges. The code solves a system of coupled non-linear equations for the area dynamics of each lumen, by computing the flux between them. A simulation starts from an initial configuration, which is defined by the user, and runs until one single lumen is left, and has siphoned all the others.

Please cite us to use this code:
[1] [Dumortier, et al., Biorxiv, 2019](https://www.biorxiv.org/content/10.1101/537985v1)

## Requirements

* Python 2.7
* Scipy library, in particular numpy
* networkx
* configparser
* Optional (only if you want to plot the networks) :
	* matplotlib.pyplot
	* matplotlib.patches

## Initialization

In order to install Python libraries, we advise you to use pip. Install pip if not already installed (see [here](https://pip.pypa.io/en/stable/installing/) for instructions). Then install the required libraries using
```pip2 install scipy networkx configparser matplotlib```

Then, make sure to enter the absolute path for file generation. Go to `outputs/example.ini` and change the `path = ABSOLUTE/PATH/git/lumen_network/outputs/config` by your actual absolute path.

## Content

* network/
	* config_lib.py : library for generation of configuration files.
	* gen_config.py : script for the generation of configuration files proper.
	* network.py 	: library for the network simulation.
	* simulation.py : script that runs a simulation.
	* submit.py		: script to submit simulation runs on a cluster, configured with the workload manager Slurm.

* patterns/
	* config_chain.pat		: template for a chain (1d) configuration.
	* config_embryo.pat		: template for the embryo-like network configuration.
	* config_embryo_noisy.pat : same but for a noisy embryo-like structure.
	* config_empty.pat		: empty template.
	* config_hexagonal.pat	: template for a hexagonal network configuration.
	* config_triangular.pat	: template for a triangular network configuration.

* outputs/
	* example.ini : configuration file for the simulation.


## Instructions

### Quick use
1. Go to the *outputs* folder (`cd lumen_network/outputs`), then :
2. Generation of a network ```../network/gen_config.py ../pattern/config_embryo.pat```
3. Launch the simulation ```cd config ; ../../network/simulation.py example.ini```

### Generating a network
Generation of network is made using the *gen_config.py* script.

#### Options
The code uses the library configparser to read template files. Template files are organized in sections, containing arguments. Below is the extensive list of sections and respective arguments. However, you can also use command lines to generate directly the configuration.

* files
	* outdir 	: *string*, path where configurations while be generated.
	* tpl 		: *string*, template for the simulations. Default is *example.ini*
	* subdir 	: *string*, name of the configuration. Default is *config*
	* nsim 		: *int*, number of simulations. Default is 1

* network
	* topology 		: *string*, topology of the network. Can be *hexagonal*, *triangular* or *chain* (1d).
	* nlayers		: *int*, number of layers of the network.
	* pbc			: *boolean*, include periodic boundary conditions if *True*. Works only for 1D *chain* configurations.
	* nbicellular 	: *int*, number of bicellular lumens.
	* seed 			: *int*, specify the seed for RNG. Default is *None*, else an integer.

* tensions
	* gamma_border	 : (positive) *float*, tension of the lumens at the border.
	* gamma_c_border : (positive) *float*, contact tension of the lumens at the border.

	**N.B.:** the tensions in the bulk are set to 1. There is a geometrical constraint, such that for ![equation](https://latex.codecogs.com/gif.latex?\gamma&space;<&space;\gamma_c&space;/&space;2), the geometry is not defined.

* noisy
	* noisy 		: *boolean*, if *True*, the configuration will be spatially noisy.
	* lumen_pos_avg : *float*, mean of the spatial noise.
	* lumen_pos_std : *float*, standard deviation of the spatial noise.

* volume
	* vol_avg : *float*, average of the initial distribution of areas.
	* vol_std : *float*, standard deviation of the initial distribution of areas.

* chimera
	* chimeras 		: *boolean*, if *True*, the generated embryo will be a chimeric one.
	* gamma_chimera : float, if chimeras, value of the tension of the mutant lumens.

* display
	* show : *boolean*, if *True*, the network is displayed.

#### Example
The generic instruction for generating a network is :

```./gen_config.py mypattern.pat```

You can either generate networks from pre-made .pat files, or simply use the command line.

* If you are in the /outputs/ folder, then write for instance
	```../network/gen_config.py ../patterns/config_embryo.pat```
* Alternatively, if you want to generate configurations from the command line directly, you may write for example
	```./network/gen_config.py outdir=outputs/ tpl=outputs/example.ini topology=hexagonal nlayers=2 ```

#### Outputs
*gen_config.py* script creates a *network/* folder containing:
* lumen.dat 			: properties of each lumen (gamma1, gamma2, gammac, initial area, lumen types).
* lumen_coord.dat 		: (x, y) coordinates of a lumen.
* lumen_lumen.dat 		: array of the lumen-lumen connections (id lumen1, id lumen2, distance).
* bridge_lumen.dat 		: array of the bridge-lumen connections (id bridge1, id lumen2, distance).
* bridge_bridge.dat 	: array of the bridge-bridge connections (id bridge1, id bridge2, distance).
* bridgeconversion.dat 	: conversion list for bridge to lumen id (bridge, corresponding empty lumen).

**N.B.1**: for symmetric lumens, we choose gamma1=gamma2.
**N.B.2**: lumen types can be ICM-multi (0), TE-multi (1), ICM-bi (2), TE-bi (3), mutant (4), TE-mutant (5).
**N.B.3**: it may happen that id bridge1 = id lumen2, but bridges and lumens have different lists of ids.

### Simulation
To run a simulation from a network configuration, use the script *simulation.py* and the library *network.py* coupled with an .ini file that stores the parameters of the simulation.

```./simulation.py example.ini```

Note that the output files will be stored in a `out/` folder.

#### Options
* network
	* path 			: *string*, path where the (*out/*) output files will be stored.
	* tube_radius 	: *float*, radius of the tube/channel in absolute units. Used for the empty lumen condition. Default is 0.01
	* friction 		: *float*, friction of the fluid. Default is 1.
* time
	* time_max 	: *float*, maximum time allowed. Default is 10000.
	* time_step : *float*, time step. Default is 0.01.
* integration
	* min_step 	 : *float*, minimum time for the integration window. Default is 0.
	* max_step 	 : *float*, max step time for the integration window. Default is 0.1
	* integrator : *string*, type of integrator. Can be **odeint**, **solve_ivp**, **RK23**. Default is odeint.
	* adaptative : *boolean*, use an adaptative time step if *True*.
	* save_area  : *boolean*, if *True*, record the areas in the *out/area.dat* file.
* swelling
	* swelling_bool : *boolean*, if *True*, add swelling to the simulation.
	* swelling_rate : *float*, value of the pumping rate. Beware that a too big value will make make the system explode.

### Ouputs
Outputs of a simulation are stored in an *out/* folder containing:
* area.dat 	: file of area dynamics at computed time steps (note that it needs to have the option *save_area=True* in the .ini to be filled).
* area.log 	: main informations of the simulation : initial and final total volumes, initial volume distribution, winning lumen id and ending time.
* error.dat : recording of the errors during simulation, such as volume loss.
* event.dat : recording of the events of the network, such as disappearance of a lumen, conversion of a lumen into a bridge, ...

### Slurm
You can launch simulations on a cluster using Slurm, using the *submit.py* script, as

```python submit.py config????```

where `config????` is the list of configurations you want to run (config0000, config0001, ...).

## Data visualisation
We provide Jupyter Notebooks for data visualisation of some results obtained with this code. Go in the *jupyter/* folder and use a notebook to visualize. You may also use a Jupyter viewer.


## About

### Authors
* Annette Mielke
* Mathieu Le Verge-Serandour
* Hervé Turlier

### Virtual Embryo
Virtual Embryo is a research team focused on the physics of morphogenesis, led by H. Turlier, at the Center for Interdisciplinary Research in Biology (CNRS UMR 7241 / INSERM U1050), in Collège de France, Paris, France.
More details here: [Virtual Embryo](https://www.virtual-embryo.com)

### Contact us
Enquiries about this project should be addressed to [Mathieu Le Verge-Serandour](mailto:mathieu.le-verge-serandour@college-de-france.fr) or [Hervé Turlier](mailto:herve.turlier@college-de-france.fr)


### License

Copyright (c) 2019 [Virtual Embryo](https://www.virtual-embryo.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

* The original paper should be cited: [Dumortier, et al., Biorxiv, 2019](https://www.biorxiv.org/content/10.1101/537985v1)

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
