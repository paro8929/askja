# askja
Hybrid Monte Carlo for SU(N) in arbitrary dimension D

Version 1.0
------------
Date: 2016-12-16

Developer: Paul Romatschke

Copy-Right: is granted, provided you cite the relevant papers

-------------
Content: 

* README: this file
* INSTALL: instructions on how to install askja
* sun.cpp: source code for HMC simulation
* measurements1d.cpp: source code for measuring observables (mostly 1d)
* measurements2d.cpp: source code for measuring observables (mostly 2d)
* measurements4d.cpp: source code for measuring observables (mostly 4d)
* sunmat.h: header file providing routines for SU(N) matrices for arbitrary N
* structs.h: header file providing routines for D-dimensional coordinate systems
* compile: bash script to compile all executables
* params_N8S10T10B20k.txt: example parameter file to start an HMC run

--------------
Description:

*'params_params_N8S10T10B20k.txt': This is an example parameter file for an HMC run. It has two columns. Left: text indicator. Right: value. The values on the rhs control the particular run. They are

NC:              Number of colors for the simulation (integer)

DIM:             Number of parent dimension (integer)

N0:              Number of lattice sites in direction of dimension 0 (integer)

N1:              Number of lattice sites in direction of dimension 1 (integer)

N2:              Number of lattice sites in direction of dimension 2 (integer)

N3:              I guess you get the point
           
BETA:            Value of lattice coupling (double, two digits after comma)

DT:              Time stepping increment (double)

STEPS:           Total Number of steps to be taken (integer)

UPDATE:          Report results in output files each time after this many steps (integer)

UPDATE2:         Check to see if current configuration can be accepted. (integer)

CONFTAKE:        Store one of of this many accepted configurations and store in file (integer)

LEVEL:           Order in Taylor expansion of matrix exp(A)=1+A+A^2/2+... (integer smaller than 9)

VERBOSE:         Can I talk now? Please????

RUNCONTINUE:     0 for a new run, 1 if pre-existing configuration should be used as initialization

RANDOM:          0 for fixed random seed, 1 for random initial seed (all other seeds will be random, though)

You can change e.g. DIM to a different value, but then you should make sure there are corresponding entries N0, N1, ... specifying the lattice sites in each of the dimensions. The minimum number of lattice sites is 1.

You can change NC as you like provided it is larger than 1. There probably is a maximum value for NC where things get numerically unstable, but I've successfully run simulations for NC=64. 

------------------------------------------------------------------

* Running askja:

Do

./askja [params_file]

to start an HMC run with the parameters specified in [params_file]

This will initialize the run and start generating a directory called output-DIM[DIM] where it will store diagnostics such with file names such as as "constraints-...", "Wilson-...", as well as gauge configurations called "U-..." with the [...] specifying the parameters for the run with which they were generated. After the number of [STEPS] has passed, the code will stop and shut down and you will have the gauge field configurations stored in the directory output-DIM[DIM]. 

* Extracting observables

To measure observables, use for instance

./measure4d output-DIM[DIM]/master-[...]

This will generate a new directory called rav-DIM[DIM] where measurements will be stored. The current measurement output are the file Nrav-[...], and the most important part is the last line of this file which has 12 columns. These are: 

           1) Configuration Number
           2) Polyakov loop ensemble average
           3) Statistical error
           4) Wilson loop ensemble average
           5) Statistical error
           6) Temporal plaquette ensemble average
           7) Statistical error
           8) Spatial plaquette ensemble average
           9) Statistical error
           10) Lattice sites in temporal (0) direction
           11) Lattice sites in spatial (1) direction
           12) BETA value
     
You can skip the first 10 configurations of a run (recommended) by e.g. doing

./measure4d output-DIM[DIM]/master-[...] 10

This will reduce sensitivity to the initial (non-thermalized) part of the HMC run.

You can also make measure be more verbose (and skip no steps) by doing

./measure4d output-DIM[DIM]/master-[...] 0 1

Finally, you can make measure measure the eigenvalue distribution of Wilson and Polyakov loops (and not skip any steps or be verbose) by doing

./measure4d output-DIM[DIM]/master-[...] 0 0 1

This will create a file called 'Hist-[...]' in rav-DIM[DIM] which is a histogramm of the eigenvalue arguments. More precisely, the eigenvalues of an SU(N) matrix will be exp(i phi) and the Hist-[...] file will tell you how phi is distributed from [-Pi;Pi]. Note that in doing that any non-trivial winding number is taken out (e.g. the sum of all eigenvalue arguments is renormalized to zero).

