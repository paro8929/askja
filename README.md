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

'params_params_N8S10T10B20k.txt': This is an example parameter file for an HMC run. It has two columns. Left: text indicator. Right: value. The values on the rhs control the particular run. They are

NC              Number of colors for the simulation (integer)
DIM             Number of parent dimension (integer)
N0              Number of lattice sites in direction of dimension 0 (integer)
N1              Number of lattice sites in direction of dimension 1 (integer)
N2              Number of lattice sites in direction of dimension 2 (integer)
N3              I guess you get the point
N4              
N5              
N6              
N7              
N8              
N9              
BETA            Value of lattice coupling (double, two digits after comma)
DT              Time stepping increment (double)
STEPS           Total Number of steps to be taken (integer)
UPDATE          Report results in output files each time after this many steps (integer)
UPDATE2         Check to see if current configuration can be accepted. (integer)
CONFTAKE        Store one of of this many accepted configurations and store in file (integer)
LEVEL           Order in Taylor expansion of matrix exp(A)=1+A+A^2/2+... (integer smaller than 9)
VERBOSE         Can I talk now? Please????
RUNCONTINUE     0 for a new run, 1 if pre-existing configuration should be used as initialization
RANDOM          0 for fixed random seed, 1 for random initial seed (all other seeds will be random, though)


