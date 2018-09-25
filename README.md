# qfnum
Numerical implementation.
Authors M.O. & C.H.E.

## Organization
In general the code has the following structure:
It is written down in .c files in which we comment how the code works. 
Variables and functions along with their parameters are predefined and fully commented in the according header.

Qfnum consists of general functions that are useful, when calculating quantum friction numerically and modules for certain geometries. 
General stuff can be found in qfhelp.c (and the according header), whereas we have (for now) two geometries in cyl.c and plate.c.
With main.c we can then calculate quantum friction. 

The aim is to provide the program with a flag for the geometry and a file in which a lenghty list of parameters is given, like this:
$ qfnum cyl cylparams.dat.
