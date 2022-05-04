# Ideal gas numerical simulation

This is an old code from my bachelor degree. Unfortunately up to now it is commented in Italian and it has some limitations in the graphical part.
In order to obtain the graphs I used an executable which calls Gnuplot with some given instructions.
This executable is not present in this directory, but essentially it gived to Gnuplot the output files of the main program.

The gas is simulated with two different boundary conditions: "periodic" and "box".

You can distinguish which file belongs to which boundary conditions from the prefix (_ can be "periodic" or "box").

What is contained in this repository:
---------

* _gas.h   : header file with all classes and functions needed for simulation. Velocity-Verlet algorithm

* _gas.cxx : main function file (In the box case header and main are not separated!)

* _gas.gif    : gif obtained with the program
