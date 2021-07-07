# channel_LES
 Large eddy simulation of turbulence in a straight channel

Author: Haining Wang

This is the C++ version of the Fortran program originally developed by Huang.

By using C++, the following advantages are achieved:

- Data and algorithms are **decoupled** by encapsulating modules into classes
- New fftw library supports arrays of **arbitrary** length

Meanwhile, several shortcomings of the original F77 implementation are avoided:

- The original implementation fails to fully utilize the unified syntax provided by the language
- The original implementation is procedural-oriented thus laking proper modularization
- The original implementation heavily depends on global variables and implicit data types throughout the code


The code includes the following modules:

- ``class Para``: manages (IO, storing, etc) case parameters for the computations
- ``class Mesh``: stores grid information, including the differential scheme coefficients and wave numbers for FFT
- ``class Field``: stores flow fields of the primitive and intermediate variables & provides tools for basic operations
  (memory copy, IO, etc)
- ``class IDM``: the actor class for the **I**mplicit **D**ecoupled **M**ethod, which is 2nd order in both time & space
- ``class Statis``: holds tools for calculating basic statistics, useful for runtime monitoring and post-processing
- ``class Martix``: holds functions for solving tri-diagonal matrix, which is the key to the high efficiency of the whole program 



Reference:

Kyoungyoun Kim, et al. (2002). "An implicit velocity decoupling procedure for the incompressible Navier–Stokes equations."
INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS 38:125–138. [[link](https://doi.org/10.1002/fld.205)]