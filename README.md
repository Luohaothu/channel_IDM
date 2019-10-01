# channel_LES
 Large eddy simulation of turbulence in a straight channel

Author: Haining Wang

This is the C++ version of the Fortran program developed by Huang.



Using C++, the following advantages are achieved:

separation of data and algorithms by encapsulating modules into classes;

new fftw library supports arrays of arbitrary length;

the old F77 codes did not make use of any edge of Fortran (such as array operations) and are totally procedure-oriented;

the old codes did not use 'implicit none', and common variables are declared in every subroutines.



The codes include the following modules:

class Para: manage (IO, storing, etc) case parameters for the computations

class Mesh: store grid information, including the difference scheme coefficients and wave numbers for fft

class Field: store flow fields of the primitive and intermediate variables, providing basic operation tools (memory copy, IO, etc)

class IDM: the implicit decoupling procedure algrithm, applying second order central difference and Crank-Nicolson scheme, dependent upon the mesh

class Statis: tools for calculating basic statistics, useful for computation state monitoring and post processing

class Martix: functions for solving tri-diagonal matrix, which is essential for the efficiency the whole program 



Reference:

Kyoungyoun Kim, et al. (2002). "An implicit velocity decoupling procedure for the incompressible Navier–Stokes equations." INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS 38:125–138.