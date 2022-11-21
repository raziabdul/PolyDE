# About PolyDE

PolyDE is not only just another library for finite element computation. 
It is rather a full FE package that can do fully automatic adaptive
solutions to multi-nature problems, i.e., coupled-physics.

This version of PolyDe is a fork from the original https://sourceforge.net/projects/polyde-fem/
It is a sort of "lagging mirror" to the most up-to-date version there. 

This fork is geared towards Linux, and a tighter integration with the latest 
Intel's oneAPI Fortran compiler and MKL library, and PETSc.

# Prerequisites to build from the source

- A rather recent Linux distro, with libc close to version 2.36
- Intel's oneAPI collection including MKL
- PETSc
- pgplot
- CMake

# Environment variable

Polyde only needs one environment variable: PolydePath

To set up, include the two lines in .bash_profile : 

PolydePath=/home/your_user_name/Polyde/

export PolydePath

# Compilation

Setup for compilation is as the standard procedure with CMake.

# Running PolyDE

The two FEMsettings.txt should be properly set up.