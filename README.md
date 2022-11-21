# About PolyDE

PolyDE is not only just another library for finite element computation. 
It is rather a full FE package that can do fully automatic adaptive
solutions to multi-nature problems, i.e., coupled-physics.

This version of PolyDe is a fork from the original https://sourceforge.net/projects/polyde-fem/
It is a sort of "lagging mirror" to the most up-to-date version there. 

This fork is geared towards Linux, and a tighter integration with the latest 
Intel's oneAPI Fortran compiler and MKL library, and PETSc.

# Available physics modes

PolyDE can solve these physical problems in 2D and 3D. The physics mode can be invoked by
specifying the following keywords: 

## 2D:

1. ACOUSTIC: Acoustic wave propagation.

2. ELECTROSTATICS: Electrostatic problem.

3. FLOW_PROFILE: Steady-state flow.

4. FLUIDINCOMPR: Incompressible fluid dynamics

5. FLUIDELINCOMPR: Coupled electrostatics - Navier Stokes equations

6. FLUIDELMASSCON: Coupled electrostatics - Navier Stokes equations with molar mass concentration

7. FLUIDELMASSTEP: Coupled electrostatics - Navier Stokes - molar mass concentration - temperature

8. HEATTR: Steady-state heat conduction in solid.

9. TRANHEATTRS: Transient heat conduction in solid.

10. LAPLACE: General Poisson's equation.

11. MAGNETICVP: Magnetics with moving conductor.

12. STATCURRENT: Stedy-state current

13. TEWAVE: TE-mode wave propogation.

14. TMWAVE: TM-mode wave propogation.

15. THERMOELECTRIC: Thermoelectric problem

16. PLANE STRAIN: Static strain problem

17. PLANE STRESS: Static stress problem

18. THERMOELASTIC: Static thermo-elasticity

19. TRANPLSTRAIN: Plain strain from transverse forces

20. GRAVITOELASTIC: Coupled gravito-elasticity

21. PIEZOELECTRIC: Piezoelectric problem

22. PIEZOPYROELEC: Coupled piezo-pyroelectric

## 3D

1. 3DHEATTR: Steady-state heat conduction in solid.

2. 3DELECTROSTATICS: Electrostatic problem.

3. 3DSOLIDMECHANICS: Static structural stress/strain.

4. 3DSEMICONDUCTOR: Electrostatic fields of semiconductors

5. 3DTHERMOELECTRIC: Thermoelectricity

6. 3DMECHTEMPSEMICOND: Electrostatic fields of piezoresistive materials

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