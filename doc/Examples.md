Examples
========

The subpages of this page contains several examples. Some of them use 
explicit, implicit or semi-implicit (IMEX) time integration methods 
implemented in PETSc (https://www.mcs.anl.gov/petsc/). To run them, 
MHYSA needs to be compiled \b with \b PETSc. Familiarity with using PETSc 
is assumed.

\b Note: 
+ In general, any example or simulation can use PETSc time-integrators (assuming
  MHYSA is compiled with PETSc) by specifying the PETSc inputs through a 
  <B>.petscrc</B> file, similar to the ones in the examples below. 
  The following file is an example of a .petscrc file (with explanatory comments). 
  - <B>mhysa/Examples/PETScInputs/.petscrc_Example</B> 
+ The PETSc example directories have a file <B>.petscrc</B> and a sym link \b petscrc
pointing to .petscrc. The file .petscrc is the actual input file; \a petscrc is needed to generate
this documentation because Doxygen does not seem to include files with names starting with a dot!
+ The inputs in the .petscrc files (any lines not starting with a #) can also be specified in the
command line, for example, 
    
    /path/to/mhysa/bin/mhysa -use-petscts -ts_type rk -ts_rk_type 4 ...

\subpage single_species_examples :
Some examples of single species flows.

\subpage multispecies_examples : 
Some examples of multispecies flows.

\page single_species_examples Single Species Examples

The following are some examples of single species flows, where the 
governing equations are the single-species Euler/Navier-Stokes equations.

\subpage sod_shock_tube_component_rec \n

\page multispecies_examples Multispecies Examples

The following are some examples of multispecies flows.

\subpage sod_shock_tube_2species_component_rec \n

\page sod_shock_tube_component_rec 1D Euler Equations (Single Species) - Sod Shock Tube (Component-Wise Reconstruction)

Description: 
-------------------

Location: \b mhysa/Examples/SingleSpecies/1D_SodShockTube_ComponentWiseRec
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + G.A. Sod, "A survey of several finite difference methods 
    for systems of nonlinear hyperbolic conservation laws," 
    J. Comput. Phys., 27, 1 (1978).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, u = 0, p = 0.1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Component-wise 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include SingleSpecies/1D_SodShockTube_ComponentWiseRec/solver.inp

\b boundary.inp
\include SingleSpecies/1D_SodShockTube_ComponentWiseRec/boundary.inp

\b physics.inp
\include SingleSpecies/1D_SodShockTube_ComponentWiseRec/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include SingleSpecies/1D_SodShockTube_ComponentWiseRec/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=0.2: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DSodShockTube_ComponentRec.png

Expected screen output:
\include SingleSpecies/1D_SodShockTube_ComponentWiseRec/output.log

\page sod_shock_tube_2species_component_rec 1D Euler Equations (Two Species) - Sod Shock Tube (Component-Wise Reconstruction)

Description: 
-------------------

Location: \b mhysa/Examples/MultiSpecies/1D_SodShockTube_2Species_ComponentWiseRec
          (This directory contains all the input files needed
          to run this case.)

This example is basically the Sod shock tube problem with two non-reacting species. In
essence, it is identical to the single-species shock tube, but tests the multispecies
code implementation.

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + G.A. Sod, "A survey of several finite difference methods 
    for systems of nonlinear hyperbolic conservation laws," 
    J. Comput. Phys., 27, 1 (1978).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, \rho_1 = m_1*\rho, \rho_2 = m_2*\rho, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, \rho_1 = m_1*\rho, \rho_2 = m_2*\rho, u = 0, p = 0.1\f$

where \f$m_1 = 0.3\f$, \f$m_2 = 0.7\f$ are the mass fractions of each species.

Numerical Method:
 + Spatial discretization (hyperbolic): Component-wise 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include MultiSpecies/1D_SodShockTube_2Species_ComponentWiseRec/solver.inp

\b boundary.inp
\include MultiSpecies/1D_SodShockTube_2Species_ComponentWiseRec/boundary.inp

\b physics.inp
\include MultiSpecies/1D_SodShockTube_2Species_ComponentWiseRec/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include MultiSpecies/1D_SodShockTube_2Species_ComponentWiseRec/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=0.2: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DSodShockTube_2Species_ComponentRec.png

Expected screen output:
\include MultiSpecies/1D_SodShockTube_2Species_ComponentWiseRec/output.log
