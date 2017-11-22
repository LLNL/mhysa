Examples
========

\subpage basic_examples :
Some basic examples that are simulated using MHYSA. They 
all use explicit time integration, and \b do \b not require MHYSA to be compiled
with PETSc. Most of them can be run on one or a small number of processors.

\subpage petsc_examples : 
Some examples that use implicit or semi-implicit (IMEX) time
integration methods implemented in PETSc. To run them, MHYSA needs to be compiled \b with \b PETSc.

\subpage ib_examples : Examples that use the immersed boundary method to simulate various geometries.

\page basic_examples Basic Examples

The following are some basic examples that are simulated using MHYSA. They 
all use explicit time integration, and \b do \b not require MHYSA to be compiled
with PETSc.

\subpage sod_shock_tube_component_rec \n
\subpage sod_shock_tube_2species_component_rec \n

\page petsc_examples PETSc Examples

The following are some examples that use explicit, implicit or semi-implicit (IMEX) time
integration methods implemented in PETSc (https://www.mcs.anl.gov/petsc/). To run them, 
MHYSA needs to be compiled \b with \b PETSc. Familiarity with using PETSc is assumed.

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



Implicit time integration:
--------------------------

Implicit-Explicit (IMEX) time integration:
------------------------------------------


\page ib_examples Immersed Boundaries Examples

The following are some examples are use the immersed boundary methodology to solve the PDE in the presence of
various geometries. To use the immersed boundary implementation in HyPar, an STL (https://en.wikipedia.org/wiki/STL_%28file_format%29) 
representation of the immersed body is necessary. Note:
+ The immersed boundary method is \b only implemented for 3-dimensional simulations (#HyPar::ndims = 3).
+ It can be used with only those physical models that define an immersed boundary implementation (#HyPar::IBFunction()), for
  example, the 3D Navier-Stokes equations (NavierStokes3DImmersedBoundary()).

3D Navier-Stokes Equations:
---------------------------


\page sod_shock_tube_component_rec 1D Euler Equations (Single Species) - Sod Shock Tube (Component-Wise Reconstruction)

Description: 
-------------------

Location: \b mhysa/Examples/1D_Euler/SodShockTube_ComponentWiseRec
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
\include 1D_Euler/SodShockTube_ComponentWiseRec/solver.inp

\b boundary.inp
\include 1D_Euler/SodShockTube_ComponentWiseRec/boundary.inp

\b physics.inp
\include 1D_Euler/SodShockTube_ComponentWiseRec/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D_Euler/SodShockTube_ComponentWiseRec/aux/init.c

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
\include 1D_Euler/SodShockTube_ComponentWiseRec/output.log

\page sod_shock_tube_2species_component_rec 1D Euler Equations (Two Species) - Sod Shock Tube (Component-Wise Reconstruction)

Description: 
-------------------

Location: \b mhysa/Examples/1D_Euler/SodShockTube_2Species_ComponentWiseRec
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
\include 1D_Euler/SodShockTube_2Species_ComponentWiseRec/solver.inp

\b boundary.inp
\include 1D_Euler/SodShockTube_2Species_ComponentWiseRec/boundary.inp

\b physics.inp
\include 1D_Euler/SodShockTube_2Species_ComponentWiseRec/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D_Euler/SodShockTube_2Species_ComponentWiseRec/aux/init.c

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
\include 1D_Euler/SodShockTube_2Species_ComponentWiseRec/output.log
