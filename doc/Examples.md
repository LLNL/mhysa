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

Basic examples
--------------

\subpage density_wave_advection_1d \n
\subpage sod_shock_tube_component_rec \n
\n
\subpage riemann_case4_component_rec \n
\subpage vortex_convection \n
\n
\subpage density_sine_wave_advection \n
\subpage isotropic_turbulence \n

Examples with immersed boundaries
---------------------------------

\subpage shock_cylinder \n

\page multispecies_examples Multispecies Examples

The following are some examples of multispecies flows.

\subpage sod_shock_tube_2species_component_rec \n
\subpage density_sine_wave_advection_2species \n

\page density_wave_advection_1d 1D Euler Equations - Density Wave Advection
Description: 
-------------------

Location: \b mhysa/Examples/SingleSpecies/1D_DensitySineWaveAdvection
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Euler equations (euler1d.h)

Domain: \f$0 \le x < 1\f$, "periodic" (#_PERIODIC_) boundaries

Initial Solution:
  \f$ \rho = \rho_\infty + \tilde{\rho} \sin(2\pi x)\f$,
  \f$ p = p_\infty, u = u_\infty\f$, where
  \f$\rho_\infty = u_\infty = 1, p_\infty = 1/\gamma, \tilde{\rho} = 0.1\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include SingleSpecies/1D_DensitySineWaveAdvection/solver.inp

\b boundary.inp
\include SingleSpecies/1D_DensitySineWaveAdvection/boundary.inp

\b physics.inp
\include SingleSpecies/1D_DensitySineWaveAdvection/physics.inp

\b lusolver.inp (optional)
\include SingleSpecies/3D_DensitySineWaveAdvection/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp (exact solution), 
compile and run the following code in the run directory:
\include SingleSpecies/1D_DensitySineWaveAdvection/aux/exact.c

Output:
-------
After running the code, there should be 21 solution output
files \b op_00000.dat, ..., \b op_00020.dat; the first one is
the initial solution, and the latter is the final solution.
All these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

The following animation shows the advection of the density wave:
@image html Solution_1DDensityWaveAdvection.gif

Since the exact solution is available at the final time, the numerical 
errors are calculated and reported on screen (see below) as well as in
the file \b errors.dat:
\include SingleSpecies/1D_DensitySineWaveAdvection/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include SingleSpecies/1D_DensitySineWaveAdvection/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

Expected screen output:
\include SingleSpecies/1D_DensitySineWaveAdvection/output.log

\page sod_shock_tube_component_rec 1D Euler Equations - Sod Shock Tube (Component-Wise Reconstruction)

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

\page density_sine_wave_advection 3D Navier-Stokes Equations - Density Sine Wave Advection

Location: \b mhysa/Examples/SingleSpecies/3D_DensitySineWaveAdvection
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)
                     By default, #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler equations 
                     are solved.)

Domain: \f$0 \le x,y,z < 1\f$, "periodic" (#_PERIODIC_) boundaries 
        everywhere.

Initial solution: 
  \f$ \rho = \rho_\infty + \tilde{\rho} \sin(2\pi x) \sin(2\pi y)\sin(2\pi z) \f$,
  \f$ p = p_\infty, u = u_\infty, v = v_\infty, w = w_\infty \f$, where
  \f$\rho_\infty = u_\infty = v_\infty = w_\infty = 1, p_\infty = 1/\gamma, \tilde{\rho} = 0.1\f$.

Other relevant parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include SingleSpecies/3D_DensitySineWaveAdvection/solver.inp

\b boundary.inp
\include SingleSpecies/3D_DensitySineWaveAdvection/boundary.inp

\b physics.inp
\include SingleSpecies/3D_DensitySineWaveAdvection/physics.inp

\b lusolver.inp (optional)
\include SingleSpecies/3D_DensitySineWaveAdvection/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp (exact solution), 
compile and run the following code in the run directory:
\include SingleSpecies/3D_DensitySineWaveAdvection/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 2 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are in the 
Tecplot 3D format.
(#HyPar::op_file_format is set to \a tecplot3d in \b solver.inp).

Since the exact solution is available at the final time, the numerical 
errors are calculated and reported on screen (see below) as well as in
the file \b errors.dat:
\include SingleSpecies/3D_DensitySineWaveAdvection/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include SingleSpecies/3D_DensitySineWaveAdvection/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

The following animation shows the advection of the density wave:
@image html Solution_3DNavStok_DensityWave.gif


Expected screen output:
\include SingleSpecies/3D_DensitySineWaveAdvection/output.log

\page density_sine_wave_advection_2species 3D Navier-Stokes Equations (Two Species) - Density Sine Wave Advection

Location: \b mhysa/Examples/MultiSpecies/3D_DensitySineWaveAdvection_2Species
          (This directory contains all the input files needed
          to run this case.)

This example is basically the 3D density wave advection problem with two 
non-reacting species. In essence, it is identical to the single-species 
density wave advection, but tests the multispecies code implementation.

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)
                     By default, #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler equations 
                     are solved.)

Domain: \f$0 \le x,y,z < 1\f$, "periodic" (#_PERIODIC_) boundaries 
        everywhere.

Initial solution: 
  \f$ \rho_1 = m_1 \rho, \rho_2 = m_2 \rho \f$, where \f$ \rho = \rho_\infty + \tilde{\rho} \sin(2\pi x) \sin(2\pi y) \sin(2\pi z) \f$, 
  and \f$ p = p_\infty, u = u_\infty, v = v_\infty, w = w_\infty \f$, where
  \f$\rho_\infty = u_\infty = v_\infty = w_\infty = 1, p_\infty = 1/\gamma, \tilde{\rho} = 0.1\f$.

where \f$m_1 = 0.3\f$, \f$m_2 = 0.7\f$ are the mass fractions of each species.

Other relevant parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/solver.inp

\b boundary.inp
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/boundary.inp

\b physics.inp
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/physics.inp

\b lusolver.inp (optional)
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp (exact solution), 
compile and run the following code in the run directory:
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 2 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are in the 
Tecplot 3D format.
(#HyPar::op_file_format is set to \a tecplot3d in \b solver.inp).

Since the exact solution is available at the final time, the numerical 
errors are calculated and reported on screen (see below) as well as in
the file \b errors.dat:
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

The following animation shows the advection of the density wave:
@image html Solution_3DNavStok_DensityWave_2Species.gif


Expected screen output:
\include MultiSpecies/3D_DensitySineWaveAdvection_2Species/output.log

\page riemann_case4_component_rec 2D Euler Equations - Riemann Problem Case 4 (Component-Wise Reconstruction)

Location: \b mhysa/Examples/SingleSpecies/2D_RiemannProblem_Case4_ComponentWiseRec
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)
                     By default, #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler equations 
                     are solved.)

Domain: \f$-0.5 \le x,y < 0.5\f$, "extrapolate" (#_EXTRAPOLATE_) boundaries along \f$x,y\f$;
        \f$0 \le z < \delta\f$, "periodic" (#_PERIODIC_) boundaries along \f$z\f$, where \f$\delta\f$ is an arbitrarily small number (this is how a 3D code is used to solve a 2D problem).

Initial solution: A 2D Riemann problem corresponding to "Case 4" in the reference below.

Other relevant parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Reference:
  + P. Lax and X.-D. Liu, Solution of two-dimensional Riemann problems of 
    gas dynamics by positive schemes, SIAM J Sci Comp 19 (1998), 319–340.

Numerical Method:
 + Spatial discretization (hyperbolic): Component-wise 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSP-RK3 (TimeRK(), #_RK_SSP3_)

Input files required:
--------------------
\b solver.inp:
\include SingleSpecies/2D_RiemannProblem_Case4_ComponentWiseRec/solver.inp

\b boundary.inp
\include SingleSpecies/2D_RiemannProblem_Case4_ComponentWiseRec/boundary.inp

\b physics.inp
\include SingleSpecies/2D_RiemannProblem_Case4_ComponentWiseRec/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include SingleSpecies/2D_RiemannProblem_Case4_ComponentWiseRec/aux/init.c

Output:
-------
After running the code, there should be 11 solution output
files \b op_00000.dat, \b op_00001.dat, ..., \b op_00010.dat; 
the first one is the initial solution, and the others are the 
solutions at subsequent times, with the last one being the final 
solution.
These files are in the Tecplot 3D format (#HyPar::op_file_format is
set to \a tecplot3d in \b solver.inp).

Final solution at t=0.25: The following figure shows the density and
is obtained by plotting \a op_00010.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_2DRiemannCase4_ComponentRec.png

Expected screen output:
\include SingleSpecies/2D_RiemannProblem_Case4_ComponentWiseRec/output.log

\page vortex_convection 2D Euler Equations - Isentropic Vortex Convection

Location: \b hypar/Examples/SingleSpecies/2D_IsentropicVortexConvection
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)
                     By default, #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler equations 
                     are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions; 
        \f$0 \le z < \delta\f$, "periodic" (#_PERIODIC_) boundaries along \f$z\f$, where \f$\delta\f$ is an arbitrarily small number (this is how a 3D code is used to solve a 2D problem).

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Other relevant parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include SingleSpecies/2D_IsentropicVortexConvection/solver.inp

\b boundary.inp
\include SingleSpecies/2D_IsentropicVortexConvection/boundary.inp

\b physics.inp
\include SingleSpecies/2D_IsentropicVortexConvection/physics.inp

\b weno.inp (optional)
\include SingleSpecies/2D_IsentropicVortexConvection/weno.inp

\b lusolver.inp (optional)
\include SingleSpecies/2D_IsentropicVortexConvection/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include SingleSpecies/2D_IsentropicVortexConvection/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      4 4 1

in \b solver.inp (i.e., 4 processors along \a x, 4 processors 
along \a y, 1 processor along \a z). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=20\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation (showing the density) shows the convection of the vortex:
@image html Solution_2DNavStokVortex.gif

Since the exact solution is available at the final time, the numerical 
errors are calculated and reported on screen (see below) as well as \b errors.dat:
\include SingleSpecies/2D_IsentropicVortexConvection/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include SingleSpecies/2D_IsentropicVortexConvection/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

Expected screen output:
\include SingleSpecies/2D_IsentropicVortexConvection/output.log

\page isotropic_turbulence 3D Navier-Stokes Equations - Isotropic Turbulence Decay

Location: \b hypar/Examples/SingleSpecies/3D_DNSIsotropicTurbulence
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 2\pi\f$, "periodic" (#_PERIODIC_) boundaries 
        everywhere.

Initial solution: Isotropic turbulent flow - The initial solution is 
specified in the Fourier space (with an energy
distribution similar to that of turbulent flow), and then transformed
to the physical space through an inverse transform.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 333.33\f$ (#NavierStokes3D::Re)
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes3D::Pr)
  + \f$M_\infty = 0.3\f$ (turbulence fluctuation Mach number) (#NavierStokes3D::Minf)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include SingleSpecies/3D_DNSIsotropicTurbulence/solver.inp

\b boundary.inp
\include SingleSpecies/3D_DNSIsotropicTurbulence/boundary.inp

\b physics.inp
\include SingleSpecies/3D_DNSIsotropicTurbulence/physics.inp

\b weno.inp (optional)
\include SingleSpecies/3D_DNSIsotropicTurbulence/weno.inp

\b lusolver.inp (optional)
\include SingleSpecies/3D_DNSIsotropicTurbulence/lusolver.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 init.c

(see the FFTW website on ways to install it).
\include SingleSpecies/3D_DNSIsotropicTurbulence/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4 4

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 64 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=5\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

To generate compute the energy spectrum from a given solution, compile 
and run the following code in the run directory. This code wants to read 
a file called \b op.bin, so make a symbolic link with that name pointing
to the solution file whose energy spectrum is to be computed. (If #HyPar::op_overwrite
is set to \a yes in \b solver.inp, then the solution file itself is called 
\b op.bin). Also note that the solution must be in binary format
(#HyPar::op_file_format must be \a binary in \b solver.inp). It will write out a 
ASCII text file \b spectrum.dat with two columns: \f$k\f$ and \f$E\left(k\right)\f$.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 fourier.c

(see the FFTW website on ways to install it). 
\include SingleSpecies/3D_DNSIsotropicTurbulence/aux/fourier.c

The following figure shows the initial and final (t=5) energy spectra:
@image html Solution_3DNavStok_IsoTurb_Spectrum.png

The following file computes the kinetic energy as a function of time
from the solution files. It writes out an ASCII text file \b energy.dat
with two colums: time and kinetic energy.
\include SingleSpecies/3D_DNSIsotropicTurbulence/aux/kineticenergy.c

The following figure shows the kinetic energy decay:
@image html Solution_3DNavStok_IsoTurb_Energy.png

The code \b mhysa/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b mhysa/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the remaining columns are the solution
components (\f$\rho, \rho u, \rho v, \rho w, e\f$).

The following figure shows the density iso-surfaces:
@image html Solution_3DNavStok_IsoTurb.png

Expected screen output:
\include SingleSpecies/3D_DNSIsotropicTurbulence/output.log

\page shock_cylinder 2D Inviscid Shock-Cylinder Interaction (Component-Wise Reconstruction)

Location: \b mhysa/Examples/SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h - by default
                     #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler
                     equations are solved.)

Domain: \f$-2.5 \le x \le 7.5\f$, \f$-5 \le y \le 5\f$
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b mhysa/Examples/STLGeometries/cylinder.stl)

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_ (with post-shock flow conditions)
  + xmax: Supersonic inflow #_SUPERSONIC_OUTFLOW_ (with pre-shock flow conditions)
  + ymin and ymax: Slip walls #_SLIP_WALL_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)

Reference:
+ O. Boiron, G. Chiavassa, R. Donat, A high-resolution penalization method for large Mach number 
  flows in the presence of obstacles, Computers & Fluids, 38 (2009), pp. 703-714, Section 5.1

Initial solution: see reference

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b mhysa/Examples/SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec/

\b solver.inp
\include SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec/solver.inp

\b boundary.inp
\include SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec/boundary.inp

\b physics.inp
\include SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b mhysa/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec/aux/init.c

Output:
-------

Note that \b iproc is set to 

      4 4 1

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 16 MPI ranks (or change \b iproc).

After running the code, there should be 81 solution files \b op_00000.bin 
(initial solution), 
\b op_00001.bin, ..., \b op_00080.bin (solution at t=2). 
Since #HyPar::op_overwrite is set to 
\a no in \b solver.inp, separate files are written for solutions at each
output time.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b mhysa/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b mhysa/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extracted slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The codes Extras/BinaryToTecplot.c and Extras/BinaryToText.c
can then be used (in the \b slices subdirectory) to convert the binary slice solution file
to Tecplot or plain text files.
\b Note that it needs the relevant \b solver.inp that can be created as follows:
+ Copy the original solver.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).

The following plot shows the density contours for the final solution (t=2):
@image html Solution_3DNavStokShockCyl_Density.png
and the following is the numerical Schlieren image (contour plot of \f$\|\nabla\rho\|_2\f$):
@image html Solution_3DNavStokShockCyl_Schlieren.png
@image html Solution_3DNavStokShockCyl_Schlieren2.png

Expected screen output:
\include SingleSpecies/2D_ShockCylinderInteraction_ComponentWiseRec/output.log

