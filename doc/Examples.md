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

\subpage sod_shock_tube  \n
\subpage lax_shock_tube \n
\subpage shu_osher \n

\subpage euler2d_riemann4 \n
\subpage euler2d_riemann6 \n
\subpage euler2d_radexp \n
\subpage euler2d_vortex

\subpage navstok2d_ldsc \n
\subpage navstok2d_flatplate

\subpage ns3d_isoturb \n

\page sod_shock_tube 1D Euler Equations - Sod Shock Tube

Description: 
-------------------

Location: \b mhysa/Examples/Euler1D/SodShockTube 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

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
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include Euler1D/SodShockTube/solver.inp

\b boundary.inp
\include Euler1D/SodShockTube/boundary.inp

\b physics.inp
\include Euler1D/SodShockTube/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include Euler1D/SodShockTube/aux/init.c

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
@image html Solution_1DSodShockTube.png

Expected screen output:
\include Euler1D/SodShockTube/output.log

\page lax_shock_tube 1D Euler Equations - Lax Shock Tube

Description: 
-------------------

Location: \b mhysa/Examples/Euler1D/LaxShockTube 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + P.D. Lax, "Weak solutions of nonlinear hyperbolic
    equations and their numerical computation," Comm.
    Pure App. Math., 7, 159 (1954).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 0.445, \rho u = 0.311, e = 8.928\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.5, \rho u = 0, e = 1.4275\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include Euler1D/LaxShockTube/solver.inp

\b boundary.inp
\include Euler1D/LaxShockTube/boundary.inp

\b physics.inp
\include Euler1D/LaxShockTube/physics.inp

\b weno.inp (optional)
\include Euler1D/LaxShockTube/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include Euler1D/LaxShockTube/aux/init.c

Output:
-------
Note that \b iproc = 2 in \b solver.inp, so run this with
2 MPI ranks (or change \b iproc to 1). After running the code, 
there should be two solution output files \b op_00000.dat and 
\b op_00001.dat; the first one is the initial solution, and the 
latter is the final solution. Both these files are ASCII text 
(#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=0.08: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DLaxShockTube.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include Euler1D/LaxShockTube/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for 
each component.

Expected screen output:
\include Euler1D/LaxShockTube/output.log


\page shu_osher 1D Euler Equations - Shu-Osher Problem

Description: 
------------

Location: \b mhysa/Examples/Euler1D/ShuOsherProblem 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + C.-W. Shu and S. Osher, "Efficient implementation of
    essentially non-oscillatory schemes ,II," J. Comput.
    Phys., 83 (1989), pp. 32–78

Domain: \f$-5 \le x \le 5\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ -5 \le x < -4\f$: \f$\rho = 27/7, u = 4\sqrt{35}/7, p = 31/3\f$
  + \f$ -4 \le x \le 5\f$: \f$\rho = 1 + 0.2\sin\left(5x\right), u = 0, p = 1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include Euler1D/ShuOsherProblem/solver.inp

\b boundary.inp
\include Euler1D/ShuOsherProblem/boundary.inp

\b physics.inp
\include Euler1D/ShuOsherProblem/physics.inp

\b weno.inp (optional)
\include Euler1D/ShuOsherProblem/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include Euler1D/ShuOsherProblem/aux/init.c

Output:
-------
After running the code, 
there should be two solution output files \b op_00000.dat and 
\b op_00001.dat; the first one is the initial solution, and the 
latter is the final solution. Both these files are ASCII text 
(#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=1.8: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DShuOsherProblem.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include Euler1D/ShuOsherProblem/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for 
each component.

Expected screen output:
\include Euler1D/ShuOsherProblem/output.log


\page euler2d_riemann4 2D Euler Equations - Riemann Problem Case 4

Location: \b mhysa/Examples/NavierStokes2D/Riemann2DCase4
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + P. Lax and X.-D. Liu, "Solution of two-dimensional Riemann
    problems of gas dynamics by positive schemes," SIAM J Sci 
    Comp 19 (1998), 319–340.

Domain: \f$-0.5 \le x,y \le 0.5\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions.

Initial solution: see \b Case \b 4 in the reference.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include NavierStokes2D/Riemann2DCase4/solver.inp

\b boundary.inp
\include NavierStokes2D/Riemann2DCase4/boundary.inp

\b physics.inp
\include NavierStokes2D/Riemann2DCase4/physics.inp

\b weno.inp (optional)
\include NavierStokes2D/Riemann2DCase4/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include NavierStokes2D/Riemann2DCase4/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=0.25\f$. Since #HyPar::op_overwrite is
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

The following plot shows the density contours at the final time t=0.25, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokRiemann4.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include NavierStokes2D/Riemann2DCase4/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for each component.

Expected screen output:
\include NavierStokes2D/Riemann2DCase4/output.log


\page euler2d_riemann6 2D Euler Equations - Riemann Problem Case 6

Location: \b mhysa/Examples/NavierStokes2D/Riemann2DCase6
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + P. Lax and X.-D. Liu, "Solution of two-dimensional Riemann
    problems of gas dynamics by positive schemes," SIAM J Sci 
    Comp 19 (1998), 319–340.

Domain: \f$-0.5 \le x,y \le 0.5\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions.

Initial solution: see \b Case \b 6 in the reference.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 3rd order MUSCL (Interp1PrimThirdOrderMUSCLChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include NavierStokes2D/Riemann2DCase6/solver.inp

\b boundary.inp
\include NavierStokes2D/Riemann2DCase6/boundary.inp

\b physics.inp
\include NavierStokes2D/Riemann2DCase6/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include NavierStokes2D/Riemann2DCase6/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=0.3\f$. Since #HyPar::op_overwrite is
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

The following plot shows the density contours at the final time t=0.3, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokRiemann6.png

Expected screen output:
\include NavierStokes2D/Riemann2DCase6/output.log


\page euler2d_radexp 2D Euler Equations - Radial Expansion Wave

Location: \b mhysa/Examples/NavierStokes2D/RadialExpansionWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: http://www.as.dlr.de/hiocfd/case_c1.5.html

Domain: \f$-0.4 \le x,y \le 0.4\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions (supersonic outflow).

Initial solution: see reference above.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 3rd order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include NavierStokes2D/RadialExpansionWave/solver.inp

\b boundary.inp
\include NavierStokes2D/RadialExpansionWave/boundary.inp

\b physics.inp
\include NavierStokes2D/RadialExpansionWave/physics.inp

\b weno.inp (optional)
\include NavierStokes2D/RadialExpansionWave/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include NavierStokes2D/RadialExpansionWave/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
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

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokRadialExpansion.png

Expected screen output:
\include NavierStokes2D/RadialExpansionWave/output.log


\page euler2d_vortex 2D Euler Equations - Isentropic Vortex Convection

Location: \b mhysa/Examples/NavierStokes2D/InviscidVortexConvection
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

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

Numerical method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include NavierStokes2D/InviscidVortexConvection/solver.inp

\b boundary.inp
\include NavierStokes2D/InviscidVortexConvection/boundary.inp

\b physics.inp
\include NavierStokes2D/InviscidVortexConvection/physics.inp

\b weno.inp (optional)
\include NavierStokes2D/InviscidVortexConvection/weno.inp

\b lusolver.inp (optional)
\include NavierStokes2D/InviscidVortexConvection/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include NavierStokes2D/InviscidVortexConvection/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

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

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokVortex.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include NavierStokes2D/InviscidVortexConvection/errors.dat
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
\include NavierStokes2D/InviscidVortexConvection/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

Expected screen output:
\include NavierStokes2D/InviscidVortexConvection/output.log


\page navstok2d_ldsc 2D Navier-Stokes Equations -  Lid-Driven Square Cavity

Location: \b mhysa/Examples/NavierStokes2D/LidDrivenCavity
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., "Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers",
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., "High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method", Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100, 1000, 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: 4th order, 4-stage Runge-Kutta (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include NavierStokes2D/LidDrivenCavity/solver.inp

\b boundary.inp
\include NavierStokes2D/LidDrivenCavity/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include NavierStokes2D/LidDrivenCavity/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes2D/LidDrivenCavity/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processor along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.dat, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
Since #HyPar::op_file_format is set to \a tecplot2d in \b solver.inp,
this file is in the ASCII Tecplot format and can be viewed in any
software that supports this format (e.g. VisIt).

Following plots show the streamlines (colored by the velocity magnitude)
for the solutions with Reynolds number 100, 1000, and 3200, respectively.
@image html Solution_2DNavStokLDSC_Re0100.png
@image html Solution_2DNavStokLDSC_Re1000.png
@image html Solution_2DNavStokLDSC_Re3200.png

Expected screen output (for Reynolds number 3200):
\include NavierStokes2D/LidDrivenCavity/output.log


\page navstok2d_flatplate 2D Navier-Stokes Equations -  Laminar Flow over Flat Plate

Location: \b mhysa/Examples/NavierStokes2D/FlatPlateLaminar/UniformGrid
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h)

Reference: 
+ Hirsch, "Numerical Computation of Internal & External Flows",
  Volume 1 (Fundamentals of Computational Fluid Dynamics), 2nd 
  Edition, Elsevier, Section 12.3.2 (page 618-625).

Domain: \f$-0.25 \le x \le 1\f$, \f$0 \le y \le 0.25\f$

Boundary conditions:
+ Symmetry BC on \f$y=0, -0.25 \le x < 0\f$ (imposed through "slip-wall" 
  #_SLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0
  wall velocity).
+ Subsonic outflow on \f$y=0.25\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.
+ Subsonic inflow on \f$x=0\f$ (#_SUBSONIC_INFLOW_) with density \f$\rho=1\f$,
  and velocity \f$(u,v) = (0.3,0)\f$.
+ Subsonic outflow on \f$x=1\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.

Initial solution: Uniform flow with \f$\rho=1, u=0.3, v=0, p=1/\gamma\f$.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100,000\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the plate length .
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.3\f$ (freestream Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

This is a steady state problem - the solution converges to a steady laminar flow over a 
flat plate, and the residuals decrease with time. The skin friction coefficient is given by
\f{equation}{
  c_f = \frac {0.664} {\sqrt{Re_x}}, Re_x = \frac {\rho u x} {\mu}
\f}
(Blasius solution).

Input files required:
---------------------

\b solver.inp
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/solver.inp

\b boundary.inp
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/boundary.inp

\b physics.inp
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/physics.inp

\b weno.inp (optional)
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 4

in \b solver.inp (i.e., 2 processors along \a x, and 4
processor along \a y). Thus, this example should be run
with 8 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b mhysa/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file (works only for 2D and 3D).
+ \b mhysa/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following plot was obtained by converting the binary file to the 
Tecplot format, and using VisIt to plot it (it shows the density 
and velocity):
@image html Solution_2DNavStokFlatPlate.png
The following plot shows a magnified view of the boundary layer:
@image html Solution_2DNavStokFlatPlateMagnified.png

The following file computes the skin friction as a function of the 
Reynolds number and writes to to a text file \b SkinFriction.dat with 4 
columns: Reynolds number, computed skin friction coefficient, exact skin 
friction coefficient (\f$0.664/\sqrt{Re_x}\f$), and normal velocity gradient 
on the plate surface (\f$\left.\partial u/\partial y\right|_{y=0}\f$). 
Compile and run it in the run directory.
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/aux/SkinFriction.c
The following figure showing the exact and computed skin friction coefficients
was obtained by plotting \b SkinFriction.dat:
@image html Solution_2DNavStokFlatPlateSkinFriction.png

Expected screen output:
\include NavierStokes2D/FlatPlateLaminar/UniformGrid/output.log


\page ns3d_isoturb 3D Navier-Stokes Equations - Isotropic Turbulence Decay

Location: \b mhysa/Examples/NavierStokes3D/DNS_IsotropicTurbulenceDecay
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

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
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/solver.inp

\b boundary.inp
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/boundary.inp

\b physics.inp
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/physics.inp

\b weno.inp (optional)
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/weno.inp

\b lusolver.inp (optional)
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/lusolver.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 init.c

(see the FFTW website on ways to install it).
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2 1

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 4 MPI ranks (or change \b iproc).

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
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/fourier.c

The following figure shows the initial and final (t=5) energy spectra:
@image html Solution_3DNavStok_IsoTurb_Spectrum.png

The following file computes the kinetic energy as a function of time
from the solution files. It writes out an ASCII text file \b energy.dat
with two colums: time and kinetic energy.
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/kineticenergy.c

The following figure shows the kinetic energy decay:
@image html Solution_3DNavStok_IsoTurb_Energy.png

The code \b mhysa/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b mhysa/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the remaining columns are the solution
components (\f$\rho, \rho u, \rho v, \rho w, e\f$).

The following figure shows the density iso-surface colored by the internal energy
(plotted in ParaView after converting the binary solution to a Tecplot file):
@image html Solution_3DNavStok_IsoTurb.png

Expected screen output:
\include NavierStokes3D/DNS_IsotropicTurbulenceDecay/output.log



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
\subpage euler2d_vortex_petsc \n
\subpage navstok2d_flatplate_petsc

Implicit-Explicit (IMEX) time integration:
------------------------------------------
\subpage euler2d_low_mach_vortex_petsc \n
\subpage ns2d_ldsc_petsc_imex \n


\page navstok2d_flatplate_petsc 2D Navier-Stokes Equations -  Laminar Flow over Flat Plate

Location: \b mhysa/Examples/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h)

Reference: 
+ Hirsch, "Numerical Computation of Internal & External Flows",
  Volume 1 (Fundamentals of Computational Fluid Dynamics), 2nd 
  Edition, Elsevier, Section 12.3.2 (page 618-625).

Domain: \f$-0.25 \le x \le 1\f$, \f$0 \le y \le 0.25\f$

Boundary conditions:
+ Symmetry BC on \f$y=0, -0.25 \le x < 0\f$ (imposed through "slip-wall" 
  #_SLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0
  wall velocity).
+ Subsonic outflow on \f$y=0.25\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.
+ Subsonic inflow on \f$x=0\f$ (#_SUBSONIC_INFLOW_) with density \f$\rho=1\f$,
  and velocity \f$(u,v) = (0.3,0)\f$.
+ Subsonic outflow on \f$x=1\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.

Initial solution: Uniform flow with \f$\rho=1, u=0.3, v=0, p=1/\gamma\f$.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100,000\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the plate length .
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.3\f$ (freestream Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX.html) -
     Although the ARK methods are semi-implicit (IMEX), here they are used in the "fully implicit" mode, i.e., the implicit 
     method is used to solve the complete equation (Note the flag \b -ts_arkimex_fully_implicit in <B>.petscrc</B>).
   - Specific method: <B>ARK2e</B> (TSARKIMEX2E - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX2E.html)

This is a steady state problem - the solution converges to a steady laminar flow over a 
flat plate, and the residuals decrease with time. The skin friction coefficient is given by
\f{equation}{
  c_f = \frac {0.664} {\sqrt{Re_x}}, Re_x = \frac {\rho u x} {\mu}
\f}
(Blasius solution).

Input files required:
---------------------

<B>.petscrc</B>
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/petscrc

\b solver.inp
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/solver.inp

\b boundary.inp
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/boundary.inp

\b physics.inp
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/physics.inp

\b weno.inp (optional)
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 3

in \b solver.inp (i.e., 2 processor along \a x, and 3
processor along \a y). Thus, this example should be run
with 6 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b mhysa/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file (works only for 2D and 3D).
+ \b mhysa/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following plot was obtained by converting the binary file to the 
Tecplot format, and using VisIt to plot it (it shows the density 
and velocity):
@image html Solution_2DNavStokFlatPlate.png
The following plot shows a magnified view of the boundary layer:
@image html Solution_2DNavStokFlatPlateMagnified.png

The following file computes the skin friction as a function of the 
Reynolds number and writes to to a text file \b SkinFriction.dat with 4 
columns: Reynolds number, computed skin friction coefficient, exact skin 
friction coefficient (\f$0.664/\sqrt{Re_x}\f$), and normal velocity gradient 
on the plate surface (\f$\left.\partial u/\partial y\right|_{y=0}\f$). 
Compile and run it in the run directory.
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/aux/SkinFriction.c
The following figure showing the exact and computed skin friction coefficients
was obtained by plotting \b SkinFriction.dat:
@image html Solution_2DNavStokFlatPlateSkinFrictionPETSc.png

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/output.log

\page ns2d_ldsc_petsc_imex 2D Navier-Stokes Equations -  Lid-Driven Square Cavity

Location: \b mhysa/Examples/NavierStokes2D/LidDrivenCavity_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., "Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers",
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., "High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method", Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX4.html)

Input files required:
---------------------

<B>.petscrc</B>
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/petscrc

\b solver.inp
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/solver.inp

\b boundary.inp
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processor along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.dat, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
Since #HyPar::op_file_format is set to \a tecplot2d in \b solver.inp,
this file is in the ASCII Tecplot format and can be viewed in any
software that supports this format (e.g. VisIt).

The following plot shows the streamlines (colored by the velocity magnitude):
@image html Solution_2DNavStokLDSC_Re3200_PETSc_IMEX.png

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output (for Reynolds number 3200):
\include NavierStokes2D/LidDrivenCavity_PETSc_IMEX/output.log




\page euler2d_vortex_petsc 2D Euler Equations - Isentropic Vortex Convection

Location: \b mhysa/Examples/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.5,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: PETSc (SolvePETSc()) - <B>Crank-Nicholson</B> (TSCN - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSCN.html)

Input files required:
---------------------

<B>.petscrc</B>
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/petscrc

\b solver.inp
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/solver.inp

\b boundary.inp
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/boundary.inp

\b physics.inp
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/physics.inp

\b lusolver.inp (optional)
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

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

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokVortexPETSc.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/errors.dat
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
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.
\b Note that the conservation error depends on the accuracy with which the 
implicit systems are solved. Since the systems are solved with relaxed tolerances
here (see \a ksp_atol, \a ksp_rtol, \a snes_atol, \a snes_rtol in <B>.petscrc</B>),
the conservation error is zero within round-off.

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/output.log


\page euler2d_low_mach_vortex_petsc 2D Euler Equations - Low-Mach Isentropic Vortex Convection

Location: \b mhysa/Examples/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
           Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
           Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
           http://dx.doi.org/10.1137/15M1044369.

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the above reference.

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

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

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX4.html)

Input files required:
---------------------

<B>.petscrc</B>
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/petscrc

\b solver.inp
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/solver.inp

\b boundary.inp
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/boundary.inp

\b physics.inp
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/physics.inp

To generate \b initial.inp (initial solution) and \b exact.inp (exact solution), 
compile and run the following code in the run 
directory:
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=100\f$. Since #HyPar::op_overwrite is
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

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokLowMachVortexPETSc.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/errors.dat
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
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.
\b Note that the conservation error depends on the accuracy with which the 
implicit systems are solved (see \a ksp_atol, \a ksp_rtol, \a snes_atol, 
\a snes_rtol in <B>.petscrc</B>).

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/output.log


\page ib_examples Immersed Boundaries Examples

The following are some examples are use the immersed boundary methodology to solve the PDE in the presence of
various geometries. To use the immersed boundary implementation in HyPar, an STL (https://en.wikipedia.org/wiki/STL_%28file_format%29) 
representation of the immersed body is necessary. Note:
+ The immersed boundary method is \b only implemented for 3-dimensional simulations (#HyPar::ndims = 3).
+ It can be used with only those physical models that define an immersed boundary implementation (#HyPar::IBFunction()), for
  example, the 3D Navier-Stokes equations (NavierStokes3DImmersedBoundary()).

3D Navier-Stokes Equations:
---------------------------

\subpage ns3d_cylinder_steady_incompressible_viscous \n
\subpage ns3d_cylinder_unsteady_incompressible_viscous

\subpage ns3d_shock_cylinder_interaction

\subpage ns3d_sphere_steady_incompressible_viscous

\page ns3d_cylinder_steady_incompressible_viscous Steady, incompressible, viscous flow around a cylinder

Location: \b mhysa/Examples/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-2,6] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b mhysa/Examples/STLGeometries/cylinder.stl)

The following images shows the grid and the cylinder:
@image html Domain3D_Cylinder.png
@image html Domain2D_Cylinder.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)

Reference:
  + Taneda, S., "Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers," Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 
  + Dennis, S. C. R., Chang, G.-Z., "Numerical solutions for steady flow past a circular
    cylinder at Reynolds numbers up to 100", Journal of Fluid Mechanics, 42 (3), 1970,
    pp. 471-489.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 10,15,20\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 20,30,40\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b mhysa/Examples/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/

\b solver.inp
\include NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/solver.inp

\b boundary.inp
\include NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 10 (corresponding to \f$Re_D\f$ of 20). To try other Reynolds 
numbers, change it here.
\include NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b mhysa/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/aux/init.c

Output:
-------

Note that \b iproc is set to 

      2 2 1

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 4 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
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

The following figure shows the flow (density and streamlines) at \f$Re_D=20\f$:
@image html Solution_3DNavStokCylinder_ReD020.png
The parameter \a Re can be changed to 20 in \a physics.inp to run this simulation for
\f$Re_D=40\f$, and following figure shows the solution:
@image html Solution_3DNavStokCylinder_ReD040.png

The following plot shows the wake length (\f$L/D\f$) as a function of the Reynolds
number (\f$Re_D\f$) for the computed solutions and experimental results reported
in the reference above:
@image html Solution_3DNavStokCylinder_WL.png

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface.dat (if #HyPar::op_overwrite
is "yes") or \a surface_nnnnn.dat (if #HyPar::op_overwrite is "no", "nnnnn" is a numerical
index) (in this example, the file \b surface.dat is written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. The following image shows the surface pressure on the cylinder (front-view):
@image html IBSurface_3DNavStokCylinder.png
Since this is a 2D simulation, the value of the surface pressure on the end-surfaces of the 
cylinder (zmin and zmax) are not physically relevant.

Expected screen output (for \f$Re_D = 20\f$):
\include NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/output.log

  
\page ns3d_cylinder_unsteady_incompressible_viscous Unsteady, incompressible, viscous flow around a cylinder (vortex shedding)

Location: \b mhysa/Examples/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-4,12] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b mhysa/Examples/STLGeometries/cylinder.stl)

The following images shows the grid and the cylinder:
@image html Domain3D_Cylinder2.png
@image html Domain2D_Cylinder2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)

Reference:
  + Taneda, S., "Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers," Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 50\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 100\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b mhysa/Examples/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/

\b solver.inp
\include NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/solver.inp

\b boundary.inp
\include NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/boundary.inp

\b physics.inp
\include NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b mhysa/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/aux/init.c

Output:
-------

Note that \b iproc is set to 

      4 2 1

in \b solver.inp (i.e., 4 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 401 solution files \b op_00000.bin, 
\b op_00001.bin, ..., \b op_00400.bin. Since #HyPar::op_overwrite is set to 
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

The following animation shows the vorticity magnitude and has been created by
plotting op_00301.bin through op_00400.bin, and shows the vortex shedding:
@image html Solution_3DNavStokCylinder_Shedding.gif

Expected screen output:
\include NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/output.log

  
\page ns3d_shock_cylinder_interaction Inviscid Shock-Cylinder Interaction 

Location: \b mhysa/Examples/NavierStokes3D/2D_Shock_Cylinder_Interaction

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
  + ymin and ymax: Slip walls $_SLIP_WALL_
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

These files are all located in: \b mhysa/Examples/NavierStokes3D/2D_Shock_Cylinder_Interaction/

\b solver.inp
\include NavierStokes3D/2D_Shock_Cylinder_Interaction/solver.inp

\b boundary.inp
\include NavierStokes3D/2D_Shock_Cylinder_Interaction/boundary.inp

\b physics.inp
\include NavierStokes3D/2D_Shock_Cylinder_Interaction/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b mhysa/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes3D/2D_Shock_Cylinder_Interaction/aux/init.c

Output:
-------

Note that \b iproc is set to 

      2 2 1

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 4 MPI ranks (or change \b iproc).

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
\include NavierStokes3D/2D_Shock_Cylinder_Interaction/output.log

  
\page ns3d_sphere_steady_incompressible_viscous Steady, incompressible, viscous flow around a sphere

Location: \b mhysa/Examples/NavierStokes3D/Sphere/Steady_Viscous_Incompressible

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the sphere defined by [-2,6] X [-2,2] X [-2,2],
        and a stretched grid beyond this zone.

Geometry: A sphere of radius 0.5 centered at (0,0)
          (\b mhysa/Examples/STLGeometries/sphere.stl)

The following image shows the sphere:
@image html Surface3D_Sphere.png

The following images shows the grid and the sphere:
@image html Domain3D_Sphere1.png
@image html Domain3D_Sphere2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_

Reference:
  + Taneda, S., “Experimental Investigation of Wake behind a Sphere at Low Reynolds Numbers,” 
    Journal of the Physical Society of Japan, Vol. 11, 1956.
  + Johnson, T.A. and Patel, V.C., “Flow Past a Sphere up to a Reynolds Number of 300,” 
    Journal of Fluid Mechanics, Vol. 378, 1999.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 100\f$ (#NavierStokes3D::Re) 
    (\b Note: since the diameter of the sphere is 1.0, the diameter-based Reynolds number 
    is the same as the specified Reynolds number \f$Re_D = Re = 100\f$).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b mhysa/Examples/NavierStokes3D/Sphere/Steady_Viscous_Incompressible/

\b solver.inp
\include NavierStokes3D/Sphere/Steady_Viscous_Incompressible/solver.inp

\b boundary.inp
\include NavierStokes3D/Sphere/Steady_Viscous_Incompressible/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 100. To try other Reynolds numbers, change it here.
\include NavierStokes3D/Sphere/Steady_Viscous_Incompressible/physics.inp

\b sphere.stl : the filename "sphere.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b mhysa/Examples/STLGeometries/sphere.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include NavierStokes3D/Sphere/Steady_Viscous_Incompressible/aux/init.c

Output:
-------

Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 2 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b mhysa/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b mhysa/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following figure shows the flow (pressure and streamlines) at \f$Re_D=100\f$:
@image html Solution_3DNavStokSphere_ReD100.png

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface.dat (if #HyPar::op_overwrite
is "yes") or \a surface_nnnnn.dat (if #HyPar::op_overwrite is "no", "nnnnn" is a numerical
index) (in this example, the file \b surface.dat is written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. The following image shows the surface pressure on the sphere (front-view):
@image html IBSurface_3DNavStokSphere.png

Expected screen output:
\include NavierStokes3D/Sphere/Steady_Viscous_Incompressible/output.log
