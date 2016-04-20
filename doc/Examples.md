Examples
========

\subpage basic_examples :
Some basic examples that are simulated using HyPar. They 
all use explicit time integration, and \b do \b not require HyPar to be compiled
with PETSc. Most of them can be run on one or a small number of processors.

\subpage petsc_examples : 
Some examples that use implicit or semi-implicit (IMEX) time
integration methods implemented in PETSc. To run them, HyPar needs to be compiled \b with \b PETSc.

\page basic_examples Basic Examples

The following are some basic examples that are simulated using HyPar. They 
all use explicit time integration, and \b do \b not require HyPar to be compiled
with PETSc.

\subpage linear_adv_sine \n
\subpage linear_adv_disc \n
\subpage linear_diff_sine 

\subpage sod_shock_tube  \n
\subpage lax_shock_tube \n
\subpage shu_osher \n
\subpage sod_shock_tube_wgrav

\subpage sw_dambreak

\subpage linear_adv_gauss \n
\subpage linear_diff_sine2d

\subpage euler2d_riemann4 \n
\subpage euler2d_riemann6 \n
\subpage euler2d_radexp \n
\subpage euler2d_vortex

\subpage euler2d_igwave \n
\subpage euler2d_rtb

\subpage navstok2d_flatplate

\subpage sw_circdambreak \n
\subpage sw_latbelt

\subpage linear_adv_3dgauss

\subpage ns3d_isoturb \n
\subpage ns3d_bubble

\subpage numa3d_bubble

\page linear_adv_sine 1D Linear Advection - Sine Wave

Location: \b hypar/Examples/1D/LinearAdvection/SineWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearAdvection/SineWave/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/SineWave/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/SineWave/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/SineWave/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/SineWave/aux/init.c

Output:
-------
After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,0.5,1: The following figure is obtained 
by plotting \a op_00000.dat (initial), \a op_00005.dat (t=0.5),
and \a op_00010.dat (final). 
@image html Solution_1DLinearAdvSine.png

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 1D/LinearAdvection/SineWave/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 1D/LinearAdvection/SineWave/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 1D/LinearAdvection/SineWave/output.log


\page linear_adv_disc 1D Linear Advection - Discontinuous Waves

Location: \b hypar/Examples/1D/LinearAdvection/DiscontinuousWaves
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$-1 \le x \le 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution:
  \f{equation}{
    u\left(x,0\right) = \left\{\begin{array}{lc} 
                          \exp\left(-\log\left(2\right)\frac{\left(x+7\right)^2}{0.0009}\right) & -0.8\le x \le -0.6 \\
                          1 & -0.4\le x \le -0.2 \\
                          1 - \left|10\left(x-0.1\right)\right| & 0\le x \le 0.2 \\
                          \sqrt{1-100\left(x-0.5\right)^2} & 0.4\le x \le 0.6 \\
                          0 & {\rm otherwise}
                        \end{array}\right.
  \f}

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearAdvection/DiscontinuousWaves/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/DiscontinuousWaves/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/DiscontinuousWaves/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/DiscontinuousWaves/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).

Final solution at t=2.0: The following figure is obtained 
by plotting \a op_00000.dat (initial) and \a op_00001.dat
(final). In both these files, the first column is grid 
index, the second column is x-coordinate, and the third 
column is the solution.
@image html Solution_1DLinearAdvDisc.png

Expected screen output:
\include 1D/LinearAdvection/DiscontinuousWaves/output.log

\page linear_diff_sine 1D Linear Diffusion - Sine Wave

Location: \b hypar/Examples/1D/LinearDiffusion/SineWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Diffusion Equation (linearadr.h)

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_) 
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical Method:
  + Spatial discretization (parabolic): 2nd order (Interp2PrimSecondOrder()),
                                        conservative (ParabolicFunctionCons1Stage())
  + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearDiffusion/SineWave/solver.inp

\b boundary.inp
\include 1D/LinearDiffusion/SineWave/boundary.inp

\b physics.inp
\include 1D/LinearDiffusion/SineWave/physics.inp

To generate \b initial.inp (initial solution) and 
\b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 1D/LinearDiffusion/SineWave/aux/exact.c

Output:
-------
After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=10\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,2,4,6,8,10: The following figure is obtained 
by plotting \a op_00000.dat (t=0, initial), \a op_00002.dat (t=2),
\a op_00004.dat (t=4), \a op_00006.dat (t=6), \a op_00008.dat 
(t=8) and \a op_00010.dat (t=10, final). 
@image html Solution_1DLinearDiffSine.png

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 1D/LinearDiffusion/SineWave/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 1D/LinearDiffusion/SineWave/output.log



\page sod_shock_tube 1D Euler Equations - Sod Shock Tube

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/SodShockTube 
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
\include 1D/Euler1D/SodShockTube/solver.inp

\b boundary.inp
\include 1D/Euler1D/SodShockTube/boundary.inp

\b physics.inp
\include 1D/Euler1D/SodShockTube/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/SodShockTube/aux/init.c

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
\include 1D/Euler1D/SodShockTube/output.log

\page lax_shock_tube 1D Euler Equations - Lax Shock Tube

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/LaxShockTube 
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
\include 1D/Euler1D/LaxShockTube/solver.inp

\b boundary.inp
\include 1D/Euler1D/LaxShockTube/boundary.inp

\b physics.inp
\include 1D/Euler1D/LaxShockTube/physics.inp

\b weno.inp (optional)
\include 1D/Euler1D/LaxShockTube/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/LaxShockTube/aux/init.c

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
\include 1D/Euler1D/LaxShockTube/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for 
each component.

Expected screen output:
\include 1D/Euler1D/LaxShockTube/output.log


\page shu_osher 1D Euler Equations - Shu-Osher Problem

Description: 
------------

Location: \b hypar/Examples/1D/Euler1D/ShuOsherProblem 
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
\include 1D/Euler1D/ShuOsherProblem/solver.inp

\b boundary.inp
\include 1D/Euler1D/ShuOsherProblem/boundary.inp

\b physics.inp
\include 1D/Euler1D/ShuOsherProblem/physics.inp

\b weno.inp (optional)
\include 1D/Euler1D/ShuOsherProblem/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/ShuOsherProblem/aux/init.c

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
\include 1D/Euler1D/ShuOsherProblem/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for 
each component.

Expected screen output:
\include 1D/Euler1D/ShuOsherProblem/output.log


\page sod_shock_tube_wgrav 1D Euler Equations - Sod Shock Tube with Gravitational Force

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/SodShockTubeWithGravity 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
    for the Gas Dynamics Equations Under Gravitational Fields",
    Journal of Scientific Computing, 54, 2013, pp. 645-662.

Domain: \f$0 \le x \le 1.0\f$, \a "slip-wall" (#_SLIP_WALL_) 
        boundary conditions (wall velocity is zero), with 
        \b uniform \b gravitational \b force \f$g=1\f$.

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, u = 0, p = 0.1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order CRWENO (Interp1PrimFifthOrderCRWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include 1D/Euler1D/SodShockTubeWithGravity/solver.inp

\b boundary.inp
\include 1D/Euler1D/SodShockTubeWithGravity/boundary.inp

\b physics.inp
\include 1D/Euler1D/SodShockTubeWithGravity/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/SodShockTubeWithGravity/aux/init.c

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
@image html Solution_1DSodShockTubeWithGravity.png

Expected screen output:
\include 1D/Euler1D/SodShockTubeWithGravity/output.log

\page sw_dambreak 1D Shallow Water Equations - Dam Breaking over Rectangular Bump

Location: \b hypar/Examples/1D/ShallowWater1D/DamBreakingRectangularBump
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Shallow Water Equations (shallowwater1d.h)

References:
  + Xing, Y., Shu, C.-W., "High order finite difference WENO
    schemes with the exact conservation property for the shallow
    water equations", Journal of Computational Physics, 208, 2005,
    pp. 206-227 (section 4.4).

Domain: \f$ 0 \le x \le 1500\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions

Initial solution:
\f{equation}{
  h\left(x\right) = \left\{\begin{array}{lc} 20 - b\left(x\right) & x <= 750 \\ 15 - b\left(x\right) & {\rm otherwise}\end{array}\right., u\left(x\right) = 0
\f}
where \f$b\left(x\right) = \left\{\begin{array}{lc} 8 & \left|x-750.0\right| \le 1500/8 \\  0 & {\rm  otherwise} \end{array}\right.\f$ is the bottom topography.

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

Note: in addition to the usual input files that HyPar needs, this
physical model needs the following input file(s):
+ \b topography.inp : file containing the bottom topography (same
  format as \b initial.inp).

\b solver.inp
\include 1D/ShallowWater1D/DamBreakingRectangularBump/solver.inp

\b boundary.inp
\include 1D/ShallowWater1D/DamBreakingRectangularBump/boundary.inp

\b physics.inp
\include 1D/ShallowWater1D/DamBreakingRectangularBump/physics.inp

\b weno.inp (optional)
\include 1D/ShallowWater1D/DamBreakingRectangularBump/weno.inp

To generate \b initial.inp and \b topography.inp, compile and run the 
following code in the run directory.
\include 1D/ShallowWater1D/DamBreakingRectangularBump/aux/init.c

Output:
-------
After running the code, there should be 5 solution output files
\b op_00000.dat, ..., \b op_00004.dat; the first one is the solution 
at \f$t=0\f$ and the final one is the solution at \f$t=60\f$. Since
#HyPar::op_overwrite is set to \a no in \b solver.inp, separate files 
are written for solutions at each output time. All the files are ASCII
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third and fourth columns are the two 
solution components.

In addition to the usual output files, the shallow water physics 
module writes out the following files:
+ \b topography_00000.dat, ..., \b topography_00004.dat: These files
  share the same format as the solution output files \b op_*.dat 
  and contains the topography \f$b\left(x\right)\f$.


Solutions at t=0,15,30,45,60: The following figure is obtained 
by plotting \a op_00000.dat (initial), \a op_00005.dat (t=0.5),
and \a op_00010.dat (final). \b Note: the figure plots 
\f$h\left(x\right)+b\left(x\right)\f$, i.e., it adds the third
column in \b op_*.dat and the third column in \b topography_*.dat.
@image html Solution_1DSWDamBreak.png

Expected screen output:
\include 1D/ShallowWater1D/DamBreakingRectangularBump/output.log


\page linear_adv_gauss 2D Linear Advection - Gaussian Pulse

Location: \b hypar/Examples/2D/LinearAdvection/GaussianPulse
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Linear Advection Equation (linearadr.h)

Domain: \f$-6 \le x < 6, -3 \le y < 3\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \exp\left[-\left(\frac{x^2}{2}+\frac{y^2}{2}\right)\right]\f$\n
Exact solution: \f$u\left(x,y,t\right) = u_0\left(x-a_xt,y-a_yt\right)\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/LinearAdvection/GaussianPulse/solver.inp

\b boundary.inp
\include 2D/LinearAdvection/GaussianPulse/boundary.inp

\b physics.inp (specifies \f$a_x\f$ and \f$a_y\f$)
\include 2D/LinearAdvection/GaussianPulse/physics.inp

\b lusolver.inp (optional)
\include 2D/LinearAdvection/GaussianPulse/lusolver.inp

\b weno.inp (optional)
\include 2D/LinearAdvection/GaussianPulse/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 2D/LinearAdvection/GaussianPulse/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 2

in \b solver.inp (i.e., 4 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 8 MPI ranks (or change \b iproc).

After running the code, there should be 21 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00020.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=12\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the solution files:
@image html Solution_2DLinearAdvGauss.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/LinearAdvection/GaussianPulse/errors.dat
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
\include 2D/LinearAdvection/GaussianPulse/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/LinearAdvection/GaussianPulse/output.log


\page linear_diff_sine2d 2D Linear Diffusion - Sine Wave

Location: \b hypar/Examples/2D/LinearDiffusion/SineWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Linear Diffusion Equation (linearadr.h)

Domain: \f$0 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \sin\left(2\pi x\right)\sin\left(2\pi y\right)\f$\n
Exact solution: \f$u\left(x,y,t\right) = \exp\left[-\pi^2 \left(4\nu_x + 4\nu_y\right) t\right] u0\left(x,y\right)\f$.

Numerical Method:
 + Spatial discretization (parabolic): 2nd order (Interp2PrimSecondOrder()), 
                                       conservative (ParabolicFunctionCons1Stage())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/LinearDiffusion/SineWave/solver.inp

\b boundary.inp
\include 2D/LinearDiffusion/SineWave/boundary.inp

\b physics.inp (specifies \f$\nu_x\f$ and \f$\nu_y\f$)
\include 2D/LinearDiffusion/SineWave/physics.inp

To generate \b initial.inp (initial solution) and 
\b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 2D/LinearDiffusion/SineWave/aux/exact.c

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
is the solution at \f$t=10\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the solution files:
@image html Solution_2DLinearDiffSine.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/LinearDiffusion/SineWave/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 2D/LinearDiffusion/SineWave/output.log


\page euler2d_riemann4 2D Euler Equations - Riemann Problem Case 4

Location: \b hypar/Examples/2D/NavierStokes2D/Riemann2DCase4
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
\include 2D/NavierStokes2D/Riemann2DCase4/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/Riemann2DCase4/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/Riemann2DCase4/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/Riemann2DCase4/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/Riemann2DCase4/aux/init.c

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
\include 2D/NavierStokes2D/Riemann2DCase4/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for each component.

Expected screen output:
\include 2D/NavierStokes2D/Riemann2DCase4/output.log


\page euler2d_riemann6 2D Euler Equations - Riemann Problem Case 6

Location: \b hypar/Examples/2D/NavierStokes2D/Riemann2DCase6
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
\include 2D/NavierStokes2D/Riemann2DCase6/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/Riemann2DCase6/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/Riemann2DCase6/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/Riemann2DCase6/aux/init.c

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
\include 2D/NavierStokes2D/Riemann2DCase6/output.log


\page euler2d_radexp 2D Euler Equations - Radial Expansion Wave

Location: \b hypar/Examples/2D/NavierStokes2D/RadialExpansionWave
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
\include 2D/NavierStokes2D/RadialExpansionWave/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RadialExpansionWave/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RadialExpansionWave/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/RadialExpansionWave/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/RadialExpansionWave/aux/init.c

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
\include 2D/NavierStokes2D/RadialExpansionWave/output.log


\page euler2d_vortex 2D Euler Equations - Isentropic Vortex Convection

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection
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
\include 2D/NavierStokes2D/InviscidVortexConvection/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection/weno.inp

\b lusolver.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include 2D/NavierStokes2D/InviscidVortexConvection/aux/exact.c

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
\include 2D/NavierStokes2D/InviscidVortexConvection/errors.dat
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
\include 2D/NavierStokes2D/InviscidVortexConvection/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection/output.log


\page euler2d_igwave 2D Euler Equations (with gravitational force) - Inertia-Gravity Waves

Location: \b hypar/Examples/2D/NavierStokes2D/InertiaGravityWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + W. C. Skamarock and J. B. Klemp, "Efficiency and accuracy of 
    the Klemp-Wilhelmson timesplitting technique", Monthly Weather 
    Review, 122 (1994), pp. 2623–2630.
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.1).

Domain: \f$0 \le x \le 300,000\,m, 0 \le y \le 10,000\,m\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions along \f$x\f$, "slip-wall" (#_SLIP_WALL_) boundary conditions along \f$y\f$.

Initial solution: See references above.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/InertiaGravityWave/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InertiaGravityWave/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InertiaGravityWave/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/InertiaGravityWave/weno.inp
\b Note: \a no_limiting is set to 1, i.e., since this is a 
smooth problem, WENO limiting is turned off, and the spatial
discretization uses a 5th order polynomial.

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/InertiaGravityWave/aux/init.c

Output:
-------
Note that \b iproc is set to 

      12 1

in \b solver.inp (i.e., 12 processors along \a x, and 1
processor along \a y). Thus, this example should be run
with 12 MPI ranks (or change \b iproc).

After running the code, there should be 26 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00025.bin; 
the first one is the solution at \f$t=0s\f$ and the final one
is the solution at \f$t=3000s\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/InertiaGravityWave/aux/PostProcess.c

The following plot shows the potential temperature perturbation
contours at the final time t=3000. It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_2DNavStokIGWave.png

If the postprocessing code above was used to write out files in 
text format, the following MATLAB script can be used to generate 
plots and visualize the solution:
\include 2D/NavierStokes2D/InertiaGravityWave/PlotSolution.m

Expected screen output:
\include 2D/NavierStokes2D/InertiaGravityWave/output.log



\page euler2d_rtb 2D Euler Equations (with gravitational force) - Rising Thermal Bubble

Location: \b hypar/Examples/2D/NavierStokes2D/RisingThermalBubble
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.2).

Domain: \f$0 \le x,y \le 1000\,m\f$, 
        "slip-wall" (#_SLIP_WALL_) boundary conditions on all sides.

Initial solution: See references above.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/RisingThermalBubble/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RisingThermalBubble/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RisingThermalBubble/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/RisingThermalBubble/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/RisingThermalBubble/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 3

in \b solver.inp (i.e., 4 processors along \a x, and 3
processor along \a y). Thus, this example should be run
with 12 MPI ranks (or change \b iproc).

After running the code, there should be 41 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00040.bin; 
the first one is the solution at \f$t=0s\f$ and the final one
is the solution at \f$t=700s\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/RisingThermalBubble/aux/PostProcess.c

The following plot shows the potential temperature perturbation
contours at the final time t=700s. It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_2DNavStokRTB.png

The following animation was created using all the transient files:
@image html Solution_2DNavStokRTB.gif

If the postprocessing code above was used to write out files in 
text format, the following MATLAB script can be used to generate 
plots and visualize the solution:
\include 2D/NavierStokes2D/RisingThermalBubble/PlotSolution.m

Expected screen output:
\include 2D/NavierStokes2D/RisingThermalBubble/output.log


\page navstok2d_flatplate 2D Navier-Stokes Equations -  Laminar Flow over Flat Plate

Location: \b hypar/Examples/2D/NavierStokes2D/FlatPlateLaminar/UniformGrid
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
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/aux/init.c

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
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file (works only for 2D and 3D).
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
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
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/aux/SkinFriction.c
The following figure showing the exact and computed skin friction coefficients
was obtained by plotting \b SkinFriction.dat:
@image html Solution_2DNavStokFlatPlateSkinFriction.png

Expected screen output:
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/output.log


\page sw_circdambreak 2D Shallow Water Equations - Circular Dam Break

Location: \b hypar/Examples/2D/ShallowWater2D/CircularDamBreak
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Shallow Water Equations (shallowwater2d.h)

References:
  + Delis, Katsaounis, "Numerical solution of the two-dimensional
    shallow water equations by the application of relaxation methods", 
    Applied Mathematical Modelling, 29 (2005), pp. 754-783 (Section 6.3).

Domain: \f$0 \le x,y \le 50 \f$, \a "extrapolate" (#_EXTRAPOLATE_) boundary
        conditions everywhere

Initial solution:
\f{equation}{
  h\left(x,y\right) = \left\{\begin{array}{cc}10 & r \le 11.0\\1 & {\rm otherwise}\end{array}\right., u\left(x,y\right) = 0.
\f}
with the bottom topography as \f$b\left(x,y\right) = 0\f$, and 
\f$r^2=\left(x-25\right)^2+\left(y-25\right)^2\f$.

Other parameters:
 + \f$g=9.8\f$ (#ShallowWater2D::g)

Numerical Method:
  + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
  + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

Note: in addition to the usual input files that HyPar needs, this
physical model needs the following input file(s):
+ \b topography.inp : file containing the bottom topography (same
  format as \b initial.inp).

However, this case has a constant (zero) bottom topography, and thus
this file can be skipped.

\b solver.inp
\include 2D/ShallowWater2D/CircularDamBreak/solver.inp

\b boundary.inp
\include 2D/ShallowWater2D/CircularDamBreak/boundary.inp

\b physics.inp
\include 2D/ShallowWater2D/CircularDamBreak/physics.inp

\b weno.inp (optional)
\include 2D/ShallowWater2D/CircularDamBreak/weno.inp

To generate \b initial.inp and \b topography.inp, compile and run the 
following code in the run directory (note: topography.inp can be skipped
since this case involves a uniform topography).
\include 2D/ShallowWater2D/CircularDamBreak/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 6 solution output files
\b op_00000.dat, ..., \b op_00005.dat; the first one is the solution 
at \f$t=0\f$ and the final one is the solution at \f$t=0.69\f$. Since
#HyPar::op_overwrite is set to \a no in \b solver.inp, separate files 
are written for solutions at each output time. 

#HyPar::op_file_format is set to \a text in \b solver.inp, and
thus, all the files are in plain text (ASCII) format. In these files, 
the data is written out as: the first two columns are grid indices, the 
next two columns are the x and y coordinates, and the remaining columns 
are the three solution components.
#HyPar::op_file_format can also be set to \a tecplot2d to get the solution
files in the Tecplot format (that can be read using any visualization 
software that supports Tecplot format).

In addition to the usual output files, the shallow water physics 
module writes out the following files:
+ \b topography_00000.dat, ..., \b topography_00005.dat: These files
  share the same format as the solution output files \b op_*.dat 
  and contains the topography \f$b\left(x\right)\f$ (for this example,
  they do not matter since the topography is flat).

The following plot shows the final solution (water height):
@image html Solution_2DShallowWater_CircDamBreak.png
It was obtained by using the following MATLAB code to plot \b op_00005.dat:
\include 2D/ShallowWater2D/CircularDamBreak/PlotSolution.m

Expected screen output:
\include 2D/ShallowWater2D/CircularDamBreak/output.log


\page sw_latbelt 2D Shallow Water Equations - Latitude Belt Flow

Location: \b hypar/Examples/2D/ShallowWater2D/LatitudeBeltFlow
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Shallow Water Equations (shallowwater2d.h)

References:
 + Zhu, Et. al., "Variational Data Assimilation with a Variable Resolution
   Finite-Element Shallow-Water Equations Model", Monthly Weather Review,
   122, 1994, pp. 946--965.\n
   Govering equations and approximation of Coriolis forces: Eqns (2.1)-(2.4)\n
   Problem definition: Eqns. (4.1)-(4.3)

Domain: \f$0 \le x \le 6,000,000\,{\rm m}, 0\le y \le 4,400,000\,{\rm m}\f$, 
        \a "periodic" (#_PERIODIC_) boundaries along \a x, 
        \a "slip wall" (#_SW_SLIP_WALL_) boundaries along \a y.

Initial solution: See equations (4.1) and (4.2) in reference,
                  constant (zero) bottom topography

Other parameters:
 + \f$g=10\f$ (#ShallowWater2D::g)
 + \f$\hat{f}=0.0001\f$ (#ShallowWater2D::fhat)
 + \f$\beta = 1.5\times 10^{-11}\f$ (#ShallowWater2D::beta)
 + \f$D=4400000\f$ (#ShallowWater2D::D)

Numerical Method:
  + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
  + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

Note: in addition to the usual input files that HyPar needs, this
physical model needs the following input file(s):
+ \b topography.inp : file containing the bottom topography (same
  format as \b initial.inp).

However, this case has a constant (zero) bottom topography, and thus
this file can be skipped.

\b solver.inp
\include 2D/ShallowWater2D/LatitudeBeltFlow/solver.inp

\b boundary.inp
\include 2D/ShallowWater2D/LatitudeBeltFlow/boundary.inp

\b physics.inp
\include 2D/ShallowWater2D/LatitudeBeltFlow/physics.inp

\b weno.inp (optional)
\include 2D/ShallowWater2D/LatitudeBeltFlow/weno.inp

To generate \b initial.inp and \b topography.inp, compile and run the 
following code in the run directory (note: topography.inp can be skipped
since this case involves a uniform topography).
\include 2D/ShallowWater2D/LatitudeBeltFlow/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 101 solution output files
\b op_00000.dat, ..., \b op_00100.dat; the first one is the solution 
at \f$t=0\f$ and the final one is the solution at \f$t=360000\f$. Since
#HyPar::op_overwrite is set to \a no in \b solver.inp, separate files 
are written for solutions at each output time. 

#HyPar::op_file_format is set to \a text in \b solver.inp, and
thus, all the files are in plain text (ASCII) format. In these files, 
the data is written out as: the first two columns are grid indices, the 
next two columns are the x and y coordinates, and the remaining columns 
are the three solution components.
#HyPar::op_file_format can also be set to \a tecplot2d to get the solution
files in the Tecplot format (that can be read using any visualization 
software that supports Tecplot format).

The following plot shows an animation of the solution: 
@image html Solution_2DShallowWater_LatitudeBelt.gif
It was obtained by using the following MATLAB code:
\include 2D/ShallowWater2D/LatitudeBeltFlow/PlotSolution.m

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/ShallowWater2D/LatitudeBeltFlow/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error of each component(#HyPar::ConservationError).
Note that that conservation error for the water height is zero to machine 
precision (the non-zero conservation "error" for the momentum components
is due to the non-zero source term - Coriolis force).

Expected screen output:
\include 2D/ShallowWater2D/LatitudeBeltFlow/output.log


\page linear_adv_3dgauss 3D Linear Advection - Gaussian Pulse

Location: \b hypar/Examples/3D/LinearAdvection/GaussianPulse
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D Linear Advection Equation (linearadr.h)

Domain: \f$-3 \le x,y,z < 3\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \exp\left[-\left(\frac{x^2}{2}+\frac{y^2}{2}+\frac{z^2}{2}\right)\right]\f$\n
Exact solution: \f$u\left(x,y,t\right) = u_0\left(x-a_xt,y-a_yt,z-a_zt\right)\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 3D/LinearAdvection/GaussianPulse/solver.inp

\b boundary.inp
\include 3D/LinearAdvection/GaussianPulse/boundary.inp

\b physics.inp (specifies \f$a_x\f$, \f$a_y\f$, and \f$a_z\f$)
\include 3D/LinearAdvection/GaussianPulse/physics.inp

\b lusolver.inp (optional)
\include 3D/LinearAdvection/GaussianPulse/lusolver.inp

\b weno.inp (optional)
\include 3D/LinearAdvection/GaussianPulse/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 3D/LinearAdvection/GaussianPulse/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, \a y, and \a z).
Thus, this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=6\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
All the files are binary (#HyPar::op_file_format is set to \a binary in \b solver.inp).
The code \b hypar/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b hypar/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the last column is the solution \a u.

The following animation was generated from the solution files 
(after converting to Tecplot format and plotting the iso-surface
in VisIt):
@image html Solution_3DLinearAdvGauss.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 3D/LinearAdvection/GaussianPulse/errors.dat
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
\include 3D/LinearAdvection/GaussianPulse/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 3D/LinearAdvection/GaussianPulse/output.log


\page ns3d_isoturb 3D Navier-Stokes Equations - Isotropic Turbulence Decay

Location: \b hypar/Examples/3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay
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
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/weno.inp

\b lusolver.inp (optional)
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/lusolver.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 init.c

(see the FFTW website on ways to install it).
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/init.c

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
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/fourier.c

The following figure shows the initial and final (t=5) energy spectra:
@image html Solution_3DNavStok_IsoTurb_Spectrum.png

The following file computes the kinetic energy as a function of time
from the solution files. It writes out an ASCII text file \b energy.dat
with two colums: time and kinetic energy.
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/kineticenergy.c

The following figure shows the kinetic energy decay:
@image html Solution_3DNavStok_IsoTurb_Energy.png

The code \b hypar/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b hypar/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the remaining columns are the solution
components (\f$\rho, \rho u, \rho v, \rho w, e\f$).

The following figure shows the density iso-surface colored by the internal energy
(plotted in ParaView after converting the binary solution to a Tecplot file):
@image html Solution_3DNavStok_IsoTurb.png

Expected screen output:
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/output.log



\page ns3d_bubble 3D Navier-Stokes Equations - Rising Thermal Bubble

Location: \b hypar/Examples/3D/NavierStokes3D/RisingThermalBubble
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "slip-wall" (#_SLIP_WALL_) boundaries 
        everywhere, with zero wall velocity.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

Initial solution: A warm bubble in cool ambient atmosphere. Note that in this example, the
gravitational forces and rising of the bubble is along the \a y-axis.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes3D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes3D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes3D::p0)
  + Hydrostatic balance type 2 (#NavierStokes3D::HB)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSP RK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 3D/NavierStokes3D/RisingThermalBubble/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/RisingThermalBubble/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/RisingThermalBubble/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/RisingThermalBubble/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/RisingThermalBubble/aux/init.c

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
is the solution at \f$t=200\,{\rm s}\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The following code (<B>Examples/3D/NavierStokes3D/RisingThermalBubble/aux/PostProcess.c</B>)
can be used to convert the binary solution file (with conserved variables 
\f$\rho,\rho u,\rho v,\rho w,e\f$) to Tecplot or plain text files with the primitive
and reference variables \f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$ where the
subscript \f$0\f$ indicates the hydrostatic mean value.
\include 3D/NavierStokes3D/RisingThermalBubble/aux/PostProcess.c

The following figure shows the density pertubation iso-surface for the initial
and final solutions (plotted in VisIt):
@image html Solution_3DNavStok_Bubble3D.png

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extract slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The following code (<B>Examples/3D/NavierStokes3D/RisingThermalBubble/aux/PostProcessSlice.c</B>)
can then be used (in the \b slices subdirectory) to convert the binary slice solution file (with conserved variables 
\f$\rho,\rho u,\rho v,\rho w,e\f$) to Tecplot or plain text files with the primitive
and reference variables \f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$.
\b Note that it needs the relevant \b solver.inp and \b physics.inp that can be created
as follows:
+ Copy the original solver.inp and physics.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).
+ In physics.inp, remove the component of \b gravity corresponding to the dimension perpendicular
  to the slice.

\include 3D/NavierStokes3D/RisingThermalBubble/aux/PostProcessSlice.c

The following figure shows the potential temperature \f$\theta\f$ along a slice at \f$z=500\,{\rm m}\f$ 
(plotted in VisIt):
@image html Solution_3DNavStok_Bubble.gif

Expected screen output:
\include 3D/NavierStokes3D/RisingThermalBubble/output.log



\page numa3d_bubble 3D NUMA (Nonhydrostatic Unified Model of the Atmosphere) Equations - Rising Thermal Bubble

Location: \b hypar/Examples/3D/NUMA/RisingThermalBubble
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D NUMA (Nonhydrostatic Unified Model of the Atmosphere) Equations (numa3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "numa-nfbc" (#_NO_FLUX_BC_) boundaries 
        everywhere.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

Initial solution: A warm bubble in cool ambient atmosphere.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#Numa3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#Numa3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ (#Numa3D::g)
  + Angular speed of earth \f$\Omega = 0\f$ (#Numa3D::Omega)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#Numa3D::Pref)
  + Reference temperature (at zero altitude) \f$T_{ref} = 300\f$ (#Numa3D::Tref)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 3D/NUMA/RisingThermalBubble/solver.inp

\b boundary.inp
\include 3D/NUMA/RisingThermalBubble/boundary.inp

\b physics.inp
\include 3D/NUMA/RisingThermalBubble/physics.inp

\b weno.inp (optional)
\include 3D/NUMA/RisingThermalBubble/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NUMA/RisingThermalBubble/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 2 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 41 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00040.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=200\,{\rm s}\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The code \b hypar/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b hypar/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the remaining columns are the solution
components (\f$\rho{'}, \rho u, \rho v, \rho w, \Theta{'}\f$).

The following figure shows the density pertubation iso-surface for the initial
and final solutions (plotted in VisIt):
@image html Solution_3DNUMA_Bubble3D.png

The code \b hypar/Extras/ExtractSlice.c can be used to extract a slice along at a specified location
along any of the dimensions (\b Note: make a subdirectory \a slices before using this).
It will write out 2D slice solutions in binary format in the \a slices subdirectory
with the same names as the original solution files. The codes \b hypar/Extras/BinaryToTecplot.c
and \b hypar/Extras/BinaryToText.c can then be used to convert these to text or Tecplot
formats for visualization. (\b Note: A copy of the original \a solver.inp is required in the \a slices subdirectory
with \b ndims reduced by 1, and \b size and \b iproc with the appropriate component removed such that
they have 2 components.)

The following figure shows the density pertubation along a slice at \f$y=500\,{\rm m}\f$ 
(plotted in VisIt):
@image html Solution_3DNUMA_Bubble.gif

Expected screen output:
\include 3D/NUMA/RisingThermalBubble/output.log



\page petsc_examples PETSc Examples

The following are some examples that use implicit or semi-implicit (IMEX) time
integration. To run them, HyPar needs to be compiled \b with \b PETSc.
\b Note that familiarity with using PETSc (https://www.mcs.anl.gov/petsc/) is assumed.
In general, any example or simulation can use PETSc time-integrators (assuming
HyPar is compiled with PETSc) by specifying the PETSc inputs through a 
<B>.petscrc</B> file. The file 
+ <B>hypar/Examples/PETScInputs/.petscrc_Example</B> 

is an example of a .petscrc file (with explanatory comments).


