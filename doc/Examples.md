
\subpage linear_adv_sine \n
\subpage linear_adv_disc \n
\subpage linear_diff_sine 

\subpage sod_shock_tube  \n
\subpage lax_shock_tube \n
\subpage shu_osher \n
\subpage sod_shock_tube_wgrav

\subpage sw_dambreak

\subpage linear_adv_gauss \n

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
After running the code, there should be two 11 output
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
After running the code, there should be two 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
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
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
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
+ \b topography_00000.dat, ..., \b topography_00001.dat: These files
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

\b physics.inp
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

After running the code, there should be two 21 output
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


