Examples
========

\section sod_shock_tube Sod Shock Tube

Description: 
-------------------
Governing equations: 1D Euler equations (euler1d.h)

Reference: G.A. Sod, "A survey of several finite difference methods 
           for systems of nonlinear hyperbolic conservation laws," 
           J. Comput. Phys., 27, 1 (1978).

Domain: 0 <= x <= 1.0, "extrapolate" boundary conditions

Initial Solution:
+ 0 <= x < 0.5: rho = 1.0, u = 0, p = 1.0
+ 0.5 <= x <= 1.0: rho = 0.125, u = 0, p = 0.1

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
Final solution at t=0.2:
@image html Solution_1DSodShockTube.png

Expected screen output:
\include 1D/Euler1D/SodShockTube/output.log
