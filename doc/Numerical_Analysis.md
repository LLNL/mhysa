Numerical Analysis
==================

The subpages of this page contains some results related to 
numerical analysis.

\subpage convergence_space_euler1d_single_species \n
\subpage convergence_time_euler1d_single_species \n
\subpage convergence_time_navstok2d_single_species \n
\n
\subpage spectral_resolution_euler1d_single_species \n

\page convergence_space_euler1d_single_species Convergence (Space) - 1D Euler (Single Species)

The spatial convergence of Mhysa is tested for two spatial reconstruction
schemes:
 + 5th order WENO (Interp1PrimFifthOrderWENO())
 + 5th order CRWENO (Interp1PrimFifthOrderCRWENO())

The 1D density wave advection problem is used since it is a smooth solution. The 
figure below shows a plot of the mesh size (dx) vs. the L2 norm of the error. Both
schemes are 5th order accurate, and shows 5th order convergence; however, the 
absolute errors for the CRWENO scheme are an order of magnitude lower than those
for the WENO scheme as predicted by Taylor series analysis.

\htmlonly <style>div.image img[src="SpatialConvergence_Euler1D.png"]{width:800px;}</style> \endhtmlonly 
@image html SpatialConvergence_Euler1D.png

The following MATLAB script (mhysa/Examples/ConvergenceTests/ConvergenceSpace_Euler1D.m) was
used to generate this plot:
\include ConvergenceTests/ConvergenceSpace_Euler1D.m

\page convergence_time_euler1d_single_species Convergence (Time) - 1D Euler (Single Species)

The time convergence of Mhysa is tested for several explicit and implicit
time integration schemes. The explicit schemes are:
 + RK2 (TimeRK(), #_RK_22_)
 + SSPRK3 (TimeRK(), #_RK_SSP3_)
 + RK4 (TimeRK(), #_RK_44_)

The implicit schemes are the 2nd, 3rd, and 4th order ARKIMEX schemes in PETSc 
(see TSARKIMEX in PETSc documentation). Although they are IMEX schemes, they 
are used in the fully-implicit mode. \b Note: Mhysa must be compiled with PETSc
to use these time integrators.

\htmlonly <style>div.image img[src="TimeConvergence_Euler1D.png"]{width:800px;}</style> \endhtmlonly 
@image html TimeConvergence_Euler1D.png

The 1D density wave advection problem is used since it is a smooth solution. The 
figure above shows a plot of the time step size (dt) vs. the L2 norm of the error.
The time integrators converge at their theoretical orders of convergence. Note the
following:
  + For this problem, dt=0.002 is unstable for #_RK_22_ and hence the error is large.
  + The 4th order methods yield round-off-level errors for the time steps considered
    here, and hence do not show a convergence.

The following MATLAB script (mhysa/Examples/ConvergenceTests/ConvergenceTime_Euler1D.m)
was used to generate the above plot:
\include ConvergenceTests/ConvergenceTime_Euler1D.m

\page convergence_time_navstok2d_single_species Convergence (Time) - 2D Euler (Single Species)

The time convergence of Mhysa is tested for several explicit and implicit
time integration schemes. The explicit schemes are:
 + RK2 (TimeRK(), #_RK_22_)
 + SSPRK3 (TimeRK(), #_RK_SSP3_)
 + RK4 (TimeRK(), #_RK_44_)

The implicit schemes are the 2nd, 3rd, and 4th order ARKIMEX schemes in PETSc 
(see TSARKIMEX in PETSc documentation). Although they are IMEX schemes, they 
are used in the fully-implicit mode. \b Note: Mhysa must be compiled with PETSc
to use these time integrators.

\htmlonly <style>div.image img[src="TimeConvergence_NS2D.png"]{width:800px;}</style> \endhtmlonly 
@image html TimeConvergence_NS2D.png

The 2D isentropic vortex convection problem is used since it is a smooth solution. The 
figure above shows a plot of the time step size (dt) vs. the L2 norm of the error.
The time integrators converge at their theoretical orders of convergence. 

The following MATLAB script 
(mhysa/Examples/ConvergenceTests/ConvergenceTime_NavierStokes2D_VortexConvection.m)
was used to generate the above plot:
\include ConvergenceTests/ConvergenceTime_NavierStokes2D_VortexConvection.m

\page spectral_resolution_euler1d_single_species Spectral Resolution - 1D Euler (Single Species)

The spectral resolution of Mhysa is studied here for the following two spatial
reconstruction schemes:
 + 5th order WENO (Interp1PrimFifthOrderWENO())
 + 5th order CRWENO (Interp1PrimFifthOrderCRWENO())

The problem consists of the advection of a density signal that is the sum of sine waves:
\f{equation}{
  \rho(x,t) = \rho_\infty + \sum_{k=1}^{N/2} \alpha\left(k\right)\cos\left[2\pi k\left(x-u_\infty t\right)+\phi\left(k\right)\right]
\f}
where \f$\alpha\left(k\right) = 0.01k^{-5/6}\f$, \f$\phi\left(k\right) = -\pi + 2\pi*{\rm random}\left[0,1\right]\f$, and \f$\rho_\infty=1, u_\infty = 1.0, p_\infty=1/\gamma\f$. 
This is simulated for 1 period over the periodic domain, and the Fourier transform of the initial and final solution is
obtained to compare the spectral resolution of the numerical method.

\htmlonly <style>div.image img[src="LinearSpectralAnalysis_Space.png"]{width:800px;}</style> \endhtmlonly 
@image html LinearSpectralAnalysis_Space.png

The figure above shows the spectral resolution for the two spatial reconstruction scheme. Though both
methods are 5th order accurate, the CRWENO schemes shows a better resolution than the WENO scheme
for the same number of grid points, and therefore, it is better suited for the simulation of turbulent
flows.

The following MATLAB script (mhysa/Examples/SpectralAnalysis/LinearSpectralAnalysis_Space_Euler1D.m) was
used to generate this plot:
\include SpectralAnalysis/LinearSpectralAnalysis_Space_Euler1D.m
