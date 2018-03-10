Scalability
===========

The subpages of this page contains some results related to scalability studies.

\subpage scalability_llnl_cab_4096\n
\subpage scalability_llnl_quartz_32768\n

\page scalability_llnl_cab_4096 Scalability (till 4096 processors) - Isotropic Turbulence Decay (3D, Single Species)

The following scaling results were obtained by simulating the isotropic turbulence decay case
(single species).

Hardware details:
  + CPU: Intel Xeon E5-2670
  + CPU Clock Speed: 2.6 GHz
  + Cores/Nodes: 16
  + Memory/Node: 32 GB
  + Interconnect: IBM QDR

Compiler details:
  + GCC version 4.7.4
  + MVAPICH2 version 2.2

Numerical methods details:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
 + Time integration: RK4 (TimeRK(), #_RK_44_)

The measured wall times are only for the solver, and does not include the wall time for
initialization, conclusion, and file I/O. These costs were negligible for the simulations
reported below.

Strong Scaling
--------------

A grid with \f$128^3\f$ (\f$\approx\f$ 2 million) points was used, and the problem
was simulated with \f$1\f$ to 4096 (\f$16^3\f$) MPI ranks (thus resulting in the number
of grid points per MPI rank varying from \f$\approx\f$ 2 million to 512). The wall time
for 100 iterations, as well as the parallel efficiency, are shown in the figures below.

\htmlonly <style>div.image img[src="mhysa_strong_scaling_walltime.png"]{width:800px;}</style> \endhtmlonly 
@image html mhysa_strong_scaling_walltime.png
\htmlonly <style>div.image img[src="mhysa_strong_scaling_efficiency.png"]{width:800px;}</style> \endhtmlonly 
@image html mhysa_strong_scaling_efficiency.png

Weak Scaling
------------

A constant number of grid points per MPI rank of 4096 (\f$16^3\f$) was used for this study, and the number
of MPI ranks varied from 1 to 4096. Consequently, the total size of the grid varied from $16^3 = 4096$ points
to \f$256^3 \approx\f$ 16.8 million points. The wall time for 1000 iterations is shown.

\htmlonly <style>div.image img[src="mhysa_weak_scaling_walltime.png"]{width:800px;}</style> \endhtmlonly 
@image html mhysa_weak_scaling_walltime.png

\page scalability_llnl_quartz_32768 Scalability (till 32768 processors) - Isotropic Turbulence Decay (3D, Single Species)

The following scaling results were obtained by simulating the isotropic turbulence decay case
(single species).

Hardware details:
  + CPU: Intel Xeon E5-2695
  + CPU Clock Speed: 2.1 GHz
  + Cores/Nodes: 36
  + Memory/Node: 128 GB
  + Interconnect: Omni-Path

Compiler details:
  + GCC version 4.9.3
  + MVAPICH2 version 2.2

Numerical methods details:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
 + Time integration: RK4 (TimeRK(), #_RK_44_)

The measured wall times are only for the solver, and does not include the wall time for
initialization, conclusion, and file I/O. These costs were negligible for the simulations
reported below.

Strong Scaling
--------------

A grid with \f$512^3\f$ (\f$\approx\f$ 134 million) points was used, and the problem
was simulated with \f$128\f$ to 32,768 MPI ranks (thus resulting in the number
of grid points per MPI rank varying from \f$\approx\f$ 1 million to 4096). The wall time
for 20 iterations, as well as the parallel efficiency, are shown in the figures below.

\htmlonly <style>div.image img[src="mhysa_strong_scaling_big_walltime.png"]{width:800px;}</style> \endhtmlonly 
@image html mhysa_strong_scaling_big_walltime.png
\htmlonly <style>div.image img[src="mhysa_strong_scaling_big_efficiency.png"]{width:800px;}</style> \endhtmlonly 
@image html mhysa_strong_scaling_big_efficiency.png
