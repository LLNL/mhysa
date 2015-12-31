/*! @file fppowersystem1bus.h
    @brief Fokker-Planck Model for 1-Bus Power System
    @author Debojyoti Ghosh

    Fokker-Planck Model for 1-Bus Power System

    \f{equation}{
      \frac {\partial p} {\partial t} 
      + \frac {\partial} {\partial x} \left[\mu\left(x,y\right)p\right]
      + \frac {\partial} {\partial y} \left[\nu\left(x,y\right)p\right]
      = D_{yx} \frac {\partial^2 p} {\partial y \partial x}
      + D_{yy} \frac {\partial^2 p} {\partial y^2}
    \f}
    where 
    \f{eqnarray}{
      \mu\left(x,y\right) &=& \omega_B\left(y-\omega_S\right) \\
      \nu\left(x,y\right) &=& \frac{\omega_S}{2H}\left[\left<P_m\right> - P_{\rm max} \sin\left(x\right) - D\left(y-\omega_S\right)\right] \\
      D_{yx} &=& \frac {\sigma^2\omega_S^2} {4H^2} \lambda^2 \omega_B \\
      D_{yy} &=& \frac {\sigma^2 \omega_S^2} {4H^2} \lambda \left( 1 - \lambda \frac {D \omega_S} {2H} \right)
    \f}
    
    Symbol                  | Name
    ----------------------- | ------------------------------------------------------------------------------
    \f$p\f$                 | probability
    \f$x\f$                 | angle between axis of generator and the magnetic field (\f$\theta\f$ in paper)
    \f$y\f$                 | generator angular speed (\f$\omega\f$ in paper)
    \f$t\f$                 | time
    \f$\omega_B\f$          | base speed
    \f$\omega_S\f$          | synchronization speed
    \f$H\f$                 | generator inertia
    \f$D\f$                 | damping factor
    \f$\left<P_m\right>\f$  | average power input
    \f$P_{\rm max}\f$       | maximum power
    \f$\sigma\f$            | square root of variance
    \f$\lambda\f$           | correlation time 

    Reference:
    + Wang, P., Barajas-Solano, D. A., Constantinescu, E. M., Abhyankar, S., 
      Ghosh, D., Smith, B. F., Huang, Z., Tartakovsky, A. M., "Probabilistic 
      Density Function Method for Stochastic ODEs of Power Systems with Uncertain 
      Power Input", SIAM/ASA Journal on Uncertainty Quantification, 3 (1), 2015, 
      pp. 873-896, http://dx.doi.org/10.1137/130940050.
*/

/*
  dp    d[mu(x,y)p]    d[nu(x,y)p]        d^2 p        d^2 p
  -- +  ----------- +  ----------- = D_yx ----- + D_yy -----
  dt        dx             dy             dy dx        dy^2

  p(x,y): probability
  x     : angle between axis of generator and the magnetic field (theta in paper)
  y     : generator angular speed (omega in paper)
  t     : time

  mu(x,y) = omegaB * (y-omegaS)
  nu(x,y) = (omegaS/2H)[Pm_avg - Pmax sin(x) - D(y-omegaS)]

             sigma^2 omegaS^2
  D_yx    =  ---------------- lambda^2 omegaB
                  4 H^2

             sigma^2 omegaS^2                      D omegaS
  D_yy    =  ---------------- lambda ( 1 - lambda ---------- )
                  4 H^2                              2 H

*/

#define _FP_POWER_SYSTEM_1BUS_  "fp-power-system-1bus"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 2
#define _MODEL_NVARS_ 1

/* define grid directions */
#define _XDIR_ 0
#define _YDIR_ 1

typedef struct fp_power_system_1bus_parameters {

  /* input parameters */
  double omegaB, omegaS, H, D, Pm_avg, Pmax, sigma, lambda;

  /* calculated, not an input */
  double pdf_integral;

} FPPowerSystem1Bus;

int FPPowerSystem1BusInitialize    (void*,void*);
int FPPowerSystem1BusCleanup       (void*);
