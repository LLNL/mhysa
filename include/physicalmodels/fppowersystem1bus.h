/*

  Fokker-Planck Model for a 1-Bus Power System

Reference: 
+ Wang P., Tartakovsky A.M., Abhyankar S., Smith B.F.,
  Huang Z., "Probabilistic Density Function Method for
  Stochastic ODEs of Power Systems with Uncertain 
  Power Input", Preprint

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
                  
  Physical Parameters:

    omegaB  : base speed
    omegaS  : synchronization speed
    H       : generator intertia
    D       : damping factor
    Pm_avg  : average power input
    Pmax    : EV/X where E is internal voltage, 
                         V is bus voltage,
                         X is total system reactance
    sigma   : square root of variance
    lambda  : correlation time

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

} FPPowerSystem1Bus;

int FPPowerSystem1BusInitialize    (void*,void*);
int FPPowerSystem1BusCleanup       (void*);
