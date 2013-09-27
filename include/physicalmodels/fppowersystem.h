/*

  Fokker-Planck Model for Power Systems

Reference: 
+ Wang P., Tartakovsky A.M., Abhyankar S., Smith B.F.,
  Huang Z., "Probabilistic Density Function Method for
  Stochastic ODEs of Power Systems with Uncertain 
  Power Input", Preprint

  dp    d[mu(x,y)p]    d[nu(x,y)p]        d^2 p
  -- +  ----------- +  ----------- = f(t) -----
  dt        dx             dy             dy^2

  mu(x,y) = O_s(y-1)
  nu(x,y) = (1/2H)[Pm - Pmax sin(x) - D(y-1)]
  f(t)    = (1/2H)^2 (l^2 q)/(lr+1) (1 - exp[-(r+1/l)t])

  where 
    r = D/2H

           [ EV/g1; t < tf
    Pmax = [ 0    ; tf < t < tcl
           [ EV/g2; tcl < t

  Physical Parameters:
    
    O_s   Synchronous Speed
    H     Inertia constant
    D     Damping constant
    l     Correlation time
    q     Noise strength
    Pm    Mean mechanical power input
    Pmax  Maximum power output
    E     Internal voltage
    V     Terminal voltage
    tf    Fault incident time
    tcl   Clearing time
    g1    Pre-fault impedance
    g2    Post-clearing impedance

*/

#define _FP_POWER_SYSTEM_  "fp-power-system"

/* define ndims and nvars for this model */
#undef _MODEL_NDIMS_
#undef _MODEL_NVARS_
#define _MODEL_NDIMS_ 2
#define _MODEL_NVARS_ 1

typedef struct fp_power_system_parameters {

  /* input parameters */
  double O_s;
  double H;
  double E;
  double V;
  double g1,g2;
  double D;
  double Pm;
  double l;
  double q;
  double tf, tcl;

  /* computed/constant parameters */
  double Pmax;

  double pdf_integral; /* not an input          */
} FPPowerSystem;

int FPPowerSystemInitialize    (void*,void*);
int FPPowerSystemCleanup       (void*);
