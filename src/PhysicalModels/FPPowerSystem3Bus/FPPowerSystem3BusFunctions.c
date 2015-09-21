/*! @file FPPowerSystem3BusFunctions.c
    @author Debojyoti Ghosh
    @brief Miscellaneous functions for the 3-bus power system model
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <math.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>

static void ComputeElectricalPower(
                                    double theta1,  /*!< Phase of generator 1 */
                                    double theta2,  /*!< Phase of generator 2 */
                                    void   *p,      /*!< Object of type #FPPowerSystem3Bus */
                                    double *Pe1,    /*!< Electrical power of generator 1 */
                                    double *Pe2,    /*!< Electrical power of generator 2 */
                                    double *Peref   /*!< Electrical power of generator 3 (reference) */
                                  )
{
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;

  *Pe1    = 0.0;
  *Pe2    = 0.0;
  *Peref  = 0.0;
}

/*! Compute the drift (advection) coefficients for the 3-bus power system */
int FPPowerSystem3BusDriftFunction(
                                    int     dir,    /*!< Spatial dimension (not used) */
                                    void    *p,     /*!< Object of type #FPPowerSystem3Bus */
                                    double  *x,     /*!< Spatial coordinates */
                                    double  t,      /*!< Current simulation time */
                                    double  *drift  /*!< Array to hold the drift velocities */
                                  )
{
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;

  double theta1 = x[0];
  double theta2 = x[1];
  double Omega1 = x[2];
  double Omega2 = x[3];
  
  double omegaB     = params->omegaB;
  double Pm1_avg    = params->Pm1_avg;
  double Pm2_avg    = params->Pm2_avg;
  double Pmref_avg  = params->Pmref_avg;
  double H1         = params->H1;
  double H2         = params->H2;
  double Href       = params->Href;
  double gamma      = params->gamma;
  
  double Pe1, Pe2, Peref;
  ComputeElectricalPower(theta1,theta2,params,&Pe1,&Pe2,&Peref);

  double F1 = Pm1_avg / (2*H1) - Pmref_avg / (2*Href);
  double F2 = Pm2_avg / (2*H2) - Pmref_avg / (2*Href);
  double S1 = Pe1 / (2*H1) - Peref / (2*Href);
  double S2 = Pe2 / (2*H2) - Peref / (2*Href);

  drift[0] = omegaB * Omega1;
  drift[1] = omegaB * Omega2;
  drift[2] = F1 - gamma*Omega1 - S1;
  drift[3] = F2 - gamma*Omega2 - S2;

  return(0);
}

/*! Compute the dissipation coefficient for the 3-bus power system */
int FPPowerSystem3BusDissipationFunction(
                                          int     dir1,   /*!< First spatial dimension for the dissipation coefficient */
                                          int     dir2,   /*!< Second spatial dimension for the dissipation coefficient */
                                          void    *p,     /*!< Object of type #FPPowerSystem3Bus */
                                          double  t,      /*!< Current simulation time */
                                          double  *dissp  /*!< Matrix of size ndims*ndims to hold the dissipation 
                                                               coefficients (row-major format)*/
                                        )
{
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;
  _ArraySetValue_(dissp,_MODEL_NDIMS_*_MODEL_NDIMS_,0.0);

  double sigma11 = params->sigma[0][0];
  double sigma12 = params->sigma[0][1];
  double sigma21 = params->sigma[1][0];
  double sigma22 = params->sigma[1][1];

  double lambda11 = params->lambda[0][0];
  double lambda12 = params->lambda[0][1];
  double lambda21 = params->lambda[1][0];
  double lambda22 = params->lambda[1][1];

  double gamma  = params->gamma;
  double omegaB = params->omegaB;

  dissp[2*_MODEL_NDIMS_+0] = sigma11*sigma11*lambda11*lambda11*omegaB;
  dissp[2*_MODEL_NDIMS_+1] = sigma12*sigma12*lambda12*lambda12*omegaB;
  dissp[3*_MODEL_NDIMS_+0] = sigma21*sigma21*lambda21*lambda21*omegaB;
  dissp[3*_MODEL_NDIMS_+1] = sigma22*sigma22*lambda22*lambda22*omegaB;

  dissp[2*_MODEL_NDIMS_+2] = sigma11*sigma11*lambda11*(1.0-gamma*lambda11);
  dissp[2*_MODEL_NDIMS_+3] = sigma12*sigma12*lambda12*(1.0-gamma*lambda12);
  dissp[3*_MODEL_NDIMS_+2] = sigma21*sigma21*lambda21*(1.0-gamma*lambda21);
  dissp[3*_MODEL_NDIMS_+3] = sigma22*sigma22*lambda22*(1.0-gamma*lambda22);

  return(0);
}
