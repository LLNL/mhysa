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

/*! Compute the electrical power of each generator, given their phases and other system parameters */
static void ComputeElectricalPower(
                                    double theta1,  /*!< Phase of generator 1 */
                                    double theta2,  /*!< Phase of generator 2 */
                                    void   *p,      /*!< Object of type #FPPowerSystem3Bus */
                                    double *Pe1,    /*!< Electrical power of generator 1 */
                                    double *Pe2,    /*!< Electrical power of generator 2 */
                                    double *Pe3     /*!< Electrical power of generator 3 */
                                  )
{
  FPPowerSystem3Bus *params = (FPPowerSystem3Bus*) p;

  double E1 = params->E1;
  double E2 = params->E2;
  double E3 = params->Eref;

  double *G = params->G;
  double *B = params->B;

  double Eph[3][2];
  Eph[0][0] = E1*cos(theta1);   Eph[0][1] = E1*sin(theta1);
  Eph[1][0] = E2*cos(theta2);   Eph[1][1] = E2*sin(theta2);
  Eph[2][0] = E3;               Eph[2][1] = 0.0;

  double Y[3][3][2];
  int i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Y[i][j][0] = G[i*3+j];
      Y[i][j][1] = B[i*3+j];
    }
  }

  double YEph[3][2];
  YEph[0][0] =   Y[0][0][0]*Eph[0][0] - Y[0][0][1]*Eph[0][1]
               + Y[0][1][0]*Eph[1][0] - Y[0][1][1]*Eph[1][1]
               + Y[0][2][0]*Eph[2][0] - Y[0][2][1]*Eph[2][1];
  YEph[0][1] =   Y[0][0][0]*Eph[0][1] + Y[0][0][1]*Eph[0][0]
               + Y[0][1][0]*Eph[1][1] + Y[0][1][1]*Eph[1][0]
               + Y[0][2][0]*Eph[2][1] + Y[0][2][1]*Eph[2][0];
  YEph[1][0] =   Y[1][0][0]*Eph[0][0] - Y[1][0][1]*Eph[0][1]
               + Y[1][1][0]*Eph[1][0] - Y[1][1][1]*Eph[1][1]
               + Y[1][2][0]*Eph[2][0] - Y[1][2][1]*Eph[2][1];
  YEph[1][1] =   Y[1][0][0]*Eph[0][1] + Y[1][0][1]*Eph[0][0]
               + Y[1][1][0]*Eph[1][1] + Y[1][1][1]*Eph[1][0]
               + Y[1][2][0]*Eph[2][1] + Y[1][2][1]*Eph[2][0];
  YEph[2][0] =   Y[2][0][0]*Eph[0][0] - Y[2][0][1]*Eph[0][1]
               + Y[2][1][0]*Eph[1][0] - Y[2][1][1]*Eph[1][1]
               + Y[2][2][0]*Eph[2][0] - Y[2][2][1]*Eph[2][1];
  YEph[2][1] =   Y[2][0][0]*Eph[0][1] + Y[2][0][1]*Eph[0][0]
               + Y[2][1][0]*Eph[1][1] + Y[2][1][1]*Eph[1][0]
               + Y[2][2][0]*Eph[2][1] + Y[2][2][1]*Eph[2][0];

  YEph[0][1] = - YEph[0][1];
  YEph[1][1] = - YEph[1][1];
  YEph[2][1] = - YEph[2][1];

  *Pe1 = Eph[0][0]*YEph[0][0] - Eph[0][1]*YEph[0][1];
  *Pe2 = Eph[1][0]*YEph[1][0] - Eph[1][1]*YEph[1][1];
  *Pe3 = Eph[2][0]*YEph[2][0] - Eph[2][1]*YEph[2][1];
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
