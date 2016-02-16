/*! @file TimeGLMGEEInitialize.c
    @brief Initialize the #_GLM_GEE_ time integrator
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <matops.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <timeintegration.h>

/*!
  Initialize the GLM-GEE (#_GLM_GEE_) time integation method: This function
  allocates the arrays to store the Butcher tableaux, and sets their 
  coefficients, as well as other parameters for the GLM-GEE methods.

  Reference:
  + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 
    2015 (http://arxiv.org/abs/1503.05166).
*/
int TimeGLMGEEInitialize(
                          char *class, /*!< Class of time integration method; must match #_GLM_GEE_ */
                          char *type,  /*!< Name of the #_GLM_GEE_ method to use */
                          void *s,     /*!< Object of type #GLMGEEParameters */
                          void *m      /*!< MPI object of type #MPIVariables */
                        )
{
  GLMGEEParameters *params = (GLMGEEParameters*) s;
  MPIVariables     *mpi    = (MPIVariables*) m;
  int i,j;

  if (!strcmp(class,_GLM_GEE_)) {

    if (!strcmp(type,_GLM_GEE_23_)) {
      params->nstages = 3;
      params->r       = 2;
      params->gamma   = 0.0;
    } else if (!strcmp(type,_GLM_GEE_24_)) {
      params->nstages = 4;
      params->r       = 2;
      params->gamma   = 0.0;
    } else if (!strcmp(type,_GLM_GEE_25I_)) {
      params->nstages = 5;
      params->r       = 2;
      params->gamma   = 0.0;
    } else if (!strcmp(type,_GLM_GEE_35_)) {
      params->nstages = 5;
      params->r       = 2;
      params->gamma   = 0.0;
    } else if (!strcmp(type,_GLM_GEE_EXRK2A_)) {
      params->nstages = 6;
      params->r       = 2;
      params->gamma   = 0.25;
    } else if (!strcmp(type,_GLM_GEE_RK32G1_)) {
      params->nstages = 8;
      params->r       = 2;
      params->gamma   = 0.0;
    } else if (!strcmp(type,_GLM_GEE_RK285EX_)) {
      params->nstages = 9;
      params->r       = 2;
      params->gamma   = 0.25;
    } else {
      fprintf(stderr,"Error in TimeGLMGEEInitialize(): %s is not a supported ",type);
      fprintf(stderr,"multi-stage time integration scheme of class %s.\n",class);
      return(1);
    }
  
    int s = params->nstages;
    int r = params->r;

    params->A_yyt = (double*) calloc (s*s,sizeof(double));
    params->B_yyt = (double*) calloc (s*r,sizeof(double));
    params->C_yyt = (double*) calloc (s*r,sizeof(double));
    params->D_yyt = (double*) calloc (r*r,sizeof(double));
    params->c_yyt = (double*) calloc (s  ,sizeof(double));
    _ArraySetValue_(params->A_yyt,s*s,0.0);
    _ArraySetValue_(params->B_yyt,s*r,0.0);
    _ArraySetValue_(params->C_yyt,s*r,0.0);
    _ArraySetValue_(params->D_yyt,r*r,0.0);
    _ArraySetValue_(params->c_yyt,s  ,0.0);
    
    params->A_yeps = (double*) calloc (s*s,sizeof(double));
    params->B_yeps = (double*) calloc (s*r,sizeof(double));
    params->C_yeps = (double*) calloc (s*r,sizeof(double));
    params->D_yeps = (double*) calloc (r*r,sizeof(double));
    params->c_yeps = (double*) calloc (s  ,sizeof(double));
    _ArraySetValue_(params->A_yeps,s*s,0.0);
    _ArraySetValue_(params->B_yeps,s*r,0.0);
    _ArraySetValue_(params->C_yeps,s*r,0.0);
    _ArraySetValue_(params->D_yeps,r*r,0.0);
    _ArraySetValue_(params->c_yeps,s  ,0.0);
    
    if (!strcmp(type,_GLM_GEE_23_)) {

      params->A_yeps[1*s+0] = 1.0;
      params->A_yeps[2*s+0] = params->A_yeps[2*s+1] = 0.25;

      params->B_yeps[0*s+0] = params->B_yeps[0*s+1] = 1.0/12.0; params->B_yeps[0*s+2] =  5.0/6.0;
      params->B_yeps[1*s+0] = params->B_yeps[1*s+1] = 1.0/12.0; params->B_yeps[1*s+2] = -1.0/6.0;

      params->C_yeps[0*r+0] = 1.0;
      params->C_yeps[1*r+0] = 1.0;  params->C_yeps[1*r+1] = 10.0;
      params->C_yeps[2*r+0] = 1.0;  params->C_yeps[2*r+1] = -1.0;

      params->D_yeps[0*r+0] = 1.0;
      params->D_yeps[1*r+1] = 1.0;

      params->A_yyt[1*s+0] = 1.0;
      params->A_yyt[2*s+0] = params->A_yyt[2*s+1] = 0.25;

      params->B_yyt[0*s+0] = params->B_yyt[0*s+1] = 1.0/12.0; params->B_yyt[0*s+2] =  5.0/6.0;
      params->B_yyt[1*s+0] = params->B_yyt[1*s+1] = 1.0/6.0;  params->B_yyt[1*s+2] =  2.0/3.0;

      params->C_yyt[0*r+0] =  1.0;
      params->C_yyt[1*r+0] = -9.0;  params->C_yyt[1*r+1] = 10.0;
      params->C_yyt[2*r+0] =  2.0;  params->C_yyt[2*r+1] = -1.0;

      params->D_yyt[0*r+0] = 1.0;
      params->D_yyt[1*r+1] = 1.0;

    } else if (!strcmp(type,_GLM_GEE_24_)) {

      params->A_yyt[1*s+0] = 0.75;
      params->A_yyt[2*s+0] = 0.25;
      params->A_yyt[2*s+1] = 29.0/60.0;
      params->A_yyt[3*s+0] = -21.0/44.0;
      params->A_yyt[3*s+1] = 145.0/44.0;
      params->A_yyt[3*s+2] = -20.0/11.0;

      params->B_yyt[0*s+0] = 109.0/275.0;
      params->B_yyt[0*s+1] = 58.0/75.0;
      params->B_yyt[0*s+2] = -37.0/110.0;
      params->B_yyt[0*s+3] = 1.0/6.0;
      params->B_yyt[1*s+0] = 3.0/11.0;
      params->B_yyt[1*s+2] = 75.0/88.0;
      params->B_yyt[1*s+3] = -1.0/8.0;

      params->C_yyt[0*r+1] = 1.0;
      params->C_yyt[1*r+0] = 75.0/58.0;
      params->C_yyt[1*r+1] = -17.0/58.0;
      params->C_yyt[2*r+1] = params->C_yyt[3*r+1] = 1.0;

      params->D_yyt[0*r+0] = 1.0;
      params->D_yyt[1*r+1] = 1.0;

      double T[r*r],Tinv[r*r];
      T[0*r+0] = 1.0;
      T[0*r+1] = 0.0;
      T[1*r+0] = 1.0;
      T[1*r+1] = 1.0-params->gamma;
      _MatrixInvert_(T,Tinv,r);

      _ArrayCopy1D_(params->A_yyt,params->A_yeps,(s*s));
      _MatrixMultiplyNonSquare_(Tinv,params->B_yyt,params->B_yeps,r,r,s);
      _MatrixMultiplyNonSquare_(params->C_yyt,T,params->C_yeps,s,r,r);
      _ArrayCopy1D_(params->D_yyt,params->D_yeps,(r*r));

    } else if (!strcmp(type,_GLM_GEE_25I_)) {

      params->A_yyt[1*s+0]=-0.94079244066783383269;
      params->A_yyt[2*s+0]= 0.64228187778301907108;
      params->A_yyt[2*s+1]= 0.10915356933958500042;
      params->A_yyt[3*s+0]=-0.51764297742287450812;
      params->A_yyt[3*s+1]= 0.74414270351096040738;
      params->A_yyt[3*s+2]=-0.71404164927824538121;
      params->A_yyt[4*s+0]=-0.44696561556825969206;
      params->A_yyt[4*s+1]=-0.76768425657590196518;
      params->A_yyt[4*s+2]= 0.20111608138142987881;
      params->A_yyt[4*s+3]= 0.93828186737840469796;

      params->B_yyt[0*s+0]=-0.029309178948150356153;
      params->B_yyt[0*s+1]=-0.49671981884013874923;
      params->B_yyt[0*s+2]= 0.34275801517650053274;
      params->B_yyt[0*s+3]= 0.32941112623949194988;
      params->B_yyt[0*s+4]= 0.85385985637229662276;
      params->B_yyt[1*s+0]= 0.78133219686062535272;
      params->B_yyt[1*s+1]= 0.074238691892675897635;
      params->B_yyt[1*s+2]= 0.57957363498384957966;
      params->B_yyt[1*s+3]=-0.24638502829674959968;
      params->B_yyt[1*s+4]=-0.18875949544040123033;

      params->C_yyt[0*r+0]= 0.16911424754448327735;
      params->C_yyt[0*r+1]= 0.83088575245551672265;
      params->C_yyt[1*r+0]= 0.53638465733199574340;
      params->C_yyt[1*r+1]= 0.46361534266800425660;
      params->C_yyt[2*r+0]= 0.39901579167169582526;
      params->C_yyt[2*r+1]= 0.60098420832830417474;
      params->C_yyt[3*r+0]= 0.87689005530618575480;
      params->C_yyt[3*r+1]= 0.12310994469381424520;
      params->C_yyt[4*r+0]= 0.99056100455550913009;
      params->C_yyt[4*r+1]= 0.0094389954444908699092;

      params->D_yyt[0*r+0] = 1.0;
      params->D_yyt[1*r+1] = 1.0;

      double T[r*r],Tinv[r*r];
      T[0*r+0] = 1.0;
      T[0*r+1] = 0.0;
      T[1*r+0] = 1.0;
      T[1*r+1] = 1.0-params->gamma;
      _MatrixInvert_(T,Tinv,r);

      _ArrayCopy1D_(params->A_yyt,params->A_yeps,(s*s));
      _MatrixMultiplyNonSquare_(Tinv,params->B_yyt,params->B_yeps,r,r,s);
      _MatrixMultiplyNonSquare_(params->C_yyt,T,params->C_yeps,s,r,r);
      _ArrayCopy1D_(params->D_yyt,params->D_yeps,(r*r));

    } else if (!strcmp(type,_GLM_GEE_35_)) {

      params->A_yyt[1*s+0] = - 2169604947363702313.0 /  24313474998937147335.0;
      params->A_yyt[2*s+0] =  46526746497697123895.0 /  94116917485856474137.0;
      params->A_yyt[2*s+1] = -10297879244026594958.0 /  49199457603717988219.0;
      params->A_yyt[3*s+0] =  23364788935845982499.0 /  87425311444725389446.0;
      params->A_yyt[3*s+1] = -79205144337496116638.0 / 148994349441340815519.0;
      params->A_yyt[3*s+2] =  40051189859317443782.0 /  36487615018004984309.0;
      params->A_yyt[4*s+0] =  42089522664062539205.0 / 124911313006412840286.0;
      params->A_yyt[4*s+1] = -15074384760342762939.0 / 137927286865289746282.0;
      params->A_yyt[4*s+2] = -62274678522253371016.0 / 125918573676298591413.0;
      params->A_yyt[4*s+3] =  13755475729852471739.0 /  79257927066651693390.0;
      
      params->B_yyt[0*s+0] =  61546696837458703723.0 /  56982519523786160813.0;
      params->B_yyt[0*s+1] = -55810892792806293355.0 / 206957624151308356511.0;
      params->B_yyt[0*s+2] =  24061048952676379087.0 / 158739347956038723465.0;
      params->B_yyt[0*s+3] =   3577972206874351339.0 /   7599733370677197135.0;
      params->B_yyt[0*s+4] = -59449832954780563947.0 / 137360038685338563670.0;
      params->B_yyt[1*s+0] = - 9738262186984159168.0 /  99299082461487742983.0;
      params->B_yyt[1*s+1] = -32797097931948613195.0 /  61521565616362163366.0;
      params->B_yyt[1*s+2] =  42895514606418420631.0 /  71714201188501437336.0;
      params->B_yyt[1*s+3] =  22608567633166065068.0 /  55371917805607957003.0;
      params->B_yyt[1*s+4] =  94655809487476459565.0 / 151517167160302729021.0;

      params->C_yyt[0*r+0] =  70820309139834661559.0 /  80863923579509469826.0;
      params->C_yyt[0*r+1] =  10043614439674808267.0 /  80863923579509469826.0;
      params->C_yyt[1*r+0] = 161694774978034105510.0 / 106187653640211060371.0;
      params->C_yyt[1*r+1] = -55507121337823045139.0 / 106187653640211060371.0;
      params->C_yyt[2*r+0] =  78486094644566264568.0 /  88171030896733822981.0;
      params->C_yyt[2*r+1] =   9684936252167558413.0 /  88171030896733822981.0;
      params->C_yyt[3*r+0] =  65394922146334854435.0 /  84570853840405479554.0;
      params->C_yyt[3*r+1] =  19175931694070625119.0 /  84570853840405479554.0;
      params->C_yyt[4*r+0] =   8607282770183754108.0 / 108658046436496925911.0;
      params->C_yyt[4*r+1] = 100050763666313171803.0 / 108658046436496925911.0;

      params->D_yyt[0*r+0] = 1.0;
      params->D_yyt[1*r+1] = 1.0;

      double T[r*r],Tinv[r*r];
      T[0*r+0] = 1.0;
      T[0*r+1] = 0.0;
      T[1*r+0] = 1.0;
      T[1*r+1] = 1.0-params->gamma;
      _MatrixInvert_(T,Tinv,r);

      _ArrayCopy1D_(params->A_yyt,params->A_yeps,(s*s));
      _MatrixMultiplyNonSquare_(Tinv,params->B_yyt,params->B_yeps,r,r,s);
      _MatrixMultiplyNonSquare_(params->C_yyt,T,params->C_yeps,s,r,r);
      _ArrayCopy1D_(params->D_yyt,params->D_yeps,(r*r));

    } else if (!strcmp(type,_GLM_GEE_EXRK2A_)) {

      params->A_yeps[1*s+0] = 1.0;
      params->A_yeps[3*s+2] = 0.5;
      params->A_yeps[4*s+2] = 0.25;
      params->A_yeps[4*s+3] = 0.25;
      params->A_yeps[5*s+2] = 0.25;
      params->A_yeps[5*s+3] = 0.25;
      params->A_yeps[5*s+4] = 0.5;

      params->B_yeps[0*s+0] = params->B_yeps[0*s+1] = 0.5;
      params->B_yeps[1*s+0] = params->B_yeps[1*s+1] = -2.0/3.0;
      params->B_yeps[1*s+2] = params->B_yeps[1*s+3] = params->B_yeps[1*s+4]
                            = params->B_yeps[1*s+5] =  1.0/3.0;

      params->C_yeps[0*r+0] = 1.0;
      params->C_yeps[1*r+0] = 1.0;
      params->C_yeps[2*r+0] = 1.0;  params->C_yeps[2*r+1] = 0.75;
      params->C_yeps[3*r+0] = 1.0;  params->C_yeps[3*r+1] = 0.75;
      params->C_yeps[4*r+0] = 1.0;  params->C_yeps[4*r+1] = 0.75;
      params->C_yeps[5*r+0] = 1.0;  params->C_yeps[5*r+1] = 0.75;

      params->D_yeps[0*r+0] = 1.0;
      params->D_yeps[1*r+1] = 1.0;

      params->A_yyt[1*s+0] = 1.0;
      params->A_yyt[3*s+2] = 0.5;
      params->A_yyt[4*s+2] = 0.25;
      params->A_yyt[4*s+3] = 0.25;
      params->A_yyt[5*s+2] = 0.25;
      params->A_yyt[5*s+3] = 0.25;
      params->A_yyt[5*s+4] = 0.5;

      params->B_yyt[0*s+0] = params->B_yyt[0*s+1] = 0.5;
      params->B_yyt[1*s+2] = params->B_yyt[1*s+3] = params->B_yyt[1*s+4]
                           = params->B_yyt[1*s+5] = 0.25;

      params->C_yyt[0*r+0] = 1.0;
      params->C_yyt[1*r+0] = 1.0;
      params->C_yyt[2*r+1] = 1.0;
      params->C_yyt[3*r+1] = 1.0;
      params->C_yyt[4*r+1] = 1.0;
      params->C_yyt[5*r+1] = 1.0;

      params->D_yyt[0*r+0] = 1.0;
      params->D_yyt[1*r+1] = 1.0;

    } else if (!strcmp(type,_GLM_GEE_RK32G1_)) {

      params->A_yeps[1*s+0] =  0.5;

      params->A_yeps[2*s+0] = -1.0;
      params->A_yeps[2*s+1] =  2.0;
        
      params->A_yeps[3*s+0] =  1.0/6.0;
      params->A_yeps[3*s+1] =  2.0/3.0;
      params->A_yeps[3*s+2] =  1.0/6.0;
        
      params->A_yeps[5*s+0] = -7.0/24.0;
      params->A_yeps[5*s+1] =  1.0/3.0;
      params->A_yeps[5*s+2] =  1.0/12.0;
      params->A_yeps[5*s+3] = -1.0/8.0;
      params->A_yeps[5*s+4] =  0.5;
        
      params->A_yeps[6*s+0] =  7.0/6.0;
      params->A_yeps[6*s+1] = -4.0/3.0;
      params->A_yeps[6*s+2] = -1.0/3.0;
      params->A_yeps[6*s+3] =  0.5;
      params->A_yeps[6*s+4] = -1.0;
      params->A_yeps[6*s+5] =  2.0;
      
      params->A_yeps[7*s+4] =  1.0/6.0;
      params->A_yeps[7*s+5] =  2.0/3.0;
      params->A_yeps[7*s+6] =  1.0/6.0;

      params->B_yeps[0*s+0] = 1.0/6.0;
      params->B_yeps[0*s+1] = 2.0/3.0;
      params->B_yeps[0*s+2] = 1.0/6.0;

      params->B_yeps[1*s+0] = -1.0/6.0;
      params->B_yeps[1*s+1] = -2.0/3.0;
      params->B_yeps[1*s+2] = -1.0/6.0;
      params->B_yeps[1*s+4] =  1.0/6.0;
      params->B_yeps[1*s+5] =  2.0/3.0;
      params->B_yeps[1*s+6] =  1.0/6.0;

      params->C_yeps[0*r+0] = 1.0;
      params->C_yeps[1*r+0] = 1.0;
      params->C_yeps[2*r+0] = 1.0;
      params->C_yeps[3*r+0] = 1.0;
      params->C_yeps[4*r+0] = 1.0;
      params->C_yeps[5*r+0] = 1.0;
      params->C_yeps[6*r+0] = 1.0;
      params->C_yeps[7*r+0] = 1.0;

      params->C_yeps[4*r+1] = 1.0;
      params->C_yeps[5*r+1] = 1.0;
      params->C_yeps[6*r+1] = 1.0;
      params->C_yeps[7*r+1] = 1.0;

      params->D_yeps[0*r+0] = 1.0;
      params->D_yeps[1*r+1] = 1.0;

      double T[r*r],Tinv[r*r];
      T[0*r+0] = 1.0;
      T[0*r+1] = 0.0;
      T[1*r+0] = 1.0;
      T[1*r+1] = 1.0-params->gamma;
      _MatrixInvert_(T,Tinv,r);

      _ArrayCopy1D_(params->A_yeps,params->A_yyt,(s*s));
      _MatrixMultiplyNonSquare_(T,params->B_yeps,params->B_yyt,r,r,s);
      _MatrixMultiplyNonSquare_(params->C_yeps,Tinv,params->C_yyt,s,r,r);
      _ArrayCopy1D_(params->D_yeps,params->D_yyt,(r*r));

    } else if (!strcmp(type,_GLM_GEE_RK285EX_)) {

      params->A_yeps[1*s+0] = 0.585786437626904966;
      params->A_yeps[2*s+0] = 0.149999999999999994;
      params->A_yeps[2*s+1] = 0.849999999999999978;
      params->A_yeps[4*s+3] = 0.292893218813452483;
      params->A_yeps[5*s+3] = 0.0749999999999999972;
      params->A_yeps[5*s+4] = 0.424999999999999989;
      params->A_yeps[6*s+3] = 0.176776695296636893;
      params->A_yeps[6*s+4] = 0.176776695296636893;
      params->A_yeps[6*s+5] = 0.146446609406726241;
      params->A_yeps[7*s+3] = 0.176776695296636893;
      params->A_yeps[7*s+4] = 0.176776695296636893;
      params->A_yeps[7*s+5] = 0.146446609406726241;
      params->A_yeps[7*s+6] = 0.292893218813452483;
      params->A_yeps[8*s+3] = 0.176776695296636893;
      params->A_yeps[8*s+4] = 0.176776695296636893;
      params->A_yeps[8*s+5] = 0.146446609406726241;
      params->A_yeps[8*s+6] = 0.0749999999999999972;
      params->A_yeps[8*s+7] = 0.424999999999999989;

      params->B_yeps[0*s+0] = 0.353553390593273786;
      params->B_yeps[0*s+1] = 0.353553390593273786;
      params->B_yeps[0*s+2] = 0.292893218813452483;

      params->B_yeps[1*s+0] =-0.471404520791031678;
      params->B_yeps[1*s+1] =-0.471404520791031678;
      params->B_yeps[1*s+2] =-0.390524291751269959;
      params->B_yeps[1*s+3] = 0.235702260395515839;
      params->B_yeps[1*s+4] = 0.235702260395515839;
      params->B_yeps[1*s+5] = 0.195262145875634979;
      params->B_yeps[1*s+6] = 0.235702260395515839;
      params->B_yeps[1*s+7] = 0.235702260395515839;
      params->B_yeps[1*s+8] = 0.195262145875634979;

      params->C_yeps[0*r+0] = 1.0;
      params->C_yeps[1*r+0] = 1.0;
      params->C_yeps[2*r+0] = 1.0;
      params->C_yeps[3*r+0] = 1.0; params->C_yeps[3*r+1] = 0.75;
      params->C_yeps[4*r+0] = 1.0; params->C_yeps[4*r+1] = 0.75;
      params->C_yeps[5*r+0] = 1.0; params->C_yeps[5*r+1] = 0.75;
      params->C_yeps[6*r+0] = 1.0; params->C_yeps[6*r+1] = 0.75;
      params->C_yeps[7*r+0] = 1.0; params->C_yeps[7*r+1] = 0.75;
      params->C_yeps[8*r+0] = 1.0; params->C_yeps[8*r+1] = 0.75;

      params->D_yeps[0*r+0] = 1.0;
      params->D_yeps[1*r+1] = 1.0;

      double T[r*r],Tinv[r*r];
      T[0*r+0] = 1.0;
      T[0*r+1] = 0.0;
      T[1*r+0] = 1.0;
      T[1*r+1] = 1.0-params->gamma;
      _MatrixInvert_(T,Tinv,r);

      _ArrayCopy1D_(params->A_yeps,params->A_yyt,(s*s));
      _MatrixMultiplyNonSquare_(T,params->B_yeps,params->B_yyt,r,r,s);
      _MatrixMultiplyNonSquare_(params->C_yeps,Tinv,params->C_yyt,s,r,r);
      _ArrayCopy1D_(params->D_yeps,params->D_yyt,(r*r));

    }

    for (i=0; i<s; i++) {
      for (j=0; j<s; j++) {
        params->c_yyt[i]  += params->A_yyt [i*s+j];
        params->c_yeps[i] += params->A_yeps[i*s+j];
      }
    }

    if (!mpi->rank) {
      FILE *in;
      int ferr;
      in = fopen("glm_gee.inp","r");
      strcat(params->ee_mode,_GLM_GEE_YEPS_);
      if (in) {
        printf("Reading GLM-GEE method parameters from glm_gee.inp.\n");
        char word[_MAX_STRING_SIZE_];
        ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
        if (!strcmp(word,"begin")) {
          while (strcmp(word,"end")) {
            ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
            if (!strcmp(word,"ee_mode")) { 
              ferr = fscanf(in,"%s",params->ee_mode); 
              if (ferr != 1) return(1); 
            } else if (strcmp(word,"end")) {
              char useless[_MAX_STRING_SIZE_];
              ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
              printf("Warning: keyword %s in file \"glm_gee.inp\" with value %s not ",word,useless);
              printf("recognized or extraneous. Ignoring.\n");
            }
          }
        } else {
          fprintf(stderr,"Error: Illegal format in file \"glm_gee.inp\".\n");
          return(1);
        }
        fclose(in);
        if (strcmp(params->ee_mode,_GLM_GEE_YEPS_) && strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
          fprintf(stderr,"Error in TimeGLMGEEInitialize(): %s is not a valid value for ",params->ee_mode);
          fprintf(stderr,"ee_mode. Acceptable inputs are %s or %s.\n",_GLM_GEE_YEPS_,_GLM_GEE_YYT_);
          strcat(params->ee_mode,_GLM_GEE_YEPS_);
        }
      }
      printf("GLM-GEE time integration error estimation mode: %s\n",params->ee_mode);
    }
    IERR MPIBroadcast_character(params->ee_mode,_MAX_STRING_SIZE_,0,&mpi->world);

    if (!strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
      params->A = params->A_yyt;
      params->B = params->B_yyt;
      params->C = params->C_yyt;
      params->D = params->D_yyt;
      params->c = params->c_yyt;
    } else {
      params->A = params->A_yeps;
      params->B = params->B_yeps;
      params->C = params->C_yeps;
      params->D = params->D_yeps;
      params->c = params->c_yeps;
    }

  } else {
    fprintf(stderr,"Error in TimeGLMGEEInitialize(): Code should not have ");
    fprintf(stderr,"reached here for %s class of time-integrators. This is a ",class);
    fprintf(stderr,"coding mistake.\n");
    return(1);
  }
  return(0);
}
