/* interpolation scheme definitions */
#define _FIRST_ORDER_UPWIND_    "1"
#define _SECOND_ORDER_CENTRAL_  "2"
#define _THIRD_ORDER_MUSCL_     "muscl3"
#define _FIFTH_ORDER_WENO_      "weno5"
#define _FIFTH_ORDER_CRWENO_    "crweno5"
#define _FIFTH_ORDER_HCWENO_    "hcweno5"

/* interpolation type definitions */
#define _CHARACTERISTIC_        "characteristic" /* characteristic-based interpolation */
#define _COMPONENTS_            "components"     /* component-wise interpolation       */

/*
  One-dimensional Interpolation Functions:-
    Functions to interpolate the primitive of a function at interfaces from its 
    cell-centered values.

  Arguments:-

    fI        double*   array of size (N+1) in the interpolation direction; will contain 
                        the interpolated interface function primitive 
                        (**needs to be allocated by the calling function)
                        (**does not have ghost points))

    fC        double*   array of size (N+2*ghosts) in the interpolation direction of the 
                        cell centered function values
                        (**does have ghost points))

    u         double*   used only by the characteristic-based interpolation schemes
                        array of size (N+2*ghosts) in the interpolation direction of the
                        cell centered conserved variable (needed to calculate the eigen-
                        decomposition at the interfaces)

    x         double*   grid point locations along the 1D line on which the interpolation
                        is being carried out
                        array of size (N+2*ghosts) 
                        Used only by non-uniform-grid interpolation schemes

    upw       int       upwind direction for non-central schemes
                        (1 -> left biased, -1 -> right biased)

    dir       int       direction/dimension along which to carry out interpolation for
                        multi-dimensional domains
                        (eg: 0 for 1D; 0 or 1 for 2D; 0,1 or 2 for 3D)

    s         void*     object of type depending on the main solver, should have at least
                        the following data:
                        + ghosts    : number of ghosts points in the domain
                        + ndims     : number of dimensions
                        + nvars     : number of variables per grid point
                                      (i.e. size of vector variable being interpolated)
                        + dim_local : local (of this process) dimensions of the domain
                                      (integer array of size ndims, with the number of 
                                      points in each dimension as elements)
    m         void*     object containing MPI domain decomposition related variables
                        (this information is used only by compact interpolation schemes)

    Notes:
      + arrangement of points along the direction of interpolation is

              |<--ghosts-->|0 1 2 ..(interior).. N-1|<--ghosts-->|

      + number of interfaces are one more than number of interior points
      + interface i lies between grid points i-1+ghosts and i+ghosts

            0   1   2   .....        ....          N-3 N-2 N-1
          | x | x | x | x | x | x | x | x | x | x | x | x | x |
          0   1   2   ...             ...    ...     N-2 N-1  N 

    Returns 0 on normal execution, non-zero on error.
*/

/* functions to interpolate the first primitive in a component-wise way
   (for conservative discretization of the 1st derivative) on a uniform grid */
/* First-order upwind */
int Interp1PrimFirstOrderUpwind           (double*,double*,double*,double*,int,int,void*,void*);
/* Third-order MUSCL scheme */
int Interp1PrimThirdOrderMUSCL            (double*,double*,double*,double*,int,int,void*,void*);
/* Fifth-order WENO scheme */
int Interp1PrimFifthOrderWENO             (double*,double*,double*,double*,int,int,void*,void*);
/* Fifth-order CRWENO scheme */
int Interp1PrimFifthOrderCRWENO           (double*,double*,double*,double*,int,int,void*,void*);
/* Fifth-order hybrid-compact WENO scheme */
int Interp1PrimFifthOrderHCWENO           (double*,double*,double*,double*,int,int,void*,void*);

/* functions to interpolate the first primitive in a characteristic-based way
   (for conservative discretization of the 1st derivative) on a uniform grid */
/* First-order upwind */
int Interp1PrimFirstOrderUpwindChar       (double*,double*,double*,double*,int,int,void*,void*);
/* Third-order MUSCL scheme */
int Interp1PrimThirdOrderMUSCLChar        (double*,double*,double*,double*,int,int,void*,void*);
/* Fifth-order WENO scheme */
int Interp1PrimFifthOrderWENOChar         (double*,double*,double*,double*,int,int,void*,void*);
/* Fifth-order CRWENO scheme */
int Interp1PrimFifthOrderCRWENOChar       (double*,double*,double*,double*,int,int,void*,void*);
/* Fifth-order hybrid-compact WENO scheme */
int Interp1PrimFifthOrderHCWENOChar       (double*,double*,double*,double*,int,int,void*,void*);

/* functions to interpolate the second primitive 
   (for conservative discretization of the 2nd derivative) */
/* Second-order central */
int Interp2PrimSecondOrder  (double*,double*,int,void*,void*);

/* Function to calculate and save the nonlinear interpolation coefficients */
int InterpSetLimiterVar(double*,double*,double*,int,void*,void*);

/* MUSCL scheme related parameters */
typedef struct paramters_muscl {
  double eps;
} MUSCLParameters;
int MUSCLInitialize(void*,void*);

/* WENO/CRWENO/HCWENO schemes related parameters and functions */
typedef struct parameters_weno {
  /* Options related to the type of WENO scheme */
  int     mapped;		    /* Use mapped weights?                                    */
  int     borges;		    /* Use Borges' implementation of weights?                 */
  int     yc;		        /* Use Yamaleev-Carpenter implementation of weights?      */
  int     no_limiting;  /* Remove limiting -> 5th order polynomial interpolation  */
  double  eps;		      /* epsilon parameter                                      */
  double	p;			      /* p parameter                                            */

  /* hybrid compact-WENO scheme related parameters 
   * References: 
   * + http://dx.doi.org/10.1006/jcph.2002.7021
   * + http://dx.doi.org/10.1016/j.jcp.2003.07.006
  */
  double  rc, xi;

  /* CRWENO scheme: sub-, main-, and super-diagonals
   * and the right-hand-side */
  double *A, *B, *C, *R;
  /* CRWENO scheme: buffer arrays for sending and
   * receiving data */
  double *sendbuf, *recvbuf;

  /* Arrays to save the WENO weights */
  double *w1, *w2, *w3;
  /* size and offset for the WENO weights arrays */
  int *offset, size;

} WENOParameters;
int WENOInitialize(void*,void*,char*,char*);
int WENOCleanup(void*);

/* define optimal weights */
#define   _WENO_OPTIMAL_WEIGHT_1_   0.1
#define   _WENO_OPTIMAL_WEIGHT_2_   0.6
#define   _WENO_OPTIMAL_WEIGHT_3_   0.3
#define   _CRWENO_OPTIMAL_WEIGHT_1_ 0.2
#define   _CRWENO_OPTIMAL_WEIGHT_2_ 0.5
#define   _CRWENO_OPTIMAL_WEIGHT_3_ 0.3

/* Macro to calculate the fifth-order WENO weights 
 *
 * Arguments:-
 *  w1,w2,w3      : the nonlinear WENO weights
 *  c1,c2,c3      : optimal coefficients
 *  m3,m2,m1,p1,p2: function value at stencil points corresponding 
 *                  to the interface i-1/2: i-3,i-2,i-1,i,i+1
 *  weno          : object of type WENOParameters containing 
 *                  parameters for the WENO method
 *
 * References:
 *  + http://dx.doi.org/10.1006/jcph.1996.0130
 *  + http://dx.doi.org/10.1016/j.jcp.2005.01.023
 *  + http://dx.doi.org/10.1016/j.jcp.2007.11.038
 *  + http://dx.doi.org/10.1016/j.jcp.2009.03.002
 *  + http://dx.doi.org/10.1007/s10915-014-9818-0
 *  
 */
#define _WENOWeights_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno) \
  { \
    /* calculate smoothness indicators */\
    double b1, b2, b3; \
    b1 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) \
         + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);  \
    b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1) \
         + one_fourth*(m2-p1)*(m2-p1);                \
    b3 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2) \
         + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);  \
    /* calculate the tau parameter for the WENO-Z and WENO-YC weights */\
    double tau; \
    if      (weno->borges) tau = absolute(b3 - b1);  \
    else if (weno->yc)     tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);  \
    else                   tau = 0;  \
    /* calculate the WENO weights */\
    double a1, a2, a3;  \
    if (weno->borges || weno->yc) { \
      a1 = c1 * (1.0 + (tau/(b1+weno->eps)) * (tau/(b1+weno->eps)) );  \
      a2 = c2 * (1.0 + (tau/(b2+weno->eps)) * (tau/(b2+weno->eps)) );  \
      a3 = c3 * (1.0 + (tau/(b3+weno->eps)) * (tau/(b3+weno->eps)) );  \
    } else {  \
      a1 = c1 / ( (b1+weno->eps) + (b1+weno->eps) );  \
      a2 = c2 / ( (b2+weno->eps) + (b2+weno->eps) );  \
      a3 = c3 / ( (b3+weno->eps) + (b3+weno->eps) );  \
    } \
    double a_sum_inv; \
    a_sum_inv = 1.0 / (a1 + a2 + a3); \
    w1 = a1 * a_sum_inv;  \
    w2 = a2 * a_sum_inv;  \
    w3 = a3 * a_sum_inv;  \
    /* apply mapping if required */\
    if (weno->mapped) { \
      a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1)); \
      a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2)); \
      a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3)); \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1 = a1 * a_sum_inv;  \
      w2 = a2 * a_sum_inv;  \
      w3 = a3 * a_sum_inv;  \
    } \
  }

#define _WENOWeights_v_JS_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators */\
    double b1[N], b2[N], b3[N]; \
    for (idx=0; idx<N; idx++) { \
      b1[idx] = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
                + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      b2[idx] = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
                + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      b3[idx] = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
                + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
    } \
    /* calculate the WENO weights */\
    double a1[N], a2[N], a3[N];  \
    for (idx=0; idx<N; idx++) { \
      a1[idx] = c1 / ( (b1[idx]+weno->eps) * (b1[idx]+weno->eps) );  \
      a2[idx] = c2 / ( (b2[idx]+weno->eps) * (b2[idx]+weno->eps) );  \
      a3[idx] = c3 / ( (b3[idx]+weno->eps) * (b3[idx]+weno->eps) );  \
    } \
    double a_sum_inv[N]; \
    for (idx=0; idx<N; idx++) a_sum_inv[idx] = 1.0 / (a1[idx] + a2[idx] + a3[idx]); \
    _ArrayMultiply1D_(w1,a1,a_sum_inv,N); \
    _ArrayMultiply1D_(w2,a2,a_sum_inv,N); \
    _ArrayMultiply1D_(w3,a3,a_sum_inv,N); \
  }

#define _WENOWeights_v_M_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators */\
    double b1[N], b2[N], b3[N]; \
    for (idx=0; idx<N; idx++) { \
      b1[idx] = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
                + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      b2[idx] = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
                + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      b3[idx] = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
                + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
    } \
    /* calculate the WENO weights */\
    double a1[N], a2[N], a3[N];  \
    for (idx=0; idx<N; idx++) { \
      a1[idx] = c1 / ( (b1[idx]+weno->eps) * (b1[idx]+weno->eps) );  \
      a2[idx] = c2 / ( (b2[idx]+weno->eps) * (b2[idx]+weno->eps) );  \
      a3[idx] = c3 / ( (b3[idx]+weno->eps) * (b3[idx]+weno->eps) );  \
    } \
    double a_sum_inv[N]; \
    for (idx=0; idx<N; idx++) a_sum_inv[idx] = 1.0 / (a1[idx] + a2[idx] + a3[idx]); \
    _ArrayMultiply1D_(w1,a1,a_sum_inv,N); \
    _ArrayMultiply1D_(w2,a2,a_sum_inv,N); \
    _ArrayMultiply1D_(w3,a3,a_sum_inv,N); \
    /* apply mapping if required */\
    for (idx=0; idx<N; idx++) { \
      a1[idx] = w1[idx] * (c1 + c1*c1 - 3*c1*w1[idx] + w1[idx]*w1[idx]) / (c1*c1 + w1[idx]*(1.0-2.0*c1)); \
      a2[idx] = w2[idx] * (c2 + c2*c2 - 3*c2*w2[idx] + w2[idx]*w2[idx]) / (c2*c2 + w2[idx]*(1.0-2.0*c2)); \
      a3[idx] = w3[idx] * (c3 + c3*c3 - 3*c3*w3[idx] + w3[idx]*w3[idx]) / (c3*c3 + w3[idx]*(1.0-2.0*c3)); \
    } \
    for (idx=0; idx<N; idx++) a_sum_inv[idx] = 1.0 / (a1[idx] + a2[idx] + a3[idx]); \
    _ArrayMultiply1D_(w1,a1,a_sum_inv,N); \
    _ArrayMultiply1D_(w2,a2,a_sum_inv,N); \
    _ArrayMultiply1D_(w3,a3,a_sum_inv,N); \
  }

#define _WENOWeights_v_Z_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators */\
    double b1[N], b2[N], b3[N]; \
    for (idx=0; idx<N; idx++) { \
      b1[idx] = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
                + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      b2[idx] = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
                + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      b3[idx] = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
                + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
    } \
    /* calculate the tau parameter for the WENO-Z and WENO-YC weights */\
    double tau[N]; \
    for (idx=0; idx<N; idx++) tau[idx] = absolute(b3[idx] - b1[idx]);  \
    /* calculate the WENO weights */\
    double a1[N], a2[N], a3[N];  \
    for (idx=0; idx<N; idx++) { \
      a1[idx] = c1 * (1.0 + (tau[idx]/(b1[idx]+weno->eps)) * (tau[idx]/(b1[idx]+weno->eps)) );  \
      a2[idx] = c2 * (1.0 + (tau[idx]/(b2[idx]+weno->eps)) * (tau[idx]/(b2[idx]+weno->eps)) );  \
      a3[idx] = c3 * (1.0 + (tau[idx]/(b3[idx]+weno->eps)) * (tau[idx]/(b3[idx]+weno->eps)) );  \
    }\
    double a_sum_inv[N]; \
    for (idx=0; idx<N; idx++) a_sum_inv[idx] = 1.0 / (a1[idx] + a2[idx] + a3[idx]); \
    _ArrayMultiply1D_(w1,a1,a_sum_inv,N); \
    _ArrayMultiply1D_(w2,a2,a_sum_inv,N); \
    _ArrayMultiply1D_(w3,a3,a_sum_inv,N); \
  }

#define _WENOWeights_v_YC_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators */\
    double b1[N], b2[N], b3[N]; \
    for (idx=0; idx<N; idx++) { \
      b1[idx] = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
                + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      b2[idx] = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
                + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      b3[idx] = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
                + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
    } \
    /* calculate the tau parameter for the WENO-Z and WENO-YC weights */\
    double tau[N]; \
    for (idx=0; idx<N; idx++) {\
      tau[idx] = (m3[idx]-4*m2[idx]+6*m1[idx]-4*p1[idx]+p2[idx])*(m3[idx]-4*m2[idx]+6*m1[idx]-4*p1[idx]+p2[idx]);  \
    } \
    /* calculate the WENO weights */\
    double a1[N], a2[N], a3[N];  \
    for (idx=0; idx<N; idx++) { \
      a1[idx] = c1 * (1.0 + (tau[idx]/(b1[idx]+weno->eps)) * (tau[idx]/(b1[idx]+weno->eps)) );  \
      a2[idx] = c2 * (1.0 + (tau[idx]/(b2[idx]+weno->eps)) * (tau[idx]/(b2[idx]+weno->eps)) );  \
      a3[idx] = c3 * (1.0 + (tau[idx]/(b3[idx]+weno->eps)) * (tau[idx]/(b3[idx]+weno->eps)) );  \
    }\
    double a_sum_inv[N]; \
    for (idx=0; idx<N; idx++) a_sum_inv[idx] = 1.0 / (a1[idx] + a2[idx] + a3[idx]); \
    _ArrayMultiply1D_(w1,a1,a_sum_inv,N); \
    _ArrayMultiply1D_(w2,a2,a_sum_inv,N); \
    _ArrayMultiply1D_(w3,a3,a_sum_inv,N); \
  }
