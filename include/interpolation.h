/* interpolation scheme definitions */
#define _FIRST_ORDER_UPWIND_    "1"
#define _SECOND_ORDER_CENTRAL_  "2"
#define _THIRD_ORDER_MUSCL_     "muscl3"
#define _FIFTH_ORDER_WENO_      "weno5"
#define _FIFTH_ORDER_CRWENO_    "crweno5"
#define _FIFTH_ORDER_HCWENO_    "hcweno5"

/* interpolation type definitions */
#define _CHARACTERISTIC_        "characteristic"
#define _COMPONENTS_            "components"  

/*
  One-dimensional Interpolation Functions
  Functions to interpolate the primitive of a function 
  at interfaces from its  cell-centered values.

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
int Interp1PrimFirstOrderUpwind           (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimThirdOrderMUSCL            (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderWENO             (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderCRWENO           (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderHCWENO           (double*,double*,double*,double*,int,int,void*,void*);

/* functions to interpolate the first primitive in a characteristic-based way
   (for conservative discretization of the 1st derivative) on a uniform grid */
int Interp1PrimFirstOrderUpwindChar       (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimThirdOrderMUSCLChar        (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderWENOChar         (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderCRWENOChar       (double*,double*,double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderHCWENOChar       (double*,double*,double*,double*,int,int,void*,void*);

/* functions to interpolate the second primitive 
   (for conservative discretization of the 2nd derivative) */
int Interp2PrimSecondOrder  (double*,double*,int,void*,void*);

/* other interpolation related functions */
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

  /* hybrid compact-WENO scheme related parameters */
  double  rc;
  double  xi;

  /* data arrays for CRWENO scheme */
  double *A, *B, *C, *R;
  double *sendbuf, *recvbuf;

  /* WENO weights */
  double *w1, *w2, *w3;
  int *offset;
  int size;

  /* function pointer to the weight calculation function */
  int (*CalculateWENOWeights) (double*,double*,double*,int,void*,void*);

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

/* macro to calculate fifth-order WENO weights */
#define _WENOWeights_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno) \
  { \
    if (weno->no_limiting) { \
      w1 = c1; \
      w2 = c2; \
      w3 = c3; \
    } else { \
      double b1, b2, b3; \
      b1 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) \
           + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);  \
      b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1) \
           + one_fourth*(m2-p1)*(m2-p1);                \
      b3 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2) \
           + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);  \
      double tau; \
      if      (weno->borges) tau = absolute(b3 - b1);  \
      else if (weno->yc)     tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);  \
      else                   tau = 0;  \
      double a1, a2, a3;  \
      if (weno->borges || weno->yc) { \
        a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));  \
        a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));  \
        a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));  \
      } else {  \
        a1 = c1 / raiseto(b1+weno->eps,weno->p);  \
        a2 = c2 / raiseto(b2+weno->eps,weno->p);  \
        a3 = c3 / raiseto(b3+weno->eps,weno->p);  \
      } \
      double a_sum_inv; \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1 = a1 * a_sum_inv;  \
      w2 = a2 * a_sum_inv;  \
      w3 = a3 * a_sum_inv;  \
      if (weno->mapped) { \
        a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1)); \
        a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2)); \
        a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3)); \
        a_sum_inv = 1.0 / (a1 + a2 + a3); \
        w1 = a1 * a_sum_inv;  \
        w2 = a2 * a_sum_inv;  \
        w3 = a3 * a_sum_inv;  \
      } \
    } \
  }
