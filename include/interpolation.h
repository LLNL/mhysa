/*! @file interpolation.h
    @brief Definitions for the functions computing the interpolated value of the primitive at the cell interfaces from the cell-centered values.
    @author Debojyoti Ghosh
*/

/*! First order upwind scheme */
#define _FIRST_ORDER_UPWIND_    "1"
/*! Second order central scheme */
#define _SECOND_ORDER_CENTRAL_  "2"
/*! Third order MUSCL scheme with Koren's limiter */
#define _THIRD_ORDER_MUSCL_     "muscl3"
/*! Fifth order Weighted Essentially Non-Oscillatory (WENO) scheme */
#define _FIFTH_ORDER_WENO_      "weno5"
/*! Fifth order Compact Reconstruction Weighted Essentially Non-Oscillatory (CRWENO) scheme */
#define _FIFTH_ORDER_CRWENO_    "crweno5"
/*! Fifth order hybrid compact-WENO scheme */
#define _FIFTH_ORDER_HCWENO_    "hcweno5"

/* interpolation type definitions */
#define _CHARACTERISTIC_        "characteristic" /*!< Characteristic-based interpolation of vectors (Physical model must define left and right eigenvectors) */
#define _COMPONENTS_            "components"     /*!< Component-wise interpolation of vectors      */

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

    uflag     int       flag to indicate where the flux function or the solution function
                        is being interpolated
                        (1 -> u; 0 -> f(u) )

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
/*! Component-wise interpolation of the first primitive at the cell interfaces using the first-order upwind scheme */
int Interp1PrimFirstOrderUpwind           (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Component-wise interpolation of the first primitive at the cell interfaces using the second-order central scheme */
int Interp1PrimSecondOrderCentral         (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Component-wise interpolation of the first primitive at the cell interfaces using the third-order MUSCL scheme */
int Interp1PrimThirdOrderMUSCL            (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Component-wise interpolation of the first primitive at the cell interfaces using the fifth-order WENO scheme */
int Interp1PrimFifthOrderWENO             (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Component-wise interpolation of the first primitive at the cell interfaces using the fifth-order CRWENO scheme */
int Interp1PrimFifthOrderCRWENO           (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Component-wise interpolation of the first primitive at the cell interfaces using the fifth-order hybrid-compact WENO scheme */
int Interp1PrimFifthOrderHCWENO           (double*,double*,double*,double*,int,int,void*,void*,int);

/* functions to interpolate the first primitive in a characteristic-based way
   (for conservative discretization of the 1st derivative) on a uniform grid */
/*! Characteristic-based interpolation of the first primitive at the cell interfaces using the first-order upwind scheme */
int Interp1PrimFirstOrderUpwindChar       (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Characteristic-based interpolation of the first primitive at the cell interfaces using the second-order central scheme */
int Interp1PrimSecondOrderCentralChar     (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Characteristic-based interpolation of the first primitive at the cell interfaces using the third-order MUSCL scheme */
int Interp1PrimThirdOrderMUSCLChar        (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Characteristic-based interpolation of the first primitive at the cell interfaces using the fifth-order WENO scheme */
int Interp1PrimFifthOrderWENOChar         (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Characteristic-based interpolation of the first primitive at the cell interfaces using the fifth-order CRWENO scheme */
int Interp1PrimFifthOrderCRWENOChar       (double*,double*,double*,double*,int,int,void*,void*,int);
/*! Characteristic-based interpolation of the first primitive at the cell interfaces using the fifth-order hybrid-compact WENO scheme */
int Interp1PrimFifthOrderHCWENOChar       (double*,double*,double*,double*,int,int,void*,void*,int);

/* functions to interpolate the second primitive 
   (for conservative discretization of the 2nd derivative) */
/*! Interpolation of the second primitive at the cell interfaces using the second-order central scheme */
int Interp2PrimSecondOrder  (double*,double*,int,void*,void*);

/*! Function to calculate and save the nonlinear interpolation coefficients for a non-linear (solution-dependent) interpolation scheme (eg. #_FIFTH_ORDER_WENO_, #_FIFTH_ORDER_CRWENO_, #_FIFTH_ORDER_HCWENO_, etc). */
int InterpSetLimiterVar(double*,double*,double*,int,void*,void*);

/*! \def MUSCLParameters
    \brief Structure of variables/parameters needed by the MUSCL scheme
 * This structure contains the variables/parameters needed by the MUSCL scheme.
*/
/*! \brief Structure of variables/parameters needed by the MUSCL scheme
 *
 * This structure contains the variables/parameters needed by the MUSCL scheme.
*/
typedef struct paramters_muscl {
  double eps; /*!< Epsilon parameter for the limiter */
} MUSCLParameters;
int MUSCLInitialize(void*,void*);

/*! \def WENOParameters
    \brief Structure of variables/parameters needed by the WENO-type scheme
 * This structure contains the variables/parameters needed by the WENO-type scheme (#_FIFTH_ORDER_WENO_, #_FIFTH_ORDER_CRWENO_, #_FIFTH_ORDER_HCWENO_).
*/
/*! \brief Structure of variables/parameters needed by the WENO-type scheme
 *
 * This structure contains the variables/parameters needed by the WENO-type scheme (#_FIFTH_ORDER_WENO_, #_FIFTH_ORDER_CRWENO_, #_FIFTH_ORDER_HCWENO_).
*/
typedef struct parameters_weno {
  int     mapped;		    /*!< Use mapped weights? (Henrick, Aslam, J. Comput. Phys., 2005) */
  int     borges;		    /*!< Use Borges' implementation of weights? (Borges, et. al, J. Comput. Phys., 2008) */
  int     yc;		        /*!< Use Yamaleev-Carpenter implementation of weights? (Yamaleev, Carpenter, J. Comput. Phys., 2009) */
  int     no_limiting;  /*!< Remove limiting -> 5th order polynomial interpolation (freeze the WENO weights to the optimal coefficients)  */
  double  eps;		      /*!< epsilon parameter */
  double	p;			      /*!< p parameter */
  double  tol;          /*!< a general tolerance parameter */

  /* hybrid compact-WENO scheme related parameters 
   * **References**: 
   * + http://dx.doi.org/10.1006/jcph.2002.7021
   * + http://dx.doi.org/10.1016/j.jcp.2003.07.006
  */
  double  rc, /*!< Parameter for the hybrid compact-WENO scheme */
          xi; /*!< Parameter for the hybrid compact-WENO scheme */

  /* CRWENO scheme: sub-, main-, and super-diagonals
   * and the right-hand-side */
  double *A, /*!< Array to save the sub-diagonal of the tridiagonal system resulting from the fifth-order CRWENO scheme */
         *B, /*!< Array to save the diagonal of the tridiagonal system resulting from the fifth-order CRWENO scheme */
         *C, /*!< Array to save the super-diagonal of the tridiagonal system resulting from the fifth-order CRWENO scheme */
         *R;/*!< Array to save the right-hand-side of the tridiagonal system resulting from the fifth-order CRWENO scheme */
  /* CRWENO scheme: buffer arrays for sending and receiving data */
  double *sendbuf, /*!< Buffer array to send data across processors */
         *recvbuf; /*!< Buffer array to receive data across processors */

  /* Arrays to save the WENO weights */
  double *w1, /*!< Array to save the first WENO weight */
         *w2, /*!< Array to save the second WENO weight */
         *w3;/*!< Array to save the third WENO weight */
  /* size and offset for the WENO weights arrays */
  int *offset /*! Array containing the offset information for the WENO weights */, 
      size /*! Size of the WENO weights array */;

} WENOParameters;
/*! Initialize the structure containing variables and parameters for WENO-type schemes */
int WENOInitialize(void*,void*,char*,char*); 
/*! Clean up the structure containing variables and parameters for WENO-type schemes */
int WENOCleanup(void*);

/* define optimal weights */
/*! Optimal value for the first fifth-order WENO weight */
#define   _WENO_OPTIMAL_WEIGHT_1_   0.1
/*! Optimal value for the second fifth-order WENO weight */
#define   _WENO_OPTIMAL_WEIGHT_2_   0.6
/*! Optimal value for the third fifth-order WENO weight */
#define   _WENO_OPTIMAL_WEIGHT_3_   0.3
/*! Optimal value for the first fifth-order CRWENO weight */
#define   _CRWENO_OPTIMAL_WEIGHT_1_ 0.2
/*! Optimal value for the second fifth-order CRWENO weight */
#define   _CRWENO_OPTIMAL_WEIGHT_2_ 0.5
/*! Optimal value for the third fifth-order CRWENO weight */
#define   _CRWENO_OPTIMAL_WEIGHT_3_ 0.3

/*! \def _WENOWeights_
 * Macro to calculate the fifth-order WENO weights
 * \n
 * Arguments:-\n
 *  + \a w1, \a w2,\a w3 are the nonlinear WENO weights.\n
 *  + \a c1, \a c2,\a c3 are optimal coefficients.\n
 *  + \a m3, \a m2,\a m1,\a p1,\a p2 are the function values at stencil points corresponding to the interface i-1/2: i-3,i-2,i-1,i,i+1\n
 *  + \a weno is an object of type #WENOParameters containing parameters for the WENO method.\n
 *\n
 * **References**:\n
 *  + http://dx.doi.org/10.1006/jcph.1996.0130\n
 *  + http://dx.doi.org/10.1016/j.jcp.2005.01.023\n
 *  + http://dx.doi.org/10.1016/j.jcp.2007.11.038\n
 *  + http://dx.doi.org/10.1016/j.jcp.2009.03.002\n
 *  + http://dx.doi.org/10.1007/s10915-014-9818-0\n
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

/*! \def _WENOWeights_v_JS_
 * Compute the WENO weights according the the Jiang & Shu formulation.
 * \n
 * Arguments:-\n
 *  + \a w1, \a w2,\a w3 are the nonlinear WENO weights.\n
 *  + \a c1, \a c2,\a c3 are optimal coefficients.\n
 *  + \a m3, \a m2,\a m1,\a p1,\a p2 are the function values at stencil points corresponding to the interface i-1/2: i-3,i-2,i-1,i,i+1\n
 *  + \a weno is an object of type #WENOParameters containing parameters for the WENO method.\n
 *  + \a N is the number of interfaces along the grid line on which this WENO-type reconstruction is happening.\n
 *\n
 * **Reference**:\n
 *  + Jiang, Shu, J. Comput. Phys., 1996. http://dx.doi.org/10.1006/jcph.1996.0130\n
 */
#define _WENOWeights_v_JS_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators and the WENO weights */\
    double b1, b2, b3, a1, a2, a3, a_sum_inv; \
    for (idx=0; idx<N; idx++) { \
      b1 = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
           + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      a1 = c1 / ( (b1+weno->eps) * (b1+weno->eps) );  \
      b2 = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
           + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      a2 = c2 / ( (b2+weno->eps) * (b2+weno->eps) );  \
      b3 = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
           + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
      a3 = c3 / ( (b3+weno->eps) * (b3+weno->eps) );  \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1[idx] = a1 * a_sum_inv; \
      w2[idx] = a2 * a_sum_inv; \
      w3[idx] = a3 * a_sum_inv; \
    } \
  }

/*! \def _WENOWeights_v_M_
 * Compute the WENO weights according the the Mapped-WENO formulation.
 * \n
 * Arguments:-\n
 *  + \a w1, \a w2,\a w3 are the nonlinear WENO weights.\n
 *  + \a c1, \a c2,\a c3 are optimal coefficients.\n
 *  + \a m3, \a m2,\a m1,\a p1,\a p2 are the function values at stencil points corresponding to the interface i-1/2: i-3,i-2,i-1,i,i+1\n
 *  + \a weno is an object of type #WENOParameters containing parameters for the WENO method.\n
 *  + \a N is the number of interfaces along the grid line on which this WENO-type reconstruction is happening.\n
 *\n
 * **Reference**:\n
 *  + Henrick, Aslam, Powers, J. Comput. Phys., 2005. http://dx.doi.org/10.1016/j.jcp.2005.01.023\n
 */
#define _WENOWeights_v_M_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators and the WENO weights */\
    double b1, b2, b3, a1, a2, a3, a_sum_inv; \
    for (idx=0; idx<N; idx++) { \
      b1 = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
           + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      a1 = c1 / ( (b1+weno->eps) * (b1+weno->eps) );  \
      b2 = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
           + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      a2 = c2 / ( (b2+weno->eps) * (b2+weno->eps) );  \
      b3 = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
           + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
      a3 = c3 / ( (b3+weno->eps) * (b3+weno->eps) );  \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1[idx] = a1 * a_sum_inv; \
      w2[idx] = a2 * a_sum_inv; \
      w3[idx] = a3 * a_sum_inv; \
      /* apply mapping if required */\
      a1 = w1[idx] * (c1 + c1*c1 - 3*c1*w1[idx] + w1[idx]*w1[idx]) / (c1*c1 + w1[idx]*(1.0-2.0*c1)); \
      a2 = w2[idx] * (c2 + c2*c2 - 3*c2*w2[idx] + w2[idx]*w2[idx]) / (c2*c2 + w2[idx]*(1.0-2.0*c2)); \
      a3 = w3[idx] * (c3 + c3*c3 - 3*c3*w3[idx] + w3[idx]*w3[idx]) / (c3*c3 + w3[idx]*(1.0-2.0*c3)); \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1[idx] = a1 * a_sum_inv; \
      w2[idx] = a2 * a_sum_inv; \
      w3[idx] = a3 * a_sum_inv; \
    } \
  }

/*! \def _WENOWeights_v_Z_
 * Compute the WENO weights according the the WENO-Z formulation.
 * \n
 * Arguments:-\n
 *  + \a w1, \a w2,\a w3 are the nonlinear WENO weights.\n
 *  + \a c1, \a c2,\a c3 are optimal coefficients.\n
 *  + \a m3, \a m2,\a m1,\a p1,\a p2 are the function values at stencil points corresponding to the interface i-1/2: i-3,i-2,i-1,i,i+1\n
 *  + \a weno is an object of type #WENOParameters containing parameters for the WENO method.\n
 *  + \a N is the number of interfaces along the grid line on which this WENO-type reconstruction is happening.\n
 *\n
 * **Reference**:\n
 *  + Borges, et. al., J. Comput. Phys., 2008. http://dx.doi.org/10.1016/j.jcp.2007.11.038\n
 */
#define _WENOWeights_v_Z_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators and the WENO weights */\
    double b1, b2, b3, a1, a2, a3, a_sum_inv, tau; \
    for (idx=0; idx<N; idx++) { \
      b1 = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
           + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      b2 = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
           + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      b3 = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
           + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
      tau = absolute(b3 - b1);  \
      a1 = c1 * (1.0 + (tau/(b1+weno->eps)) * (tau/(b1+weno->eps)) );  \
      a2 = c2 * (1.0 + (tau/(b2+weno->eps)) * (tau/(b2+weno->eps)) );  \
      a3 = c3 * (1.0 + (tau/(b3+weno->eps)) * (tau/(b3+weno->eps)) );  \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1[idx] = a1 * a_sum_inv; \
      w2[idx] = a2 * a_sum_inv; \
      w3[idx] = a3 * a_sum_inv; \
    } \
  }

/*! \def _WENOWeights_v_YC_
 * Compute the WENO weights according the the ESWENO formulation.
 * \n
 * Arguments:-\n
 *  + \a w1, \a w2,\a w3 are the nonlinear WENO weights.\n
 *  + \a c1, \a c2,\a c3 are optimal coefficients.\n
 *  + \a m3, \a m2,\a m1,\a p1,\a p2 are the function values at stencil points corresponding to the interface i-1/2: i-3,i-2,i-1,i,i+1\n
 *  + \a weno is an object of type #WENOParameters containing parameters for the WENO method.\n
 *  + \a N is the number of interfaces along the grid line on which this WENO-type reconstruction is happening.\n
 *\n
 * **Reference**:\n
 *  + Yamaleev, Carpenter, J. Comput. Phys., 2009. http://dx.doi.org/10.1016/j.jcp.2009.03.002\n
 */
#define _WENOWeights_v_YC_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno,N) \
  { \
    int idx; \
    /* calculate smoothness indicators and the WENO weights */\
    double b1, b2, b3, a1, a2, a3, a_sum_inv, tau; \
    for (idx=0; idx<N; idx++) { \
      b1 = thirteen_by_twelve*(m3[idx]-2*m2[idx]+m1[idx])*(m3[idx]-2*m2[idx]+m1[idx]) \
           + one_fourth*(m3[idx]-4*m2[idx]+3*m1[idx])*(m3[idx]-4*m2[idx]+3*m1[idx]);  \
      b2 = thirteen_by_twelve*(m2[idx]-2*m1[idx]+p1[idx])*(m2[idx]-2*m1[idx]+p1[idx]) \
           + one_fourth*(m2[idx]-p1[idx])*(m2[idx]-p1[idx]);                \
      b3 = thirteen_by_twelve*(m1[idx]-2*p1[idx]+p2[idx])*(m1[idx]-2*p1[idx]+p2[idx]) \
           + one_fourth*(3*m1[idx]-4*p1[idx]+p2[idx])*(3*m1[idx]-4*p1[idx]+p2[idx]);  \
      tau = (m3[idx]-4*m2[idx]+6*m1[idx]-4*p1[idx]+p2[idx])*(m3[idx]-4*m2[idx]+6*m1[idx]-4*p1[idx]+p2[idx]);  \
      a1 = c1 * (1.0 + (tau/(b1+weno->eps)) * (tau/(b1+weno->eps)) );  \
      a2 = c2 * (1.0 + (tau/(b2+weno->eps)) * (tau/(b2+weno->eps)) );  \
      a3 = c3 * (1.0 + (tau/(b3+weno->eps)) * (tau/(b3+weno->eps)) );  \
      a_sum_inv = 1.0 / (a1 + a2 + a3); \
      w1[idx] = a1 * a_sum_inv; \
      w2[idx] = a2 * a_sum_inv; \
      w3[idx] = a3 * a_sum_inv; \
    } \
  }
