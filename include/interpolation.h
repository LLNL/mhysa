/* definitions */
#define _FIRST_ORDER_UPWIND_    "1"
#define _FIFTH_ORDER_WENO_      "weno5"

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

/* functions to interpolate the first primitive (for 1st derivative) */
int Interp1PrimFirstOrderUpwind (double*,double*,int,int,void*,void*);
int Interp1PrimFifthOrderWENO   (double*,double*,int,int,void*,void*);

/* WENO scheme related parameters and functions */
typedef struct parameters_weno {
  /* Options related to the type of WENO scheme */
  int     mapped;		    /* Use mapped weights?                                    */
  int     borges;		    /* Use Borges' implementation of weights?                 */
  int     yc;		        /* Use Yamaleev-Carpenter implementation of weights?      */
  int     no_limiting;  /* Remove limiting -> 5th order polynomial interpolation  */
  double  eps;		      /* epsilon parameter                                      */
  double	p;			      /* p parameter                                            */
} WENOParameters;
int WENOInitialize(void*,void*);

