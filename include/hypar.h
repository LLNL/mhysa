#include <basic.h>

typedef struct main_parameters {

  /* Solver parameters */
  int     ndims;                      /* number of dimensions                             */
  int     nvars;                      /* number of components of the state vector         */
  int     *dim_global;                /* global dimensions                                */
  int     *dim_local;                 /* local dimensions                                 */
  int     npoints_global;             /* total number of points (= product of dim_global) */
  int     npoints_local;              /* total number of points (= product of dim_local ) */
  int     ghosts;                     /* number of ghost points                           */
  int     n_iter;                     /* number of time iterations                        */
  double  dt;                         /* time step size                                   */

  char    time_scheme[_MAX_STRING_SIZE_];       /* time-integration scheme class (eg. RK) */
  char    time_scheme_type[_MAX_STRING_SIZE_];  /* specific time-integration scheme type  */
  int     spatial_scheme_hyp;         /* spatial discretization scheme for hyperbolic term*/
  int     spatial_scheme_par;         /* spatial discretization scheme for parabolic  term*/

  /* Data arrays */
  int    *index;                      /* ndims-dimensional variable index                 */
  double *x;                          /* coordinate vector                                */
  double *dxinv;                      /* 1/dx                                             */
  double *u;                          /* state vector                                     */
  double *hyp;                        /* array to hold the hyperbolic terms               */
  double *par;                        /* array to hold the parabolic terms                */
  double *source;                     /* array to hold the source    terms                */

  /* Boundary conditions */
  int   nBoundaryZones;               /* number of boundary zones                         */
  void  *boundary;                    /* pointer to boundary zones                        */    

  /* I/O parameters */
  int screen_op_iter;                     /* frequency of screen output                   */
  int file_op_iter;                       /* frequency of file output                     */
  int write_residual;                     /* write residual to file                       */
  char op_file_format[_MAX_STRING_SIZE_]; /* output file format                           */
  char op_overwrite  [_MAX_STRING_SIZE_]; /* overwrite output file?                       */
  char op_filename   [_MAX_STRING_SIZE_]; /* output filename                              */

  /* Functions */
  int (*WriteOutput)              (int,int,int*,double*,double*,char*,int*);  
  int (*ApplyBoundaryConditions)  (void*,void*,double*);                     
  int (*TimeIntegrate)            (void*);                                  
  int (*InterpolateInterfacesHyp) (double*,double*,int,int,void*,void*);
  int (*InterpolateInterfacesPar) ();
  int (*HyperbolicFunction)       (double*,double*,void*,void*);
  int (*ParabolicFunction)        (double*,double*,void*,void*);
  int (*SourceFunction)           (double*,double*,void*,void*);

  /* Physics  */
  char model[_MAX_STRING_SIZE_];          /* name of model, ie, linear advection, euler...*/
  void *physics;                          /* object containing the physics                */
  /* Physical model specific functions                                                    */
  /* These functions are mandatory; if not required in a particular model, 
      they should be set to NULL                                                          */
  double (*ComputeCFL)         (void*,void*,double);
  double (*ComputeDiffNumber)  (void*,void*,double);
  int    (*FFunction)          (double*,double*,int,void*);
  int    (*GFunction)          ();
  int    (*SFunction)          ();
  int    (*Upwind)             (double*,double*,double*,double*,int,void*);

  /* Other parameters */
  void *interp;                         /* Interpolation-related parameters         */
  void *msti;                           /* Multi-stage time-integration parameters  */

} HyPar;

/* Functions */
int Cleanup                 (void*,void*);
int Initialize              (void*,void*);
int InitializeBoundaries    (void*,void*);
int InitializePhysics       (void*,void*);
int InitializeSolvers       (void*,void*);
int InitialSolution         (void*,void*);
int OutputSolution          (void*,void*);
int ReadInputs              (void*,void*);
int Solve                   (void*,void*);
