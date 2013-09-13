#include <basic.h>

typedef struct solver_parameters {

  /* Solver parameters */
  int     ndims;                      /* number of dimensions                             */
  int     nvars;                      /* number of components of the state vector         */
  int     *dim_global;                /* global dimensions                                */
  int     *dim_local;                 /* local dimensions                                 */
  int     npoints_global;             /* total number of points (= product of dim_global) */
  int     npoints_local;              /* total number of points (= product of dim_local ) */
  int     ghosts;                     /* number of ghost points                           */
  int     n_iter;                     /* number of time iterations                        */
  int     hyp_space_scheme;           /* spatial discretization term for hyperbolic term  */
  int     par_space_scheme;           /* spatial discretization term for parabolic term   */
  int     time_scheme;                /* time-integration scheme                          */
  double  dt;                         /* time step size                                   */

  /* Data arrays */
  int    *index;                          /* ndims-dimensional variable index             */
  double *x;                              /* coordinate vector                            */
  double *u;                              /* state vector                                 */

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
  int (*WriteOutput)  (int,int,int*,double*,double*,char*,int*);  /* write data to file    */
  int (*TimeIntegrate)(void*);                                    /* time integration      */

  /* Physics  */
  void *physics;                          /* object containing the physics of the case    */


} HyPar;

/* Functions */
int ReadInputs            (void*,void*);
int Initialize            (void*,void*);
int InitialSolution       (void*,void*);
int InitializeBoundaries  (void*,void*);
int InitializeSolvers     (void*,void*);
int Solve                 (void*,void*);
int OutputSolution        (void*,void*);
int Cleanup               (void*,void*);
