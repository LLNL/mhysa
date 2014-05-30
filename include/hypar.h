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
  int     restart_iter;               /* restart from this iteration number               */
  double  dt;                         /* time step size                                   */

  char    time_scheme         [_MAX_STRING_SIZE_];/* time-integration scheme class (eg. RK)           */
  char    time_scheme_type    [_MAX_STRING_SIZE_];/* specific time-integration scheme type            */
  char    spatial_scheme_hyp  [_MAX_STRING_SIZE_];/* spatial discretization scheme for hyperbolic term*/
  char    spatial_type_par    [_MAX_STRING_SIZE_];/* spatial discretization type for parabolic  term  */
  char    spatial_scheme_par  [_MAX_STRING_SIZE_];/* spatial discretization scheme for parabolic  term*/
  char    interp_type         [_MAX_STRING_SIZE_];/* type of interpolation - characteristic 
                                                     or component-wise                                */
  char    SplitHyperbolicFlux [_MAX_STRING_SIZE_];/* split the hyperbolic flux into two terms         */

  /* Data arrays */
  int    *index;                      /* ndims-dimensional variable index                 */
  double *x;                          /* coordinate vector                                */
  double *dxinv;                      /* 1/dx                                             */
  double *u;                          /* state vector                                     */
  double *hyp;                        /* array to hold the hyperbolic terms               */
  double *par;                        /* array to hold the parabolic terms                */
  double *source;                     /* array to hold the source    terms                */
  /* arrays to hold temporary data during computations */
  double *fluxC, *fluxI, *Deriv1, *Deriv2;
  double *uL, *uR, *fL, *fR;

  /* Boundary conditions */
  int   nBoundaryZones;               /* number of boundary zones                         */
  void  *boundary;                    /* pointer to boundary zones                        */    

  /* I/O parameters */
  int screen_op_iter;                     /* frequency of screen output                   */
  int file_op_iter;                       /* frequency of file output                     */
  int write_residual;                     /* write residual to file                       */
  char ip_file_type  [_MAX_STRING_SIZE_]; /* whether initial solution file is ascii or bin*/
  char input_mode    [_MAX_STRING_SIZE_]; /* initial solution read in serial or parallel  */
  char op_file_format[_MAX_STRING_SIZE_]; /* output file format                           */
  char op_overwrite  [_MAX_STRING_SIZE_]; /* overwrite output file?                       */
  char op_filename   [_MAX_STRING_SIZE_]; /* output filename                              */

  /* Functions */
  int (*WriteOutput)              (int,int,int*,double*,double*,char*,int*);  
  int (*ApplyBoundaryConditions)  (void*,void*,double*,double*,int,double);                     
  int (*TimeIntegrate)            (void*);                                  
  int (*InterpolateInterfacesHyp) (double*,double*,double*,double*,int,int,void*,void*);
  int (*InterpolateInterfacesPar) (double*,double*,int,void*,void*);
  int (*SetInterpLimiterVar)      (double*,double*,double*,int,void*,void*);
  int (*FirstDerivativePar)       (double*,double*,int,void*,void*);
  int (*SecondDerivativePar)      (double*,double*,int,void*,void*);

  /* right hand side functions */
  int (*HyperbolicFunction) (double*,double*,void*,void*,double,int,
                             int(*)(double*,double*,int,void*,double)); /* hyperbolic terms */
  int (*ParabolicFunction)  (double*,double*,void*,void*,double);       /* parabolic terms  */
  int (*SourceFunction)     (double*,double*,void*,void*,double);       /* source terms     */

  /* Physics  */
  char model[_MAX_STRING_SIZE_];          /* name of model, ie, linear advection, euler...*/
  void *physics;                          /* object containing the physics                */
  /* Physical model specific functions                                                    */
  /* These functions are mandatory; if not required in a particular model, 
      they should be set to NULL                                                          */
  double (*ComputeCFL)         (void*,void*,double,double);
  double (*ComputeDiffNumber)  (void*,void*,double,double);
  int    (*FFunction)          (double*,double*,int,void*,double);
  int    (*dFFunction)         (double*,double*,int,void*,double);
  int    (*GFunction)          (double*,double*,int,void*,double);
  int    (*HFunction)          (double*,double*,int,int,void*,double);
  int    (*SFunction)          (double*,double*,void*,double);
  int    (*Upwind)             (double*,double*,double*,double*,double*,double*,
                                int,void*,double);
  /* physics-specific pre/post-time-step/stage functions */
  int    (*PreStage)           (int,double**,void*,void*,double);
  int    (*PostStage)          (int,double**,void*,void*,double);
  int    (*PreStep)            (double*,void*,void*,double);
  int    (*PostStep)           (double*,void*,void*,double);
  int    (*PrintStep)          (void*,void*,double);

  /* Physics-specific functions for characteristic-based interpolation */
  int   (*AveragingFunction)   (double*,double*,double*,void*);
  int   (*GetLeftEigenvectors) (double*,double*,void*,int);
  int   (*GetRightEigenvectors)(double*,double*,void*,int);

  /* Other parameters */
  void *interp;       /* Interpolation-related parameters         */
  void *msti;         /* Multi-stage time-integration parameters  */
  void *lusolver;     /* Tridiagonal LU solver parameters         */

  /* Errors */
  double error[3];                /* L1,L2,Linf errors, if calculated                             */
  char   ConservationCheck[_MAX_STRING_SIZE_];
  double *VolumeIntegral;         /* volume integral of the solution variable                     */
  double *VolumeIntegralInitial;  /* initial volume integral of the solution variable             */
  double *StageBoundaryIntegral;  /* boundary integral of flux for a time-integration stage       */
  double *StepBoundaryIntegral;   /* boundary integral of flux for a time-integration step        */
  double *TotalBoundaryIntegral;  /* total boundary integral of flux                              */
  int    (*VolumeIntegralFunction)    (double*,double*,void*,void*);
  int    (*BoundaryIntegralFunction)  (void*,void*);

#ifdef with_petsc
  /* PETSc */
  int     use_petscTS;  /* Use PETSc time-integration? */
  double  *uref;        /* copy of state vector        */
  double  *rhsref;      /* copy of the RHS vector      */
  double  *rhs;         /* RHS vector                  */
#endif

} HyPar;

/* Functions */
int CalculateError          (void*,void*);
int Cleanup                 (void*,void*);
int Initialize              (void*,void*);
int InitializeBoundaries    (void*,void*);
int InitializePhysics       (void*,void*);
int InitializeSolvers       (void*,void*);
int InitialSolution         (void*,void*);
int OutputSolution          (void*,void*);
int ReadInputs              (void*,void*);
int Solve                   (void*,void*);
#ifdef with_petsc
int SolvePETSc            (void*,void*);
#endif


/* Some definitions - types of discretizations available 
   for the parabolic (2nd derivative) term                */
#define _NC_1STAGE_   "nonconservative-1stage"/* Non-conservative, direct evaluation of the 2nd deriv  */
#define _NC_2STAGE_   "nonconservative-2stage"/* Non-conservative, two-stage evaluation of the 2nd deriv  */
#define _CONS_1STAGE_ "conservative-1stage"   /* Conservative, direct evaluation of the 2nd deriv      */

