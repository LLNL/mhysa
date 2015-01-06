#include <basic.h>

/* structure containing all solver-specific variables and functions */
typedef struct main_parameters {

  /* Number of dimensions */
  int     ndims;

  /* Number of variables or DoFs at a grid point */
  int     nvars;

  /* Global dimensions: array of size ndims containing the global grid size in each dimension */
  int     *dim_global;
  
  /* Local dimensions: array of size ndims containing the local grid size in each dimension */
  int     *dim_local;

  /* Global number of grid points 
   * (product of all the elements of dim_global) */
  int     npoints_global;

  /* Local number of grid points 
   * (product of all the elements of dim_local) */
  int     npoints_local, npoints_local_wghosts;

  /* Number of ghost points at the boundary - it's the same along all dimensions */
  int     ghosts;

  /* Number of time steps */
  int     n_iter;

  /* If restart run, time step iteration at which to restart. 0 -> not a restart run */
  int     restart_iter;

  /* time step size */
  double  dt;

  /*  choice of time integration class (eg RK) */
  char    time_scheme         [_MAX_STRING_SIZE_];

  /* specific time-integration scheme in that class (eg. rk44, ssprk3) */
  char    time_scheme_type    [_MAX_STRING_SIZE_];

  /* choice of spatial discretization scheme for the hyperbolic terms (eg: weno5, crweno5, muscl3) */
  char    spatial_scheme_hyp  [_MAX_STRING_SIZE_];

  /* type of reconstruction for spatial discretization of hyperbolic term 
   * (characteristic or component-wise) */
  char    interp_type         [_MAX_STRING_SIZE_];

  /* split the hyperbolic flux into two terms - for implicit-explicit time-integration or
   * for any other purpose */
  char    SplitHyperbolicFlux [_MAX_STRING_SIZE_];

  /* type of spatial discretization for the parabolic term 
   * conservative-1stage, nonconservative-1stage, or nonconservative-2stage */
  char    spatial_type_par    [_MAX_STRING_SIZE_];

  /* choice of spatial discretization scheme for the parabolic term */
  char    spatial_scheme_par  [_MAX_STRING_SIZE_];

  /* an ndims-dimensional integer array used to reference grid points */
  int    *index;

  /* the coordinate vector: one 1D array containing the grid coordinates along each dimension,
   * one dimension after the other.
   * Use _GetCoordinate_ in basic.h to access the grid coordinate at a specific grid point */
  double *x;                          /* coordinate vector                                */

  /* array containing (1.0/dx): layout same as that of x */
  double *dxinv;

  /* Solution vector: the ndims-dimensional solution vector with nvars components at each 
   * grid point is stored as a 1D array. **Include ghost points**
   * Use _ArrayIndex1D_ (in arrayfunctions.h) to calculate the index in the 1D array corres-
   * ponding to an ndims-dimensional index (i_0, i_1, i_2, ..., i_{ndims-1}) */
  double *u;

  /* Array to hold the discretized hyperbolic term: Same layout as u */
  double *hyp;

  /* Array to hold the discretized parabolic term: Same layout as u */
  double *par;
  
  /* Array to hold the source term: Same layout as u */
  double *source;

  /* Array to hold the cell-centered hyperbolic flux: Same layout as u */
  double *fluxC;
  
  /* Array to hold the cell-centered modified solution: Same layout as u */
  double *uC;
  
  /* Array to hold the interface hyperbolic flux for a conservative finite-difference method
   * Since number of interfaces is one more than the number of cell-centers, the dimensions 
   * of fluxI is accordingly increased during allocation. **Does not have ghost points.** */
  double *fluxI;
  
  /* arrays to hold the left and right biased recontructed solution and hyperbolic flux at
   * the interfaces. Same layout as fluxI, **No ghost points** */
  double *uL, *uR, *fL, *fR;

  /* arrays to hold approximations to the first and second derivatives to a given function
   * at grid points. Layout is same as u, hyp, par, source. **Includes ghost points** */
  double *Deriv1, *Deriv2;

  /* object to hold a sparse banded matrix - used for storing the Jacobian */
  void  *Jac;

  /* Boundary conditions: Number of boundary zones  */
  int   nBoundaryZones;
  /* Pointer to the boundary zones: boundary zone type is defined in boundaryconditions.h */
  void  *boundary;

  /* frequency (iterations) of writing iteration information (dt,CFL,norm,etc) to screen */
  int screen_op_iter;
  
  /* frequency (iterations) of writing solution to file */
  int file_op_iter;

  /* flag to control if residual is written to file */
  int write_residual;

  /* mode of reading in initial solution: serial, parallel or mpi-io */
  char input_mode    [_MAX_STRING_SIZE_];

  /* type of initial solution file: ascii or binary */
  char ip_file_type  [_MAX_STRING_SIZE_];

  /* mode of writing solution to file: serial or parallel */
  char output_mode   [_MAX_STRING_SIZE_];

  /* solution output file format: binary, text, tecplot2d, tecplot 3d  */
  char op_file_format[_MAX_STRING_SIZE_];

  /* overwrite solution file when writing new one? ("yes" - for steady solutions or if interested
   * only in the final solution; "no" - to keep the solutions written at intermediate time steps). */
  char op_overwrite  [_MAX_STRING_SIZE_];

  /* filename index for files written every few iterations */
  char *filename_index;
  int  index_length;
  /* solution filename extension */
  char solnfilename_extn[_MAX_STRING_SIZE_];

  /* Function to write the solution to file */
  int (*WriteOutput)              (int,int,int*,double*,double*,char*,int*);  

  /* Function to apply the physical boundary conditions to the solution */
  int (*ApplyBoundaryConditions)  (void*,void*,double*,double*,int,double);

  /* Function to integrate the solution in time */
  int (*TimeIntegrate)            (void*);

  /* Function to interpolate a function at the grid interfaces from the cell-centered values
   * for the hyperbolic flux */
  int (*InterpolateInterfacesHyp) (double*,double*,double*,double*,int,int,void*,void*,int);

  /* Function to pre-calculate the nonlinear interpolation coefficients for the hyperbolic
   * flux interpolation */
  int (*NonlinearInterp)    (double*,void*,void*,double,
                             int(*)(double*,double*,int,void*,double));

  /* function to calculate the non-linear interpolation coefficients of the given scheme */
  int (*SetInterpLimiterVar)      (double*,double*,double*,int,void*,void*);

  /* Function to interpolate a function at the grid interfaces from the cell-centered values
   * for the parabolic flux (needed for a conservative 1-stage discretization) */
  int (*InterpolateInterfacesPar) (double*,double*,int,void*,void*);

  /* Function to calculate the cell-centered first derivative of a given function */
  int (*FirstDerivativePar)       (double*,double*,int,int,void*,void*);

  /* Function to calculate the cell-centered second derivative of a given function */
  int (*SecondDerivativePar)      (double*,double*,int,void*,void*);

  /* Function to calculate the hyperbolic term */
  int (*HyperbolicFunction) (double*,double*,void*,void*,double,int,
                             int(*)(double*,double*,int,void*,double),
                             int(*)(double*,double*,double*,double*,double*,
                                    double*,int,void*,double));

  /* Function to calculate the parabolic term */
  int (*ParabolicFunction)  (double*,double*,void*,void*,double);

  /* Function to calculate the source term */
  int (*SourceFunction)     (double*,double*,void*,void*,double);

  /* name of physical model (defined in the header files in folder physicalmodels) */
  char model[_MAX_STRING_SIZE_];

  /* object providing the physics of the PDE being solved */
  void *physics;

  /* Function to calculate the CFL number */
  double (*ComputeCFL)         (void*,void*,double,double);

  /* Function to calculate the diffusion number */
  double (*ComputeDiffNumber)  (void*,void*,double,double);

  /* Function to calculate the hyperbolic flux function */
  int    (*FFunction)          (double*,double*,int,void*,double);
  /* If hyperbolic flux is split, function to calculate the split part of
   * the hyperbolic flux function. The splitting is
   * F(u) = [F(u) - dF(u)] + dF(u) */ 
  int    (*dFFunction)         (double*,double*,int,void*,double);

  /* Function to calculate the upwind interface flux, given the left- and right-biased fluxes */
  int    (*Upwind)             (double*,double*,double*,double*,double*,double*,int,void*,double);
  /* Function to calculate the upwind interface split flux, given the left- and right-biased fluxes */
  int    (*UpwinddF)           (double*,double*,double*,double*,double*,double*,int,void*,double);

  /* Function to calculate the parabolic function */
  int    (*GFunction)          (double*,double*,int,void*,double);

  /* Function to calculate the parabolic function in the 2-stage treatment */
  int    (*HFunction)          (double*,double*,int,int,void*,double);

  /* Function to calculate the source function */
  int    (*SFunction)          (double*,double*,void*,void*,double);

  /* Function to calculate the modified solution for upwinding */
  int    (*UFunction)          (double*,double*,int,void*,void*,double);

  /* Function to calculate the Jacobian for implicit time-integration */
  int    (*JFunction)          (void*,double*,void*,void*,double,double);
  /* Function to calculate the preconditioning matrix for implicit time-integration */
  int    (*PFunction)          (void*,double*,void*,void*,double,double);

  /* Function to do some pre-time-integration-stage computations, if required */
  int    (*PreStage)           (int,double**,void*,void*,double);
  /* Function to do some post-time-integration-stage computations, if required */
  int    (*PostStage)          (int,double**,void*,void*,double);
  /* Function to do some pre-time-integration-step computations, if required */
  int    (*PreStep)            (double*,void*,void*,double);
  /* Function to do some post-time-integration-step computations, if required */
  int    (*PostStep)           (double*,void*,void*,double);
  /* Function to do print some physics-specific time-integration-step information,
   * if required */
  int    (*PrintStep)          (void*,void*,double);

  /* Function to calculate the averaged solution at the interface (provided by the physics) */
  int   (*AveragingFunction)   (double*,double*,double*,void*);
  
  /* Function to calculate the left eigenvectors at the interface (provided by the physics) */
  int   (*GetLeftEigenvectors) (double*,double*,void*,int);
  /* Function to calculate the right eigenvectors at the interface (provided by the physics) */
  int   (*GetRightEigenvectors)(double*,double*,void*,int);

  /* object containing interpolation-related parameters (of hyperbolic flux reconstruction */
  void *interp;
  /* object containing multi-stage time-integration (RK-type) related parameters */
  void *msti;
  /* object containing parameters for the tridiagonal solver */
  void *lusolver;

  /* Errors - L1, L2 and L_inf; calculated only if an exact solution file is supplied */
  double error[3];
  /* Conservation error in the solution -- does not indicate anything if parabolic and source terms 
   * are non-zero */
  double *ConservationError;
  /* check for conservation error? */
  char   ConservationCheck[_MAX_STRING_SIZE_];
  /* Volume integral of the solution over the global domain */
  double *VolumeIntegral;
  /* Volume integral of the initial solution over the global domain */
  double *VolumeIntegralInitial;
  /* Surface integral of the flux at the boundary for each time-integration stage */
  double *StageBoundaryIntegral;
  /* Surface integral of the flux at the boundary for a time-integration step */
  double *StepBoundaryIntegral;
  /* Total surface integral of the flux over the global domain boundary */
  double *TotalBoundaryIntegral;
  /* Function to calculate the volume integral of a given function */
  int    (*VolumeIntegralFunction)    (double*,double*,void*,void*);
  /* Function to calculate the boundary integral of the flux */
  int    (*BoundaryIntegralFunction)  (void*,void*);
  /* Function to calculate the conservation error */
  int    (*CalculateConservationError)(void*,void*);

#ifdef with_petsc
  /* PETSc */
  int     use_petscTS;  /* Use PETSc time-integration? */
  double  *uref;        /* copy of solution vector        */
  double  *rhsref;      /* copy of the RHS vector      */
  double  *rhs;         /* RHS vector                  */
#endif

  /* "hack-y" parameters */
  int flag_nonlinearinterp; /* flag to globally switch on/off non-linear interpolation */

  /* strides along each dimension for arrays with and without ghost padding */
  int *stride_with_ghosts;
  int *stride_without_ghosts;

  /* function call counts */
  int count_hyp,
      count_par,
      count_sou;
#ifdef with_petsc
  int count_RHSFunction,
      count_IFunction,
      count_IJacobian,
      count_IJacFunction;
#endif

} HyPar;

/* The following functions are called by main() */
int CalculateError          (void*,void*);  /* Calculate the error in the final solution */
int Cleanup                 (void*,void*);  /* Clean up: deallocate all arrays and objects */
int Initialize              (void*,void*);  /* Initialize the solver */
int InitializeBoundaries    (void*,void*);  /* Initialize the boundary conditions */
int InitializePhysics       (void*,void*);  /* Initialize the physics */
int InitializeSolvers       (void*,void*);  /* Initialize the solvers */
int InitialSolution         (void*,void*);  /* Read the initial solution */
int OutputSolution          (void*,void*);  /* Write solution to file */
int ReadInputs              (void*,void*);  /* Read the input parameters */
int Solve                   (void*,void*);  /* Solve the PDE - time-integration */
#ifdef with_petsc
int SolvePETSc            (void*,void*);    /* Solve the PDE using PETSc TS */
#endif


/* Some definitions - types of discretizations available 
   for the parabolic (2nd derivative) term                */
#define _NC_1STAGE_   "nonconservative-1stage"    /* Non-conservative, direct evaluation of the 2nd deriv               */
#define _NC_2STAGE_   "nonconservative-2stage"    /* Non-conservative, two-stage evaluation of the 2nd deriv            */
#define _NC_1_5STAGE_ "nonconservative-1.5stage"  /* Non-conservative, "1.5"-stage evaluation of the 2nd deriv  
                                                     + direct, 1-stage evaluation of the 2nd derivative in one variable
                                                     + 2-stage evaluation of cross-derivatives                          */
#define _CONS_1STAGE_ "conservative-1stage"       /* Conservative, direct evaluation of the 2nd deriv                   */

