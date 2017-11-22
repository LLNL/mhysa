/*! @file hypar.h
    @brief Contains structure definition for hypar.
    @author Debojyoti Ghosh
 */

#include <basic.h>

/*! \def HyPar
    \brief Structure containing all solver-specific variables and functions 
 * This structure contains all the variables and function pointers for the 
 * main PDE solver.
 *
*/

/*! \brief Structure containing all solver-specific variables and functions 
 *
 * This structure contains all the variables and function pointers for the 
 * main PDE solver.
*/
typedef struct main_parameters {

  /*! Number of spatial/coordinate dimensions (input - \b solver.inp ) */
  int     ndims;

  /*! Number of variables or DoFs at a grid point (input - \b solver.inp ) */
  int     nvars;

  /*! Global dimensions: array of size #HyPar::ndims containing the global grid size in each spatial/coordinate dimension  (input - \b solver.inp )*/
  int     *dim_global;
  
  /*! Local dimensions: array of size #HyPar::ndims containing the local grid size in each spatial/coordinate dimension (computed, based on the number of processors) */
  int     *dim_local;

  /*! Global number of grid points (product of all the elements of dim_global) */
  int     npoints_global;

  int     npoints_local,          /*!< Local number of grid points (product of all the elements of dim_local) */
          npoints_local_wghosts;  /*!< Local number of grid points with ghost points (product of the [elements of dim_local + 2*#HyPar::ghosts]) */

  /*! Number of ghost points at the boundary - it's the same along all dimensions (input - \b solver.inp ) */
  int     ghosts;

  /*! Number of time steps (input - \b solver.inp ) */
  int     n_iter;

  /*! If restart run, time step iteration at which to restart. 0 -> not a restart run (input - \b solver.inp ) */
  int     restart_iter;

  /*! time step size (input - \b solver.inp ) */
  double  dt;

  /*!  choice of time integration class (eg RK) (input - \b solver.inp ) */
  char    time_scheme         [_MAX_STRING_SIZE_];

  /*! specific time-integration scheme in that class (eg. rk44, ssprk3) (input - \b solver.inp ) */
  char    time_scheme_type    [_MAX_STRING_SIZE_];

  /*! choice of spatial discretization scheme for the hyperbolic terms (eg: weno5, crweno5, muscl3) (input - \b solver.inp ) */
  char    spatial_scheme_hyp  [_MAX_STRING_SIZE_];

  /*! type of reconstruction for spatial discretization of hyperbolic term 
   * (characteristic or component-wise) (input - \b solver.inp ) */
  char    interp_type         [_MAX_STRING_SIZE_];

  /*! split the hyperbolic flux into two terms - for implicit-explicit time-integration or
   * for any other purpose (input - \b solver.inp ) */
  char    SplitHyperbolicFlux [_MAX_STRING_SIZE_];

  /*! type of spatial discretization for the parabolic term 
   * conservative-1stage, nonconservative-1stage, or nonconservative-2stage (input - \b solver.inp )*/
  char    spatial_type_par    [_MAX_STRING_SIZE_];

  /*! choice of spatial discretization scheme for the parabolic term (input - \b solver.inp ) */
  char    spatial_scheme_par  [_MAX_STRING_SIZE_];

  /*! a #HyPar::ndims-dimensional integer array used to reference grid points */
  int    *index;

  /*! the coordinate vector: one 1D array containing the spatial coordinates along each dimension
   * of the grid points, one dimension after the other.
   * Use #_GetCoordinate_ to access the spatial coordinate at a specific grid point */
  double *x;                          

  /*! array containing (1.0/dx): layout same as that of x */
  double *dxinv;

  /*! Solution vector: the #HyPar::ndims-dimensional solution vector with nvars components at each 
   * grid point is stored as a 1D array. **Includes ghost points**
   * Use #_ArrayIndex1D_ to calculate the index in the 1D array 
   * corresponding to a #HyPar::ndims-dimensional index (i_0, i_1, i_2, ..., i_{ndims-1}) */
  double *u;

  /*! Array to hold the discretized hyperbolic term: Same layout as u */
  double *hyp;

  /*! Array to hold the discretized parabolic term: Same layout as u */
  double *par;
  
  /*! Array to hold the source term: Same layout as u */
  double *source;

  /*! Array to hold the cell-centered hyperbolic flux: Same layout as u */
  double *fluxC;
  
  /*! Array to hold the cell-centered modified solution: Same layout as u */
  double *uC;
  
  /*! Array to hold the interface hyperbolic flux for a conservative finite-difference method
   * Since number of interfaces is one more than the number of cell-centers, the dimensions 
   * of fluxI is accordingly increased during allocation. **Does not have ghost points.** */
  double *fluxI;
  
  double *uL, /*!< Array to hold the left-biased reconstructed solution 
                   at the interface. Same layout as fluxI. **No ghost points** */
         *uR, /*!< Array to hold the right-biased reconstructed solution 
                   at the interface. Same layout as fluxI. **No ghost points** */
         *fL, /*!< Array to hold the left-biased reconstructed hyperbolic flux 
                   at the interface. Same layout as fluxI. **No ghost points** */
         *fR; /*!< Array to hold the right-biased reconstructed hyperbolic flux 
                   at the interface. Same layout as fluxI. **No ghost points** */

  /*! arrays to hold approximations to the first and second derivatives to a given function
   * at grid points. Layout is same as u, hyp, par, source. **Includes ghost points** */
  double  *Deriv1, /*!< Array to hold approximations to the first derivative to 
                        a given function at grid points. Layout is same as u, hyp, 
                        par, source. **Includes ghost points** */
          *Deriv2; /*!< Array to hold approximations to the second derivative to 
                        a given function at grid points. Layout is same as u, hyp, 
                        par, source. **Includes ghost points** */

  /*! Boundary conditions: Number of boundary zones  */
  int   nBoundaryZones;
  /*! Pointer to the boundary zones: boundary zone type is defined in boundaryconditions.h */
  void  *boundary;
  /*! Pointer to array of size #HyPar::ndims: each element is 1 if the domain is periodic along
      that spatial dimension; zero otherwise. */
  int   *isPeriodic;

  /*! pointer to the time-integration object */
  void *time_integrator;

  /*! frequency (iterations) of writing iteration information (dt,CFL,norm,etc) to screen (input - \b solver.inp )*/
  int screen_op_iter;
  
  /*! frequency (iterations) of writing solution to file (input - \b solver.inp )*/
  int file_op_iter;

  /*! flag to control if residual is written to file (input - \b solver.inp )*/
  int write_residual;

  /*! mode of reading in initial solution: serial, parallel or mpi-io (input - \b solver.inp )*/
  char input_mode    [_MAX_STRING_SIZE_];

  /*! type of initial solution file: ascii or binary (input - \b solver.inp )*/
  char ip_file_type  [_MAX_STRING_SIZE_];

  /*! mode of writing solution to file: serial or parallel (input - \b solver.inp )*/
  char output_mode   [_MAX_STRING_SIZE_];

  /*! solution output file format: binary, text, tecplot2d, tecplot 3d  (input - \b solver.inp )*/
  char op_file_format[_MAX_STRING_SIZE_];

  /*! overwrite solution file when writing new one? ("yes" - for steady solutions or if interested
   * only in the final solution; "no" - to keep the solutions written at intermediate time steps). 
     (input - \b solver.inp )*/
  char op_overwrite  [_MAX_STRING_SIZE_];

  /*! filename index for files written every few iterations */
  char *filename_index;
  /*! length of filename_index - should be sufficient for the number of files expected to be written */
  int  index_length;
  /*! solution filename extension */
  char solnfilename_extn[_MAX_STRING_SIZE_];

  /*! Pointer to the function to write the solution to file, assigned in InitializeSolvers() */
  int (*WriteOutput)              (int,int,int*,double*,double*,char*,int*);  

  /*! Pointer to the function to apply the physical boundary conditions to the solution, assigned in InitializeSolvers() */
  int (*ApplyBoundaryConditions)  (void*,void*,double*,double*,double);

  /*! Pointer to the function to apply the immersed boundary conditions to the solution, assigned in InitializeSolvers() */
  int (*ApplyIBConditions)  (void*,void*,double*,double);

  /*! Pointer to the function to integrate the solution in time, assigned in InitializeSolvers() */
  int (*TimeIntegrate)            (void*);

  /*! Pointer to the function to interpolate a function at the grid interfaces from the cell-centered values
   * for the hyperbolic flux, assigned in InitializeSolvers() */
  int (*InterpolateInterfacesHyp) (double*,double*,double*,double*,int,int,void*,void*,int);

  /*! Pointer to the function to pre-calculate the nonlinear interpolation coefficients for the hyperbolic
   * flux interpolation, assigned in InitializeSolvers() */
  int (*NonlinearInterp)    (double*,void*,void*,double,
                             int(*)(double*,double*,int,void*,double));

  /*! Pointer to the function to calculate the non-linear interpolation coefficients of the given scheme, 
      assigned by the initialization function of the non-linear interpolation method 
      (eg. WENOInitialize() )*/
  int (*SetInterpLimiterVar)      (double*,double*,double*,int,void*,void*);

  /*! Pointer to the function to interpolate a function at the grid interfaces from the cell-centered values
   * for the parabolic flux (needed for a conservative 1-stage discretization),
     assigned in InitializeSolvers()*/
  int (*InterpolateInterfacesPar) (double*,double*,int,void*,void*);

  /*! Pointer to the function to calculate the cell-centered first derivative of a given function, for the
      evaluation of the parabolic term; assigned in InitializeSolvers()*/
  int (*FirstDerivativePar)       (double*,double*,int,int,void*,void*);

  /*! Pointer to the function to calculate the cell-centered second derivative of a given function, for the
      evaluation of the parabolic term; assigned in InitializeSolvers()*/
  int (*SecondDerivativePar)      (double*,double*,int,void*,void*);

  /*! Pointer to the function to calculate the hyperbolic term (assigned in InitializeSolvers()) */
  int (*HyperbolicFunction) (double*,double*,void*,void*,double,int,
                             int(*)(double*,double*,int,void*,double),
                             int(*)(double*,double*,double*,double*,double*,
                                    double*,int,void*,double));

  /*! Pointer to the function to calculate the parabolic term (assigned in InitializeSolvers())*/
  int (*ParabolicFunction)  (double*,double*,void*,void*,double);

  /*! Pointer to the function to calculate the source term (assigned in InitializeSolvers())*/
  int (*SourceFunction)     (double*,double*,void*,void*,double);

  /*! name of physical model (defined in the header files in folder physicalmodels) 
      (input - \b solver.inp ) */
  char model[_MAX_STRING_SIZE_];

  /*! object providing the physics of the PDE being solved */
  void *physics;

  /*! Pointer to the function to calculate the CFL number (assigned in the physical model initialization called from InitializePhysics()) */
  double (*ComputeCFL)         (void*,void*,double,double);

  /*! Pointer to the function to calculate the diffusion number (assigned in the physical model initialization called from InitializePhysics()) */
  double (*ComputeDiffNumber)  (void*,void*,double,double);

  /*! Pointer to the function to calculate the hyperbolic flux function (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*FFunction)          (double*,double*,int,void*,double);
  /*! If hyperbolic flux is split as \f$f\left(u\right) = \left[f\left(u\right)-df\left(u\right)\right] + df\left(u\right)\f$
      (see #HyPar::SplitHyperbolicFlux), 
      function to calculate \f$df\left(u\right)\f$ (assigned in the physical model initialization called from InitializePhysics()) */ 
  int    (*dFFunction)         (double*,double*,int,void*,double);
  /*! If hyperbolic flux is split as \f$f\left(u\right) = \left[f\left(u\right)-df\left(u\right)\right] + df\left(u\right)\f$, 
      (see #HyPar::SplitHyperbolicFlux), 
      function to calculate \f$\left[f-df\right]\left(u\right)\f$(assigned in the physical model initialization called from InitializePhysics()) .\n
   Specifying this is optional; if the physical model does not explicitly specify this, then it is computed by subtracting 
   HyPar::dFFunction() from HyPar::FFunction(). \sa #HyPar::flag_fdf_specified */ 
  int    (*FdFFunction)         (double*,double*,int,void*,double);
  /*! Flag indicating whether the physical model has explicitly specified the function to compute 
      \f$\left[f-df\right]\left(u\right)\f$ (HyPar::FdFFunction()) or not. Relevant if the hyperbolic flux
      is being partitioned. \sa #HyPar::SplitHyperbolicFlux, HyPar::dFFunction(), HyPar::FdFFunction() */
  int flag_fdf_specified;

  /*! Pointer to the function to calculate the upwind interface flux, given the left- and right-biased fluxes (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*Upwind)             (double*,double*,double*,double*,double*,double*,int,void*,double);
  /*! Pointer to the function to calculate the upwind interface split flux \f$df\left(u\right)\f$, given the left- and right-biased 
      approximations (assigned in the physical model initialization called from InitializePhysics()). Relevant only if the hyperbolic flux is being partitioned as 
      \f$f\left(u\right) = \left[f\left(u\right)-df\left(u\right)\right] + df\left(u\right)\f$
      \sa #HyPar::SplitHyperbolicFlux, HyPar::dFFunction()*/
  int    (*UpwinddF)           (double*,double*,double*,double*,double*,double*,int,void*,double);
  /*! Pointer to the function to calculate the upwind interface split flux \f$\left[f-df\right]\left(u\right)\f$, given its left- 
      and right-biased approximations (assigned in the physical model initialization called from InitializePhysics()). Relevant only if the hyperbolic flux is being partitioned as 
      \f$f\left(u\right) = \left[f\left(u\right)-df\left(u\right)\right] + df\left(u\right)\f$ (see
      #HyPar::SplitHyperbolicFlux, HyPar::dFFunction()), and if the split part
      \f$\left[f-df\right]\left(u\right)\f$ is being specified explicitly (see HyPar::FdFFunction(),
      #HyPar::flag_fdf_specified). */
  int    (*UpwindFdF)          (double*,double*,double*,double*,double*,double*,int,void*,double);

  /*! Pointer to the function to calculate the parabolic function with no cross-derivatives (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*GFunction)          (double*,double*,int,void*,double);

  /*! Pointer to the function to calculate the parabolic function with cross-derivatives (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*HFunction)          (double*,double*,int,int,void*,double);

  /*! Pointer to the function to calculate the source function (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*SFunction)          (double*,double*,void*,void*,double);

  /*! Pointer to the function to calculate the modified solution for upwinding (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*UFunction)          (double*,double*,int,void*,void*,double);

  /*! Pointer to the function to calculate the flux Jacobian for a given solution state (assigned in the physical model initialization called from InitializePhysics()). The
      flux Jacobian is the Jacobian of the analytical (*not* spatially discretized) flux for a given solution state (for example, at a grid point). The size is (#HyPar::nvars)^2 
      and the matrix is stored as  a 1D array in row-major format. */
  int    (*JFunction)          (double*,double*,void*,int,int);

  /*! Pointer to the function to do some pre-time-integration-stage computations, if required (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*PreStage)           (int,double**,void*,void*,double);
  /*! Pointer to the function to do some post-time-integration-stage computations, if required (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*PostStage)          (double*,void*,void*,double);
  /*! Pointer to the function to do some pre-time-integration-step computations, if required (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*PreStep)            (double*,void*,void*,double);
  /*! Pointer to the function to do some post-time-integration-step computations, if required (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*PostStep)           (double*,void*,void*,double,int);
  /*! Pointer to the function to print some physics-specific time-integration-step information,
   * if required (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*PrintStep)          (void*,void*,double);
  /*! Pointer to the function to write (to file) physics-related data that may not be 
      a part of the solution (assigned in the physical model initialization called from InitializePhysics()) */
  int    (*PhysicsOutput)      (void*,void*);

  /*! Pointer to the function to calculate the averaged solution state, given two solution states (assigned in the physical model initialization called from InitializePhysics()) */
  int   (*AveragingFunction)   (double*,double*,double*,void*);
  
  /*! Pointer to the function to calculate the left eigenvectors of the flux Jacobian for a given solution state (assigned in the physical model initialization called from InitializePhysics()) */
  int   (*GetLeftEigenvectors) (double*,double*,void*,int);
  /*! Pointer to the function to calculate the right eigenvectors of the flux Jacobian for a given solution state (assigned in the physical model initialization called from InitializePhysics()) */
  int   (*GetRightEigenvectors)(double*,double*,void*,int);

  /*! object containing interpolation-related parameters (of hyperbolic flux reconstruction */
  void *interp;
  /*! object containing arrays needed for compact finite-difference methods */
  void *compact;
  /*! object containing multi-stage time-integration (RK-type) related parameters */
  void *msti;
  /*! object containing parameters for the tridiagonal solver */
  void *lusolver;

  /*! Errors - L1, L2 and L_inf; calculated only if an exact solution file is supplied */
  double error[3];
  /*! Conservation error in the solution -- does not indicate anything if parabolic and source terms 
   * are non-zero */
  double *ConservationError;
  /*! check for conservation error? */
  char   ConservationCheck[_MAX_STRING_SIZE_];
  /*! Volume integral of the solution over the global domain */
  double *VolumeIntegral;
  /*! Volume integral of the initial solution over the global domain */
  double *VolumeIntegralInitial;
  /*! Surface integral of the flux at the boundary for each time-integration stage */
  double *StageBoundaryIntegral;
  /*! Surface integral of the flux at the boundary for a time-integration step */
  double *StepBoundaryIntegral;
  /*! Total surface integral of the flux over the global domain boundary */
  double *TotalBoundaryIntegral;
  /*! Pointer to the function to calculate the volume integral of a given function */
  int    (*VolumeIntegralFunction)    (double*,double*,void*,void*);
  /*! Pointer to the function to calculate the boundary integral of the flux */
  int    (*BoundaryIntegralFunction)  (void*,void*);
  /*! Pointer to the function to calculate the conservation error */
  int    (*CalculateConservationError)(void*,void*);

#ifdef with_petsc
  int     use_petscTS;  /*!< Use PETSc time-integration? */
  double  *u0;          /*!< copy of solution vector        */
  double  *uref;        /*!< copy of solution vector        */
  double  *rhsref;      /*!< copy of the RHS vector      */
  double  *rhs;         /*!< RHS vector                  */
#endif

  /*! flag to globally switch on/off non-linear interpolation */
  int flag_nonlinearinterp; 

  /*! strides along each dimension for an array with ghost points */
  int *stride_with_ghosts;
  /*! strides along each dimension for an array without ghost points */
  int *stride_without_ghosts;

  int count_hyp, /*!< number of times the hyperbolic function is called */
      count_par, /*!< number of times the parabolic function is called */
      count_sou; /*!< number of times the source function is called */
#ifdef with_petsc
  int count_RHSFunction,    /*!< number of times the RHSFunction is called */
      count_IFunction,      /*!< number of times the IFunction is called */
      count_IJacobian,      /*!< number of times the IJacobian is called */
      count_RHSJacobian,    /*!< number of times the RHSJacobian is called */
      count_IJacFunction,   /*!< number of times the IJacFunction is called */
      count_RHSJacFunction; /*!< number of times the RHSJacFunction is called */
#endif

  /*! blanking array: of same size and layout as #HyPar::u (but with 1 component
      per grid point, it has a value 1 for all valid grid points and 0 for grid 
      points that are blanked out. It is essentially an integer array, but the
      declared as a \a double type to use functions defined for the \a double 
      data type.
  */        
  double *iblank;          

  /*! Name of immersed body STL file (input - \b solver.inp )*/
  char ib_filename[_MAX_STRING_SIZE_];
  /*! Flag to indicate if immersed boundaries are in use */
  int flag_ib;
  /*! Immersed boundary object */
  void *ib;

  /*! Physics-specific immersed boundary treatment function (assigned in the physical model initialization called from InitializePhysics()) */
  int (*IBFunction) (void*,void*,double*,double);

} HyPar;

/* The following functions are called by main() */
int CalculateError                (void*,void*);/*!< Calculate the error in the final solution */
int Cleanup                       (void*,void*);/*!< Clean up: deallocate all arrays and objects */
int Initialize                    (void*,void*);/*!< Initialize the solver */
int InitializeBoundaries          (void*,void*);/*!< Initialize the boundary conditions */
int InitializeImmersedBoundaries  (void*,void*);/*!< Initialize the immersed boundary conditions */
int InitializePhysics             (void*,void*);/*!< Initialize the physics */
int InitializeSolvers             (void*,void*);/*!< Initialize the solvers */
int InitialSolution               (void*,void*);/*!< Read the initial solution */
int OutputSolution                (void*,void*);/*!< Write solution to file */
int ReadInputs                    (void*,void*);/*!< Read the input parameters */
int Solve                         (void*,void*);/*!< Solve the PDE - time-integration */
#ifdef with_petsc
int SolvePETSc                    (void*,void*);  /*!< Solve the PDE using PETSc TS */
#endif

/* Some definitions - types of discretizations available 
   for the parabolic (2nd derivative) term  */
/*! Non-conservative, direct evaluation of the 2nd deriv */
#define _NC_1STAGE_   "nonconservative-1stage"    
/*! Non-conservative, two-stage evaluation of the 2nd deriv */
#define _NC_2STAGE_   "nonconservative-2stage"    
/*! Non-conservative, "1.5"-stage evaluation of the 2nd deriv  - 
 * direct, 1-stage evaluation of the 2nd derivative in one variable
 * 2-stage evaluation of cross-derivatives */
#define _NC_1_5STAGE_ "nonconservative-1.5stage"  
/*! Conservative, direct evaluation of the 2nd deriv */
#define _CONS_1STAGE_ "conservative-1stage"       
