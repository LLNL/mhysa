/*! @file timeintegration.h
    @author Debojyoti Ghosh
    @brief Contains structures and function declarations for time integration
*/

#include <basic.h>

/* definitions */
/*! Forward Euler time integration */
#define _FORWARD_EULER_ "euler"
/*! Runge-Kutta time integration method */
#define _RK_            "rk"
/*! General Linear Methods with Global Error Estimators */
#define _GLM_GEE_       "glm-gee"

/*! \def TimeIntegration
    \brief Structure of variables/parameters and function pointers for time integration
*/
/*! \brief Structure of variables/parameters and function pointers for time integration
 *
 * This structure contains all the variables, parameters, and function pointers 
 * required for integrating the spatially-discretized semi-discrete ordinary 
 * differential equation in time
*/
typedef struct time_integration_variables {
  /*! Current iteration number */
  int     iter;         
  /*! Total number of iterations */
  int     n_iter;       
  /*! Restart iteration number (0 for a non-restart simulation) */
  int     restart_iter; 
  /*! Current solution time */
  double  waqt;         
  /*! Time step size */
  double  dt;           
  /*! Norm of the change in the solution at a time step */
  double  norm;         
  /*! Maximum CFL at a time step */
  double  max_cfl;      
  /*! Maximum diffusion number at a time step */
  double  max_diff;     

  /*! Solver object of type #HyPar */
  void    *solver;      
  /*! MPI object of type #MPIVariables */
  void    *mpi;         
  /*! Array to store the current solution */
  double  *u;           

  /*! Array to store the right-hand side */ 
  double  *rhs;         

  /*! Arrays to store stage values for a multi-stage time-integration method */
  double  **U; 
  /*! Arrays to store stage right-hand-sides for a multi-stage time-integration method */
  double  **Udot;

  /*! Array to store the flux integral at the physical boundary at each stage of 
      a multi-stage time-integration method (to compute conservation errors) */
  double **BoundaryFlux;

  /*! Pointer to file to write residual history if required */
  void *ResidualFile; 

  /*! Pointer to the function that takes one time step using the desired method */
  int (*TimeIntegrate) (void*);
  /*! Pointer to the function that computes the right-hand-side */
  int (*RHSFunction)   (double*,double*,void*,void*,double);
} TimeIntegration;

/* Explicit Runge-Kutta Methods */
#define _RK_1FE_        "1fe"     /*!< Forward Euler                        */
#define _RK_22_         "22"      /*!< 2 stage, 2nd order                   */
#define _RK_33_         "33"      /*!< 3 stage, 3rd order                   */
#define _RK_44_         "44"      /*!< 4 stage, 4th order                   */
#define _RK_SSP3_       "ssprk3"  /*!< 3 stage, 3rd order SSP               */
#define _RK_TVD3_       "tvdrk3"  /*!< Same as ssprk3                       */
/*! \def ExplicitRKParameters
    \brief Structure containing the parameters for an explicit Runge-Kutta method
*/
/*! \brief Structure containing the parameters for an explicit Runge-Kutta method

    Contains the parameters defining an explicit Runge Kutta time integration 
    method, expressed as follows:
    \f{align}{
      {\bf U}^{\left(i\right)} &= {\bf u}_n + \delta t \sum_{j=1}^{i-1} A_{ij} \dot{\bf U}^{\left(j\right)}\ {\rm (Stage\ values)},\\
      {\bf u}_{n+1} &= {\bf u}_n + \delta t \sum_{i=1}^{n} b_i \dot{\bf U}^{\left(i\right)}\ {\rm (Step\ completion)},
    \f}
    where \f${\bf u}_n\f$ is the current solution, \f${\bf u}_{n+1}\f$ is the solution
    at the next time step, \f$A_{ij}\f$ and \f$b_i\f$ are the stage computation and 
    step completion coefficients in the Butcher tableau form.\n\n
    For explanation about the Butcher tableaux, see:
    + Butcher, J., Numerical Methods for Ordinary Differential Equations, Wiley, 2003.
*/
typedef struct _explicit_rungekutta_time_integration_ {
  int nstages; /*!< number of stages */
  double *A, /*!< Stage computation coefficients (Butcher tableau form),
                  saved as a 1D-array in row-major form */
         *b, /*!< Step completion coefficients (Butcher tableau form) */
         *c; /*!< Stage time coefficients (Butcher tableau form) */
} ExplicitRKParameters;
/*! Initialize the explicit Runge-Kutta time-integration method */
int TimeExplicitRKInitialize(char*,char*,void*,void*);
/*! Clean up variables related to the explicit Runge-Kutta time-integration method */
int TimeExplicitRKCleanup   (void*);

/* General Linear Methods with Global Error Estimate */
#define _GLM_GEE_YYT_     "yyt" /*!< \f$y-\tilde{y}\f$ form */
#define _GLM_GEE_YEPS_    "yeps" /*!< \f$y-\epsilon\f$ form */
#define _GLM_GEE_23_      "23"  /*!< A 3-stage, 2nd order method */
#define _GLM_GEE_24_      "24"  /*!< A 4-stage, 2nd order method */
#define _GLM_GEE_25I_     "25i" /*!< A 5-stage, 2nd order method, with good imaginary stability */
#define _GLM_GEE_35_      "35"  /*!< A 5-stage, 3rd order method */
#define _GLM_GEE_EXRK2A_  "exrk2a"  /*!< RK-2a with an error estimator */
#define _GLM_GEE_RK32G1_  "rk32g1"  /*!< A 3rd order method */
#define _GLM_GEE_RK285EX_ "rk285ex" /*!< A 2nd order RK method with an error estimator */
/*! \def GLMGEEParameters
    \brief Structure containing the parameters for the General Linear Methods with Global Error Estimators (GLM-GEE)
*/
/*! \brief Structure containing the parameters for the General Linear Methods with Global Error Estimators (GLM-GEE)

    Contains the parameters defining a general linear time-integration method with global error estimators.
    \n\n
    See:
    + Constantinescu, E.M., "Estimating Global Errors in Time Stepping", SIAM Journal on Numerical Analysis,
      (In press).
*/
typedef struct _glm_gee_time_integration_ {
  int nstages, /*!< Number of stages */
      r;       /*!< Number of auxiliary solutions propagated with solution */
  char ee_mode[_MAX_STRING_SIZE_]; /*!< Error estimation mode (#_GLM_GEE_YYT_ or #_GLM_GEE_YEPS_) */
  double  *A_yyt,   /*!< Stage computation coefficients (\f$y-\tilde{y}\f$ error-estimation mode) (row-major) */
          *B_yyt ,  /*!< Step completion coefficients (\f$y-\tilde{y}\f$ error-estimation mode) (row-major) */
          *C_yyt ,  /*!< Stage computation coefficients (\f$y-\tilde{y}\f$ error-estimation mode) (row-major) */
          *D_yyt ,  /*!< Step completion coefficients (\f$y-\tilde{y}\f$ error-estimation mode) (row-major) */
          *c_yyt;   /*!< Stage time coefficients (\f$y-\tilde{y}\f$ error-estimation mode) (row-major) */
  double  *A_yeps,  /*!< Stage computation coefficients (\f$y-\epsilon\f$ error-estimation mode) (row-major) */
          *B_yeps,  /*!< Step completion coefficients (\f$y-\epsilon\f$ error-estimation mode) (row-major) */
          *C_yeps,  /*!< Stage computation coefficients (\f$y-\epsilon\f$ error-estimation mode) (row-major) */
          *D_yeps,  /*!< Step completion coefficients (\f$y-\epsilon\f$ error-estimation mode) (row-major) */
          *c_yeps;  /*!< Stage time coefficients (\f$y-\epsilon\f$ error-estimation mode) (row-major) */
  double *A, /*!< Pointer to the stage computation coefficients */
         *B, /*!< Pointer to the step completion coefficients */
         *C, /*!< Pointer to the stage computation coefficients */
         *D, /*!< Pointer to the step completion coefficients */
         *c; /*!< Pointer to stage time coefficients */
  double gamma; /*!< Gamma parameter */
} GLMGEEParameters;
/*! Initialize the General Linear Methods with Global Error Estimators (GLM-GEE) */
int TimeGLMGEEInitialize(char*,char*,void*,void*);
/*! Clean up variables relted to General Linear Methods with Global Error Estimators (GLM-GEE) */
int TimeGLMGEECleanup   (void*);

/*! Initialize the time integration */
int TimeInitialize      (void*,void*,void*);
/*! Clean up variables related to time integration */
int TimeCleanup         (void*);
/*! Function called at the beginning of a time step */
int TimePreStep         (void*);
/*! Take one step in time */
int TimeStep            (void*);
/*! Function called at the end of a time step */
int TimePostStep        (void*);
/*! Print time integration related information */
int TimePrintStep       (void*);
/*! Compute/estimate error in solution */
int TimeError           (void*,void*,double*);
/*! Function to get auxiliary solutions if available (for example, in GLM-GEE methods) */
int TimeGetAuxSolutions (int*,double**,void*,int);

/*! Take a step in time using the Forward Euler method */
int TimeForwardEuler  (void*);
/*! Take a step in time using the explicit Runge-Kutta method */
int TimeRK            (void*);
/*! Take a step in time using the General Linear Methods with Global Error Estimators */
int TimeGLMGEE        (void*);
