/*! @file timeintegration.h
    @author Debojyoti Ghosh
    @brief Contains structures and function declarations for time integration
*/

#include <sys/time.h>
#include <basic.h>

/* definitions */
/*! Forward Euler time integration \sa TimeForwardEuler() */
#define _FORWARD_EULER_ "euler"
/*! Runge-Kutta time integration method \sa TimeRK() */
#define _RK_            "rk"
/*! General Linear Methods with Global Error Estimators \sa TimeGLMGEE() */
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
  /*! Final simulation time */
  double  t_final;
  /*! Restart iteration number (0 for a non-restart simulation) */
  int     restart_iter; 
  /*! Restart simulation time (0 for a non-restart simulation) */
  double  restart_time; 
  /*! Current solution time */
  double  waqt;         
  /*! Time step size */
  double  dt;           
  /*! CFL (for CFL-based time steps) */
  double  cfl;           
  /*! Norm of the change in the solution at a time step */
  double  norm;         
  /*! Maximum CFL at a time step */
  double  max_cfl;      

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

  /*! iteration start time */
  struct timeval iter_start_time;
  /*! iteration end time */
  struct timeval iter_end_time;
  /*! iteration wallclock time (in seconds) */
  double iter_wctime;

} TimeIntegration;

/* Explicit Runge-Kutta Methods */
/*!
  Forward Euler, written as a Runge Kutta method (1 stage, 1st order)

  \f$\begin{array}{c|c} 0 & 0 \\ \hline & 1 \end{array}\f$

  \sa TimeExplicitRKInitialize()
*/
#define _RK_1FE_        "1fe"     
/*!
  2-stage, 2nd order Runge Kutta method, often known as "RK2a"

  \f$\begin{array}{c|cc} 0 & & \\ 1 & 1 & \\ \hline & \frac{1}{2} & \frac{1}{2}\end{array}\f$

  \sa TimeExplicitRKInitialize()
*/
#define _RK_22_         "22"      
/*!
  3-stage, 3rd order Runge Kutta method

  \f$\begin{array}{c|ccc} 0 & & & \\ \frac{2}{3} & \frac{2}{3} & & \\ \frac{2}{3} & \frac{5}{12} & \frac{1}{4} & \\ \hline & \frac{1}{4} & -\frac{1}{4} & 1\end{array}\f$

  \sa TimeExplicitRKInitialize()
*/
#define _RK_33_         "33"
/*!
  The classical 4-stage, 4th order Runge Kutta method

  \f$ \begin{array}{c|cccc} 0 & & & & \\ \frac{1}{2} & \frac{1}{2} & & & \\ \frac{1}{2} & & \frac{1}{2} & & \\ 1 & & & 1 & \\ \hline & \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6} \end{array} \f$

  \sa TimeExplicitRKInitialize()
*/
#define _RK_44_         "44"      
/*!
  Strong-Stability-Preserving (SSP) 3-stage, 3rd order Runge Kutta method

  \f$ \begin{array}{c|ccc} 0 & & & \\ 1 & 1 & & \\ \frac{1}{2} & \frac{1}{4} & \frac{1}{4} & \\ \hline & \frac{1}{6} & \frac{1}{6} & \frac{2}{3} \end{array}\f$

  \sa TimeExplicitRKInitialize()

  Reference:
  + Gottlieb, S., Ketcheson, D. I., and Shu, C.-W., High Order Strong Stability 
    Preserving Time Discretizations, J. Sci. Comput., 38 (3), 2009, pp. 251-289,
    http://dx.doi.org/10.1007/s10915-008-9239-z.
*/
#define _RK_SSP3_       "ssprk3"  
#define _RK_TVD3_       "tvdrk3"  /*!< Same as #_RK_SSP3_                   */
/*! \def ExplicitRKParameters
    \brief Structure containing the parameters for an explicit Runge-Kutta method
*/
/*! \brief Structure containing the parameters for an explicit Runge-Kutta method

    Contains the parameters defining an explicit Runge Kutta time integration 
    method.

    \sa TimeRK()
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
/*! \f$y-\tilde{y}\f$ form 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_YYT_     "yyt" 
/*! \f$y-\epsilon\f$ form 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_YEPS_    "yeps" 
/*! A 3-stage, 2nd order method 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_23_      "23"  
/*! A 4-stage, 2nd order method 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_24_      "24"  
/*! A 5-stage, 2nd order method, with good imaginary stability 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_25I_     "25i" 
/*! A 5-stage, 3rd order method 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_35_      "35"  
/*! RK-2a with an error estimator 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_EXRK2A_  "exrk2a"  
/*! A 3rd order method 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_RK32G1_  "rk32g1"  
/*! A 2nd order RK method with an error estimator 

    \sa TimeGLMGEEInitialize()

    Reference:
    + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
#define _GLM_GEE_RK285EX_ "rk285ex" 
/*! \def GLMGEEParameters
    \brief Structure containing the parameters for the General Linear Methods with Global Error Estimators (GLM-GEE)
*/
/*! \brief Structure containing the parameters for the General Linear Methods with Global Error Estimators (GLM-GEE)

    Contains the parameters defining a General Linear time-integration method with global error estimators.

    \sa TimeGLMGEE()
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
int TimePrintStep       (void*,int);
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
