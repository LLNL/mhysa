#include <basic.h>

/* definitions */
#define _FORWARD_EULER_ "euler"
#define _RK_            "rk"
#define _GLM_GEE_       "glm-gee"

typedef struct time_integration_variables {
  int     iter;         /* iteration number                     */
  int     n_iter;       /* Total number of iterations           */
  int     restart_iter; /* Restart iteration number             */
  double  waqt;         /* time                                 */
  double  dt;           /* time step                            */
  double  norm;         /* norm of the solution                 */
  double  max_cfl;      /* max CFL for a time step              */
  double  max_diff;     /* max diffusion number for a time step */

  void    *solver;      /* solver object                        */
  void    *mpi;         /* mpi    object                        */
  double  *u;           /* array to store solution              */

  double  *rhs;         /* right-hand side array                */ 

  /* arrays for multi-stage schemes */
  double **U,**Udot;    /* stage values and RHS                 */
  double **BoundaryFlux;/* Boundary integral for each stage     */

  void *ResidualFile; /* file to write residual             */

  /* Functions */
  int (*TimeIntegrate) (void*);
  int (*RHSFunction)   (double*,double*,void*,void*,double);
} TimeIntegration;

/* Explicit Runge-Kutta Methods */
#define _RK_1FE_        "1fe"     /* Forward Euler                        */
#define _RK_22_         "22"      /* 2 stage, 2nd order                   */
#define _RK_33_         "33"      /* 3 stage, 3rd order                   */
#define _RK_44_         "44"      /* 4 stage, 4th order                   */
#define _RK_SSP3_       "ssprk3"  /* 3 stage, 3rd order SSP               */
#define _RK_TVD3_       "tvdrk3"  /* Same as ssprk3                       */
typedef struct _explicit_rungekutta_time_integration_ {
  int nstages;    /* number of stages */
  double *A,*b,*c;/* Butcher tableaux */
} ExplicitRKParameters;
int TimeExplicitRKInitialize(char*,char*,void*,void*);
int TimeExplicitRKCleanup   (void*);

/* General Linear Methods with Global Error Estimate */
#define _GLM_GEE_YYT_    "yyt"
#define _GLM_GEE_YEPS_   "yeps"
#define _GLM_GEE_23_     "23"
#define _GLM_GEE_24_     "24"
#define _GLM_GEE_35_     "35"
#define _GLM_GEE_EXRK2A_ "exrk2a"
#define _GLM_GEE_RK32G1_ "rk32g1"
typedef struct _glm_gee_time_integration_ {
  int nstages,    /* number of stages */
      r;
  char ee_mode[_MAX_STRING_SIZE_];
  double *A_yyt ,*B_yyt ,*C_yyt ,*D_yyt , *c_yyt;
  double *A_yeps,*B_yeps,*C_yeps,*D_yeps, *c_yeps;
  double *A, *B, *C, *D, *c;
  double gamma;
} GLMGEEParameters;
int TimeGLMGEEInitialize(char*,char*,void*,void*);
int TimeGLMGEECleanup   (void*);

/* functions */
int TimeInitialize      (void*,void*,void*);
int TimeCleanup         (void*);
int TimePreStep         (void*);
int TimeStep            (void*);
int TimePostStep        (void*);
int TimePrintStep       (void*);
int TimeError           (void*,void*,double*);
int TimeGetAuxSolutions (int*,double**,void*,int);

/* Time Integration Functions */
int TimeForwardEuler  (void*);
int TimeRK            (void*);
int TimeGLMGEE        (void*);
