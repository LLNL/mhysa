/* definitions */
#define _FORWARD_EULER_ "euler"
#define _RK_            "rk"

#define _RK_1FE_        "1fe"     /* Forward Euler                        */
#define _RK_22_         "22"      /* 2 stage, 2nd order                   */
#define _RK_33_         "33"      /* 3 stage, 3rd order                   */
#define _RK_44_         "44"      /* 4 stage, 4th order                   */
#define _RK_SSP3_       "ssprk3"  /* 3 stage, 3rd order SSP               */
#define _RK_TVD3_       "tvdrk3"  /* Same as ssprk3                       */

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

typedef struct _explicit_rungekutta_time_integration_ {
  int nstages;    /* number of stages */
  double *A,*b,*c;/* Butcher tableaux */
} ExplicitRKParameters;

/* functions */
int TimeExplicitRKInitialize(char*,char*,void*);
int TimeExplicitRKCleanup   (void*);

int TimeInitialize    (void*,void*,void*);
int TimeCleanup       (void*);
int TimePreStep       (void*);
int TimeStep          (void*);
int TimePostStep      (void*);
int TimePrintStep     (void*);

/* Time Integration Functions */
int TimeForwardEuler  (void*);
int TimeRK            (void*);
