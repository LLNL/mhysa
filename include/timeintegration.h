/* definitions */
#define _FORWARD_EULER_ "euler"
#define _RK_            "rk"

#define _RK_1FE_        "1fe"     /* Forward Euler                        */
#define _RK_22_         "22"      /* 2 stage, 2nd order                   */
#define _RK_33_         "33"      /* 3 stage, 3rd order                   */
#define _RK_44_         "44"      /* 4 stage, 4th order                   */
#define _RK_SSP3_       "ssprk3"  /* 3 stage, 3rd order SSP               */

typedef struct time_integration_variables {
  int     iter;     /* iteration number                     */
  int     n_iter;   /* Total number of iterations           */
  double  waqt;     /* time                                 */
  double  dt;       /* time step                            */
  double  norm;     /* norm of the solution                 */
  double  max_cfl;  /* max CFL for a time step              */
  double  max_diff; /* max diffusion number for a time step */

  void    *solver;  /* solver object                        */
  void    *mpi;     /* mpi    object                        */
  double  *u;       /* array to store solution              */

  double  *rhs;     /* right-hand side array                */ 

  /* arrays for multi-stage schemes */
  double **U,**Udot; /* stage values and RHS                */

  void *ResidualFile; /* file to write residual             */

  /* Functions */
  int (*TimeIntegrate) (void*);
  int (*RHSFunction)   (double*,double*,void*,void*);
} TimeIntegration;

typedef struct _multistage_time_integration_ {
  int nstages;    /* number of stages */
  double *A,*b;   /* Butcher tableaux */
} MSTIParameters;

/* functions */
int TimeMSTIInitialize(char*,char*,void*);
int TimeMSTICleanup   (void*);
int TimeInitialize    (void*,void*,void*);
int TimeCleanup       (void*);
int TimePreStep       (void*);
int TimeStep          (void*);
int TimePostStep      (void*);
int TimePrintStep     (void*);

/* Time Integration Functions */
int TimeForwardEuler  (void*);
int TimeRK            (void*);
