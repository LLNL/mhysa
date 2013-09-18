/* definitions */
#define _FORWARD_EULER_ "euler"
#define _RK_            "rk"

#define _RK_1FE_        "1fe"

typedef struct time_integration_variables {
  int     iter;     /* iteration number           */
  int     n_iter;   /* Total number of iterations */
  double  waqt;     /* time                       */
  double  dt;       /* time step                  */
  double  norm;     /* norm of the solution       */
  double  max_cfl;  /* max CFL for a time step    */

  void    *solver;  /* solver object              */
  void    *mpi;     /* mpi    object              */
  double  *u;       /* array to store solution    */

  double  *rhs;     /* right-hand side array      */ 

  /* arrays for multi-stage schemes */
  double **U,**Udot; /* stage values and RHS      */

  void *ResidualFile; /* file to write residual   */
  int (*TimeIntegrate)(void*);/* time integration */
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
