/* definitions */
#define _FORWARD_EULER_ 1

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

  void*   *ResidualFile;
  int (*TimeIntegrate)(void*);/* time integration */
} TimeIntegration;

/* functions */
int TimeInitialize (void*,void*,void*);
int TimeCleanup    (void*);
int TimePreStep    (void*);
int TimeStep       (void*);
int TimePostStep   (void*);
int TimePrintStep  (void*);

/* Time Integration Functions */
int TimeForwardEuler  (void*);
