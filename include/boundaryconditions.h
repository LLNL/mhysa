#define _MAX_STRING_SIZE_ 500

#define _PERIODIC_    "periodic"
#define _EXTRAPOLATE_ "extrapolate"
#define _DIRICHLET_   "dirichlet"

typedef struct domain_boundaries {
  char    bctype [_MAX_STRING_SIZE_]; /* Type of boundary condition                           */
  int     var;                        /* variable to apply this BC on                         */
  int     dim;                        /* dimension along which this BC applies                */
  int     face;                       /* 1 -> left/min, -1 -> right/max                       */
  double  *xmin,*xmax;                /* extent of this boundary condition                    */

  int on_this_proc;   /* flag indicating if this BC is applicable on this process             */
  int *is, *ie;       /* Index range on which to apply this BC on this process                */

  int (*BCFunction)(void*,void*,int,int,int*,int,double*); /* the boundary condition function */

  double *DirichletValue;   /* specified value for steady Dirichlet BC */
} DomainBoundary;

/* Functions */
int BCInitialize(void*);
int BCCleanUp   (void*);

int BCPeriodic    (void*,void*,int,int,int*,int,double*);    /* Periodic boundary conditions    */
int BCExtrapolate (void*,void*,int,int,int*,int,double*);    /* extrapolate boundary conditions */
int BCDirichlet   (void*,void*,int,int,int*,int,double*);    /* extrapolate boundary conditions */
