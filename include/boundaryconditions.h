#include <basic.h>

#define _PERIODIC_    "periodic"
#define _EXTRAPOLATE_ "extrapolate"
#define _DIRICHLET_   "dirichlet"
#define _REFLECT_     "reflect"

/* some BC types unique to the euler/navier-stokes systems */
#define _NOSLIP_WALL_       "noslip-wall"
#define _SLIP_WALL_         "slip-wall"
#define _SUBSONIC_INFLOW_   "subsonic-inflow"
#define _SUBSONIC_OUTFLOW_  "subsonic-outflow"
/* note: supersonic inflow/outflow can be enforced by dirichlet/extrapolate bcs */

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

  /* variables specific to Navier-Stokes/Euler equations BCs */
  double gamma,                                   /* ratio of specific heats  */
         FlowDensity,*FlowVelocity,FlowPressure;  /* boundary flow conditions */

} DomainBoundary;

/* Functions */
int BCInitialize(void*);
int BCCleanup   (void*);

int BCPeriodic    (void*,void*,int,int,int*,int,double*);    /* Periodic boundary conditions    */
int BCExtrapolate (void*,void*,int,int,int*,int,double*);    /* extrapolate boundary conditions */
int BCDirichlet   (void*,void*,int,int,int*,int,double*);    /* Dirichlet boundary conditions   */
int BCReflect     (void*,void*,int,int,int*,int,double*);    /* Reflection boundary conditions  */

/* BCs specified for Euler/Navier-Stokes equations */
int BCNoslipWall      (void*,void*,int,int,int*,int,double*);    /* No-slip wall (viscous) boundary conditions  */
int BCSlipWall        (void*,void*,int,int,int*,int,double*);    /* Slip (inviscid) wall   boundary conditions  */
int BCSubsonicInflow  (void*,void*,int,int,int*,int,double*);    /* Subsonic inflow        boundary conditions  */
int BCSubsonicOutflow (void*,void*,int,int,int*,int,double*);    /* Subsonic outflow       boundary conditions  */
