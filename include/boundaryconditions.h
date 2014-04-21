#include <basic.h>

#define _PERIODIC_            "periodic"
#define _EXTRAPOLATE_         "extrapolate"
#define _DIRICHLET_           "dirichlet"
#define _REFLECT_             "reflect"

/* some BC types unique to the euler/navier-stokes systems */
#define _NOSLIP_WALL_         "noslip-wall"
#define _SLIP_WALL_           "slip-wall"
#define _SUBSONIC_INFLOW_     "subsonic-inflow"
#define _SUBSONIC_OUTFLOW_    "subsonic-outflow"
#define _SUPERSONIC_INFLOW_   "supersonic-inflow"
#define _SUPERSONIC_OUTFLOW_  "supersonic-outflow"
/* note: supersonic inflow/outflow can be enforced by dirichlet/extrapolate bcs */

typedef struct domain_boundaries {
  char    bctype [_MAX_STRING_SIZE_]; /* Type of boundary condition                           */
  int     dim;                        /* dimension along which this BC applies                */
  int     face;                       /* 1 -> left/min, -1 -> right/max                       */
  double  *xmin,*xmax;                /* extent of this boundary condition                    */

  int on_this_proc;   /* flag indicating if this BC is applicable on this process             */
  int *is, *ie;       /* Index range on which to apply this BC on this process                */

  /* the boundary condition function for the solution vector U */
  int (*BCFunctionU) (void*,void*,int,int,int*,int,double*);
  /* the boundary condition function for the vector \Delta U (needed for implicit time-integration */
  int (*BCFunctionDU)(void*,void*,int,int,int*,int,double*,double*);

  double *DirichletValue;   /* specified value for steady Dirichlet BC */

  /* variables specific to Navier-Stokes/Euler equations BCs */
  double gamma,                                   /* ratio of specific heats  */
         FlowDensity,*FlowVelocity,FlowPressure;  /* boundary flow conditions */

} DomainBoundary;

/* Functions */
int BCInitialize(void*);
int BCCleanup   (void*);

/* Boundary condition implementations for the solution vector U */
int BCPeriodicU           (void*,void*,int,int,int*,int,double*);    /* Periodic boundary conditions    */
int BCExtrapolateU        (void*,void*,int,int,int*,int,double*);    /* extrapolate boundary conditions */
int BCDirichletU          (void*,void*,int,int,int*,int,double*);    /* Dirichlet boundary conditions   */
int BCReflectU            (void*,void*,int,int,int*,int,double*);    /* Reflection boundary conditions  */
int BCNoslipWallU         (void*,void*,int,int,int*,int,double*);    /* No-slip wall (viscous) boundary conditions  */
int BCSlipWallU           (void*,void*,int,int,int*,int,double*);    /* Slip (inviscid) wall   boundary conditions  */
int BCSubsonicInflowU     (void*,void*,int,int,int*,int,double*);    /* Subsonic inflow        boundary conditions  */
int BCSubsonicOutflowU    (void*,void*,int,int,int*,int,double*);    /* Subsonic outflow       boundary conditions  */
int BCSupersonicInflowU   (void*,void*,int,int,int*,int,double*);    /* Supersonic inflow      boundary conditions  */
int BCSupersonicOutflowU  (void*,void*,int,int,int*,int,double*);    /* Supersonic outflow     boundary conditions  */

/* Boundary condition implementations for the (\Delta U) */
int BCPeriodicDU          (void*,void*,int,int,int*,int,double*,double*);    /* Periodic boundary conditions    */
int BCExtrapolateDU       (void*,void*,int,int,int*,int,double*,double*);    /* extrapolate boundary conditions */
int BCDirichletDU         (void*,void*,int,int,int*,int,double*,double*);    /* Dirichlet boundary conditions   */
int BCReflectDU           (void*,void*,int,int,int*,int,double*,double*);    /* Reflection boundary conditions  */
int BCNoslipWallDU        (void*,void*,int,int,int*,int,double*,double*);    /* No-slip wall (viscous) boundary conditions  */
int BCSlipWallDU          (void*,void*,int,int,int*,int,double*,double*);    /* Slip (inviscid) wall   boundary conditions  */
int BCSubsonicInflowDU    (void*,void*,int,int,int*,int,double*,double*);    /* Subsonic inflow        boundary conditions  */
int BCSubsonicOutflowDU   (void*,void*,int,int,int*,int,double*,double*);    /* Subsonic outflow       boundary conditions  */
int BCSupersonicInflowDU  (void*,void*,int,int,int*,int,double*,double*);    /* Supersonic inflow      boundary conditions  */
int BCSupersonicOutflowDU (void*,void*,int,int,int*,int,double*,double*);    /* Supersonic outflow     boundary conditions  */
