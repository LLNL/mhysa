/*! @file boundaryconditions.h
    @brief Containts the structures and definitions for boundary condition implementation
    @author Debojyoti Ghosh
*/

#include <basic.h>

/*! Periodic boundary conditions \sa #BCPeriodicU */
#define _PERIODIC_                      "periodic"
/*! Extrapolative boundary conditions (values at ghost cells are copied from interior) \sa #BCExtrapolateU */
#define _EXTRAPOLATE_                   "extrapolate"
/*! Dirichlet boundary conditions (values at ghost cells specified through input) \sa #BCDirichletU */
#define _DIRICHLET_                     "dirichlet"
/*! Reflective boundary conditions (values at ghost cells negative of interior) \sa #BCReflectU */
#define _REFLECT_                       "reflect"
/*! Sponge boundary conditions \sa #BCSpongeSource, #BCSpongeUDummy */
#define _SPONGE_                        "sponge"

/* some BC types unique to the euler/navier-stokes systems */
/*! Viscous wall boundary condition (specific to Navier-Stokes) \sa #BCNoslipWallU */
#define _NOSLIP_WALL_                   "noslip-wall"
/*! Inviscid wall boundary condition (specific to Euler/Navier-Stokes) \sa BCSlipWallU */
#define _SLIP_WALL_                     "slip-wall"
/*! Subsonic inflow boundary condition: density and velocity are specified in the input, pressure is extrapolated from the interior (specific to Euler/Navier-Stokes) \sa #BCSubsonicInflowU*/
#define _SUBSONIC_INFLOW_               "subsonic-inflow"
/*! Subsonic outflow boundary condition: pressure is specified in the input, density and velocity are extrapolated from the interior (specific to Euler/Navier-Stokes) \sa #BCSubsonicOutflowU */
#define _SUBSONIC_OUTFLOW_              "subsonic-outflow"
/*! Subsonic "ambivalent" boundary condition: (specific to Euler/Navier-Stokes) the velocity at the boundary is extrapolated from the interior, and based on that, either subsonic outflow
    or inflow boundary conditions are applied. Specify all flow quantities (density, velocity, and pressure) in the input; depending on whether it is inflow or outflow, the appropriate
    quantities will be used. \sa #BCSubsonicAmbivalentU 
*/
#define _SUBSONIC_AMBIVALENT_              "subsonic-ambivalent"
/*! Supersonic inflow boundary condition: density, velocity, and pressure are specified in the input (specific to Euler/Navier-Stokes) \sa #BCSupersonicInflowU */
#define _SUPERSONIC_INFLOW_             "supersonic-inflow"
/*! Supersonic outflow boundary condition: all flow quantities are extrapolated from the interior (specific to Euler/Navier-Stokes) \sa BCSupersonicOutflowU */
#define _SUPERSONIC_OUTFLOW_            "supersonic-outflow"
/*! Turbulent, supersonic inflow boundary condition: density, velocity, and pressure are specified in the input, along with turbulent fluctuations (specific to Euler/Navier-Stokes) \sa #BCTurbulentSupersonicInflowU, #BCReadTurbulentInflowData */
#define _TURBULENT_SUPERSONIC_INFLOW_   "turbulent-supersonic-inflow"

/*! \def DomainBoundary
    \brief Structure containing the variables and function pointers defining a boundary
 * This structure contains all the variables and function pointers needed to specify
 * a boundary zone.
*/

/*!\brief Structure containing the variables and function pointers defining a boundary
 *
 * This structure contains all the variables and function pointers needed to specify
 * a boundary zone.
*/
typedef struct domain_boundaries {
  /*! Type of boundary (#_PERIODIC_, #_EXTRAPOLATE_, #_DIRICHLET_, etc) */
  char    bctype [_MAX_STRING_SIZE_]; 
  /*! Dimension along which this BC applies (For an \a n -dimensional problem, dimensions are indexed \a 0 to \a n-1 ) */
  int     dim;                        
  /*! Face on which this BC applies (1 -> left/min, or -1 -> right/max) */
  int     face;                       
  /*! Spatial extent of this boundary condition: \a xmin is an array of size \a n for a \a n -dimensional problem containing the starting spatial coordinates of the zone where this boundary applies */
  double  *xmin;                
  /*! Spatial extent of this boundary condition: \a xmax is an array of size \a n for a \a n -dimensional problem containing the ending spatial coordinates of the zone where this boundary applies */
  double  *xmax;                

  int on_this_proc;   /*!< Flag indicating if this BC is applicable on this process  (not an input) */
  int *is, *ie;       /*!< Index range on which to apply this BC on this process (not an input) */

  /*! Pointer to the specific boundary condition function for the solution vector U */
  int (*BCFunctionU) (void*,void*,int,int,int*,int,double*,double);

  double *DirichletValue;   /*!< Specified value for steady Dirichlet BC */
  double *SpongeValue;      /*!< Specified value for steady Sponge    BC */

  int    *UnsteadyDirichletSize; /*!< Size of array to hold unsteady Dirichlet data */
  double *UnsteadyDirichletData; /*!< Array to hold unsteady Dirichlet data         */
  /*! Filename to read in unsteady Dirichlet data from */
  char    UnsteadyDirichletFilename[_MAX_STRING_SIZE_]; 

  /* variables specific to Navier-Stokes/Euler equations BCs */
  double *FlowDensity,      /*!< Boundary flow densities of each species (specific to Euler/Navier-Stokes) */
         *FlowVelocity,     /*!< Boundary flow velocity (specific to Euler/Navier-Stokes) */
         FlowPressure;      /*!< Boundary flow pressure (specific to Euler/Navier-Stokes) */

  /* Physics object that allows physics-specific BC implementations */
  void **physics;

} DomainBoundary;

/* Functions */
int BCInitialize(void*); /*!< Function to initialize the boundary conditions */
int BCCleanup   (void*); /*!< Function to clean up boundary conditions-related variables and arrays */

/* Boundary condition implementations for the solution vector U */
/*! Periodic boundary conditions for the solution vector U */
int BCPeriodicU                     (void*,void*,int,int,int*,int,double*,double);    
/*! extrapolate boundary conditions for the solution vector U */
int BCExtrapolateU                  (void*,void*,int,int,int*,int,double*,double);    
/*! Dirichlet boundary conditions for the solution vector U */
int BCDirichletU                    (void*,void*,int,int,int*,int,double*,double);    
/*! Reflection boundary conditions for the solution vector U */
int BCReflectU                      (void*,void*,int,int,int*,int,double*,double);    
/*! No-slip wall (viscous) boundary conditions for the solution vector U */
int BCNoslipWallU                   (void*,void*,int,int,int*,int,double*,double);    
/*! Slip (inviscid) wall boundary conditions for the solution vector U */
int BCSlipWallU                     (void*,void*,int,int,int*,int,double*,double);    
/*! Subsonic inflow boundary conditions for the solution vector U */
int BCSubsonicInflowU               (void*,void*,int,int,int*,int,double*,double);    
/*! Subsonic outflow boundary conditions for the solution vector U */
int BCSubsonicOutflowU              (void*,void*,int,int,int*,int,double*,double);    
/*! Subsonic "ambivalent" boundary conditions for the solution vector U */
int BCSubsonicAmbivalentU           (void*,void*,int,int,int*,int,double*,double);    
/*! Supersonic inflow boundary conditions for the solution vector U */
int BCSupersonicInflowU             (void*,void*,int,int,int*,int,double*,double);    
/*! Supersonic outflow boundary conditions for the solution vector U */
int BCSupersonicOutflowU            (void*,void*,int,int,int*,int,double*,double);    
/*! Turbulent Supersonic inflow boundary conditions for the solution vector U */
int BCTurbulentSupersonicInflowU    (void*,void*,int,int,int*,int,double*,double);    

/*! a special BC enforcement - an absorbent sponge - enforced through a source term */
int BCSpongeSource        (void*,int,int,int,int*,double*,double*,double*);
/*! dummy function that get called during applying BCs - they don't do anything */
int BCSpongeUDummy        (void*,void*,int,int,int*,int,double*,double);

/*! Function to read in unsteady boundary data for turbulent inflow */
int BCReadTurbulentInflowData(void*,void*,int,int,int*);
