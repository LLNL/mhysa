/*! @file BCInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize boundary-conditions-related function pointers.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boundaryconditions.h>

/*! Assign the function pointers for boundary condition application depending on the boundary
    type, for a given boundary object */
int BCInitialize(void *b /*!< Boundary object of type #DomainBoundary*/)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  if      (!strcmp(boundary->bctype,_PERIODIC_                    )) boundary->BCFunctionU = BCPeriodicU;
  else if (!strcmp(boundary->bctype,_EXTRAPOLATE_                 )) boundary->BCFunctionU = BCExtrapolateU;
  else if (!strcmp(boundary->bctype,_DIRICHLET_                   )) boundary->BCFunctionU = BCDirichletU;  
  else if (!strcmp(boundary->bctype,_REFLECT_                     )) boundary->BCFunctionU = BCReflectU;    
  else if (!strcmp(boundary->bctype,_SPONGE_                      )) boundary->BCFunctionU = BCSpongeUDummy;    
  else if (!strcmp(boundary->bctype,_NOSLIP_WALL_                 )) boundary->BCFunctionU = BCNoslipWallU;    
  else if (!strcmp(boundary->bctype,_SLIP_WALL_                   )) boundary->BCFunctionU = BCSlipWallU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_OUTFLOW_            )) boundary->BCFunctionU = BCSubsonicOutflowU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_INFLOW_             )) boundary->BCFunctionU = BCSubsonicInflowU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_AMBIVALENT_         )) boundary->BCFunctionU = BCSubsonicAmbivalentU;    
  else if (!strcmp(boundary->bctype,_SUPERSONIC_OUTFLOW_          )) boundary->BCFunctionU = BCSupersonicOutflowU;    
  else if (!strcmp(boundary->bctype,_SUPERSONIC_INFLOW_           )) boundary->BCFunctionU = BCSupersonicInflowU;    
  else if (!strcmp(boundary->bctype,_TURBULENT_SUPERSONIC_INFLOW_ )) boundary->BCFunctionU = BCTurbulentSupersonicInflowU;    
  else {
    fprintf(stderr,"Error in BCInitialize(): \"%s\" is not a supported boundary condition.\n",
            boundary->bctype);
    return(1);
  }

  return(0);
}
