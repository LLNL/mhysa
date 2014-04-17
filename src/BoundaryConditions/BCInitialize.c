#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boundaryconditions.h>

int BCInitialize(void *b)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

  if      (!strcmp(boundary->bctype,_PERIODIC_        )) boundary->BCFunctionU = BCPeriodicU;
  else if (!strcmp(boundary->bctype,_EXTRAPOLATE_     )) boundary->BCFunctionU = BCExtrapolateU;
  else if (!strcmp(boundary->bctype,_DIRICHLET_       )) boundary->BCFunctionU = BCDirichletU;  
  else if (!strcmp(boundary->bctype,_REFLECT_         )) boundary->BCFunctionU = BCReflectU;    
  else if (!strcmp(boundary->bctype,_NOSLIP_WALL_     )) boundary->BCFunctionU = BCNoslipWallU;    
  else if (!strcmp(boundary->bctype,_SLIP_WALL_       )) boundary->BCFunctionU = BCSlipWallU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_OUTFLOW_)) boundary->BCFunctionU = BCSubsonicOutflowU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_INFLOW_ )) boundary->BCFunctionU = BCSubsonicInflowU;    
  else {
    fprintf(stderr,"Error in BCInitialize(): \"%s\" is not a supported boundary condition.\n",
            boundary->bctype);
    return(1);
  }

  if      (!strcmp(boundary->bctype,_PERIODIC_        )) boundary->BCFunctionDU = BCPeriodicDU;
  else if (!strcmp(boundary->bctype,_EXTRAPOLATE_     )) boundary->BCFunctionDU = BCExtrapolateDU;
  else if (!strcmp(boundary->bctype,_DIRICHLET_       )) boundary->BCFunctionDU = BCDirichletDU;  
  else if (!strcmp(boundary->bctype,_REFLECT_         )) boundary->BCFunctionDU = BCReflectDU;    
  else if (!strcmp(boundary->bctype,_NOSLIP_WALL_     )) boundary->BCFunctionDU = BCNoslipWallDU;    
  else if (!strcmp(boundary->bctype,_SLIP_WALL_       )) boundary->BCFunctionDU = BCSlipWallDU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_OUTFLOW_)) boundary->BCFunctionDU = BCSubsonicOutflowDU;    
  else if (!strcmp(boundary->bctype,_SUBSONIC_INFLOW_ )) boundary->BCFunctionDU = BCSubsonicInflowDU;    
  else {
    fprintf(stderr,"Error in BCInitialize(): \"%s\" is not a supported boundary condition.\n",
            boundary->bctype);
    return(1);
  }

  return(0);
}
