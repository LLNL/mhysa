#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boundaryconditions.h>

int BCInitialize(void *b)
{
  DomainBoundary *boundary = (DomainBoundary*) b;

//  if      (!strcmp(boundary->bctype,_PERIODIC_   )) boundary->BCFunction = BCPeriodic;
//  else if (!strcmp(boundary->bctype,_EXTRAPOLATE_)) boundary->BCFunction = BCExtrapolate;
//  else {
//    fprintf(stderr,"Error in BCInitialize(): \"%s\" is not a supported boundary condition.\n",
//            boundary->bctype);
//    return(1);
//  }

  return(0);
}
