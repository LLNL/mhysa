#include <stdlib.h>
#include <string.h>
#include <timeintegration.h>

int TimeExplicitRKCleanup(void *s)
{
  ExplicitRKParameters *params = (ExplicitRKParameters*) s;
  if (params->A) free(params->A);
  if (params->b) free(params->b);
  if (params->c) free(params->c);
  return(0);
}
