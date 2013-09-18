#include <stdlib.h>
#include <string.h>
#include <timeintegration.h>

int TimeMSTICleanup(void *s)
{
  MSTIParameters *params = (MSTIParameters*) s;
  if (params->A) free(params->A);
  if (params->b) free(params->b);
  return(0);
}
