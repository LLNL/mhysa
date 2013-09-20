#include <stdlib.h>
#include <linearadr.h>

int LinearADRCleanup(void *s)
{
  LinearADR *physics = (LinearADR*) s;
  if (physics->a) free(physics->a);
  if (physics->d) free(physics->d);

  return(0);
}
