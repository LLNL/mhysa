#include <stdlib.h>
#include <linearadr.h>

int LinearADRCleanup(void *s)
{
  LinearADR *physics = (LinearADR*) s;
  if (physics->a) free(physics->a);

  return(0);
}
