#include <stdlib.h>
#include <bandedmatrix.h>

int BandedMatrixDestroy(void *A)
{
  BandedMatrix *B = (BandedMatrix*) A;

  free(B->ncol);
  free(B->nrow);
  free(B->data);

  return(0);
}
