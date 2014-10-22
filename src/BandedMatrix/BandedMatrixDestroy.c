#include <stdlib.h>
#include <bandedmatrix.h>

int BandedMatrixDestroy(void *A)
{
  BandedMatrix *B = (BandedMatrix*) A;

  if (B->ncol) free(B->ncol);
  if (B->nrow) free(B->nrow);
  if (B->data) free(B->data);

  return(0);
}
