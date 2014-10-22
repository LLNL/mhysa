#include <stdlib.h>
#include <bandedmatrix.h>

int BandedMatrixInitialize(void *A)
{
  BandedMatrix *B = (BandedMatrix*) A;

  B->nbands       = 0;
  B->nrows_local  = 0;
  B->BlockSize    = 0;

  B->ncol = NULL;
  B->nrow = NULL;
  B->data = NULL;

  return(0);
}
