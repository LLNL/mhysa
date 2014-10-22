#include <stdlib.h>
#include <bandedmatrix.h>

int BandedMatrixPreallocate(void *A,int nbands,int nrows_local,int BlockSize)
{
  BandedMatrix *B = (BandedMatrix*) A;

  B->nbands       = nbands;
  B->nrows_local  = nrows_local;
  B->BlockSize    = BlockSize;

  B->ncol = (int*) calloc (nrows_local*nbands, sizeof(int));
  B->nrow = (int*) calloc (nrows_local, sizeof(int));
  B->data = (double*) calloc (nrows_local*nbands*BlockSize*BlockSize, sizeof(double));

  return(0);
}
