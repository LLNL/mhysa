/*! @file BandedMatrixPreallocate.c
    @author Debojyoti Ghosh
    @brief Preallocate memory for a banded matrix object.
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Preallocate memory for a banded matrix object. */
int BandedMatrixPreallocate(
                              void  *A,           /*!< Banded matrix object of the type BandedMatrix */
                              int   nbands,       /*!< Number of bands */
                              int   nrows_local,  /*!< Local number of rows */
                              int   BlockSize     /*!< Block size */
                           )
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
