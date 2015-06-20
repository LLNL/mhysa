/*! @file BandedMatrixInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize a banded matrix object
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Initialize a newly-created banded matrix object. */
int BandedMatrixInitialize(void *A /*!< Banded matrix object of type BandedMatrix*/)
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
