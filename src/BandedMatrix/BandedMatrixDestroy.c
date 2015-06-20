/*! @file BandedMatrixDestroy.c
    @author Debojyoti Ghosh
    @brief Destroy a banded matrix object
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Free up allocations inside a banded matrix object */
int BandedMatrixDestroy(void *A /*!< Banded matrix object of type #BandedMatrix */)
{
  BandedMatrix *B = (BandedMatrix*) A;

  if (B->ncol) free(B->ncol);
  if (B->nrow) free(B->nrow);
  if (B->data) free(B->data);

  return(0);
}
