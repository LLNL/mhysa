/*! @file bandedmatrix.h
    @brief Data structure and some function declarations for banded blocked matrices
    @author Debojyoti Ghosh
*/

/*! \def BandedMatix
    \brief Structure of variables defining a banded matrix
 * This structure contains all the variables for defining a banded matrix.
*/

/*! \brief Structure for definined a banded matrix
 *
 * This structure contains all the variables for defining a banded matrix.
*/
typedef struct banded_matrix {
  int nbands;                     /*! number of block bands                  */
  int nrows_local;                /*! number of block rows (local)           */
  int BlockSize;                  /*! block size                             */
  int *ncol;                      /*! global column numbers for each block   */
  int *nrow;                      /*! global row numbers for each block      */
  double *data;                   /*! array containing the matrix elements   */
} BandedMatrix;

int BandedMatrixDestroy     (void*);  /*!< Destroy a banded matrix object */
int BandedMatrixInitialize  (void*);  /*!< Initialize a banded matrix object */
int BandedMatrixPreallocate (void*,int,int,int);  /*!< Preallocate memory for a banded matrix */
