/* 
  Data structure and some function declarations
  for banded blocked matrices
*/

typedef struct banded_matrix {
  int nbands;                     /* number of block bands                  */
  int nrows_local;                /* number of block rows (local)           */
  int BlockSize;                  /* block size                             */
  int *ncol;                      /* global column numbers for each block   */
  int *nrow;                      /* global row numbers for each block      */
  double *data;
} BandedMatrix;

int BandedMatrixDestroy     (void*);
int BandedMatrixInitialize  (void*);
int BandedMatrixPreallocate (void*,int,int,int);
