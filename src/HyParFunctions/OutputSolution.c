#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/* Function declarations */
static void IncrementFilename (char*);

int OutputSolution(void *s, void *m)
{
  HyPar         *solver = (HyPar*)       s;
  MPIVariables  *mpi    = (MPIVariables*)m;
  int           d;
  _DECLARE_IERR_;

  /* if WriteOutput() is NULL, then return */
  if (!solver->WriteOutput) return(0);

  /* root process: allocate global output arrays */
  double *ug, *xg;
  if (!mpi->rank) {
    int size_global;

    size_global = 1;
    for (d=0; d<solver->ndims; d++) size_global *= solver->dim_global[d];
    ug = (double*) calloc (size_global*solver->nvars,sizeof(double));
    _ArraySetValue_(ug,size_global*solver->nvars,0.0);

    size_global = 0;
    for (d=0; d<solver->ndims; d++) size_global += solver->dim_global[d];
    xg = (double*) calloc (size_global,sizeof(double));
    _ArraySetValue_(xg,size_global,0.0); CHECKERR(ierr);

  } else {

    /* null pointers on non-root processes */
    ug = xg = NULL;

  }

  /* Assemble the local output arrays into the global output arrays */
  IERR MPIGatherArraynD(solver->ndims,mpi,ug,solver->u,solver->dim_global,
                          solver->dim_local,solver->ghosts,solver->nvars);  CHECKERR(ierr);
  int offset_global, offset_local;
  offset_global = offset_local = 0;
  for (d=0; d<solver->ndims; d++) {
    IERR MPIGatherArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),
                            &solver->x[offset_local+solver->ghosts],
                            mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    offset_global += solver->dim_global[d];
    offset_local  += solver->dim_local [d] + 2*solver->ghosts;
  }

  if (!mpi->rank) {
    /* write output file to disk */
    char filename[_MAX_STRING_SIZE_] = "";
    strcat(filename,solver->op_filename);
    IERR solver->WriteOutput(solver->ndims,solver->nvars,solver->dim_global,xg,ug,
                               filename,solver->index); CHECKERR(ierr);
    if (!strcmp(solver->op_overwrite,"no"))  IncrementFilename(solver->op_filename);

    /* Clean up output arrays */
    free(xg);
    free(ug);
  }

  return(0);
}

void IncrementFilename(char *f)
{
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
}
