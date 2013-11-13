#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int InitialSolution(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i,d, ferr;
  int           offset_global, offset_local;
  _DECLARE_IERR_;

  /* Only root process reads in initial solution file */
  double *ug,*xg,*dxinvg; /* global solution vector and grid arrays */
  if (!mpi->rank) {
    int size,offset;
    /* allocate global solution array */
    size = solver->npoints_global*solver->nvars;
    ug = (double*) calloc(size,sizeof(double));
    size = 0;
    for (d=0; d<solver->ndims; d++) size += solver->dim_global[d];
    xg      = (double*) calloc(size,sizeof(double));
    dxinvg  = (double*) calloc(size,sizeof(double));

    /* Reading grid and initial solution */
    printf("Reading grid and initial conditions from file \"initial.inp\".\n");
    FILE *in; in = fopen("initial.inp","r");
    if (!in) {
      fprintf(stderr,"Error: initial solution file \"initial.inp\" not found.\n");
      return(1);
    }

    /* read grid and calculate dxinv*/
    offset = 0;
    for (d = 0; d < solver->ndims; d++) {
      for (i = 0; i < solver->dim_global[d]; i++) ferr = fscanf(in,"%lf",&xg[i+offset]);
      for (i = 0; i < solver->dim_global[d]; i++) {
        if      (i == 0)                        dxinvg[i+offset] = 1.0/(xg[i+1+offset]-xg[i  +offset]);
        else if (i == solver->dim_global[d]-1)  dxinvg[i+offset] = 1.0/(xg[i  +offset]-xg[i-1+offset]);
        else                                    dxinvg[i+offset] = 2.0/(xg[i+1+offset]-xg[i-1+offset]);
      }
      offset += solver->dim_global[d];
    }

    /* read solution */
    for (i = 0; i < solver->nvars; i++) {
      int *index = solver->index;
      int done = 0; _ArraySetValue_(index,solver->ndims,0);
      while (!done) {
        int p; _ArrayIndex1D_(solver->ndims,solver->dim_global,index,0,p);
        ferr = fscanf(in,"%lf",&ug[p*solver->nvars+i]);
        if (ferr != 1) return(1);
        _ArrayIncrementIndex_(solver->ndims,solver->dim_global,index,done);
      }
    }

    fclose(in);

  } else {
    ug      = NULL;
    xg      = NULL;
    dxinvg  = NULL;
  }

  /* partition initial solution across the processes */
  IERR MPIPartitionArraynD(solver->ndims,mpi,(mpi->rank?NULL:ug),solver->u,
                             solver->dim_global,solver->dim_local,
                             solver->ghosts,solver->nvars); CHECKERR(ierr);

  /* partition x vector across the processes */
  offset_global = offset_local = 0;
  for (d=0; d<solver->ndims; d++) {
    IERR MPIPartitionArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),
                                    &solver->x[offset_local+solver->ghosts],
                                    mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    IERR MPIPartitionArray1D(mpi,(mpi->rank?NULL:&dxinvg[offset_global]),
                                    &solver->dxinv[offset_local+solver->ghosts],
                                    mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    offset_global += solver->dim_global[d];
    offset_local  += solver->dim_local [d] + 2*solver->ghosts;
  }

  /* exchange MPI-boundary values of x between processors */
  offset_local = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->x[offset_local],solver->dim_local[d],
                                   solver->ghosts,d,solver->ndims); CHECKERR(ierr);
    IERR MPIExchangeBoundaries1D(mpi,&solver->dxinv[offset_local],solver->dim_local[d],
                                   solver->ghosts,d,solver->ndims); CHECKERR(ierr);
    offset_local  += solver->dim_local [d] + 2*solver->ghosts;
  }

  /* fill in ghost values of x and dxinv at physical boundaries by extrapolation */
  offset_local = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *x     = &solver->x    [offset_local];
    double *dxinv = &solver->dxinv[offset_local];
    int    ghosts = solver->ghosts;
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) {
        int delta = ghosts - i;
        dxinv[i] = dxinv[ghosts];
        x[i] = x[ghosts] + ((double) delta) * (x[ghosts]-x[ghosts+1]);
      }
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
        int delta = i - (dim[d]+ghosts-1);
        dxinv[i] = dxinv[dim[d]+ghosts-1];
        x[i] =  x[dim[d]+ghosts-1] 
              + ((double) delta) * (x[dim[d]+ghosts-1]-x[dim[d]+ghosts-2]);
      }
    }
    offset_local  += dim[d] + 2*ghosts;
  }
  

  if (!mpi->rank) {
    free(ug);
    free(xg);
    free(dxinvg);
  }
  return(0);
}
