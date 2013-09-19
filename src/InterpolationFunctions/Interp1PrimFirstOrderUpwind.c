#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  First order upwind interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 1

int Interp1PrimFirstOrderUpwind(double *fI,double *fC,int upw,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0, d;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  if ((!fI) || (!fC)) {
    fprintf(stderr, "Error in FirstOrderUpwind(): input arrays not allocated.\n");
    return(1);
  }
  if (ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in FirstOrderUpwind(): insufficient number of ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int *indexC       = (int*) calloc (ndims,sizeof(int));
  int *indexI       = (int*) calloc (ndims,sizeof(int));
  int *index_outer  = (int*) calloc (ndims,sizeof(int));
  int *bounds_outer = (int*) calloc (ndims,sizeof(int));
  int *bounds_inter = (int*) calloc (ndims,sizeof(int));
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  int done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,indexC,ndims); CHECKERR(ierr);
    ierr = ArrayCopy1D_int(index_outer,indexI,ndims); CHECKERR(ierr);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      indexC[dir] = (upw > 0 ? indexI[dir]-1 : indexI[dir]);
      int p = ArrayIndex1D(ndims,bounds_inter,indexI,NULL,0     );
      int q = ArrayIndex1D(ndims,dim         ,indexC,NULL,ghosts);
      int v; for (v=0; v<nvars; v++)  fI[p*nvars+v] = fC[q*nvars+v];
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  free(indexC);
  free(indexI);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);
  
  return(0);
}
