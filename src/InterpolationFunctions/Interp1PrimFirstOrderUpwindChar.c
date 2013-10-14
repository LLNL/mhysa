#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  First order upwind characteristic-based interpolation (uniform grid )
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 1

int Interp1PrimFirstOrderUpwindChar(double *fI,double *fC,double *u,int upw,int dir,void *s,void *m)
{
  HyPar         *solver = (HyPar*) s;
  int           ierr    = 0, k, v;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  /* checks */
  if ((!fI) || (!fC) || (!u)) {
    fprintf(stderr, "Error in Interp1PrimFirstOrderUpwindChar(): input arrays not allocated.\n");
    return(1);
  }
  if (ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in Interp1PrimFirstOrderUpwindChar(): insufficient number of ghosts.\n");
    return(1);
  }
  if (nvars == 1) {
    fprintf(stderr, "Error in Interp1PrimFirstOrderUpwindChar(): number of variables = 1.\n");
    fprintf(stderr, "Code shouldn't have reached this function.\n");
    return(1);
  }
  if (      (!solver->AveragingFunction) 
        ||  (!solver->GetLeftEigenvectors)
        ||  (!solver->GetRightEigenvectors) ) {
    fprintf(stderr, "Error in Interp1PrimFirstOrderUpwindChar(): One of the required functions undefined.\n");
    fprintf(stderr, "AveragingFunction(), GetLeftEigenvectors() or GetRightEigenvectors().\n");
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

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double **R, **L, *uavg, *fchar;
  R     = (double**) calloc (nvars,sizeof(double*));
  L     = (double**) calloc (nvars,sizeof(double*));
  for (k = 0; k < nvars; k++){
    R[k]    = (double*) calloc (nvars,sizeof(double));
    L[k]    = (double*) calloc (nvars,sizeof(double));
  }
  uavg  = (double*)  calloc (nvars,sizeof(double ));
  fchar = (double*)  calloc (nvars,sizeof(double ));

  int done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {

    ierr = ArrayCopy1D_int(index_outer,indexC,ndims); CHECKERR(ierr);
    ierr = ArrayCopy1D_int(index_outer,indexI,ndims); CHECKERR(ierr);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      int pL, pR; /* 1D index of the left and right cells */
      indexC[dir] = indexI[dir]-1; pL = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      indexC[dir] = indexI[dir]  ; pR = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      int p; /* 1D index of the interface */
      p = ArrayIndex1D(ndims,bounds_inter,indexI,NULL,0);

      /* find averaged state at this interface */
      ierr = solver->AveragingFunction(uavg,&u[nvars*pL],&u[nvars*pR],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      ierr = solver->GetLeftEigenvectors  (uavg,L,solver->physics); CHECKERR(ierr);
      ierr = solver->GetRightEigenvectors (uavg,R,solver->physics); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double fcL = 0, fcR = 0;
        for (k = 0; k < nvars; k++) {
          fcL += L[v][k] * fC[pL*nvars+k];
          fcR += L[v][k] * fC[pR*nvars+k];
        }

        /* first order upwind approximation of the characteristic flux */
        fchar[v] = (upw > 0 ? fcL : fcR);

      }

      /* calculate the interface u from the characteristic u */
      ierr = MatVecMult(nvars,&fI[nvars*p],R,fchar); CHECKERR(ierr);

    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  for (k = 0; k < nvars; k++) {
    free(R[k]);
    free(L[k]);
  }
  free(R);
  free(L);
  free(uavg);
  free(fchar);

  free(indexC);
  free(indexI);
  free(index_outer);
  free(bounds_outer);
  free(bounds_inter);
  
  return(0);
}
