#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Fifth order WENO weights calculation for uniform grid
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

static int WENOFifthOrderCalculateWeightsJS(double*,double*,double*,int,void*,void*);
static int WENOFifthOrderCalculateWeightsM (double*,double*,double*,int,void*,void*);
static int WENOFifthOrderCalculateWeightsZ (double*,double*,double*,int,void*,void*);
static int WENOFifthOrderCalculateWeightsYC(double*,double*,double*,int,void*,void*);

int WENOFifthOrderCalculateWeights(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;

  if (weno->yc)           return(WENOFifthOrderCalculateWeightsYC (fC,uC,x,dir,solver,mpi));
  else if (weno->borges)  return(WENOFifthOrderCalculateWeightsZ  (fC,uC,x,dir,solver,mpi));
  else if (weno->mapped)  return(WENOFifthOrderCalculateWeightsM  (fC,uC,x,dir,solver,mpi));
  else                    return(WENOFifthOrderCalculateWeightsJS (fC,uC,x,dir,solver,mpi));
}

int WENOFifthOrderCalculateWeightsJS(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  int             upw, i;
  double          *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
  double          *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int *stride= solver->stride_with_ghosts;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  /* calculate dimension offset */
  int offset = weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  ww1LF = weno->w1 + offset;
  ww2LF = weno->w2 + offset;
  ww3LF = weno->w3 + offset;
  ww1RF = weno->w1 + 2*weno->size + offset;
  ww2RF = weno->w2 + 2*weno->size + offset;
  ww3RF = weno->w3 + 2*weno->size + offset;
  ww1LU = weno->w1 + weno->size + offset;
  ww2LU = weno->w2 + weno->size + offset;
  ww3LU = weno->w3 + weno->size + offset;
  ww1RU = weno->w1 + 2*weno->size + weno->size + offset;
  ww2RU = weno->w2 + 2*weno->size + weno->size + offset;
  ww3RU = weno->w3 + 2*weno->size + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1L,qm2L,qm3L,qp1L,qp2L,p,qm1R,qm2R,qm3R,qp1R,qp2R;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride[dir];
      qm2L = qm1L -   stride[dir];
      qp1L = qm1L +   stride[dir];
      qp2L = qm1L + 2*stride[dir];

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride[dir];
      qm2R = qm1R +   stride[dir];
      qp1R = qm1R -   stride[dir];
      qp2R = qm1R - 2*stride[dir];
        
      /* Defining stencil points */
      double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
      m3LF = (fC+qm3L*nvars);
      m2LF = (fC+qm2L*nvars);
      m1LF = (fC+qm1L*nvars);
      p1LF = (fC+qp1L*nvars);
      p2LF = (fC+qp2L*nvars);
      double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
      m3RF = (fC+qm3R*nvars);
      m2RF = (fC+qm2R*nvars);
      m1RF = (fC+qm1R*nvars);
      p1RF = (fC+qp1R*nvars);
      p2RF = (fC+qp2R*nvars);
      double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
      m3LU = (uC+qm3L*nvars);
      m2LU = (uC+qm2L*nvars);
      m1LU = (uC+qm1L*nvars);
      p1LU = (uC+qp1L*nvars);
      p2LU = (uC+qp2L*nvars);
      double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
      m3RU = (uC+qm3R*nvars);
      m2RU = (uC+qm2R*nvars);
      m1RU = (uC+qm1R*nvars);
      p1RU = (uC+qp1R*nvars);
      p2RU = (uC+qp2R*nvars);

      /* optimal weights*/
      double c1, c2, c3;
      if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        } else {
          /* CRWENO5 at the interior points */
          c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
          c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
          c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
        }
      } else {
        /* WENO5 and HCWENO5 */
        c1 = _WENO_OPTIMAL_WEIGHT_1_;
        c2 = _WENO_OPTIMAL_WEIGHT_2_;
        c3 = _WENO_OPTIMAL_WEIGHT_3_;
      }

      /* calculate WENO weights */
      upw =  1; _WENOWeights_v_JS_((ww1LF+p*nvars),(ww2LF+p*nvars),(ww3LF+p*nvars),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno,nvars);
      upw = -1; _WENOWeights_v_JS_((ww1RF+p*nvars),(ww2RF+p*nvars),(ww3RF+p*nvars),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno,nvars);
      upw =  1; _WENOWeights_v_JS_((ww1LU+p*nvars),(ww2LU+p*nvars),(ww3LU+p*nvars),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno,nvars);
      upw = -1; _WENOWeights_v_JS_((ww1RU+p*nvars),(ww2RU+p*nvars),(ww3RU+p*nvars),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno,nvars);
    }
  }

  return(0);
}

int WENOFifthOrderCalculateWeightsM(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  int             upw, i;
  double          *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
  double          *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int *stride= solver->stride_with_ghosts;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  /* calculate dimension offset */
  int offset = weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  ww1LF = weno->w1 + offset;
  ww2LF = weno->w2 + offset;
  ww3LF = weno->w3 + offset;
  ww1RF = weno->w1 + 2*weno->size + offset;
  ww2RF = weno->w2 + 2*weno->size + offset;
  ww3RF = weno->w3 + 2*weno->size + offset;
  ww1LU = weno->w1 + weno->size + offset;
  ww2LU = weno->w2 + weno->size + offset;
  ww3LU = weno->w3 + weno->size + offset;
  ww1RU = weno->w1 + 2*weno->size + weno->size + offset;
  ww2RU = weno->w2 + 2*weno->size + weno->size + offset;
  ww3RU = weno->w3 + 2*weno->size + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1L,qm2L,qm3L,qp1L,qp2L,p,qm1R,qm2R,qm3R,qp1R,qp2R;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride[dir];
      qm2L = qm1L -   stride[dir];
      qp1L = qm1L +   stride[dir];
      qp2L = qm1L + 2*stride[dir];

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride[dir];
      qm2R = qm1R +   stride[dir];
      qp1R = qm1R -   stride[dir];
      qp2R = qm1R - 2*stride[dir];
        
      /* Defining stencil points */
      double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
      m3LF = (fC+qm3L*nvars);
      m2LF = (fC+qm2L*nvars);
      m1LF = (fC+qm1L*nvars);
      p1LF = (fC+qp1L*nvars);
      p2LF = (fC+qp2L*nvars);
      double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
      m3RF = (fC+qm3R*nvars);
      m2RF = (fC+qm2R*nvars);
      m1RF = (fC+qm1R*nvars);
      p1RF = (fC+qp1R*nvars);
      p2RF = (fC+qp2R*nvars);
      double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
      m3LU = (uC+qm3L*nvars);
      m2LU = (uC+qm2L*nvars);
      m1LU = (uC+qm1L*nvars);
      p1LU = (uC+qp1L*nvars);
      p2LU = (uC+qp2L*nvars);
      double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
      m3RU = (uC+qm3R*nvars);
      m2RU = (uC+qm2R*nvars);
      m1RU = (uC+qm1R*nvars);
      p1RU = (uC+qp1R*nvars);
      p2RU = (uC+qp2R*nvars);

      /* optimal weights*/
      double c1, c2, c3;
      if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        } else {
          /* CRWENO5 at the interior points */
          c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
          c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
          c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
        }
      } else {
        /* WENO5 and HCWENO5 */
        c1 = _WENO_OPTIMAL_WEIGHT_1_;
        c2 = _WENO_OPTIMAL_WEIGHT_2_;
        c3 = _WENO_OPTIMAL_WEIGHT_3_;
      }

      /* calculate WENO weights */
      upw =  1; _WENOWeights_v_M_((ww1LF+p*nvars),(ww2LF+p*nvars),(ww3LF+p*nvars),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno,nvars);
      upw = -1; _WENOWeights_v_M_((ww1RF+p*nvars),(ww2RF+p*nvars),(ww3RF+p*nvars),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno,nvars);
      upw =  1; _WENOWeights_v_M_((ww1LU+p*nvars),(ww2LU+p*nvars),(ww3LU+p*nvars),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno,nvars);
      upw = -1; _WENOWeights_v_M_((ww1RU+p*nvars),(ww2RU+p*nvars),(ww3RU+p*nvars),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno,nvars);
    }
  }

  return(0);
}

int WENOFifthOrderCalculateWeightsZ(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  int             upw, i;
  double          *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
  double          *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int *stride= solver->stride_with_ghosts;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  /* calculate dimension offset */
  int offset = weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  ww1LF = weno->w1 + offset;
  ww2LF = weno->w2 + offset;
  ww3LF = weno->w3 + offset;
  ww1RF = weno->w1 + 2*weno->size + offset;
  ww2RF = weno->w2 + 2*weno->size + offset;
  ww3RF = weno->w3 + 2*weno->size + offset;
  ww1LU = weno->w1 + weno->size + offset;
  ww2LU = weno->w2 + weno->size + offset;
  ww3LU = weno->w3 + weno->size + offset;
  ww1RU = weno->w1 + 2*weno->size + weno->size + offset;
  ww2RU = weno->w2 + 2*weno->size + weno->size + offset;
  ww3RU = weno->w3 + 2*weno->size + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1L,qm2L,qm3L,qp1L,qp2L,p,qm1R,qm2R,qm3R,qp1R,qp2R;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride[dir];
      qm2L = qm1L -   stride[dir];
      qp1L = qm1L +   stride[dir];
      qp2L = qm1L + 2*stride[dir];

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride[dir];
      qm2R = qm1R +   stride[dir];
      qp1R = qm1R -   stride[dir];
      qp2R = qm1R - 2*stride[dir];
        
      /* Defining stencil points */
      double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
      m3LF = (fC+qm3L*nvars);
      m2LF = (fC+qm2L*nvars);
      m1LF = (fC+qm1L*nvars);
      p1LF = (fC+qp1L*nvars);
      p2LF = (fC+qp2L*nvars);
      double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
      m3RF = (fC+qm3R*nvars);
      m2RF = (fC+qm2R*nvars);
      m1RF = (fC+qm1R*nvars);
      p1RF = (fC+qp1R*nvars);
      p2RF = (fC+qp2R*nvars);
      double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
      m3LU = (uC+qm3L*nvars);
      m2LU = (uC+qm2L*nvars);
      m1LU = (uC+qm1L*nvars);
      p1LU = (uC+qp1L*nvars);
      p2LU = (uC+qp2L*nvars);
      double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
      m3RU = (uC+qm3R*nvars);
      m2RU = (uC+qm2R*nvars);
      m1RU = (uC+qm1R*nvars);
      p1RU = (uC+qp1R*nvars);
      p2RU = (uC+qp2R*nvars);

      /* optimal weights*/
      double c1, c2, c3;
      if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        } else {
          /* CRWENO5 at the interior points */
          c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
          c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
          c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
        }
      } else {
        /* WENO5 and HCWENO5 */
        c1 = _WENO_OPTIMAL_WEIGHT_1_;
        c2 = _WENO_OPTIMAL_WEIGHT_2_;
        c3 = _WENO_OPTIMAL_WEIGHT_3_;
      }

      /* calculate WENO weights */
      upw =  1; _WENOWeights_v_Z_((ww1LF+p*nvars),(ww2LF+p*nvars),(ww3LF+p*nvars),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno,nvars);
      upw = -1; _WENOWeights_v_Z_((ww1RF+p*nvars),(ww2RF+p*nvars),(ww3RF+p*nvars),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno,nvars);
      upw =  1; _WENOWeights_v_Z_((ww1LU+p*nvars),(ww2LU+p*nvars),(ww3LU+p*nvars),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno,nvars);
      upw = -1; _WENOWeights_v_Z_((ww1RU+p*nvars),(ww2RU+p*nvars),(ww3RU+p*nvars),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno,nvars);
    }
  }

  return(0);
}

int WENOFifthOrderCalculateWeightsYC(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  int             upw, i;
  double          *ww1LF, *ww2LF, *ww3LF, *ww1RF, *ww2RF, *ww3RF;
  double          *ww1LU, *ww2LU, *ww3LU, *ww1RU, *ww2RU, *ww3RU;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int *stride= solver->stride_with_ghosts;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  /* calculate dimension offset */
  int offset = weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  ww1LF = weno->w1 + offset;
  ww2LF = weno->w2 + offset;
  ww3LF = weno->w3 + offset;
  ww1RF = weno->w1 + 2*weno->size + offset;
  ww2RF = weno->w2 + 2*weno->size + offset;
  ww3RF = weno->w3 + 2*weno->size + offset;
  ww1LU = weno->w1 + weno->size + offset;
  ww2LU = weno->w2 + weno->size + offset;
  ww3LU = weno->w3 + weno->size + offset;
  ww1RU = weno->w1 + 2*weno->size + weno->size + offset;
  ww2RU = weno->w2 + 2*weno->size + weno->size + offset;
  ww3RU = weno->w3 + 2*weno->size + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1L,qm2L,qm3L,qp1L,qp2L,p,qm1R,qm2R,qm3R,qp1R,qp2R;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1L);
      qm3L = qm1L - 2*stride[dir];
      qm2L = qm1L -   stride[dir];
      qp1L = qm1L +   stride[dir];
      qp2L = qm1L + 2*stride[dir];

      indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1R);
      qm3R = qm1R + 2*stride[dir];
      qm2R = qm1R +   stride[dir];
      qp1R = qm1R -   stride[dir];
      qp2R = qm1R - 2*stride[dir];
        
      /* Defining stencil points */
      double *m3LF, *m2LF, *m1LF, *p1LF, *p2LF;
      m3LF = (fC+qm3L*nvars);
      m2LF = (fC+qm2L*nvars);
      m1LF = (fC+qm1L*nvars);
      p1LF = (fC+qp1L*nvars);
      p2LF = (fC+qp2L*nvars);
      double *m3RF, *m2RF, *m1RF, *p1RF, *p2RF;
      m3RF = (fC+qm3R*nvars);
      m2RF = (fC+qm2R*nvars);
      m1RF = (fC+qm1R*nvars);
      p1RF = (fC+qp1R*nvars);
      p2RF = (fC+qp2R*nvars);
      double *m3LU, *m2LU, *m1LU, *p1LU, *p2LU;
      m3LU = (uC+qm3L*nvars);
      m2LU = (uC+qm2L*nvars);
      m1LU = (uC+qm1L*nvars);
      p1LU = (uC+qp1L*nvars);
      p2LU = (uC+qp2L*nvars);
      double *m3RU, *m2RU, *m1RU, *p1RU, *p2RU;
      m3RU = (uC+qm3R*nvars);
      m2RU = (uC+qm2R*nvars);
      m1RU = (uC+qm1R*nvars);
      p1RU = (uC+qp1R*nvars);
      p2RU = (uC+qp2R*nvars);

      /* optimal weights*/
      double c1, c2, c3;
      if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
        if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
            || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
          /* Use WENO5 at the physical boundaries */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        } else {
          /* CRWENO5 at the interior points */
          c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
          c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
          c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
        }
      } else {
        /* WENO5 and HCWENO5 */
        c1 = _WENO_OPTIMAL_WEIGHT_1_;
        c2 = _WENO_OPTIMAL_WEIGHT_2_;
        c3 = _WENO_OPTIMAL_WEIGHT_3_;
      }

      /* calculate WENO weights */
      upw =  1; _WENOWeights_v_YC_((ww1LF+p*nvars),(ww2LF+p*nvars),(ww3LF+p*nvars),c1,c2,c3,m3LF,m2LF,m1LF,p1LF,p2LF,weno,nvars);
      upw = -1; _WENOWeights_v_YC_((ww1RF+p*nvars),(ww2RF+p*nvars),(ww3RF+p*nvars),c1,c2,c3,m3RF,m2RF,m1RF,p1RF,p2RF,weno,nvars);
      upw =  1; _WENOWeights_v_YC_((ww1LU+p*nvars),(ww2LU+p*nvars),(ww3LU+p*nvars),c1,c2,c3,m3LU,m2LU,m1LU,p1LU,p2LU,weno,nvars);
      upw = -1; _WENOWeights_v_YC_((ww1RU+p*nvars),(ww2RU+p*nvars),(ww3RU+p*nvars),c1,c2,c3,m3RU,m2RU,m1RU,p1RU,p2RU,weno,nvars);
    }
  }

  return(0);
}

int WENOFifthOrderCalculateWeightsChar(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*) s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             i, k, v, upw;
  double          *ww1, *ww2, *ww3;

  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  /* calculate dimension offset */
  int offset = weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double R[nvars*nvars], L[nvars*nvars], uavg[nvars], fchar[nvars];

  /* calculate weights for a left-biased interpolation */
  upw = 1;

  ww1 = weno->w1 + offset;
  ww2 = weno->w2 + offset;
  ww3 = weno->w3 + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      /* 1D indices of the stencil grid points */
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }

      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&uC[nvars*qm1],&uC[nvars*qp1],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double m3, m2, m1, p1, p2;
        m3 = m2 = m1 = p1 = p2 = 0;
        for (k = 0; k < nvars; k++) {
          m3 += L[v*nvars+k] * fC[qm3*nvars+k];
          m2 += L[v*nvars+k] * fC[qm2*nvars+k];
          m1 += L[v*nvars+k] * fC[qm1*nvars+k];
          p1 += L[v*nvars+k] * fC[qp1*nvars+k];
          p2 += L[v*nvars+k] * fC[qp2*nvars+k];
        }

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* calculate WENO weights */
        double w1,w2,w3;
        _WENOWeights_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno);

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }
    }
  }

  ww1 = weno->w1 + weno->size + offset;
  ww2 = weno->w2 + weno->size + offset;
  ww3 = weno->w3 + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      /* 1D indices of the stencil grid points */
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }

      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&uC[nvars*qm1],&uC[nvars*qp1],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double m3, m2, m1, p1, p2;
        m3 = m2 = m1 = p1 = p2 = 0;
        for (k = 0; k < nvars; k++) {
          m3 += L[v*nvars+k] * uC[qm3*nvars+k];
          m2 += L[v*nvars+k] * uC[qm2*nvars+k];
          m1 += L[v*nvars+k] * uC[qm1*nvars+k];
          p1 += L[v*nvars+k] * uC[qp1*nvars+k];
          p2 += L[v*nvars+k] * uC[qp2*nvars+k];
        }

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* calculate WENO weights */
        double w1,w2,w3;
        _WENOWeights_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno);

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }
    }
  }

  /* calculate weights for a right-biased interpolation */
  upw = -1;

  ww1 = weno->w1 + 2*weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      /* 1D indices of the stencil grid points */
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }

      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&uC[nvars*qm1],&uC[nvars*qp1],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double m3, m2, m1, p1, p2;
        m3 = m2 = m1 = p1 = p2 = 0;
        for (k = 0; k < nvars; k++) {
          m3 += L[v*nvars+k] * fC[qm3*nvars+k];
          m2 += L[v*nvars+k] * fC[qm2*nvars+k];
          m1 += L[v*nvars+k] * fC[qm1*nvars+k];
          p1 += L[v*nvars+k] * fC[qp1*nvars+k];
          p2 += L[v*nvars+k] * fC[qp2*nvars+k];
        }

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* calculate WENO weights */
        double w1,w2,w3;
        _WENOWeights_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno);

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }
    }
  }

  ww1 = weno->w1 + 2*weno->size + weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,k,v,R,L,uavg,fchar,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);

    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {

      /* 1D indices of the stencil grid points */
      int qm1,qm2,qm3,qp1,qp2;
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      } else {
        indexC[dir] = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      }

      int p; /* 1D index of the interface */
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);

      /* find averaged state at this interface */
      IERR solver->AveragingFunction(uavg,&uC[nvars*qm1],&uC[nvars*qp1],solver->physics); CHECKERR(ierr);

      /* Get the left and right eigenvectors */
      IERR solver->GetLeftEigenvectors  (uavg,L,solver->physics,dir); CHECKERR(ierr);
      IERR solver->GetRightEigenvectors (uavg,R,solver->physics,dir); CHECKERR(ierr);

      /* For each characteristic field */
      for (v = 0; v < nvars; v++) {

        /* calculate the characteristic flux components along this characteristic */
        double m3, m2, m1, p1, p2;
        m3 = m2 = m1 = p1 = p2 = 0;
        for (k = 0; k < nvars; k++) {
          m3 += L[v*nvars+k] * uC[qm3*nvars+k];
          m2 += L[v*nvars+k] * uC[qm2*nvars+k];
          m1 += L[v*nvars+k] * uC[qm1*nvars+k];
          p1 += L[v*nvars+k] * uC[qp1*nvars+k];
          p2 += L[v*nvars+k] * uC[qp2*nvars+k];
        }

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* calculate WENO weights */
        double w1,w2,w3;
        _WENOWeights_(w1,w2,w3,c1,c2,c3,m3,m2,m1,p1,p2,weno);

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }
    }
  }

  return(0);
}
