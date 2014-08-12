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

int WENOFifthOrderCalculateWeights(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  int             upw, i;
  double          *ww1, *ww2, *ww3;

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

  /* calculate weights for a left-biased interpolation */
  upw = 1;

  ww1 = weno->w1 + offset;
  ww2 = weno->w2 + offset;
  ww3 = weno->w3 + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
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
        
      /* Defining stencil points */
      double m3[nvars], m2[nvars], m1[nvars], p1[nvars], p2[nvars];
      _ArrayCopy1D_((fC+qm3*nvars),m3,nvars);
      _ArrayCopy1D_((fC+qm2*nvars),m2,nvars);
      _ArrayCopy1D_((fC+qm1*nvars),m1,nvars);
      _ArrayCopy1D_((fC+qp1*nvars),p1,nvars);
      _ArrayCopy1D_((fC+qp2*nvars),p2,nvars);

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
      double w1[nvars],w2[nvars],w3[nvars];
      _WENOWeights_v_((ww1+p*nvars),(ww2+p*nvars),(ww3+p*nvars),c1,c2,c3,m3,m2,m1,p1,p2,weno,nvars);
    }
  }

  ww1 = weno->w1 + weno->size + offset;
  ww2 = weno->w2 + weno->size + offset;
  ww3 = weno->w3 + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
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

      /* Defining stencil points */
      double m3[nvars], m2[nvars], m1[nvars], p1[nvars], p2[nvars];
      _ArrayCopy1D_((uC+qm3*nvars),m3,nvars);
      _ArrayCopy1D_((uC+qm2*nvars),m2,nvars);
      _ArrayCopy1D_((uC+qm1*nvars),m1,nvars);
      _ArrayCopy1D_((uC+qp1*nvars),p1,nvars);
      _ArrayCopy1D_((uC+qp2*nvars),p2,nvars);

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
      double w1[nvars],w2[nvars],w3[nvars];
      _WENOWeights_v_((ww1+p*nvars),(ww2+p*nvars),(ww3+p*nvars),c1,c2,c3,m3,m2,m1,p1,p2,weno,nvars);
    }
  }

  /* calculate weights for a right-biased interpolation */
  upw = -1;

  ww1 = weno->w1 + 2*weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
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

      /* Defining stencil points */
      double m3[nvars], m2[nvars], m1[nvars], p1[nvars], p2[nvars];
      _ArrayCopy1D_((fC+qm3*nvars),m3,nvars);
      _ArrayCopy1D_((fC+qm2*nvars),m2,nvars);
      _ArrayCopy1D_((fC+qm1*nvars),m1,nvars);
      _ArrayCopy1D_((fC+qp1*nvars),p1,nvars);
      _ArrayCopy1D_((fC+qp2*nvars),p2,nvars);

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
      double w1[nvars],w2[nvars],w3[nvars];
      _WENOWeights_v_((ww1+p*nvars),(ww2+p*nvars),(ww3+p*nvars),c1,c2,c3,m3,m2,m1,p1,p2,weno,nvars);
    }
  }

  ww1 = weno->w1 + 2*weno->size + weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + weno->size + offset;
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
  for (i=0; i<N_outer; i++) {
    _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
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

      /* Defining stencil points */
      double m3[nvars], m2[nvars], m1[nvars], p1[nvars], p2[nvars];
      _ArrayCopy1D_((uC+qm3*nvars),m3,nvars);
      _ArrayCopy1D_((uC+qm2*nvars),m2,nvars);
      _ArrayCopy1D_((uC+qm1*nvars),m1,nvars);
      _ArrayCopy1D_((uC+qp1*nvars),p1,nvars);
      _ArrayCopy1D_((uC+qp2*nvars),p2,nvars);

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
      double w1[nvars],w2[nvars],w3[nvars];
      _WENOWeights_v_((ww1+p*nvars),(ww2+p*nvars),(ww3+p*nvars),c1,c2,c3,m3,m2,m1,p1,p2,weno,nvars);
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
