#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Fifth order WENO characteristic-based weight calculation for uniform grid 
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int WENOFifthOrderCalculateWeightsChar(double *fC,double *uC,double *x,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*) s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             k, v, upw, done;
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

  /* allocate arrays for the averaged state, eigenvectors and characteristic interpolated f */
  double R[nvars*nvars], L[nvars*nvars], uavg[nvars], fchar[nvars];

  /* calculate weights for a left-biased interpolation */
  upw = 1;

  ww1 = weno->w1 + offset;
  ww2 = weno->w2 + offset;
  ww3 = weno->w3 + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {

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
            c1 = 0.1;
            c2 = 0.6;
            c3 = 0.3;
          } else {
            /* CRWENO5 at the interior points */
            c1 = 0.2;
            c2 = 0.5;
            c3 = 0.3;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = 0.2;
          c2 = 0.5;
          c3 = 0.3;
        }

        /* calculate WENO weights */
        double w1,w2,w3;

        if (weno->no_limiting) {
          /* fifth order polynomial interpolation */
          w1 = c1;
          w2 = c2;
          w3 = c3;
        } else {
          /* Smoothness indicators */
          double b1, b2, b3;
          b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
               + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
          b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
               + one_fourth*(m2-p1)*(m2-p1);
          b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
               + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
          /* This parameter is needed for Borges' or Yamaleev-Carpenter
             implementation of non-linear weights  */
          double tau;
          if (weno->borges) {
            tau = absolute(b3 - b1);
          } else if (weno->yc) {
            tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
          } else {
            tau = 0;
          }
  
          /* Defining the non-linear weights */
          double a1, a2, a3;
          if (weno->borges || weno->yc) {
            a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));
            a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));
            a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));
          } else {
            a1 = c1 / raiseto(b1+weno->eps,weno->p);
            a2 = c2 / raiseto(b2+weno->eps,weno->p);
            a3 = c3 / raiseto(b3+weno->eps,weno->p);
          }

          /* Convexity */
          double a_sum_inv;
          a_sum_inv = 1.0 / (a1 + a2 + a3);
          w1 = a1 * a_sum_inv;
          w2 = a2 * a_sum_inv;
          w3 = a3 * a_sum_inv;
  
          /* Mapping the weights, if needed */
          if (weno->mapped) {
            a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
            a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
            a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
            a_sum_inv = 1.0 / (a1 + a2 + a3);
            w1 = a1 * a_sum_inv;
            w2 = a2 * a_sum_inv;
            w3 = a3 * a_sum_inv;
          }
  
        }

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  ww1 = weno->w1 + weno->size + offset;
  ww2 = weno->w2 + weno->size + offset;
  ww3 = weno->w3 + weno->size + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {

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
            c1 = 0.1;
            c2 = 0.6;
            c3 = 0.3;
          } else {
            /* CRWENO5 at the interior points */
            c1 = 0.2;
            c2 = 0.5;
            c3 = 0.3;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = 0.2;
          c2 = 0.5;
          c3 = 0.3;
        }

        /* calculate WENO weights */
        double w1,w2,w3;

        if (weno->no_limiting) {
          /* fifth order polynomial interpolation */
          w1 = c1;
          w2 = c2;
          w3 = c3;
        } else {
          /* Smoothness indicators */
          double b1, b2, b3;
          b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
               + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
          b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
               + one_fourth*(m2-p1)*(m2-p1);
          b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
               + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
          /* This parameter is needed for Borges' or Yamaleev-Carpenter
             implementation of non-linear weights  */
          double tau;
          if (weno->borges) {
            tau = absolute(b3 - b1);
          } else if (weno->yc) {
            tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
          } else {
            tau = 0;
          }
  
          /* Defining the non-linear weights */
          double a1, a2, a3;
          if (weno->borges || weno->yc) {
            a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));
            a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));
            a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));
          } else {
            a1 = c1 / raiseto(b1+weno->eps,weno->p);
            a2 = c2 / raiseto(b2+weno->eps,weno->p);
            a3 = c3 / raiseto(b3+weno->eps,weno->p);
          }

          /* Convexity */
          double a_sum_inv;
          a_sum_inv = 1.0 / (a1 + a2 + a3);
          w1 = a1 * a_sum_inv;
          w2 = a2 * a_sum_inv;
          w3 = a3 * a_sum_inv;
  
          /* Mapping the weights, if needed */
          if (weno->mapped) {
            a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
            a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
            a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
            a_sum_inv = 1.0 / (a1 + a2 + a3);
            w1 = a1 * a_sum_inv;
            w2 = a2 * a_sum_inv;
            w3 = a3 * a_sum_inv;
          }
  
        }

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  /* calculate weights for a left-biased interpolation */
  upw = -1;

  ww1 = weno->w1 + 2*weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {

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
            c1 = 0.1;
            c2 = 0.6;
            c3 = 0.3;
          } else {
            /* CRWENO5 at the interior points */
            c1 = 0.2;
            c2 = 0.5;
            c3 = 0.3;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = 0.2;
          c2 = 0.5;
          c3 = 0.3;
        }

        /* calculate WENO weights */
        double w1,w2,w3;

        if (weno->no_limiting) {
          /* fifth order polynomial interpolation */
          w1 = c1;
          w2 = c2;
          w3 = c3;
        } else {
          /* Smoothness indicators */
          double b1, b2, b3;
          b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
               + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
          b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
               + one_fourth*(m2-p1)*(m2-p1);
          b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
               + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
          /* This parameter is needed for Borges' or Yamaleev-Carpenter
             implementation of non-linear weights  */
          double tau;
          if (weno->borges) {
            tau = absolute(b3 - b1);
          } else if (weno->yc) {
            tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
          } else {
            tau = 0;
          }
  
          /* Defining the non-linear weights */
          double a1, a2, a3;
          if (weno->borges || weno->yc) {
            a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));
            a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));
            a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));
          } else {
            a1 = c1 / raiseto(b1+weno->eps,weno->p);
            a2 = c2 / raiseto(b2+weno->eps,weno->p);
            a3 = c3 / raiseto(b3+weno->eps,weno->p);
          }

          /* Convexity */
          double a_sum_inv;
          a_sum_inv = 1.0 / (a1 + a2 + a3);
          w1 = a1 * a_sum_inv;
          w2 = a2 * a_sum_inv;
          w3 = a3 * a_sum_inv;
  
          /* Mapping the weights, if needed */
          if (weno->mapped) {
            a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
            a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
            a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
            a_sum_inv = 1.0 / (a1 + a2 + a3);
            w1 = a1 * a_sum_inv;
            w2 = a2 * a_sum_inv;
            w3 = a3 * a_sum_inv;
          }
  
        }

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  ww1 = weno->w1 + 2*weno->size + weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + weno->size + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {

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
            c1 = 0.1;
            c2 = 0.6;
            c3 = 0.3;
          } else {
            /* CRWENO5 at the interior points */
            c1 = 0.2;
            c2 = 0.5;
            c3 = 0.3;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = 0.2;
          c2 = 0.5;
          c3 = 0.3;
        }

        /* calculate WENO weights */
        double w1,w2,w3;

        if (weno->no_limiting) {
          /* fifth order polynomial interpolation */
          w1 = c1;
          w2 = c2;
          w3 = c3;
        } else {
          /* Smoothness indicators */
          double b1, b2, b3;
          b1 = thirteen_by_twelve*(m1-2*p1+p2)*(m1-2*p1+p2)
               + one_fourth*(3*m1-4*p1+p2)*(3*m1-4*p1+p2);
          b2 = thirteen_by_twelve*(m2-2*m1+p1)*(m2-2*m1+p1)
               + one_fourth*(m2-p1)*(m2-p1);
          b3 = thirteen_by_twelve*(m3-2*m2+m1)*(m3-2*m2+m1) 
               + one_fourth*(m3-4*m2+3*m1)*(m3-4*m2+3*m1);
  
          /* This parameter is needed for Borges' or Yamaleev-Carpenter
             implementation of non-linear weights  */
          double tau;
          if (weno->borges) {
            tau = absolute(b3 - b1);
          } else if (weno->yc) {
            tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
          } else {
            tau = 0;
          }
  
          /* Defining the non-linear weights */
          double a1, a2, a3;
          if (weno->borges || weno->yc) {
            a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));
            a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));
            a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));
          } else {
            a1 = c1 / raiseto(b1+weno->eps,weno->p);
            a2 = c2 / raiseto(b2+weno->eps,weno->p);
            a3 = c3 / raiseto(b3+weno->eps,weno->p);
          }

          /* Convexity */
          double a_sum_inv;
          a_sum_inv = 1.0 / (a1 + a2 + a3);
          w1 = a1 * a_sum_inv;
          w2 = a2 * a_sum_inv;
          w3 = a3 * a_sum_inv;
  
          /* Mapping the weights, if needed */
          if (weno->mapped) {
            a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
            a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
            a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
            a_sum_inv = 1.0 / (a1 + a2 + a3);
            w1 = a1 * a_sum_inv;
            w2 = a2 * a_sum_inv;
            w3 = a3 * a_sum_inv;
          }
  
        }

        /* save the weights */
        *(ww1+p*nvars+v) = w1;
        *(ww2+p*nvars+v) = w2;
        *(ww3+p*nvars+v) = w3;

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
