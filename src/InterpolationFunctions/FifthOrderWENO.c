#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Fifth order WENO interpolation (uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

int FifthOrderWENO(double *fI,double *fC,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  int             ierr    = 0, d;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  double one_sixth          = 1.0/6.0;
  double thirteen_by_twelve = 13.0/12.0;
  double one_fourth         = 1.0/4.0;

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
      int qm1,qm2,qm3,qp1,qp2,p;
      p = ArrayIndex1D(ndims,bounds_inter,indexI,NULL,0);
      if (upw > 0) {
        indexC[dir] = indexI[dir]-3; qm3 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-2; qm2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-1; qm1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]  ; qp1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]+1; qp2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      } else {
        indexC[dir] = indexI[dir]+2; qm3 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]+1; qm2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]  ; qm1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-1; qp1 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
        indexC[dir] = indexI[dir]-2; qp2 = ArrayIndex1D(ndims,dim,indexC,NULL,ghosts);
      }
      int v; 
      for (v=0; v<nvars; v++)  {
        /* Defining stencil points */
        double m3, m2, m1, p1, p2;
        m3 = fC[qm3*nvars+v];
        m2 = fC[qm2*nvars+v];
        m1 = fC[qm1*nvars+v];
        p1 = fC[qp1*nvars+v];
        p2 = fC[qp2*nvars+v];

        /* Candidate stencils and their optimal weights*/
        double f1, f2, f3, c1, c2, c3;
        f1 = (2*one_sixth)*m1 + (5*one_sixth)*p1 - (one_sixth)*p2;        c1 = 0.3;
        f2 = (-one_sixth)*m2 + (5.0*one_sixth)*m1 + (2*one_sixth)*p1;     c2 = 0.6;
        f3 = (2*one_sixth)*m3 - (7.0*one_sixth)*m2 + (11.0*one_sixth)*m1; c3 = 0.1;

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
        fI[p*nvars+v] = w1*f1 + w2*f2 + w3*f3;
      }
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
