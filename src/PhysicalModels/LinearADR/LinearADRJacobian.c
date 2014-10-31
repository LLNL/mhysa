#include <basic.h>
#include <arrayfunctions.h>
#include <bandedmatrix.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

int LinearADRJacobian(void *J,double *u,void *s,void *m,double t,double alpha)
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  LinearADR     *physics = (LinearADR*)     solver->physics;
  BandedMatrix  *Jac     = (BandedMatrix*)  J;

  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int ndims   = solver->ndims;
  int nvars   = solver->nvars, bs2 = nvars*nvars;
  int nbands  = Jac->nbands;
  int index[ndims], d;

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,0     ,p);

    /* global row and column number of this grid point */
    int qc_global; _ArrayIndex1DWO_(ndims,solver->dim_global,index,mpi->is,0,qc_global);
    Jac->nrow[p] = qc_global;

    /* set the diagonal element of the matrix */
    Jac->ncol[nbands*p+ndims] = qc_global;
    Jac->data[nbands*bs2*p+ndims] = alpha;

    /* for each dimension, calculate the corresponding off-diagonal elements */
    for (d=0; d<ndims; d++) {
      int neighborL[ndims], neighborR[ndims];
      _ArrayCopy1D_(index,neighborL,ndims); neighborL[d]--;
      _ArrayCopy1D_(index,neighborR,ndims); neighborR[d]++;
      int qL_global; _ArrayIndex1DWO_(ndims,solver->dim_global,neighborL,mpi->is,0,qL_global);
      int qR_global; _ArrayIndex1DWO_(ndims,solver->dim_global,neighborR,mpi->is,0,qR_global);
      double dxinvL, dxinvC, dxinvR, dxinv;
      _GetCoordinate_(d,neighborL[d],dim,ghosts,solver->dxinv,dxinvL);
      _GetCoordinate_(d,index[d]    ,dim,ghosts,solver->dxinv,dxinvC);
      _GetCoordinate_(d,neighborR[d],dim,ghosts,solver->dxinv,dxinvR);

      if ((index[d] == 0) && (mpi->ip[d] == 0)) {
        /* left boundary along this dimension - forward  difference */
        Jac->ncol[nbands*p+ndims-(d+1)] = qR_global+1;
        Jac->ncol[nbands*p+ndims+(d+1)] = qR_global;
        dxinv = 1.0 / (0.5 * (1.0/dxinvR + 1.0/dxinvC)) ;
        Jac->data[nbands*bs2*p+(ndims-(d+1))]  =   0.0;
        Jac->data[nbands*bs2*p+(ndims      )] +=   physics->a[d] * dxinv;
        Jac->data[nbands*bs2*p+(ndims+(d+1))]  = - physics->a[d] * dxinv;
      } else if ((index[d] == dim[d]-1) && (mpi->ip[d] == mpi->iproc[d]-1)) {
        /* right boundary along this dimension -  backward difference */
        Jac->ncol[nbands*p+ndims-(d+1)] = qL_global;
        Jac->ncol[nbands*p+ndims+(d+1)] = qL_global-1;
        dxinv = 1.0 / (0.5 * (1.0/dxinvL + 1.0/dxinvC)) ;
        Jac->data[nbands*bs2*p+(ndims-(d+1))]  = physics->a[d] * dxinv;
        Jac->data[nbands*bs2*p+(ndims      )] -= physics->a[d] * dxinv;
        Jac->data[nbands*bs2*p+(ndims+(d+1))]  = 0.0;
      } else {
        /* interior */
        Jac->ncol[nbands*p+ndims-(d+1)] = qL_global;
        Jac->ncol[nbands*p+ndims+(d+1)] = qR_global;
        if (physics->a[d] >= 0) {
          /* backward difference */
          dxinv = 1.0 / (0.5 * (1.0/dxinvL + 1.0/dxinvC)) ;
          Jac->data[nbands*bs2*p+(ndims-(d+1))]  = physics->a[d] * dxinv;
          Jac->data[nbands*bs2*p+(ndims      )] -= physics->a[d] * dxinv;
          Jac->data[nbands*bs2*p+(ndims+(d+1))]  = 0.0;
        } else {
          /* forward  difference */
          dxinv = 1.0 / (0.5 * (1.0/dxinvR + 1.0/dxinvC)) ;
          Jac->data[nbands*bs2*p+(ndims-(d+1))]  =   0.0;
          Jac->data[nbands*bs2*p+(ndims      )] +=   physics->a[d] * dxinv;
          Jac->data[nbands*bs2*p+(ndims+(d+1))]  = - physics->a[d] * dxinv;
        }
      }
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(0);
}
