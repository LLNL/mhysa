#include <physicalmodels/fpdoublewell.h>
#include <mpivars.h>
#include <hypar.h>

double FPDoubleWellComputeDiffNumber(void *s,void *m,double dt,double t)
{
  HyPar         *solver = (HyPar*)        s;
  FPDoubleWell  *params = (FPDoubleWell*) solver->physics;
  int           d, i, v;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;

  int     offset  = 0;
  double  max_diffno = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double local_diffno =   0.5 * params->q * dt 
                              * dxinv[offset+ghosts+i] 
                              * dxinv[offset+ghosts+i];
        if (local_diffno > max_diffno) max_diffno = local_diffno;
      }
    }
    offset += dim[d] + 2*ghosts;
  }

  return(max_diffno);
}
