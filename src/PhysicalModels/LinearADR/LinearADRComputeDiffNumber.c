#include <linearadr.h>
#include <mpivars.h>
#include <hypar.h>

double LinearADRComputeDiffNumber(void *s,void *m,double dt)
{
  HyPar         *solver = (HyPar*)        s;
  LinearADR     *params = (LinearADR*)    solver->physics;
  int           d, i, v;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;

  int     offset  = 0;
  double  max_diffno = 0;
  for (d = 0; d < ndims; d++) {
    for (i = 0; i < dim[d]; i++) {
      for (v = 0; v < nvars; v++) {
        double local_diffno =   params->d[nvars*d+v] * dt 
                              * dxinv[offset+i] * dxinv[offset+i];
        if (local_diffno > max_diffno) max_diffno = local_diffno;
      }
    }
    offset += dim[d];
  }

  return(max_diffno);
}
