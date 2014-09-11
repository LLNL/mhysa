#include <basic.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

int NavierStokes2DGravityField(void*,void*,double*);

int NavierStokes2DPreStage(int stage,double **U,void *s,void *m,double waqt)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  IERR NavierStokes2DGravityField(solver,mpi,U[stage]); CHECKERR(ierr);

  return(0);
}
