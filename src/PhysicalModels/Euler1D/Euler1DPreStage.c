#include <basic.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

int Euler1DGravityField(void*,void*,double*);

int Euler1DPreStage(int stage,double **U,void *s,void *m,double waqt)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  IERR Euler1DGravityField(solver,mpi,U[stage]); CHECKERR(ierr);

  return(0);
}
