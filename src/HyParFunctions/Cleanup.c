#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <tridiagLU.h>
#include <boundaryconditions.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <hypar.h>

/* include header files for each physical model */
#include <physicalmodels/linearadr.h>
#include <physicalmodels/fpdoublewell.h>
#include <physicalmodels/fppowersystem.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <physicalmodels/euler1d.h>
#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes2d.h>
#include <physicalmodels/navierstokes3d.h>
#include <physicalmodels/numa2d.h>
#include <physicalmodels/numa3d.h>

int Cleanup(void *s,void *m)
{
  HyPar           *solver   = (HyPar*)          s;
  MPIVariables    *mpi      = (MPIVariables*)   m;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  int             i;
  _DECLARE_IERR_;

  if (!mpi->rank) printf("Deallocating arrays.\n");

  /* Clean up boundary zones */
  for (i = 0; i < solver->nBoundaryZones; i++) {
    IERR BCCleanup(&boundary[i]); CHECKERR(ierr);
  }
  free(solver->boundary);

  /* Clean up any allocations in physical model */
  if (!strcmp(solver->model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {
    IERR LinearADRCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_FP_DOUBLE_WELL_)) {
    IERR FPDoubleWellCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_FP_POWER_SYSTEM_)) {
    IERR FPPowerSystemCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_FP_POWER_SYSTEM_1BUS_)) {
    IERR FPPowerSystem1BusCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_FP_POWER_SYSTEM_3BUS_)) {
    IERR FPPowerSystem3BusCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_EULER_1D_)) {
    IERR Euler1DCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_EULER_2D_)) {
    IERR Euler2DCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_NAVIER_STOKES_2D_)) {
    IERR NavierStokes2DCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_NAVIER_STOKES_3D_)) {
    IERR NavierStokes3DCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_NUMA2D_)) {
    IERR Numa2DCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_NUMA3D_)) {
    IERR Numa3DCleanup(solver->physics); CHECKERR(ierr);
  }
  free(solver->physics);

  /* Clean up any allocations from time-integration */
  if (solver->msti) {
    IERR TimeMSTICleanup(solver->msti); CHECKERR(ierr);
    free(solver->msti);
  }

  /* Clean up any spatial reconstruction related allocations */
  if (   (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_WENO_  )) 
      || (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_))
      || (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_)) ) {
    IERR WENOCleanup(solver->interp); CHECKERR(ierr);
  }
  if (solver->interp)   free(solver->interp);
  if (solver->lusolver) free(solver->lusolver);

  /* Free the communicators created */
  IERR MPIFreeCommunicators(solver->ndims,mpi); CHECKERR(ierr);

  /* These variables are allocated in Initialize.c */
  free(solver->dim_global);
  free(solver->dim_local);
  free(solver->index);
  free(solver->u);
#ifdef with_petsc
  if (solver->uref)   free(solver->uref);
  if (solver->rhsref) free(solver->rhsref);
  if (solver->rhs)    free(solver->rhs);   
#endif
  free(solver->hyp);
  free(solver->par);
  free(solver->source);
  free(solver->fluxC);
  free(solver->uC);
  free(solver->Deriv1);
  free(solver->Deriv2);
  free(solver->fluxI);
  free(solver->uL);
  free(solver->uR);
  free(solver->fL);
  free(solver->fR);
  free(solver->x);
  free(solver->dxinv);
  free(mpi->iproc);
  free(mpi->ip);
  free(mpi->is);
  free(mpi->ie);
  free(mpi->bcperiodic);
  free(mpi->sendbuf);
  free(mpi->recvbuf);
  free(solver->VolumeIntegral);
  free(solver->VolumeIntegralInitial);
  free(solver->StageBoundaryIntegral);
  free(solver->StepBoundaryIntegral);
  free(solver->TotalBoundaryIntegral);
  free(solver->ConservationError);
  free(solver->stride_with_ghosts);
  free(solver->stride_without_ghosts);

  return(0);
}
