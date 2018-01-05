/*! @file Cleanup.c
    @author Debojyoti Ghosh
    @brief Clean up and free memory after simulation is complete.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <bandedmatrix.h>
#include <mpivars.h>
#include <tridiagLU.h>
#include <boundaryconditions.h>
#include <immersedboundaries.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <hypar.h>

/* include header files for each physical model */
#include <physicalmodels/euler1d.h>
#include <physicalmodels/navierstokes3d.h>

/*! Cleans up and frees the memory after the completion of the simulation. */
int Cleanup(
              void *s, /*!< Solver object of type #HyPar */
              void *m  /*!< MPI object of type #MPIVariables */
           )
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

  /* Clean up immersed boundaries */
  if (solver->flag_ib) {
    IERR IBCleanup(solver->ib);
    free(solver->ib);
  }

  /* Clean up any allocations in physical model */
  if (!strcmp(solver->model,_EULER_1D_)) {
    IERR Euler1DCleanup(solver->physics); CHECKERR(ierr);
  } else if (!strcmp(solver->model,_NAVIER_STOKES_3D_)) {
    IERR NavierStokes3DCleanup(solver->physics); CHECKERR(ierr);
  }
  free(solver->physics);

  /* Clean up any allocations from time-integration */
#ifdef with_petsc
  if (!solver->use_petscTS) {
    if (!strcmp(solver->time_scheme,_RK_)) {
      IERR TimeExplicitRKCleanup(solver->msti); CHECKERR(ierr);
      free(solver->msti);
    } else if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
      IERR TimeGLMGEECleanup(solver->msti); CHECKERR(ierr);
      free(solver->msti);
    }
  }
#else
  if (!strcmp(solver->time_scheme,_RK_)) {
    IERR TimeExplicitRKCleanup(solver->msti); CHECKERR(ierr);
    free(solver->msti);
  } else if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
    IERR TimeGLMGEECleanup(solver->msti); CHECKERR(ierr);
    free(solver->msti);
  }
#endif

  /* Clean up any spatial reconstruction related allocations */
  if (   (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_WENO_  )) 
      || (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_))
      || (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_)) ) {
    IERR WENOCleanup(solver->interp); CHECKERR(ierr);
  }
  if (solver->interp)   free(solver->interp);
  if (   (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_COMPACT_UPWIND_ )) 
      || (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_         ))
      || (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_         )) ) {
    IERR CompactSchemeCleanup(solver->compact); CHECKERR(ierr);
  }
  if (solver->compact)  free(solver->compact);
  if (solver->lusolver) free(solver->lusolver);

  /* Free the communicators created */
  IERR MPIFreeCommunicators(solver->ndims,mpi); CHECKERR(ierr);

  /* These variables are allocated in Initialize.c */
  free(solver->dim_global);
  free(solver->dim_local);
  free(solver->index);
  free(solver->u);
#ifdef with_petsc
  if (solver->u0)     free(solver->u0);
  if (solver->uref)   free(solver->uref);
  if (solver->rhsref) free(solver->rhsref);
  if (solver->rhs)    free(solver->rhs);   
#endif
  free(solver->iblank);
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
  free(solver->isPeriodic);
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

  if (solver->filename_index) free(solver->filename_index);

  return(0);
}
