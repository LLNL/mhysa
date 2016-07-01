/*! @file PetscCreatePointList.c
    @brief Create list of computational points.
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <stdio.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

/*!
  Create a list of computational points: This is a list of all the 
  grid points on which the PDE is solved. Thus, it is the total
  number of grid points minus the ghost points and blanked out 
  points.

  Note: this list is local, not global.
*/
int PetscCreatePointList(void *obj /*!< Object of type #PETScContext */)
{
  PETScContext  *ctxt    = (PETScContext*) obj;
  HyPar         *solver  = (HyPar*) ctxt->solver;
  MPIVariables  *mpi     = (MPIVariables*) ctxt->mpi;


  int *dim  = solver->dim_local,
      ndims = solver->ndims,
      ghosts= solver->ghosts,
      done, index[ndims], npoints;

  /* count the number of computational points */
  npoints = 0; 
  done    = 0; 
  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    if (1) npoints++;
    //if (solver->iblank[p] == 1) npoints++;
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  ctxt->npoints = npoints;

  int nv = ndims+1;
  ctxt->points = (int*) calloc (ctxt->npoints*nv,sizeof(int));
  npoints = 0; 
  done    = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    //if (solver->iblank[p] == 1) {
    if (1) {
      _ArrayCopy1D_(index,(ctxt->points+npoints*nv),ndims);
      (ctxt->points+npoints*nv)[ndims] = p;
      npoints++;
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  if (npoints != ctxt->npoints) {
    fprintf(stderr,"Error in PetscCreatePointList() on rank %d: Inconsistency in point count - %d, %d.\n",
            mpi->rank, ctxt->npoints, npoints);
    return(1);
  }

  int global_npoints;
  MPISum_integer(&global_npoints,&npoints,1,&mpi->world);
  if (!mpi->rank) {
    printf("PETSc: total number of computational points is %d.\n",global_npoints);
  }

  return(0);
}

#endif
