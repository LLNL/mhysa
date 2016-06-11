/*! @file PetscGlobalDOF.c
    @author Debojyoti Ghosh
    @brief Compute the global DOF index for all the grid points
*/

#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

static int ApplyPeriodicity(
                              int     dir,    /*!< Spatial dimension along which to apply periodicity */
                              int     ndims,  /*!< Number of spatial dimensions */
                              int     *size,  /*!< Integer array with the number of grid points in 
                                                   each spatial dimension */
                              int     ghosts, /*!< Number of ghost points */
                              double  *phi    /*!< The array on which to apply the boundary condition */
                           )
{
  int bounds[ndims], index1[ndims], index2[ndims], offset[ndims], 
      done, p1 = 0, p2 = 0;
  _ArrayCopy1D_(size,bounds,ndims); bounds[dir] = ghosts;

  done = 0; _ArraySetValue_(index1,ndims,0);
  while (!done) {
    _ArraySetValue_(offset,ndims,0); offset[dir] = -ghosts;
    _ArrayIndex1DWO_(ndims,size,index1,offset,ghosts,p1);
    _ArrayCopy1D_(index1,index2,ndims);
    index2[dir] = index1[dir] + size[dir]-ghosts;
    _ArrayIndex1D_(ndims,size,index2,ghosts,p2);

    phi[p1] = phi[p2];
    _ArrayIncrementIndex_(ndims,bounds,index1,done);
  }

  done = 0; _ArraySetValue_(index1,ndims,0);
  while (!done) {
    _ArraySetValue_(offset,ndims,0); offset[dir] = size[dir];
    _ArrayIndex1DWO_(ndims,size,index1,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,size,index1,ghosts,p2);

    phi[p1] = phi[p2];
    _ArrayIncrementIndex_(ndims,bounds,index1,done);
  }
  return(0);
}

/*! Compute the global DOF index for all the grid points: The "global DOF index"
    is the component number (or block component number for #HyPar::nvars > 1) of
    a grid point in the global solution vector. It is also the row number (or
    block row number) of the grid point in the global matrix representing, for 
    example, the Jacobian of the right-hand-side.

    #PETScContext::globalDOF is an integer array with the same layout as the solution
    array #HyPar::u (but with one component) containing the global DOF index for the
    corresponding grid points. It has the same number of ghost points as #HyPar::u.
    + This array is initialized to -1.
    + The global DOF indices are computed for all non-ghost grid points.
    + If any boundaries are periodic, periodic boundary conditions are applied to fill
      the appropriate ghost points.
    + Ghost points corresponding to internal (MPI) boundaries are filled using 
      MPIExchangeBoundariesnD().
    + Thus, ghost points corresponding to physical, non-periodic boundaries retain the 
      initial value of -1.
*/
int PetscGlobalDOF(void *c /*!< Object of type #PETScContext*/)
{
  PETScContext  *ctxt   = (PETScContext*) c;
  HyPar         *solver = (HyPar*) ctxt->solver;
  MPIVariables  *mpi    = (MPIVariables*) ctxt->mpi;
  _DECLARE_IERR_;

  int   *dim      = solver->dim_local,
        ndims     = solver->ndims,
        ghosts    = solver->ghosts,
        size_wg   = solver->npoints_local_wghosts,
        rank      = mpi->rank,
        nproc     = mpi->nproc,
        nv        = ndims + 1, i;

  /* if globalDOF already allocated, free it */
  if (ctxt->globalDOF) free(ctxt->globalDOF);
  ctxt->globalDOF = (double*) calloc(size_wg,sizeof(double));
  _ArraySetValue_(ctxt->globalDOF,size_wg,-1.0);

  int local_sizes[nproc];
  _ArraySetValue_(local_sizes,nproc,0);
  local_sizes[rank] = ctxt->npoints;
  MPIMax_integer(local_sizes,local_sizes,nproc,&mpi->world);

  int myOffset = 0;
  for (i=0; i<rank; i++) myOffset += local_sizes[i];

  for (i=0; i<ctxt->npoints; i++) {
    int p = (ctxt->points+i*nv)[ndims];
    ctxt->globalDOF[p] = (double) (i + myOffset);
  }

  for (i=0; i<ndims; i++) {
    if (solver->isPeriodic[i]) {
      IERR ApplyPeriodicity(i,ndims,dim,ghosts,ctxt->globalDOF); CHECKERR(ierr);
    }
  }
  IERR MPIExchangeBoundariesnD(ndims,1,dim,ghosts,mpi,ctxt->globalDOF);

  return(0);
}

#endif
