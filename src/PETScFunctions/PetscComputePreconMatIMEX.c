/*! @file PetscComputePreconMatIMEX.c
    @brief Contains the function to assemble the preconditioning matrix
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscComputePreconMatIMEX"
/*!
  Function to compute and assemble the preconditioning matrix. It is an approximation to
  the actual Jacobian of the left-hand-side \a dy/dt - \a g(y) for the implicit-explicit 
  time integration of the ODE \a dy/dt - \a g(y) = \a f(y). This can then be used with
  a suitable preconditioner in PETSc. See 
  http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/index.html
  for more information on PETSc preconditioners.\n\n
  The approximation is constructed by spatially discretizing the flux Jacobian with the
  first order upwind scheme.
*/
int PetscComputePreconMatIMEX(
                              Mat Pmat,   /*!< Preconditioning matrix to construct */
                              Vec Y,      /*!< Solution vector */
                              void *ctxt  /*!< Application context */
                             )
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = context->solver;
  MPIVariables    *mpi     = context->mpi;
  PetscErrorCode  ierr;
  int             ndims   = solver->ndims,
                  nvars   = solver->nvars,
                  ghosts  = solver->ghosts,
                  *dim    = solver->dim_local,
                  *dim_g  = solver->dim_global,
                  index[ndims],indexL[ndims],indexR[ndims],
                  v,dir,done,rows[nvars],cols[nvars];
  double          *u   = solver->u, dxinv, values[nvars*nvars];

  PetscFunctionBegin;
  /* copy solution from PETSc vector */
  ierr = TransferVecFromPETSc(u,Y,context); CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,context->waqt); CHECKERR(ierr);
  ierr = MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u); CHECKERR(ierr);
  /* initialize preconditioning matrix to zero */
  ierr = MatZeroEntries(Pmat); CHKERRQ(ierr);

  /* loop through all grid points */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    /* compute local and global 1D index of this grid point */
    int p;  _ArrayIndex1D_(ndims,dim,index,ghosts,p); /* local - for accessing u */
    int pg; _ArrayIndex1DWO_(ndims,dim_g,index,mpi->is,0,pg); /* global - row number in Pmat */
    /* compute the contributions from the flux derivatives along each dimension */
    for (dir = 0; dir < ndims; dir++) {
      /* compute indices and global 1D indices for left and right neighbors */
      _ArrayCopy1D_(index,indexL,ndims); indexL[dir]--;
      _ArrayCopy1D_(index,indexR,ndims); indexR[dir]++;
      int pgL; _ArrayIndex1DWO_(ndims,dim_g,indexL,mpi->is,0,pgL);
      int pgR; _ArrayIndex1DWO_(ndims,dim_g,indexR,mpi->is,0,pgR);
      int pL;  _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      int pR;  _ArrayIndex1D_(ndims,dim,indexR,ghosts,pL);
      /* Retrieve 1/delta-x at this grid point */
      _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

      /* diagonal element */
      for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pg + v; }
      ierr = solver->JFunction(values,(u+nvars*p),solver->physics,dir,0);
      _ArrayScale1D_(values,dxinv,(nvars*nvars));
      for (v=0; v<nvars; v++) values[nvars*v+v] += context->shift;
      ierr = MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES); CHKERRQ(ierr);

      /* left neighbor */
      if ((indexL[dir]+mpi->is[dir]) >= 0) {
        for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgL + v; }
        ierr = solver->JFunction(values,(u+nvars*pL),solver->physics,dir,1);
        _ArrayScale1D_(values,-dxinv,(nvars*nvars));
        ierr = MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES); CHKERRQ(ierr);
      }
      
      /* right neighbor */
      if ((indexR[dir]+mpi->is[dir]) < dim_g[dir]) {
        for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgR + v; }
        ierr = solver->JFunction(values,(u+nvars*pR),solver->physics,dir,-1);
        _ArrayScale1D_(values,-dxinv,(nvars*nvars));
        ierr = MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES); CHKERRQ(ierr);
      }
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  ierr = MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (Pmat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif
