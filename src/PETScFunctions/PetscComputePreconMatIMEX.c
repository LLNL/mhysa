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
  Compute and assemble the preconditioning matrix for the implicit-explicit (IMEX) time integration
  of the governing equations: The ODE, obtained after discretizing the governing PDE in space,
  is expressed as follows (for the purpose of IMEX time integration):
  \f{eqnarray}{
    \frac {d{\bf U}}{dt} &=& {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right), \\
    \Rightarrow \frac {d{\bf U}}{dt} - {\bf G}\left({\bf U}\right) &=& {\bf F}\left({\bf U}\right), 
  \f}
  where \f${\bf F}\f$ is non-stiff and integrated in time explicitly, and \f${\bf G}\f$
  is stiff and integrated in time implicitly, and \f${\bf U}\f$ represents the entire
  solution vector.

  The Jacobian of the implicit part is thus given by:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \frac {\partial {\bf G}} {\partial {\bf U}} \right]
  \f}
  where \f$\alpha\f$ is the shift coefficient (#PETScContext::shift) of the time integration method.
  
  Currently, this function is implemented for the case where the implicitly-treated \f${\bf G}\f$
  is the spatially discretized part or whole of a hyperbolic term, i.e.,
  \f{equation}{
    {\bf G}\left({\bf U}\right) = \mathcal{D}\left\{{\bf g}\left({\bf u}\right)\right\} \approx \nabla \cdot {\bf g}\left({\bf u}\right),
  \f}
  with the governing PDE as
  \f{equation}{
    \frac {\partial {\bf u}} {\partial t} + \cdots + \nabla \cdot {\bf g}\left({\bf u}\right) + \cdots = \cdots,
  \f}
  and \f$\mathcal{D}\f$ representing the spatial discretization method. Thus, the Jacobian can be written as
  follows:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \mathcal{D}\left\{\frac {\partial {\bf g}} {\partial {\bf u}}\right\} \right]
  \f}
  The preconditioning matrix is usually a close approximation of the actual Jacobian matrix, where the actual
  Jacobian may be too expensive to evaluate and assemble. In this function, the preconditioning matrix is 
  the following approximation of the actual Jacobian:
  \f{equation}{
    {\bf J}_p = \left[\alpha{\bf I} - \mathcal{D}^{\left(1\right)}\left\{\frac {\partial {\bf g}} {\partial {\bf u}}\right\} \right] \approx {\bf J},
  \f}
  where \f$\mathcal{D}^{\left(1\right)}\f$ represents a 1st order upwind discretization operator. The matrix \f${\bf J}_p\f$
  is provided to the preconditioner. Note that #HyPar::JFunction is defined by the specific physics being solved, and computes 
  \f$\partial {\bf g}/ \partial {\bf u}\f$ at a grid point.

  + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/index.html for more information on PETSc preconditioners.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
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
      int pR;  _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
      /* Retrieve 1/delta-x at this grid point */
      _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

      /* diagonal element */
      for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pg + v; }
      ierr = solver->JFunction(values,(u+nvars*p),solver->physics,dir,0);
      _ArrayScale1D_(values,dxinv,(nvars*nvars));
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
  
  ierr = MatShift(Pmat,context->shift); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif
