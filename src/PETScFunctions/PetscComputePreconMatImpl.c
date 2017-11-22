/*! @file PetscComputePreconMatImpl.c
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
#define __FUNCT__ "PetscComputePreconMatImpl"
/*!
  Compute and assemble the preconditioning matrix for the implicit time integration
  of the governing equations: The ODE, obtained after discretizing the governing PDE in space,
  is expressed as follows:
  \f{equation}{
    \frac {d{\bf U}}{dt} = {\bf F}\left({\bf U}\right) 
    \Rightarrow \frac {d{\bf U}}{dt} - {\bf F}\left({\bf U}\right) = 0, 
  \f}
  where \f${\bf F}\f$ is the spatially discretized right-hand-side, and \f${\bf U}\f$ 
  represents the entire solution vector.

  The Jacobian is thus given by:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \frac {\partial {\bf F}} {\partial {\bf U}} \right]
  \f}
  where \f$\alpha\f$ is the shift coefficient (#PETScContext::shift) of the time integration method.

  \f$\mathcal{D}\f$ represents the spatial discretization method. Thus, the Jacobian can be written as
  follows:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \mathcal{D}\left\{\frac {\partial {\bf f}} {\partial {\bf u}}\right\} \right]
  \f}
  The preconditioning matrix is usually a close approximation of the actual Jacobian matrix, where the actual
  Jacobian may be too expensive to evaluate and assemble. In this function, the preconditioning matrix is 
  the following approximation of the actual Jacobian:
  \f{equation}{
    {\bf J}_p = \left[\alpha{\bf I} - \mathcal{D}^{\left(1\right)}\left\{\frac {\partial {\bf f}} {\partial {\bf u}}\right\} \right] \approx {\bf J},
  \f}
  where \f$\mathcal{D}^{\left(1\right)}\f$ represents a 1st order upwind discretization operator. The matrix \f${\bf J}_p\f$
  is provided to the preconditioner. Note that #HyPar::JFunction is defined by the specific physics being solved, and computes 
  \f$\partial {\bf f}/ \partial {\bf u}\f$ at a grid point.

  \b Note: Currently, this function is implemented only for the case where the implicitly-treated \f${\bf F}\f$
  is the spatially discretized hyperbolic term, i.e.,
  \f{equation}{
    {\bf F}\left({\bf U}\right) = \mathcal{D}\left\{{\bf f}\left({\bf u}\right)\right\} \approx \nabla \cdot {\bf f}\left({\bf u}\right).
  \f}
  If the governing PDE has a parabolic and/or source term in addition to the hyperbolic term, i.e.,
  \f{equation}{
    \frac {\partial {\bf u}} {\partial t} + \nabla \cdot {\bf f}\left({\bf u}\right) + \cdots\ \left({\rm parabolic\ term}\right) = \cdots\ \left({\rm source\ term}\right),
  \f}
  this function computes the approximation of the Jacobian of only the hyperbolic term. 

  + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/index.html for more information on PETSc preconditioners.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
int PetscComputePreconMatImpl(
                              Mat Pmat,   /*!< Preconditioning matrix to construct */
                              Vec Y,      /*!< Solution vector */
                              void *ctxt  /*!< Application context */
                             )
{
  PETScContext    *context = (PETScContext*) ctxt;
  HyPar           *solver  = context->solver;
  MPIVariables    *mpi     = context->mpi;
  PetscErrorCode  ierr;
  int             ndims       = solver->ndims,
                  nvars       = solver->nvars,
                  npoints     = context->npoints,
                  ghosts      = solver->ghosts,
                  *dim        = solver->dim_local,
                  *isPeriodic = solver->isPeriodic,
                  *points     = context->points,
                  index[ndims],indexL[ndims],indexR[ndims],
                  v,n,dir,done,rows[nvars],cols[nvars];
  double          *u      = solver->u, 
                  *iblank = solver->iblank,
                  dxinv, values[nvars*nvars];

  PetscFunctionBegin;
  /* copy solution from PETSc vector */
  ierr = TransferVecFromPETSc(u,Y,context); CHECKERR(ierr);
  /* apply boundary conditions and exchange data over MPI interfaces */
  ierr = solver->ApplyBoundaryConditions(solver,mpi,u,NULL,context->waqt); CHECKERR(ierr);
  ierr = MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u); CHECKERR(ierr);
  /* initialize preconditioning matrix to zero */
  ierr = MatZeroEntries(Pmat); CHKERRQ(ierr);

  /* loop through all computational points */
  for (n = 0; n < npoints; n++) {
    int *this_point = points + n*(ndims+1);
    int p = this_point[ndims];
    int index[ndims]; _ArrayCopy1D_(this_point,index,ndims);

    double iblank = solver->iblank[p];

    /* compute the contributions from the flux derivatives along each dimension */
    for (dir = 0; dir < ndims; dir++) {

      /* compute indices and global 1D indices for left and right neighbors */
      _ArrayCopy1D_(index,indexL,ndims); indexL[dir]--;
      int pL;  _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);

      _ArrayCopy1D_(index,indexR,ndims); indexR[dir]++;
      int pR;  _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

      int pg, pgL, pgR;
      pg  = (int) context->globalDOF[p];
      pgL = (int) context->globalDOF[pL];
      pgR = (int) context->globalDOF[pR];

      /* Retrieve 1/delta-x at this grid point */
      _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

      /* diagonal element */
      for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pg + v; }
      ierr = solver->JFunction(values,(u+nvars*p),solver->physics,dir,0);
      _ArrayScale1D_(values,(dxinv*iblank),(nvars*nvars));
      ierr = MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES); CHKERRQ(ierr);

      /* left neighbor */
      if (pgL >= 0) {
        for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgL + v; }
        ierr = solver->JFunction(values,(u+nvars*pL),solver->physics,dir,1);
        _ArrayScale1D_(values,(-dxinv*iblank),(nvars*nvars));
        ierr = MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES); CHKERRQ(ierr);
      }
      
      /* right neighbor */
      if (pgR >= 0) {
        for (v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgR + v; }
        ierr = solver->JFunction(values,(u+nvars*pR),solver->physics,dir,-1);
        _ArrayScale1D_(values,(-dxinv*iblank),(nvars*nvars));
        ierr = MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (Pmat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  ierr = MatShift(Pmat,context->shift); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif
