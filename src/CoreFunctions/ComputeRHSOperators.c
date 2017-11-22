/*! @file ComputeRHSOperators.c
    @author Debojyoti Ghosh
    @brief Compute the matrices representing the right-hand-side Jacobians.
*/
#ifdef compute_rhs_operators

#include <basic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Computes the matrices representing the Jacobians of the hyperbolic
    (#HyPar::HyperbolicFunction), parabolic (#HyPar::ParabolicFunction), 
    and the source (#HyPar::SourceFunction) terms. \n\n
    Each element is compute through finite-differences, by perturbing 
    one DOF of the solution at a time. The matrices are not stored in the 
    memory, instead the non-zero are written to file as soon as they are 
    computed. The filenames for the matrices are "Mat_*Function_nnnnn.dat", 
    where "nnnnn" is a time-dependent index. 
    + The format (ASCII text) of the file is as follows:\n
      \a n \n
      \a i \a j \a val \n
      \a i \a j \a val \n
      ...\n
      \a i \a j \a val \n
      where \a n is the size of the matrix, followed by all the non-zero 
      elements (each line contains the indices \a i,\a j and the value \a val 
      of one non-zero element).
    + This function is called after #HyPar::file_op_iter iterations (i.e.
      the same frequency at which solution files are written).
    + If a splitting for the hyperbolic flux is defined, then the Jacobians
      of the complete hyperbolic term, as well as the split terms are computed.
    + The eigenvalues of these matrices can be computed and plotted in MATLAB 
      using the scripts Examples/Matlab/ComputeEvals.m and Examples/Matlab/PlotEvals.m
      respectively.
    + To use this function, the code must be compiled with the flag 
      \b -Dcompute_rhs_operators, and run on a single processor. This 
      is very slow for larger domains, and should be used only for 
      the numerical analysis of test cases.
    + The current solution stored in #HyPar::u is used as the reference state
      at which the Jacobians are computed (this is relevant to know for 
      non-linear right-hand-sides functions).

   Note: Parabolic terms have not yet been included.
*/
int ComputeRHSOperators(
                          void *s,  /*!< Solver object of type #HyPar */
                          void *m,  /*!< MPI object of type MPIVariables */
                          double t  /*!< Current simulation time */
                       )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i, j, size, ndof;
  double        *f, *f0, *u, *u0, *rhs, *drhs;
  FILE          *fout;
  char          filename[_MAX_STRING_SIZE_];

  if (mpi->nproc > 1) return(0);

  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int ghosts = solver->ghosts;
  int *dim   = solver->dim_local;
  int index[ndims];

  double epsilon = 1e-6;
  double tolerance = 1e-15;

  /* allocate arrays */
  size = solver->npoints_local_wghosts * nvars;
  u    = (double*) calloc (size,sizeof(double));
  u0   = (double*) calloc (size,sizeof(double));
  rhs  = (double*) calloc (size,sizeof(double));
  drhs = (double*) calloc (size,sizeof(double));
  ndof = solver->npoints_local * nvars;
  f    = (double*) calloc (ndof,sizeof(double));
  f0   = (double*) calloc (ndof,sizeof(double));

  /* copy the current solution to u0 */
  _ArrayCopy1D_(solver->u,u0,size); 
  /* apply boundary conditions to the solution u0 */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u0,NULL,t);CHECKERR(ierr);
  IERR solver->ApplyIBConditions(solver,mpi,u0,t);CHECKERR(ierr);
  IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u0); CHECKERR(ierr);

  /* compute linearized matrix for the hyperbolic FFunction */
  if (solver->FFunction) {
    strcpy(filename,"Mat_FFunction_");
    strcat(filename,solver->filename_index);
    strcat(filename,".dat");
    printf("ComputeRHSOperators(): Computing linearized matrix operator for FFunction. ndof=%d.\n",ndof);
    printf("ComputeRHSOperators(): Writing to sparse matrix file %s.\n",filename);
    fout = fopen(filename,"w");
    fprintf(fout,"%d\n",ndof);
    /* compute the FFunction of u0 */
    _ArraySetValue_(f0,ndof,0.0);
    IERR solver->HyperbolicFunction(rhs,u0,solver,mpi,t,1,solver->FFunction,
                                    solver->Upwind); CHECKERR(ierr);
    _ArrayScale1D_(rhs,-1.0,size);
    /* transfer to an array f0 which as no ghost points */
    IERR ArrayCopynD(ndims,rhs,f0,dim,ghosts,0,index,nvars); CHECKERR(ierr);
    for (i=0; i<ndof; i++) {
      /* copy the solution u0 to u */
      _ArrayCopy1D_(u0,u,size);
      /* find the 1D index p in an array with ghosts points (u),
       * corresponding to index i which assumes no ghost points */
      int ii, p, v;
      ii = i / nvars;
      v  = i - ii*nvars;
      _ArrayIndexnD_(ndims,ii,dim,index,0); 
      _ArrayIndex1D_(ndims,dim,index,ghosts,p);
      /* add a perturbation */
      u[nvars*p+v] += epsilon;
      /* apply boundary conditions to the perturbed solution u */
      IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);CHECKERR(ierr);
      IERR solver->ApplyIBConditions(solver,mpi,u,t);CHECKERR(ierr);
      IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u); CHECKERR(ierr);
      /* compute the FFunction of u */
      IERR solver->HyperbolicFunction(rhs,u,solver,mpi,t,1,solver->FFunction,
                                      solver->Upwind); CHECKERR(ierr);
      _ArrayScale1D_(rhs,-1.0,size);
      /* transfer to an array f which as no ghost points */
      _ArraySetValue_(f,ndof,0.0);
      IERR ArrayCopynD(ndims,rhs,f,dim,ghosts,0,index,nvars); CHECKERR(ierr);
      /* subtract f0 from f */
      _ArrayAXPY_(f0,-1.0,f,ndof);
      /* f is now espilon*column of the matrix */
      for (j=0; j<ndof; j++) {
        double mat_elem = f[j] / epsilon;
        /* write to file if element is non-zero */
        if (absolute(mat_elem) > tolerance) fprintf(fout,"%5d %5d %+1.16e\n",j+1,i+1,mat_elem);
      }
    }
    fclose(fout);
  }
  /* compute linearized matrix for the split hyperbolic (FFunction-dFFunction) and dFFunction
     functions, if dFFunction is available */
  if (solver->dFFunction) {
    strcpy(filename,"Mat_FdFFunction_");
    strcat(filename,solver->filename_index);
    strcat(filename,".dat");
    printf("ComputeRHSOperators(): Computing linearized matrix operator for (FFunction-dFFunction). ndof=%d.\n",ndof);
    printf("ComputeRHSOperators(): Writing to sparse matrix file %s.\n",filename);
    fout = fopen(filename,"w");
    fprintf(fout,"%d\n",ndof);
    /* compute the FFunction of u0 */
    _ArraySetValue_(f0,ndof,0.0);
    if (solver->flag_fdf_specified) {
      IERR solver->HyperbolicFunction(rhs,u0,solver,mpi,t,1,solver->FdFFunction,
                                      solver->UpwindFdF); CHECKERR(ierr);
      _ArrayScale1D_(rhs,-1.0,size);
    } else {
      IERR solver->HyperbolicFunction(rhs,u0,solver,mpi,t,1,solver->FFunction,
                                      solver->Upwind); CHECKERR(ierr);
      IERR solver->HyperbolicFunction(drhs,u0,solver,mpi,t,0,solver->dFFunction,
                                      solver->UpwinddF); CHECKERR(ierr);
      _ArrayScale1D_(rhs,-1.0,size);
      _ArrayAXPY_(drhs,1.0,rhs,size);
    }
    /* transfer to an array f0 which as no ghost points */
    IERR ArrayCopynD(ndims,rhs,f0,dim,ghosts,0,index,nvars); CHECKERR(ierr);
    for (i=0; i<ndof; i++) {
      /* copy the solution u0 to u */
      _ArrayCopy1D_(u0,u,size);
      /* find the 1D index p in an array with ghosts points (u),
       * corresponding to index i which assumes no ghost points */
      int ii, p, v;
      ii = i / nvars;
      v  = i - ii*nvars;
      _ArrayIndexnD_(ndims,ii,dim,index,0); 
      _ArrayIndex1D_(ndims,dim,index,ghosts,p);
      /* add a perturbation */
      u[nvars*p+v] += epsilon;
      /* apply boundary conditions to the perturbed solution u */
      IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);CHECKERR(ierr);
      IERR solver->ApplyIBConditions(solver,mpi,u,t);CHECKERR(ierr);
      IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u); CHECKERR(ierr);
      /* compute the FFunction of u */
      if (solver->flag_fdf_specified) {
        IERR solver->HyperbolicFunction(rhs,u,solver,mpi,t,1,solver->FdFFunction,
                                        solver->UpwindFdF); CHECKERR(ierr);
        _ArrayScale1D_(rhs,-1.0,size);
      } else {
        IERR solver->HyperbolicFunction(rhs,u,solver,mpi,t,1,solver->FFunction,
                                        solver->Upwind); CHECKERR(ierr);
        IERR solver->HyperbolicFunction(drhs,u,solver,mpi,t,0,solver->dFFunction,
                                        solver->UpwinddF); CHECKERR(ierr);
        _ArrayScale1D_(rhs,-1.0,size);
        _ArrayAXPY_(drhs,1.0,rhs,size);
      }
      /* transfer to an array f which as no ghost points */
      _ArraySetValue_(f,ndof,0.0);
      IERR ArrayCopynD(ndims,rhs,f,dim,ghosts,0,index,nvars); CHECKERR(ierr);
      /* subtract f0 from f */
      _ArrayAXPY_(f0,-1.0,f,ndof);
      /* f is now espilon*column of the matrix */
      for (j=0; j<ndof; j++) {
        double mat_elem = f[j] / epsilon;
        /* write to file if element is non-zero */
        if (absolute(mat_elem) > tolerance) fprintf(fout,"%5d %5d %+1.16e\n",j+1,i+1,mat_elem);
      }
    }
    fclose(fout);

    strcpy(filename,"Mat_dFFunction_");
    strcat(filename,solver->filename_index);
    strcat(filename,".dat");
    printf("ComputeRHSOperators(): Computing linearized matrix operator for dFFunction. ndof=%d.\n",ndof);
    printf("ComputeRHSOperators(): Writing to sparse matrix file %s.\n",filename);
    fout = fopen(filename,"w");
    fprintf(fout,"%d\n",ndof);
    /* compute the FFunction of u0 */
    _ArraySetValue_(f0,ndof,0.0);
    IERR solver->HyperbolicFunction(rhs,u0,solver,mpi,t,0,solver->dFFunction,
                                    solver->UpwinddF); CHECKERR(ierr);
    _ArrayScale1D_(rhs,-1.0,size);
    /* transfer to an array f0 which as no ghost points */
    IERR ArrayCopynD(ndims,rhs,f0,dim,ghosts,0,index,nvars); CHECKERR(ierr);
    for (i=0; i<ndof; i++) {
      /* copy the solution u0 to u */
      _ArrayCopy1D_(u0,u,size);
      /* find the 1D index p in an array with ghosts points (u),
       * corresponding to index i which assumes no ghost points */
      int ii, p, v;
      ii = i / nvars;
      v  = i - ii*nvars;
      _ArrayIndexnD_(ndims,ii,dim,index,0); 
      _ArrayIndex1D_(ndims,dim,index,ghosts,p);
      /* add a perturbation */
      u[nvars*p+v] += epsilon;
      /* apply boundary conditions to the perturbed solution u */
      IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);CHECKERR(ierr);
      IERR solver->ApplyIBConditions(solver,mpi,u,t);CHECKERR(ierr);
      IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u); CHECKERR(ierr);
      /* compute the FFunction of u */
      IERR solver->HyperbolicFunction(rhs,u,solver,mpi,t,0,solver->dFFunction,
                                      solver->UpwinddF); CHECKERR(ierr);
      _ArrayScale1D_(rhs,-1.0,size);
      /* transfer to an array f which as no ghost points */
      _ArraySetValue_(f,ndof,0.0);
      IERR ArrayCopynD(ndims,rhs,f,dim,ghosts,0,index,nvars); CHECKERR(ierr);
      /* subtract f0 from f */
      _ArrayAXPY_(f0,-1.0,f,ndof);
      /* f is now espilon*column of the matrix */
      for (j=0; j<ndof; j++) {
        double mat_elem = f[j] / epsilon;
        /* write to file if element is non-zero */
        if (absolute(mat_elem) > tolerance) fprintf(fout,"%5d %5d %+1.16e\n",j+1,i+1,mat_elem);
      }
    }
    fclose(fout);
  }

  /* compute linearized matrix for the source SFunction */
  if (solver->SFunction) {
    strcpy(filename,"Mat_SFunction_");
    strcat(filename,solver->filename_index);
    strcat(filename,".dat");
    printf("ComputeRHSOperators(): Computing linearized matrix operator for SFunction. ndof=%d.\n",ndof);
    printf("ComputeRHSOperators(): Writing to sparse matrix file %s.\n",filename);
    fout = fopen(filename,"w");
    fprintf(fout,"%d\n",ndof);
    /* compute the FFunction of u0 */
    _ArraySetValue_(f0,ndof,0.0);
    IERR solver->SourceFunction(rhs,u0,solver,mpi,t); CHECKERR(ierr);
    /* transfer to an array f0 which as no ghost points */
    IERR ArrayCopynD(ndims,rhs,f0,dim,ghosts,0,index,nvars); CHECKERR(ierr);
    for (i=0; i<ndof; i++) {
      /* copy the solution u0 to u */
      _ArrayCopy1D_(u0,u,size);
      /* find the 1D index p in an array with ghosts points (u),
       * corresponding to index i which assumes no ghost points */
      int ii, p, v;
      ii = i / nvars;
      v  = i - ii*nvars;
      _ArrayIndexnD_(ndims,ii,dim,index,0); 
      _ArrayIndex1D_(ndims,dim,index,ghosts,p);
      /* add a perturbation */
      u[nvars*p+v] += epsilon;
      /* apply boundary conditions to the perturbed solution u */
      IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);CHECKERR(ierr);
      IERR solver->ApplyIBConditions(solver,mpi,u,t);CHECKERR(ierr);
      IERR MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u); CHECKERR(ierr);
      /* compute the SFunction of u */
      IERR solver->SourceFunction(rhs,u,solver,mpi,t); CHECKERR(ierr);
      /* transfer to an array f which as no ghost points */
      _ArraySetValue_(f,ndof,0.0);
      IERR ArrayCopynD(ndims,rhs,f,dim,ghosts,0,index,nvars); CHECKERR(ierr);
      /* subtract f0 from f */
      _ArrayAXPY_(f0,-1.0,f,ndof);
      /* f is now espilon*column of the matrix */
      for (j=0; j<ndof; j++) {
        double mat_elem = f[j] / epsilon;
        /* write to file if element is non-zero */
        if (absolute(mat_elem) > tolerance) fprintf(fout,"%5d %5d %+1.16e\n",j+1,i+1,mat_elem);
      }
    }
    fclose(fout);
  }

  /* clean up */
  free(u);
  free(u0);
  free(rhs);
  free(drhs);
  free(f);
  free(f0);
  return(0);
}

#endif
