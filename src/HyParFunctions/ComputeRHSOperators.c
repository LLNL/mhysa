/*
  Function to compute and print (to file) matrix representations
  of the right-hand side hyperbolic, parabolic and source functions.

  Implemented only for single-processor simulations.
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

int ComputeRHSOperators(void *s,void *m,double t)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i, j, size, ndof;
  double        *f, *f0, *u, *u0, *rhs;
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
  ndof = solver->npoints_local * nvars;
  f    = (double*) calloc (ndof,sizeof(double));
  f0   = (double*) calloc (ndof,sizeof(double));

  /* copy the current solution to u0 */
  _ArrayCopy1D_(solver->u,u0,size); 
  /* apply boundary conditions to the solution u0 */
  IERR solver->ApplyBoundaryConditions(solver,mpi,u0,NULL,0,t);CHECKERR(ierr);
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
      IERR solver->ApplyBoundaryConditions(solver,mpi,u,NULL,0,t);CHECKERR(ierr);
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

  /* clean up */
  free(u);
  free(u0);
  free(rhs);
  free(f);
  free(f0);
  return(0);
}

#endif
