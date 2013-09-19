#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int ParabolicFunctionCons1Stage(double *par,double *u,void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0, d, v, i, done;
  double        *FluxI  = NULL; /* interface flux     array */
  double        *Func   = NULL; /* diffusion function array */

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  double  *dxinv = solver->dxinv;

  int *index          = (int*) calloc (ndims,sizeof(int));
  int *index1         = (int*) calloc (ndims,sizeof(int));
  int *index2         = (int*) calloc (ndims,sizeof(int));
  int *dim_interface  = (int*) calloc (ndims,sizeof(int));

  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  ierr = ArraySetValue_double(par,size*nvars,0.0); CHECKERR(ierr);
  if (!solver->FFunction) return(0); /* zero parabolic term */

  int offset = 0;
  for (d = 0; d < ndims; d++) {

    /* allocate array for the diffusion function */
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    Func = (double*) calloc (size_cellcenter*nvars,sizeof(double));
    /* evaluate cell-centered flux */
    ierr = solver->GFunction(Func,u,d,solver); CHECKERR(ierr);

    /* calculate interface flux array dimensions */
    ierr = ArrayCopy1D_int(dim,dim_interface,ndims); CHECKERR(ierr); dim_interface[d]++;
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];
    /* allocate interface array for conservative discretization */
    FluxI = (double*) calloc (size_interface*nvars,sizeof(double));
    /* compute interface fluxes */
    ierr = solver->InterpolateInterfacesPar(FluxI,Func,d,solver,mpi); CHECKERR(ierr);

    /* calculate the second derivative */
    done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
    while (!done) {
      ierr = ArrayCopy1D_int(index,index1,ndims); CHECKERR(ierr);
      ierr = ArrayCopy1D_int(index,index2,ndims); CHECKERR(ierr); index2[d]++;
      int p  = ArrayIndex1D(ndims,dim          ,index ,NULL,ghosts);
      int p1 = ArrayIndex1D(ndims,dim_interface,index1,NULL,0     );
      int p2 = ArrayIndex1D(ndims,dim_interface,index2,NULL,0     );
      for (v=0; v<nvars; v++) 
        par[nvars*p+v] =  (dxinv[offset+index[d]] * dxinv[offset+index[d]]) 
                        * (FluxI[nvars*p2+v] - FluxI[nvars*p1+v]);
      done = ArrayIncrementIndex(ndims,dim,index);
    }

    /* free interface array */
    free(FluxI);
    free(Func);

    offset += dim[d];
  }

  free(index );
  free(index1);
  free(index2);

  return(0);
}
