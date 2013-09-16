#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <advectiondiffusionreaction.h>
#include <mpivars.h>
#include <hypar.h>

static int AdvectionFlux(double*,double*,int*,int,void*,int,int*,int*,int*,int,int);

int LinearADRAdvection(void *s,void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           ierr    = 0, d, v, i, done;
  double        *FluxI  = NULL; /* interface flux     */
  double        *FluxC  = NULL; /* cell centered flux */

  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int ghosts = solver->ghosts;
  int *dim   = solver->dim_local;

  int *index          = (int*) calloc (ndims,sizeof(int));
  int *index1         = (int*) calloc (ndims,sizeof(int));
  int *index2         = (int*) calloc (ndims,sizeof(int));
  int *dim_interface  = (int*) calloc (ndims,sizeof(int));

  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);

  ierr = ArraySetValue_double(solver->hyp,size*nvars,0.0); CHECKERR(ierr);

  for (d = 0; d < ndims; d++) {

    /* allocate array for cell-centered flux */
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    FluxC = (double*) calloc (size_cellcenter*nvars,sizeof(double));
    /* evaluate cell-centered flux */
    ierr = AdvectionFlux(solver->u,FluxC,dim,ghosts,solver->physics,d,
                         index,index1,index2,ndims,nvars); CHECKERR(ierr);

    /* calculate interface flux array dimensions */
    ierr = ArrayCopy1D_int(dim,dim_interface,ndims); CHECKERR(ierr); dim_interface[d]++;
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];
    /* allocate interface array for conservative discretization */
    FluxI = (double*) calloc (size_interface*nvars,sizeof(double));
    /* compute interface fluxes */

    /* calculate the first derivative */
    done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
    while (!done) {
      ierr = ArrayCopy1D_int(index,index1,ndims); CHECKERR(ierr);
      ierr = ArrayCopy1D_int(index,index2,ndims); CHECKERR(ierr); index2[d]++;
      int p  = ArrayIndex1D(ndims,dim          ,index ,NULL,ghosts);
      int p1 = ArrayIndex1D(ndims,dim_interface,index1,NULL,0     );
      int p2 = ArrayIndex1D(ndims,dim_interface,index2,NULL,0     );
      for (v=0; v<nvars; v++) solver->hyp[nvars*p+v] = dxinv * (FluxI[nvars*p2+v] - FluxI[nvars*p1+v]);
      done = ArrayIncrementIndex(ndims,dim,index);
    }

    /* free interface array */
    free(FluxI);
  }

  free(index );
  free(index1);
  free(index2);

  return(0);
}

int AdvectionFlux(double *u,double *f,int *dim,int ghosts,void *prm,int dir,
                  int *index,int *bounds,int *offset,int ndims,int nvars)
{
  LinearADR *param = (LinearADR*) prm;
  int       ierr   = 0, i, v;

  /* set bounds for array index to include ghost points */
  ierr = ArrayCopy1D_int(dim,bounds,ndims); CHECKERR(ierr);
  for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  ierr = ArraySetValue_int(offset,ndims,-ghosts); CHECKERR(ierr);

  int done = 0; ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,offset,ghosts);
    for (v = 0; v < nvars; v++) f[nvars*p+v] = param->a[nvars*dir+v] * u[nvars*p+v];
    done = ArrayIncrementIndex(ndims,bounds,index);
  }
  return(0);
}
