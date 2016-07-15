/*! @file InitializeImmersedBoundaries.c
    @author Debojyoti Ghosh
    @brief Initialize the immersed boundary implementation
*/

#include <stdio.h>
#include <stdlib.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>
#include <io.h>
#include <hypar.h>

/*! Initialize the immersed boundaries, if present.
    + Read in immersed body from STL file.
    + Allocate and set up #ImmersedBoundary object.
    + Identify blanked-out grid points based on immersed body geometry.
    + Identify and make a list of immersed boundary points on each rank.
    + For each immersed boundary point, find the "nearest" facet.
*/
int InitializeImmersedBoundaries(
                                  void *s,  /*!< Solver object of type #HyPar */
                                  void *m   /*!< MPI object of type #MPIVariables */
                                )
{
  HyPar             *solver   = (HyPar*)        s;
  MPIVariables      *mpi      = (MPIVariables*) m;
  ImmersedBoundary  *ib       = NULL;
  Body3D            *body     = NULL;
  int               stat, d, ndims = solver->ndims;

  if ((!solver->flag_ib) || (ndims != _IB_NDIMS_)) {
    solver->ib = NULL;
    return(0);
  }

  /* Read in immersed body from file */
  IBReadBodySTL(&body,solver->ib_filename,mpi,&stat);
  if (stat) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in InitializeImmersedBoundaries(): Unable to ");
      fprintf(stderr,"read immersed body from file %s.\n",solver->ib_filename);
    }
    solver->flag_ib = 0;
    solver->ib = NULL;
    return(1);
  }
  IBComputeBoundingBox(body);

  /* allocate immersed boundary object and set it up */
  ib = (ImmersedBoundary*) calloc (1, sizeof(ImmersedBoundary));
  ib->tolerance = 1e-12;
  ib->delta     = 1e-6;
  ib->itr_max   = 500;
  ib->body      = body;
  solver->ib    = ib;

  int     offset_global, offset_local,
          *dim_local  = solver->dim_local,
          *dim_global = solver->dim_global,
          ghosts      = solver->ghosts,
          size        = dim_global[0] + dim_global[1] + dim_global[2],
          count       = 0;
  double  *Xg         = (double*) calloc(size,sizeof(double));

  /* assemble the global grid on rank 0 */
  offset_global = offset_local = 0;
  for (d=0; d<ndims; d++) {
    IERR MPIGatherArray1D(mpi,(mpi->rank?NULL:&Xg[offset_global]),
                          &solver->x[offset_local+ghosts],
                          mpi->is[d],mpi->ie[d],dim_local[d],0); CHECKERR(ierr);
    offset_global += dim_global[d];
    offset_local  += dim_local [d] + 2*ghosts;
  }
  /* send the global grid to other ranks */
  MPIBroadcast_double(Xg,size,0,&mpi->world);

  /* identify whether this is a 3D or "pseudo-2D" simulation */
  IBIdentifyMode(Xg,dim_global,solver->ib);

  /* identify grid points inside the immersed body */
  int count_inside_body = 0;
  count = IBIdentifyBody(solver->ib,dim_global,dim_local,ghosts,mpi,Xg,solver->iblank);
  MPISum_integer(&count_inside_body,&count,1,&mpi->world);
  free(Xg);

  /* At ghost points corresponding to the physical boundary, extrapolate from the interior 
     (this should also work for bodies that are adjacent to physical boundaries). At interior
     (MPI) boundaries, exchange iblank across MPI ranks.
  */
  int indexb[ndims], indexi[ndims], bounds[ndims], offset[ndims];
  for (d = 0; d < ndims; d++) {
    /* left boundary */
    if (!mpi->ip[d]) {
      _ArrayCopy1D_(dim_local,bounds,ndims); bounds[d] = ghosts;
      _ArraySetValue_(offset,ndims,0); offset[d] = -ghosts;
      int done = 0; _ArraySetValue_(indexb,ndims,0);
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,ndims); indexi[d] = ghosts-1-indexb[d];
        int p1; _ArrayIndex1DWO_(ndims,dim_local,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (ndims,dim_local,indexi,ghosts,p2);
        solver->iblank[p1] = solver->iblank[p2];
        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }
    /* right boundary */
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      _ArrayCopy1D_(dim_local,bounds,ndims); bounds[d] = ghosts;
      _ArraySetValue_(offset,ndims,0); offset[d] = dim_local[d];
      int done = 0; _ArraySetValue_(indexb,ndims,0);
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,ndims); indexi[d] = dim_local[d]-1-indexb[d];
        int p1; _ArrayIndex1DWO_(ndims,dim_local,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (ndims,dim_local,indexi,ghosts,p2);
        solver->iblank[p1] = solver->iblank[p2];
        _ArrayIncrementIndex_(ndims,bounds,indexb,done);
      }
    }
  }
  MPIExchangeBoundariesnD(ndims,1,dim_local,ghosts,mpi,solver->iblank);

  /* identify and create a list of immersed boundary points on each rank */
  int count_boundary_points = 0;
  count = IBIdentifyBoundary(solver->ib,mpi,dim_local,ghosts,solver->iblank);
  MPISum_integer(&count_boundary_points,&count,1,&mpi->world);

  /* find the nearest facet for each immersed boundary point */
  double ld = 0, xmin, xmax, ymin, ymax, zmin, zmax;
  _GetCoordinate_(0,0             ,dim_local,ghosts,solver->x,xmin);
  _GetCoordinate_(0,dim_local[0]-1,dim_local,ghosts,solver->x,xmax);
  _GetCoordinate_(1,0             ,dim_local,ghosts,solver->x,ymin);
  _GetCoordinate_(1,dim_local[1]-1,dim_local,ghosts,solver->x,ymax);
  _GetCoordinate_(2,0             ,dim_local,ghosts,solver->x,zmin);
  _GetCoordinate_(2,dim_local[2]-1,dim_local,ghosts,solver->x,zmax);
  double xlen = xmax - xmin;
  double ylen = ymax - ymin;
  double zlen = zmax - zmin;
  ld = max3(xlen,ylen,zlen);
  count = IBNearestFacetNormal(solver->ib,mpi,solver->x,ld,dim_local,ghosts);
  if (count) return(count);

  /* For the immersed boundary points, find the interior points for extrapolation,
     and compute their interpolation coefficients */
  count = IBInterpCoeffs(solver->ib,mpi,solver->x,dim_local,ghosts,solver->iblank);
  if (count) return(count);

  /* Create facet mapping */;
  count = IBCreateFacetMapping(ib,mpi,solver->x,dim_local,ghosts);
  if (count) return(count);
  
  /* Done */
  if (!mpi->rank) {
    double percentage;
    printf("Immersed body read from %s:\n",solver->ib_filename);
    printf("    Number of facets: %d\n    Bounding box: [%3.1f,%3.1lf] X [%3.1f,%3.1lf] X [%3.1f,%3.1lf]\n",
           body->nfacets,body->xmin,body->xmax,body->ymin,body->ymax,body->zmin,body->zmax);
    percentage = ((double)count_inside_body)/((double)solver->npoints_global)*100.0;
    printf("    Number of grid points inside immersed body: %d (%4.1f%%).\n",count_inside_body,percentage);
    percentage = ((double)count_boundary_points)/((double)solver->npoints_global)*100.0;
    printf("    Number of immersed boundary points        : %d (%4.1f%%).\n",count_boundary_points,percentage);
    printf("    Immersed body simulation mode             : %s.\n", ib->mode);
  }

  return(0);
}
