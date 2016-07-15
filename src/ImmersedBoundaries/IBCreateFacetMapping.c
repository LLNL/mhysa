/*! @file IBCreateFacetMapping.c
    @brief Create facet mapping
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <string.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>

/*! is x inside the interval [a,b]? */
static inline int isInside(
                            double x, /*!< the value to check for */
                            double a, /*!< small end of the interval */
                            double b  /*!< big end of the interval */
                          )
{
  return ((x >= a) && (x <= b));
}

/*! Given a point in the 3D space (xc, yc, zc), this function finds the
    indices of the 8 grid points that define the grid cell the given 
    point is in, as well as the trilinear interpolation coefficients
    for each of the surrounding grid points.
*/
static int interpNodesCoeffs(
                              void    *m,       /*!< MPI object of type #MPIVariables */
                              double  xc,       /*!< x-coordinate of the point */
                              double  yc,       /*!< y-coordinate of the point */
                              double  zc,       /*!< z-coordinate of the point */
                              double  *x,       /*!< array of x-coordinates of the grid */
                              double  *y,       /*!< array of y-coordinates of the grid */
                              double  *z,       /*!< array of z-coordinates of the grid */
                              int     *dim,     /*!< local dimensions of the grid */
                              int     ghosts,   /*!< number of ghost points */
                              char    *mode,    /*!< "mode", i.e., #ImmersedBoundary::mode */
                              int     *ii,      /*!< i-index of the surrounding node at the high end
                                                     (i.e. smallest i such that x[i] > xc) */
                              int     *jj,      /*!< j-index of the surrounding node at the high end
                                                     (i.e. smallest j such that y[j] > yc) */
                              int     *kk,      /*!< k-index of the surrounding node at the high end
                                                     (i.e. smallest k such that z[k] > zc) */
                              int     *inodes,  /*!< array to store the indices of the surrounding nodes */
                              double  *icoeffs  /*!< array to store the interpolation coefficients of the surrounding nodes */
                            )
{
  MPIVariables *mpi = (MPIVariables*) m;

  int i, j, k, ic, jc, kc;
  ic = jc = kc = -1;
  
  double  xmin = 0.5 * (x[ghosts-1]         + x[ghosts]),
          xmax = 0.5 * (x[dim[0]+ghosts-1]  + x[dim[0]+ghosts]),
          ymin = 0.5 * (y[ghosts-1]         + y[ghosts]),
          ymax = 0.5 * (y[dim[1]+ghosts-1]  + y[dim[1]+ghosts]),
          zmin = 0.5 * (z[ghosts-1]         + z[ghosts]),
          zmax = 0.5 * (z[dim[2]+ghosts-1]  + z[dim[2]+ghosts]);

  for (i = 0; i < dim[0]+2*ghosts-1; i++) {
    if (isInside(xc,x[i],x[i+1])) {
      ic = i;
      break;
    }
  }
  if      (ic <= ghosts-1)        ic = ghosts;
  else if (ic >= dim[0]+ghosts-1) ic = dim[0]+ghosts-2;
  
  for (j = 0; j < dim[1]+2*ghosts-1; j++) {
    if (isInside(yc,y[j],y[j+1])) {
      jc = j;
      break;
    }
  }
  if      (jc <= ghosts-1)        jc = ghosts;
  else if (jc >= dim[1]+ghosts-1) jc = dim[1]+ghosts-2;

  for (k = 0; k < dim[2]+2*ghosts-1; k++) {
    if (isInside(zc,z[k],z[k+1])) {
      kc = k;
      break;
    }
  }
  if      (kc <= ghosts-1)        kc = ghosts;
  else if (kc >= dim[2]+ghosts-1) kc = dim[2]+ghosts-2;
  
  if      (!strcmp(mode,_IB_XY_))  { kc = ghosts; zc = 0.5*(zmin+zmax); }
  else if (!strcmp(mode,_IB_XZ_))  { jc = ghosts; yc = 0.5*(ymin+ymax); }
  else if (!strcmp(mode,_IB_YZ_))  { ic = ghosts; xc = 0.5*(xmin+xmax); }
  
  if (ic == -1) {
    fprintf(stderr,"Error in interpNodesCoeffs() (in ImmersedBoundaries/IBCreateFacetMapping.c) on rank %d: ic = -1.\n", mpi->rank);
    return(1);
  }
  if (jc == -1) {
    fprintf(stderr,"Error in interpNodesCoeffs() (in ImmersedBoundaries/IBCreateFacetMapping.c) on rank %d: jc = -1.\n", mpi->rank);
    return(1);
  }
  if (kc == -1) {
    fprintf(stderr,"Error in interpNodesCoeffs() (in ImmersedBoundaries/IBCreateFacetMapping.c) on rank %d: kc = -1.\n", mpi->rank);
    return(1);
  }
  ic++;
  jc++;
  kc++;

  if (ii) *ii = ic;
  if (jj) *jj = jc;
  if (kk) *kk = kc;
  
  int pc[_IB_NNODES_], index[_IB_NDIMS_];
  index[0]=ic-1-ghosts; index[1]=jc-1-ghosts; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[0]);
  index[0]=ic-ghosts  ; index[1]=jc-1-ghosts; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[1]);
  index[0]=ic-1-ghosts; index[1]=jc-ghosts  ; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[2]);
  index[0]=ic-ghosts  ; index[1]=jc-ghosts  ; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[3]);
  index[0]=ic-1-ghosts; index[1]=jc-1-ghosts; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[4]);
  index[0]=ic-ghosts  ; index[1]=jc-1-ghosts; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[5]);
  index[0]=ic-1-ghosts; index[1]=jc-ghosts  ; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[6]);
  index[0]=ic-ghosts  ; index[1]=jc-ghosts  ; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[7]);
  _ArrayCopy1D_(pc,inodes,_IB_NNODES_);
  
  double coeffs[_IB_NNODES_];
  TrilinearInterpCoeffs(x[ic-1],x[ic],y[jc-1],y[jc],z[kc-1],z[kc],xc,yc,zc,&coeffs[0]);
  _ArrayCopy1D_(coeffs,icoeffs,_IB_NNODES_);

  return(0);
}

/*!
  This function creates a "facet map", i.e., on each MPI rank, it does the following:
  + Makes a list of facets (defining the immersed body surface) that lie within the
    local computational domain of this MPI rank ("local facets").
  + For each local facet, finds and stores the indices of the grid points that 
    surround it, as well as the trilinear interpolation coefficients.
  + For each local facet, finds a "near-surface" point, i.e., a point near the surface
    ("near" in terms of the local grid spacing) along the outward surface normal (i.e.,
    outside the body), and finds and stores the indices of the grid points that
    surround it, as well as the trilinear interpolation coefficients.

Note: each MPI rank has a copy of the entire immersed body, i.e., all the facets.
*/
int IBCreateFacetMapping(
                          void    *ib,    /*!< Immersed boundary object of type #ImmersedBoundary */
                          void    *m,     /*!< MPI object of type #MPIVariables */
                          double  *X,     /*!< Array of local spatial coordinates */
                          int     *dim,   /*!< Local dimensions */
                          int     ghosts  /*!< Number of ghost points */
                        )
{
  ImmersedBoundary  *IB     = (ImmersedBoundary*) ib;
  MPIVariables      *mpi    = (MPIVariables*) m;
  Body3D            *body   = IB->body;
  int               nfacets = body->nfacets, n, count, ierr;
  Facet3D           *facets = body->surface;

  double  *x = X, 
          *y = (x + dim[0] + 2*ghosts), 
          *z = (y + dim[1] + 2*ghosts);

  double  xmin = 0.5 * (x[ghosts-1]         + x[ghosts]),
          xmax = 0.5 * (x[dim[0]+ghosts-1]  + x[dim[0]+ghosts]),
          ymin = 0.5 * (y[ghosts-1]         + y[ghosts]),
          ymax = 0.5 * (y[dim[1]+ghosts-1]  + y[dim[1]+ghosts]),
          zmin = 0.5 * (z[ghosts-1]         + z[ghosts]),
          zmax = 0.5 * (z[dim[2]+ghosts-1]  + z[dim[2]+ghosts]);

  count = 0;
  for (n = 0; n < nfacets; n++) {

    /* find facet centroid */
    double xc, yc, zc;
    xc = (facets[n].x1 + facets[n].x2 + facets[n].x3) / 3.0;
    yc = (facets[n].y1 + facets[n].y2 + facets[n].y3) / 3.0;
    zc = (facets[n].z1 + facets[n].z2 + facets[n].z3) / 3.0;

    if (!strcmp(IB->mode,_IB_3D_)) {
      if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) count++;
    } else if (!strcmp(IB->mode,_IB_XY_)) {
      if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax)) count++;
    } else if (!strcmp(IB->mode,_IB_XZ_)) {
      if (isInside(xc,xmin,xmax) && isInside(zc,zmin,zmax)) count++;
    } else if (!strcmp(IB->mode,_IB_YZ_)) {
      if (isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) count++;
    }
  }

  int nfacets_local = count;
  int nfacets_global;
  MPISum_integer(&nfacets_global,&nfacets_local,1,&mpi->world);
  if (nfacets_global != nfacets) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in IBCreateFacetMapping(): nfacets_global = %d, ", nfacets_global);
      fprintf(stderr,"but nfacets = %d.\n", nfacets);
    }
    return(1);
  }

  if (nfacets_local > 0) {

    FacetMap *fmap = (FacetMap*) calloc (nfacets_local, sizeof(FacetMap));
    count = 0;

    for (n = 0; n < nfacets; n++) {

      double xc, yc, zc;
      xc = (facets[n].x1 + facets[n].x2 + facets[n].x3) / 3.0;
      yc = (facets[n].y1 + facets[n].y2 + facets[n].y3) / 3.0;
      zc = (facets[n].z1 + facets[n].z2 + facets[n].z3) / 3.0;

      int flag = 0;
      if (!strcmp(IB->mode,_IB_3D_)) {
        if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) flag = 1;
      } else if (!strcmp(IB->mode,_IB_XY_)) {
        if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax)) flag = 1;
      } else if (!strcmp(IB->mode,_IB_XZ_)) {
        if (isInside(xc,xmin,xmax) && isInside(zc,zmin,zmax)) flag = 1;
      } else if (!strcmp(IB->mode,_IB_YZ_)) {
        if (isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) flag = 1;
      }

      if (flag == 1) {
        fmap[count].facet = facets + n;
        fmap[count].index = n;
        
        fmap[count].xc = xc;
        fmap[count].yc = yc;
        fmap[count].zc = zc;

        int ic,jc, kc;

        ierr = interpNodesCoeffs(mpi,xc,yc,zc,x,y,z,dim,ghosts,IB->mode,&ic,&jc,&kc,fmap[count].interp_nodes,fmap[count].interp_coeffs);
        if (ierr) return(ierr);

        double dx = x[ic] - x[ic-1];
        double dy = y[jc] - y[jc-1];
        double dz = z[kc] - z[kc-1];
        double ds;
        if      (!strcmp(IB->mode,_IB_XY_)) ds = min(dx,dy);
        else if (!strcmp(IB->mode,_IB_XZ_)) ds = min(dx,dz);
        else if (!strcmp(IB->mode,_IB_YZ_)) ds = min(dy,dz);
        else                                ds = min3(dx,dy,dz);

        double nx = fmap[count].facet->nx;
        double ny = fmap[count].facet->ny;
        double nz = fmap[count].facet->nz;

        if (nx == 0.0) nx += IB->delta*ds;
        if (ny == 0.0) ny += IB->delta*ds;
        if (nz == 0.0) nz += IB->delta*ds;

        double xns = xc + sign(nx)*ds;
        double yns = yc + sign(ny)*ds;
        double zns = zc + sign(nz)*ds;

        fmap[count].dx = xns - xc;
        fmap[count].dy = yns - yc;
        fmap[count].dz = zns - zc;

        ierr = interpNodesCoeffs(mpi,xns,yns,zns,x,y,z,dim,ghosts,IB->mode,NULL,NULL,NULL,fmap[count].interp_nodes_ns,fmap[count].interp_coeffs_ns);
        if (ierr) return(ierr);

        count++;
      }
    }

    IB->nfacets_local = nfacets_local;
    IB->fmap = fmap;  

  } else {

    IB->nfacets_local = 0;
    IB->fmap = NULL;

  }

  return(0);
}
