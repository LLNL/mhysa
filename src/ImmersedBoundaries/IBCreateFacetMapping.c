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

/*!
  This function creates a "facet map", i.e., on each MPI rank, it does the following:
  + Makes a list of facets (defining the immersed body surface) that lie within the
    local computational domain of this MPI rank ("local facets").
  + For each local facets, finds and stores the indices of the grid points that 
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
  int               nfacets = body->nfacets, n, count;
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
        int i, j, k, ic, jc, kc;

        ic = jc = kc = -1;
        for (i = ghosts-1; i < dim[0]+ghosts; i++) {
          if (isInside(xc,x[i],x[i+1])) {
            ic = i;
            break;
          }
        }
        for (j = ghosts-1; j < dim[1]+ghosts; j++) {
          if (isInside(yc,y[j],y[j+1])) {
            jc = j;
            break;
          }
        }
        for (k = ghosts-1; k < dim[2]+ghosts; k++) {
          if (isInside(zc,z[k],z[k+1])) {
            kc = k;
            break;
          }
        }

        if      (!strcmp(IB->mode,_IB_XY_))  { kc = ghosts; zc = 0.5*(zmin+zmax); }
        else if (!strcmp(IB->mode,_IB_XZ_))  { jc = ghosts; yc = 0.5*(ymin+ymax); }
        else if (!strcmp(IB->mode,_IB_YZ_))  { ic = ghosts; xc = 0.5*(xmin+xmax); }

        if (ic == -1) {
          fprintf(stderr,"Error in IBCreateFacetMapping() on rank %d: ic = -1.\n", mpi->rank);
          return(1);
        }
        if (jc == -1) {
          fprintf(stderr,"Error in IBCreateFacetMapping() on rank %d: jc = -1.\n", mpi->rank);
          return(1);
        }
        if (kc == -1) {
          fprintf(stderr,"Error in IBCreateFacetMapping() on rank %d: kc = -1.\n", mpi->rank);
          return(1);
        }
        ic++;
        jc++;
        kc++;

        int pc[_IB_NNODES_], index[_IB_NDIMS_];
        index[0]=ic-1-ghosts; index[1]=jc-1-ghosts; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[0]);
        index[0]=ic-ghosts  ; index[1]=jc-1-ghosts; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[1]);
        index[0]=ic-1-ghosts; index[1]=jc-ghosts  ; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[2]);
        index[0]=ic-ghosts  ; index[1]=jc-ghosts  ; index[2]=kc-1-ghosts; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[3]);
        index[0]=ic-1-ghosts; index[1]=jc-1-ghosts; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[4]);
        index[0]=ic-ghosts  ; index[1]=jc-1-ghosts; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[5]);
        index[0]=ic-1-ghosts; index[1]=jc-ghosts  ; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[6]);
        index[0]=ic-ghosts  ; index[1]=jc-ghosts  ; index[2]=kc-ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,dim,index,ghosts,pc[7]);
        _ArrayCopy1D_(pc,fmap[count].interp_nodes,_IB_NNODES_);

        double coeffs[_IB_NNODES_];
        TrilinearInterpCoeffs(x[ic-1],x[ic],y[jc-1],y[jc],z[kc-1],z[kc],xc,yc,zc,&coeffs[0]);
        _ArrayCopy1D_(coeffs,fmap[count].interp_coeffs,_IB_NNODES_);

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
