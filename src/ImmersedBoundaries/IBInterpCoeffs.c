/*! @file IBInterpCoeffs.c
    @brief Compute interpolation nodes and coefficients for immersed boundary points.
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <immersedboundaries.h>

/*! 
  Compute the interpolation nodes and coefficients for immersed boundary points: For each
  immersed boundary point, do the following:
  + From the immersed boundary point, extend a probe in the direction defined by the outward
    normal of the "nearest" facet (computed in IBNearestFacetNormal()), till the probe tip
    reaches a point in space such that all surrounding (#_IB_NNODES_) grid points are "interior" points,
    i.e., outside the immersed body (they are "interior" to the computational domain).
  + Store the indices of the surrounding grid points, as well as the trilinear interpolation
    coefficients to interpolate a variable from the surrounding points to the probe tip.
*/
int IBInterpCoeffs(
                    void    *ib,    /*!< Immersed boundary object of type #ImmersedBoundary */
                    void    *m,     /*!< MPI object of type #MPIVariables */
                    double  *X,     /*!< Array of (local) spatial coordinates */
                    int     *dim_l, /*!< Integer array of local grid size in each spatial dimension */
                    int     ghosts, /*!< Number of ghost points */
                    double  *blank  /*!< Blanking array: for grid points within the
                                         body, this value will be set to 0 */
                  )
{
  ImmersedBoundary  *IB       = (ImmersedBoundary*) ib;
  MPIVariables      *mpi      = (MPIVariables*) m;
  IBNode            *boundary = IB->boundary;

  double  eps         = IB->tolerance;
  int     maxiter     = IB->itr_max,
          n_boundary  = IB->n_boundary_nodes;

  int imax        = dim_l[0],
      jmax        = dim_l[1],
      kmax        = dim_l[2];

  int is = mpi->is[0],
      js = mpi->is[1],
      ks = mpi->is[2];

  static int index[_IB_NDIMS_];

	int dg;
	for (dg = 0; dg < n_boundary; dg++) {
		int    i, j, k, p;
		double xb, yb, zb;
		double nx, ny, nz;
		double xx, yy, zz;
		double dx, dy, dz;
    double ds, dist;
		double xtip, ytip, ztip;

		index[0] = i = boundary[dg].i;
		index[1] = j = boundary[dg].j;
		index[2] = k = boundary[dg].k;
    _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,p);

		xb = boundary[dg].x;
		yb = boundary[dg].y;
		zb = boundary[dg].z;

		nx = boundary[dg].face->nx;
		ny = boundary[dg].face->ny;
		nz = boundary[dg].face->nz;
		xx = boundary[dg].face->x1;
		yy = boundary[dg].face->y1;
		zz = boundary[dg].face->z1;

		dist = nx*(xx-xb) + ny*(yy-yb) + nz*(zz-zb);

    double x1, x2, y1, y2, z1, z2;
    _GetCoordinate_(0,(i+1),dim_l,ghosts,X,x1);
    _GetCoordinate_(0,(i-1),dim_l,ghosts,X,x2);
    _GetCoordinate_(1,(j+1),dim_l,ghosts,X,y1);
    _GetCoordinate_(1,(j-1),dim_l,ghosts,X,y2);
    _GetCoordinate_(2,(k+1),dim_l,ghosts,X,z1);
    _GetCoordinate_(2,(k-1),dim_l,ghosts,X,z2);
    dx = 0.5 * (x1 - x2);
    dy = 0.5 * (y1 - y2);
    dz = 0.5 * (z1 - z2);
		ds = min3(dx, dy, dz);

		xtip = xb + dist*nx;
		ytip = yb + dist*ny;
		ztip = zb + dist*nz;

		int is_it_in = 0;
		int iter = 0;
		int itip, jtip, ktip;
		while(!is_it_in && (iter < maxiter)) {
			iter++;
			itip = i;
			jtip = j;
			ktip = k;
		  
			if (xtip > xb)  {
        double xx;
        _GetCoordinate_(0,itip,dim_l,ghosts,X,xx);
        while ((xx < xtip) && (itip < imax+ghosts-1)) {
          itip++;
          _GetCoordinate_(0,itip,dim_l,ghosts,X,xx);
        }
      }	else {
        double xx;
        _GetCoordinate_(0,(itip-1),dim_l,ghosts,X,xx);
        while ((xx > xtip) && (itip > -ghosts)) {
          itip--;
          _GetCoordinate_(0,(itip-1),dim_l,ghosts,X,xx);
        }
      }

			if (ytip > yb) {
        double yy;
        _GetCoordinate_(1,jtip,dim_l,ghosts,X,yy);
        while ((yy < ytip) && (jtip < jmax+ghosts-1)) {
          jtip++;
          _GetCoordinate_(1,jtip,dim_l,ghosts,X,yy);
        }
      } else {
        double yy;
        _GetCoordinate_(1,(jtip-1),dim_l,ghosts,X,yy);
        while ((yy > ytip) && (jtip > -ghosts)) {
          jtip--;
          _GetCoordinate_(1,(jtip-1),dim_l,ghosts,X,yy);
        }
      }

			if (ztip > zb) {
        double zz;
        _GetCoordinate_(2,ktip,dim_l,ghosts,X,zz);
        while ((zz < ztip) && (ktip < kmax+ghosts-1))	{
          ktip++;
          _GetCoordinate_(2,ktip,dim_l,ghosts,X,zz);
        }
      } else {
        double zz;
        _GetCoordinate_(2,(ktip-1),dim_l,ghosts,X,zz);
        while ((zz > ztip) && (ktip > -ghosts))	{
          ktip--;
          _GetCoordinate_(2,(ktip-1),dim_l,ghosts,X,zz);
        }
      }

      int ptip[_IB_NNODES_];
      index[0] = itip  ; index[1] = jtip  ; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[0]);
      index[0] = itip-1; index[1] = jtip  ; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[1]);
      index[0] = itip  ; index[1] = jtip-1; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[2]);
      index[0] = itip  ; index[1] = jtip  ; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[3]);
      index[0] = itip-1; index[1] = jtip-1; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[4]);
      index[0] = itip  ; index[1] = jtip-1; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[5]);
      index[0] = itip-1; index[1] = jtip  ; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[6]);
      index[0] = itip-1; index[1] = jtip-1; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[7]);

			int nflow = 0;
			nflow += blank[ptip[0]];
			nflow += blank[ptip[1]];
			nflow += blank[ptip[2]];
			nflow += blank[ptip[3]];
			nflow += blank[ptip[4]];
			nflow += blank[ptip[5]];
			nflow += blank[ptip[6]];
			nflow += blank[ptip[7]];
			if (nflow == _IB_NNODES_) {
				is_it_in = 1;
			} else if (nflow < _IB_NNODES_) {
				is_it_in = 0;
				xtip += nx*absolute(ds);
				ytip += ny*absolute(ds);
				ztip += nz*absolute(ds);
			} else {
				fprintf(stderr,"Error in IBInterpCoeffs() (Bug in code) - counting interior points surrounding probe tip \n");
        fprintf(stderr,"on rank %d.\n", mpi->rank);
        fprintf(stderr,"Value of nflow is %d but can only be positive and <= %d.\n",nflow,_IB_NNODES_);
        return(1);
			}
		}

    if (!is_it_in) {
      fprintf(stderr,"Error in IBInterpCoeffs() on rank %d - interior point not found for immersed boundary point (%d,%d,%d)!\n",
              mpi->rank, i, j, k);
      return(1);
    }    
	
		double tlx[2],tly[2],tlz[2];
    _GetCoordinate_(0,(itip-1),dim_l,ghosts,X,tlx[0]);
    _GetCoordinate_(0,(itip  ),dim_l,ghosts,X,tlx[1]);
    _GetCoordinate_(1,(jtip-1),dim_l,ghosts,X,tly[0]);
    _GetCoordinate_(1,(jtip  ),dim_l,ghosts,X,tly[1]);
    _GetCoordinate_(2,(ktip-1),dim_l,ghosts,X,tlz[0]);
    _GetCoordinate_(2,(ktip  ),dim_l,ghosts,X,tlz[1]);

    int ptip[_IB_NNODES_];
    index[0]=itip-1; index[1]=jtip-1; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[0]);
    index[0]=itip  ; index[1]=jtip-1; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[1]);
    index[0]=itip-1; index[1]=jtip  ; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[2]);
    index[0]=itip  ; index[1]=jtip  ; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[3]);
    index[0]=itip-1; index[1]=jtip-1; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[4]);
    index[0]=itip  ; index[1]=jtip-1; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[5]);
    index[0]=itip-1; index[1]=jtip  ; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[6]);
    index[0]=itip  ; index[1]=jtip  ; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,ptip[7]);
    _ArrayCopy1D_(ptip,boundary[dg].interp_nodes,_IB_NNODES_);

    double coeffs[_IB_NNODES_];
    TrilinearInterpCoeffs(tlx[0],tlx[1],tly[0],tly[1],tlz[0],tlz[1],xtip,ytip,ztip,&coeffs[0]);
    _ArrayCopy1D_(coeffs,boundary[dg].interp_coeffs,_IB_NNODES_);

		double tipdist = absolute(nx*(xx-xtip) + ny*(yy-ytip) + nz*(zz-ztip));
    boundary[dg].interp_node_distance = tipdist;
    boundary[dg].surface_distance = absolute(dist);
		if (tipdist < eps) {
			fprintf(stderr,"Warning in IBInterpCoeffs() on rank %d - how can probe tip be on surface? Tipdist = %e\n",
              mpi->rank,tipdist);
		}

	}
  return(0);
}
