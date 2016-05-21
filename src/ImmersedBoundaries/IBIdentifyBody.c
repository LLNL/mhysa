/*! @file IBIdentifyBody.c
    @author Debojyoti Ghosh
    @brief Identify grid points inside immersed body
*/

#include <stdio.h>
#include <stdlib.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>

/*!
  Identify the grid points of the given grid that are 
  inside a given body whose surface is defined as an
  unstructured triangulation. This function uses the
  ray-tracing method and is a copy of a FORTRAN function
  originally written by Dr. Jay Sitaraman.
*/
int IBIdentifyBody(
                    void   *ib,     /*!< Immersed boundary object of type #ImmersedBoundary */
                    int    *dim_g,  /*!< global dimensions */
                    int    *dim_l,  /*!< local dimensions */
                    int    ghosts,  /*!< number of ghost points */
                    void   *m,      /*!< MPI object of type #MPIVariables */
                    double *X,      /*!< Array of global spatial coordinates */
                    double *blank   /*!< Blanking array: for grid points within the
                                         body, this value will be set to 0 */
                  )
{
  ImmersedBoundary  *IB     = (ImmersedBoundary*) ib;
  MPIVariables      *mpi    = (MPIVariables*) m;
  Body3D            *body   = IB->body;
  double            *x      = X, 
                    *y      = X+dim_g[0], 
                    *z      = X+dim_g[0]+dim_g[1],
                    eps     = IB->tolerance;
  int               itr_max = IB->itr_max,
	                  i, j, k, n, v;

	double xmax = body->xmax, 
         xmin = body->xmin, 
         ymax = body->ymax, 
         ymin = body->ymin, 
         zmax = body->zmax, 
         zmin = body->zmin;
  
  double fac = 1.5;
	double Lx, Ly, Lz;
	Lx = (xmax - xmin);
	Ly = (ymax - ymin);
	Lz = (zmax - zmin);
	double xc, yc, zc;
	xc = 0.5 * (xmin + xmax);
	yc = 0.5 * (ymin + ymax);
	zc = 0.5 * (zmin + zmax);
	xmax = xc + fac * Lx/2;
	xmin = xc - fac * Lx/2;
	ymax = yc + fac * Ly/2;
	ymin = yc - fac * Ly/2;
	zmax = zc + fac * Lz/2;
	zmin = zc - fac * Lz/2;

	int imin, imax, jmin, jmax, kmin, kmax;
	imin = dim_g[0]-1;	imax = 0;
	jmin = dim_g[1]-1;	jmax = 0;
	kmin = dim_g[2]-1;	kmax = 0;
	for (i = 0; i < dim_g[0]; i++) {
		for (j = 0; j < dim_g[1]; j++) {
			for (k = 0; k < dim_g[2]; k++) {
				if (   ((x[i]-xmin)*(x[i]-xmax) < 0) 
            && ((y[j]-ymin)*(y[j]-ymax) < 0) 
            && ((z[k]-zmin)*(z[k]-zmax) < 0)) {
					imin = min(i, imin);
					imax = max(i, imax);
					jmin = min(j, jmin);
					jmax = max(j, jmax);
					kmin = min(k, kmin);
					kmax = max(k, kmax);
				}
			}
		}
	}

	double *cof[5], *xd, *dist;
  xd    = (double*) calloc (itr_max,sizeof(double));
  dist  = (double*) calloc (itr_max,sizeof(double));
	for (v = 0; v < 5; v++) cof[v] = (double*) calloc(body->nfacets, sizeof(double));

	for (n = 0; n < body->nfacets; n++) {
		cof[0][n] = body->surface[n].y1 - body->surface[n].y3;
		cof[1][n] = body->surface[n].y2 - body->surface[n].y3;
    cof[2][n] = body->surface[n].z1 - body->surface[n].z3;
    cof[3][n] = body->surface[n].z2 - body->surface[n].z3;
		double den = cof[0][n]*cof[3][n] - cof[1][n]*cof[2][n];
		if (absolute(den) > eps)	cof[4][n] = 1.0 / den;
		else			                cof[4][n] = 0;		
	}

  int count = 0;
	for (j = jmin; j <= jmax; j++) {
		for (k = kmin; k <= kmax; k++) {
			int itr = 0;
			for (n = 0; n < body->nfacets; n++) {
				if (cof[4][n] != 0) {
					double yy, zz;
					yy = body->surface[n].y3 - y[j];
					zz = body->surface[n].z3 - z[k];
					double l1, l2, l3;
					l1 = (cof[1][n]*zz - cof[3][n]*yy) * cof[4][n];
					l2 = (cof[2][n]*yy - cof[0][n]*zz) * cof[4][n];
					l3 = 1 - l1 - l2;
					if ((l1 > -eps) && (l2 > -eps) && (l3 > -eps)){
						xd[itr]   = l1*body->surface[n].x1 + l2*body->surface[n].x2 + l3*body->surface[n].x3;
						dist[itr] = abs(x[imin]-xd[itr]);
						itr++;
					}
				}
			}
			if (itr > itr_max) {
        if (!mpi->rank) {
				  fprintf(stderr,"Error: In IBIdentyBody() - itr > %d. Recompilation of code needed.\n",itr_max);
				  fprintf(stderr,"Increase the value of \"itr_max\" in IBIdentyBody().\n");
        }
        return(1);
			}
	
			if (itr > 0) {
				int ii, jj;
	
				for (ii = 0; ii < itr; ii++) {
					for (jj = ii+1; jj < itr; jj++) {
						if (dist[jj] < dist[ii]) {
							double temp;
							temp     = dist[ii];
							dist[ii] = dist[jj];
							dist[jj] = temp;
							temp     = xd[ii];
							xd[ii]   = xd[jj];
							xd[jj]   = temp;
						}
					}
				}

				for (ii = 1; ii < itr-1; ii++) {
					if (abs(xd[ii]-xd[ii-1]) < eps) {
						for (jj = ii+1; jj < itr; jj++) {
							xd[jj-1]    = xd[jj];
							dist[jj-1]  = dist[jj];
						}
						itr--;
						ii--;
					}
				}

				ii = 0;
				int inside = 0;
				for (i = imin; i <= imax; i++) {
					if ((x[i]-xd[ii])*(x[i]-xd[ii+1]) < eps) {
						inside = 1;
            /* this point is inside                      */
            /* check if (i,j,k) lies within this process */
            /* if so, set blank                          */
            if (   ((i-mpi->is[0])*(i-mpi->ie[0]) <= 0) 
                && ((j-mpi->is[1])*(j-mpi->ie[1]) <= 0) 
                && ((k-mpi->is[2])*(k-mpi->ie[2]) <= 0)) {
              /* calculate local indices */
              int index[_IB_NDIMS_];
              index[0] = i-mpi->is[0];
              index[1] = j-mpi->is[1];
              index[2] = k-mpi->is[2];
              int p; _ArrayIndex1D_(_IB_NDIMS_,dim_l,index,ghosts,p);
              blank[p] = 0;
              count++;
            }
					} else {
						if (inside) {
							if (ii+2 < itr-1)	ii += 2;
							inside = 0;
						}
					}
				}

			}
		}
	}

  free(xd);
  free(dist);
  return(count);
}
