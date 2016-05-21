/*! @file IBIdentifyBoundary.c
    @brief Identiy the boundary nodes of the immersed body
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <immersedboundaries.h>

/*! Count the number of immersed boundary points: boundary points are those
    grid points inside the immersed body that are within stencil-width-distance of
    a grid point outside the body.
*/
static int CountBoundaryPoints(
                                int     imax,   /*!< Number of grid points in x */
                                int     jmax,   /*!< Number of grid points in y */
                                int     kmax,   /*!< Number of grid points in z */
                                int     ghosts, /*!< Number of ghost points */
                                double  *blank  /*!< blanking array where entries are zero
                                                     for grid points inside, and one for
                                                     grid points outside. */
                              )
{
  static int dim[_IB_NDIMS_], indexC[_IB_NDIMS_], indexN[_IB_NDIMS_], 
             i, j, k, p, q, count;
  dim[0] = imax;
  dim[1] = jmax;
  dim[2] = kmax;

  count = 0;
	for (i = 0; i < imax; i++) {
		for (j = 0; j < jmax; j++) {
			for (k = 0; k < kmax; k++) {
        indexC[0] = i;
        indexC[1] = j;
        indexC[2] = k;
				_ArrayIndex1D_(_IB_NDIMS_,dim,indexC,ghosts,p);
        /* if this point is inside the body (0), find out if any */
        /* of the neighboring points are outside (1)              */
				if (!blank[p]){

          int g, flag = 0;
					for (g = 1; g <= ghosts; g++){

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

					}
          if (flag)   count++;
				}
			}
		}
	}
  return(count);
}

/*! Set the indices of the immersed boundary points.*/
static int SetBoundaryPoints(
                                int     imax,     /*!< Number of grid points in x */
                                int     jmax,     /*!< Number of grid points in y */
                                int     kmax,     /*!< Number of grid points in z */
                                int     ghosts,   /*!< Number of ghost points */
                                double  *blank,   /*!< blanking array where entries are zero
                                                     for grid points inside, and one for
                                                     grid points outside. */
                                void    *b        /*!< Array of immersed boundary points of type #IBNode */
                            )
{
  IBNode *boundary = (IBNode*) b;
  static int dim[_IB_NDIMS_], indexC[_IB_NDIMS_], indexN[_IB_NDIMS_], 
             i, j, k, p, q, count;
  dim[0] = imax;
  dim[1] = jmax;
  dim[2] = kmax;

  count = 0;
	for (i = 0; i < imax; i++) {
		for (j = 0; j < jmax; j++) {
			for (k = 0; k < kmax; k++) {
        indexC[0] = i;
        indexC[1] = j;
        indexC[2] = k;
				_ArrayIndex1D_(_IB_NDIMS_,dim,indexC,ghosts,p);
        /* if this point is inside the body (0), find out if any */
        /* of the neighboring points are outside (1)              */
				if (!blank[p]){

          int g, flag = 0;
					for (g = 1; g <= ghosts; g++){

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,ghosts,q);
						if (blank[q])	flag = 1;

					}
          if (flag) {
            boundary[count].i = i;
            boundary[count].j = j;
            boundary[count].k = k;
            boundary[count].p = p;
            count++;
          }
				}
			}
		}
	}
  return(count);
}

/*! Identify the immersed boundary points: an immersed boundary point is any grid point
    inside the immersed body that is within stencil-width-distance of a grid point outside
    the immersed body. This function does the following:
    + count the number of immersed boundary points.
    + allocate the array of immersed boundary points and set their indices.
*/
int IBIdentifyBoundary(
                        void   *ib,     /*!< Immersed boundary object of type #ImmersedBoundary */
                        void   *m,      /*!< MPI object of type #MPIVariables */
                        int    *dim_l,  /*!< local dimensions */
                        int    ghosts,  /*!< number of ghost points */
                        double *blank   /*!< Blanking array: for grid points within the
                                             immersed body, this value will be set to 0 */
                      )
{
  ImmersedBoundary  *IB     = (ImmersedBoundary*) ib;
  MPIVariables      *mpi    = (MPIVariables*) m;
  Body3D            *body   = IB->body;

  int imax = dim_l[0],
      jmax = dim_l[1],
      kmax = dim_l[2];

  int n_boundary_nodes = CountBoundaryPoints(imax,jmax,kmax,ghosts,blank);
  IB->n_boundary_nodes = n_boundary_nodes;
  if (n_boundary_nodes == 0) IB->boundary = NULL;
  else {
    IB->boundary = (IBNode*) calloc (n_boundary_nodes, sizeof(IBNode));
    int check = SetBoundaryPoints(imax,jmax,kmax,ghosts,blank,IB->boundary);
    if (check != n_boundary_nodes) {
      fprintf(stderr,"Error in IBIdentifyBoundary(): Inconsistency encountered when setting boundary indices. ");
      fprintf(stderr,"on rank %d.\n",mpi->rank);
      fprintf(stderr,"SetBoundaryPoints() returned %d, while n_boundary_nodes is %d.\n",check,n_boundary_nodes);
    }
  }

  return(n_boundary_nodes);
}
