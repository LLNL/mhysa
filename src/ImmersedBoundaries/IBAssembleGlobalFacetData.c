/*! @file IBAssembleGlobalFacetData.c
    @brief Assemble a local array with some facet var to a global one
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <arrayfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>

/*! Assemble any local data computed at immersed body surface (facet centroids)
 *  into a global array.
 *
 *  The local array *must* be of length 
 *  (#ImmersedBoundary::nfacets_local X nvars).
 *
 *  The global array must be NULL at input; it will have the size
 *  (#Body3D::nfacets X nvars) at output on rank 0;
 *  it will remain NULL on other ranks.
*/
int IBAssembleGlobalFacetData(void*               m,          /*!< MPI object of type #MPIVariables */
                              void*               ib,         /*!< Immersed boundary object of type #ImmersedBoundary */
                              const double* const local_var,  /*!< Local array */
                              double** const      global_var, /*!< Array to store the gradient; must be NULL at input */
                              int                 nvars       /*!< Number of components in var */
                             )
{
  MPIVariables     *mpi = (MPIVariables*) m; 
  ImmersedBoundary *IB  = (ImmersedBoundary*) ib;

  if ((*global_var) != NULL) {
    fprintf(stderr,"Error in IBAssembleGlobalFacetData()\n");
    fprintf(stderr," global_var is not NULL on rank %d\n",
            mpi->rank );
    return 1;
  }

  int nfacets_local = IB->nfacets_local;
  FacetMap *fmap = IB->fmap;

  if ((nfacets_local == 0) && (local_var != NULL)) {
    fprintf(stderr, "Error in IBAssembleGlobalFacetData()\n");
    fprintf(stderr, " nfacets_local is 0 but local_var is not NULL!\n");
    return 1;
  }
  if ((nfacets_local > 0) && (local_var == NULL)) {
    fprintf(stderr, "Error in IBAssembleGlobalFacetData()\n");
    fprintf(stderr, " nfacets_local > 0 but local_var is NULL!\n");
    return 1;
  }

#ifndef serial
  MPI_Status status;
#endif

	if (!mpi->rank) {

    int nfacets_global = IB->body->nfacets;

    /* allocate arrays for whole domain */
		*global_var = (double*) calloc (nfacets_global*nvars, sizeof(double));

    /* check array - to avoid double counting of facets */
		int *check = (int*) calloc (nfacets_global, sizeof(int));
		for (int n = 0; n < nfacets_global; n++) check[n] = 0;
		int check_total_facets = 0;

    /* local data */
		for (int n = 0; n < nfacets_local; n++) {

			if (!check[fmap[n].index]) {

        _ArrayCopy1D_((local_var+n), ((*global_var)+fmap[n].index), nvars);
				check[fmap[n].index] = 1;

			} else {

        fprintf(stderr,"Error in IBAssembleGlobalFacetData(): ");
        fprintf(stderr,"Facet %d has already been assigned a value. Double counting of facet!\n",fmap[n].index);
        return 1;

			}

		}
		check_total_facets += nfacets_local;

#ifndef serial
		for (int proc = 1; proc < mpi->nproc; proc++) {

			int nf_incoming;
			MPI_Recv(&nf_incoming, 1, MPI_INT, proc, 98927, MPI_COMM_WORLD, &status);
			check_total_facets += nf_incoming;

      if (nf_incoming > 0) {

			  int    *indices_incoming = (int*)    calloc(nf_incoming, sizeof(int));
			  double *var_incoming     = (double*) calloc(nf_incoming*nvars, sizeof(double));

			  MPI_Recv(indices_incoming, nf_incoming, MPI_INT, proc, 98928, mpi->world, &status);
			  MPI_Recv(var_incoming, nf_incoming*nvars, MPI_DOUBLE, proc, 98929, mpi->world, &status);

			  for (int n = 0; n < nf_incoming; n++) {
				  if (!check[indices_incoming[n]]) {
            _ArrayCopy1D_((var_incoming+n), ((*global_var)+indices_incoming[n]), nvars);
					  check[indices_incoming[n]] = 1;
				  } else {
            fprintf(stderr,"Error in IBAssembleGlobalFacetData(): ");
            fprintf(stderr,"Facet %d has already been assigned a value. Double counting of facet!\n",
                    indices_incoming[n]);
            return 1;
				  }
			  }
        
        free(indices_incoming);
        free(var_incoming);
      }
		}
#endif

		if (check_total_facets != nfacets_global)	{
      fprintf(stderr,"Error in IBAssembleGlobalFacetData(): mismatch in total facet count.\n");
    }

	} else {

#ifndef serial

		MPI_Send(&nfacets_local, 1, MPI_INT, 0, 98927, MPI_COMM_WORLD);

    if (nfacets_local > 0) {

      int i, *indices = (int*) calloc (nfacets_local, sizeof(int));
      for (i = 0; i < nfacets_local; i++) indices[i] = fmap[i].index;

		  MPI_Send(indices, nfacets_local, MPI_INT, 0, 98928, mpi->world);
		  MPI_Send(local_var, nfacets_local*nvars, MPI_DOUBLE, 0, 98929, mpi->world);

      free(indices);
    }
#endif

	}

  return(0);
}
