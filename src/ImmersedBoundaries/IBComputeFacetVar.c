/*! @file IBComputeFacetVar.c
    @brief Compute the interpolated value of a provided variable on immersed body surface
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <arrayfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>
#include <hypar.h>

/*! Calculate the interpolated value of a provided variables on the
 *  immersed body surface: for each "local facet", i.e., facets 
 *  (of the immersed body surface) that lie within the local 
 *  computational domain of this MPI rank, compute the interpolated
 *  value.
 *
 *  The variable should be a grid variable with size/layout the
 *  same as the solution variable (#HyPar::u) with the appropriate
 *  number of ghost points.
 *
 *  If the incoming variable has multiple components, this
 *  function will calculate interpolated value of each component.
 *
 *  The grad var array should be NULL at input; at output, it 
 *  will point to an array of size (#ImmersedBoundary::nfacets_local X nvars)
 *  that contains the interpolated value at each local facet. If there
 *  are no local facets, it will remain NULL.
 *
 *  The interpolation is bi/tri-linear (second-order).
*/
int IBComputeFacetVar(void*               s,       /*!< Solver object of type #HyPar */
                      void*               m,       /*!< MPI object of type #MPIVariables */
                      const double* const var,     /*!< Variable to compute the interpolated value of */
                      int                 nvars,   /*!< Number of components in var */
                      double** const      face_var /*!< Array to store the interpolated value; 
                                                        must be NULL at input */
                     )
{
  HyPar             *solver  = (HyPar*)          s;
  MPIVariables      *mpi     = (MPIVariables*)   m; 
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->ib;

  if (!solver->flag_ib) return(0);

  if ((*face_var) != NULL) {
    fprintf(stderr,"Error in IBComputeFacetVar()\n");
    fprintf(stderr," face_var is not NULL on rank %d\n",
            mpi->rank );
    return 1;
  }

  int nfacets_local = IB->nfacets_local;
  FacetMap *fmap = IB->fmap;

  if (nfacets_local > 0) {
    (*face_var) = (double*) calloc (nvars*nfacets_local, sizeof(double));

    for (int n = 0; n < nfacets_local; n++) {
  
      double *alpha;
      int    *nodes, j, k;
  
      double *v_c = (*face_var) + n*nvars;
      alpha = &(fmap[n].interp_coeffs[0]);
      nodes = &(fmap[n].interp_nodes[0]);
      _ArraySetValue_(v_c,nvars,0.0);
      for (j=0; j<_IB_NNODES_; j++) {
        for (k=0; k<nvars; k++) {
          v_c[k] += ( alpha[j] * var[nvars*nodes[j]+k] );
        }
      }
    }

  }

  return(0);
}
