/*! @file IBComputeNormalGradient.c
    @brief Compute normal gradient of a provided variable on immersed body surface
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <arrayfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>
#include <hypar.h>

/*! Calculate the normal gradient of a provided variables on the
 *  immersed body surface: for each "local facet", i.e., facets 
 *  (of the immersed body surface) that lie within the local 
 *  computational domain of this MPI rank, compute the normal 
 *  gradient.
 *
 *  The variable should be a grid variable with size/layout the
 *  same as the solution variable (#HyPar::u) with the appropriate
 *  number of ghost points.
 *
 *  If the incoming variable has multiple components, this
 *  function will calculate normal gradient of each component.
 *
 *  The grad var array should be NULL at input; at output, it 
 *  will point to an array of size (#ImmersedBoundary::nfacets_local X nvars)
 *  that contains the normal gradient at each local facet. If there
 *  are no local facets, it will remain NULL.
 *
 *  The gradient calculation is currently first-order.
*/
int IBComputeNormalGradient(void*               s,       /*!< Solver object of type #HyPar */
                            void*               m,       /*!< MPI object of type #MPIVariables */
                            const double* const var,     /*!< Variable to compute the gradient of */
                            int                 nvars,   /*!< Number of components in var */
                            double** const      grad_var /*!< Array to store the gradient; must be NULL at input */
                          )
{
  HyPar             *solver  = (HyPar*)          s;
  MPIVariables      *mpi     = (MPIVariables*)   m; 
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->ib;

  if (!solver->flag_ib) return(0);

  if ((*grad_var) != NULL) {
    fprintf(stderr,"Error in IBComputeNormalGradient()\n");
    fprintf(stderr," grad_var is not NULL on rank %d\n",
            mpi->rank );
    return 1;
  }

  int nfacets_local = IB->nfacets_local;
  FacetMap *fmap = IB->fmap;

  if (nfacets_local > 0) {
    (*grad_var) = (double*) calloc (nvars*nfacets_local, sizeof(double));

    for (int n = 0; n < nfacets_local; n++) {
  
      double *alpha;
      int    *nodes, j, k;
  
      double v_c[nvars];
      alpha = &(fmap[n].interp_coeffs[0]);
      nodes = &(fmap[n].interp_nodes[0]);
      _ArraySetValue_(v_c,nvars,0.0);
      for (j=0; j<_IB_NNODES_; j++) {
        for (k=0; k<nvars; k++) {
          v_c[k] += ( alpha[j] * var[nvars*nodes[j]+k] );
        }
      }
  
      double v_ns[nvars];
      alpha = &(fmap[n].interp_coeffs_ns[0]);
      nodes = &(fmap[n].interp_nodes_ns[0]);
      _ArraySetValue_(v_ns,nvars,0.0);
      for (j=0; j<_IB_NNODES_; j++) {
        for (k=0; k<nvars; k++) {
          v_ns[k] += ( alpha[j] * var[nvars*nodes[j]+k] );
        }
      }
  
     double nx = fmap[n].facet->nx;
     double ny = fmap[n].facet->ny;
     double nz = fmap[n].facet->nz;
  
      for (k=0; k<nvars; k++) {
        double dv_dx = (v_ns[k] - v_c[k]) / fmap[n].dx;
        double dv_dy = (v_ns[k] - v_c[k]) / fmap[n].dy;
        double dv_dz = (v_ns[k] - v_c[k]) / fmap[n].dz;
        (*grad_var)[n*nvars+k] = dv_dx*nx + dv_dy*ny + dv_dz*nz;
      }
  
    }

  }

  return(0);
}
