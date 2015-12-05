/*! @file WENOFifthOrderInitializeWeights.c
    @brief Initialize the arrays for nonlinear weights of WENO-type schemes
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
  This function initializes the arrays to store the nonlinear weights for fifth order WENO-type schemes
  to their optimal values (i.e. values of the weights for a perfectly smooth solution). 5th order methods
  need three weights at each grid interface.

  \b Notes:
  + The three weights at all the grid interfaces along all the dimensions are stored in three big arrays
    (#WENOParameters::w1, #WENOParameters::w2, #WENOParameters::w3), one for each weight. 
  + This function initializes the values for the spatial dimension given by the input \a dir.
  + \a offset indicates the location in #WENOParameters::w1,#WENOParameters::w2, and #WENOParameters::w3
    from where the chunk of memory with the weights for the spatial dimension \a dir starts.
  + This chunk of memory (1D array) represents a multidimensional array of the following size: 
    4 X (number of interfaces) X (number of solution components #HyPar::nvars). The factor 4 comes from the
    fact that the array stores weights for the interpolation of left-biased flux function, left-biased 
    solution, right-biased flux function, and right-biased solution, one after the other, in this order.
  + The weights are initialized to their optimal values so that, if no limiting is specified 
    (#WENOParameters::no_limiting = 1), then no further computation of weights are required.
*/
int WENOFifthOrderInitializeWeights(
                                      int   dir,  /*!< Spatial dimension */
                                      void  *s,   /*!< Solver object of type #HyPar */
                                      void  *m    /*!< MPI object of type #MPIVariables */
                                   )
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  int             done;
  double          *ww1, *ww2, *ww3;


  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* calculate dimension offset */
  int offset = weno->offset[dir];

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  /* calculate weights for a left-biased interpolation */
  ww1 = weno->w1 + offset;
  ww2 = weno->w2 + offset;
  ww3 = weno->w3 + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p, v;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      for (v=0; v<nvars; v++)  {
        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* save the weights */
        *(ww1+p*nvars+v) = c1;
        *(ww2+p*nvars+v) = c2;
        *(ww3+p*nvars+v) = c3;
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  ww1 = weno->w1 + weno->size + offset;
  ww2 = weno->w2 + weno->size + offset;
  ww3 = weno->w3 + weno->size + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p, v;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      for (v=0; v<nvars; v++)  {

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* save the weights */
        *(ww1+p*nvars+v) = c1;
        *(ww2+p*nvars+v) = c2;
        *(ww3+p*nvars+v) = c3;
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  /* calculate weights for a right-biased interpolation */
  ww1 = weno->w1 + 2*weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p, v;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      for (v=0; v<nvars; v++)  {

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* save the weights */
        *(ww1+p*nvars+v) = c1;
        *(ww2+p*nvars+v) = c2;
        *(ww3+p*nvars+v) = c3;
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  ww1 = weno->w1 + 2*weno->size + weno->size + offset;
  ww2 = weno->w2 + 2*weno->size + weno->size + offset;
  ww3 = weno->w3 + 2*weno->size + weno->size + offset;
  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int p, v;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      for (v=0; v<nvars; v++)  {

        /* optimal weights*/
        double c1, c2, c3;
        if (!strcmp(solver->spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_)) {
          if (   ((mpi->ip[dir] == 0                ) && (indexI[dir] == 0       ))
              || ((mpi->ip[dir] == mpi->iproc[dir]-1) && (indexI[dir] == dim[dir])) ) {
            /* Use WENO5 at the physical boundaries */
            c1 = _WENO_OPTIMAL_WEIGHT_1_;
            c2 = _WENO_OPTIMAL_WEIGHT_2_;
            c3 = _WENO_OPTIMAL_WEIGHT_3_;
          } else {
            /* CRWENO5 at the interior points */
            c1 = _CRWENO_OPTIMAL_WEIGHT_1_;
            c2 = _CRWENO_OPTIMAL_WEIGHT_2_;
            c3 = _CRWENO_OPTIMAL_WEIGHT_3_;
          }
        } else {
          /* WENO5 and HCWENO5 */
          c1 = _WENO_OPTIMAL_WEIGHT_1_;
          c2 = _WENO_OPTIMAL_WEIGHT_2_;
          c3 = _WENO_OPTIMAL_WEIGHT_3_;
        }

        /* save the weights */
        *(ww1+p*nvars+v) = c1;
        *(ww2+p*nvars+v) = c2;
        *(ww3+p*nvars+v) = c3;
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
