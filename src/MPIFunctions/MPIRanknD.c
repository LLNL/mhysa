/*! @file MPIRanknD.c
    @brief Compute the n-dimensional rank
    @author Debojyoti Ghosh
*/
#include <mpivars.h>

/*!
  This function computes the n-dimensional rank, given the 1D rank and the total
  number of MPI ranks along each spatial dimension.

  \b 1D \b Rank: This is the rank of the process in the communicator.\n
  \b n-Dimensional \b Rank: This represents an integer array, where each element 
  is the rank of the process along a spatial dimension.

  Consider a 2D simulation running with 21 MPI ranks - 7 along the \a x direction, 
  and 3 along the \a y direction, as shown in the following figure:
  @image html nd_ranks.png
  @image latex nd_ranks.eps width=0.9\textwidth
  The boldface number in the parentheses is the n-dimensional rank (n=2), while the 
  number below it in normal typeface is the 1D rank, corresponding to the rank in 
  the MPI communicator.

  Note:
  + the array to hold the computed n-dimensional rank must be preallocated and the size
    must be equal to the number of spatial dimensions.
*/
int MPIRanknD(
                int ndims,    /*!< Number of spatial dimensions */
                int rank,     /*!< 1D rank*/
                int *iproc,   /*!< Integer array whose elements are the number of MPI ranks along each dimension */
                int *ip       /*!< Integer array whose elements are the rank of this process along each dimesion */
             )
{
  int i,term    = 1;
  for (i=0; i<ndims; i++) term *= iproc[i];
  for (i=ndims-1; i>=0; i--) {
    term /= iproc[i];
    ip[i] = rank/term;
    rank -= ip[i]*term;
  }
  return(0);
}
