/*! @file FindInterval.c
    @author Debojyoti Ghosh
    @brief Find grid point indices corresponding to a spatial interval.
*/

#include <mathfunctions.h>

/*! Given an interval \f$\left[a,b\right], a\leq b\f$, find grid indices \a imin 
    and \a imax, such that
    \f{align}{
      imin &= \min\ i\ {\rm satisfying}\ x_i \geq a\\
      imax &= \max\ i\ {\rm satisfying}\  x_i \leq b
    \f}
    where \f$\left\{x_i; 0\leq i < N , x_i < x_{i+1} \forall i \right\}\f$ 
    represents a 1-dimensional grid.
    \n\n
    Note: This function handles 1-dimensional intervals and grids only.
*/
void FindInterval(
                  double  a,      /*!< Lower bound of interval */
                  double  b,      /*!< Upper bound of interval */
                  double  *x,     /*!< Array of spatial coordinates representing a grid */
                  int     N,      /*!< Number of grid points / size of x */
                  int     *imin,  /*!< Lowest grid index within [a,b] */
                  int     *imax   /*!< Highest grid index within [a,b] */
                 )
{
  int i;
  *imax = -1;
  *imin =  N;

  double min_dx = x[1] - x[0];
  for (i = 2; i < N; i++) {
    double dx = x[i] - x[i-1];
    if (dx < min_dx) min_dx = dx;
  }
  double tol = 1e-10 * min_dx;

  for (i = 0; i < N; i++) {
    if (x[i] <= (b+tol)) *imax = i+1;
  }
  for (i = N-1; i > -1; i--) {
    if (x[i] >= (a-tol)) *imin = i;
  }

  return;
}
