/*! @file TrilinearInterpolation.c
    @brief Compute coefficients for trilinear interpolation
    @author Debojyoti Ghosh
*/

#include <mathfunctions.h>

/*!
  This function computes the coefficients for a trilinear interpolation at a given
  point (x,y,z) inside a cube defined by [xmin,xmax] X [ymin,ymax] X [zmin,zmax].
  The coefficients are stored in an array of size 8 with each element corresponding
  to a corner of the cube in the following order:\n
  coeffs[0] => xmin,ymin,zmin\n
  coeffs[1] => xmax,ymin,zmin\n
  coeffs[2] => xmin,ymax,zmin\n
  coeffs[3] => xmax,ymax,zmin\n
  coeffs[4] => xmin,ymin,zmax\n
  coeffs[5] => xmax,ymin,zmax\n
  coeffs[6] => xmin,ymax,zmax\n
  coeffs[7] => xmax,ymax,zmax
*/
void TrilinearInterpCoeffs(
                            double xmin,  /*!< x-coordinate of the lower-end */
                            double xmax,  /*!< x-coordinate of the higher-end */  
                            double ymin,  /*!< y-coordinate of the lower-end */
                            double ymax,  /*!< y-coordinate of the higher-end */
                            double zmin,  /*!< z-coordinate of the lower-end */
                            double zmax,  /*!< z-coordinate of the higher-end */
                            double x,     /*!< x-coordinate of the point to interpolate at */
                            double y,     /*!< y-coordinate of the point to interpolate at */
                            double z,     /*!< z-coordinate of the point to interpolate at */
                            double *coeffs/*!< array of size 8 (pre-allocated) to store the coefficients in */
                          )
{
	double vol_inv = 1 / ((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
	double tldx1 = x - xmin;
	double tldx2 = xmax - x;
	double tldy1 = y - ymin;
	double tldy2 = ymax - y;
	double tldz1 = z - zmin;
	double tldz2 = zmax - z;

	coeffs[0] = tldz2 * tldy2 * tldx2 * vol_inv;
	coeffs[1] = tldz2 * tldy2 * tldx1 * vol_inv;
	coeffs[2] = tldz2 * tldy1 * tldx2 * vol_inv;
	coeffs[3] = tldz2 * tldy1 * tldx1 * vol_inv;
	coeffs[4] = tldz1 * tldy2 * tldx2 * vol_inv;
	coeffs[5] = tldz1 * tldy2 * tldx1 * vol_inv;
	coeffs[6] = tldz1 * tldy1 * tldx2 * vol_inv;
	coeffs[7] = tldz1 * tldy1 * tldx1 * vol_inv;

	return;
}
