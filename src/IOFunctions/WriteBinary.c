/*! @file WriteBinary.c
    @author Debojyoti Ghosh
    @brief Write a vector field and its grid to a binary file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

/*! Write a vector field and its grid to a binary file. The data is written
    in the following format:\n\n
    {\n
      ndims\n
      nvars\n
      dim[0] dim[1] dim[2] ... dim[ndims-1]\n
      x0_i (0 <= i < dim[0])\n
      x1_i (0 <= i < dim[1])\n
      ...\n
      x{ndims-1}_i (0 <= i < dim[ndims-1])\n
      [u0,u1,...,u{nvars-1}]_p (0 <= p < N) (with no commas)\n
    }\n\n
    where \n
    x0, x1, ..., x{ndims-1} represent the spatial dimensions (for a 3D problem, x0 = x, x1 = y, x2 = z),\n
    u0, u1, ..., u{nvars-1} are each component of the vector u at a grid point,\n
    N = dim[0]*dim[1]*...*dim[ndims-1] is the total number of points,\n
    and p = i0 + dim[0]*( i1 + dim[1]*( i2 + dim[2]*( ... + dim[ndims-2] * i{ndims-1} ))) (see #_ArrayIndexnD_)\n
    with i0, i1, i2, etc representing grid indices along each spatial dimension, i.e.,\n
    0 <= i0 < dim[0]-1\n
    0 <= i1 < dim[1]-1\n
    ...\n
    0 <= i{ndims-1} < dim[ndims=1]-1\n
*/
int WriteBinary(
                  int ndims,  /*!< Number of spatial dimensions */
                  int nvars,  /*!< Number of variables at each grid point */
                  int *dim,   /*!< Integer array with the number of grid points 
                                   in each spatial dimension as its entries */
                  double *x,  /*!< Array of spatial coordinates representing a
                                   Cartesian grid (no ghost points) */
                  double *u,  /*!< Array containing the vector field to write 
                                   (no ghost points) */
                  char *f,    /*!< Filename */
                  int *index  /*!< Preallocated integer array of size ndims */
               )
{
  int size, d;
  size_t bytes;
  FILE *out;
  out = fopen(f,"wb");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* write ndims, nvars */
  bytes = fwrite(&ndims,sizeof(int),1,out);
  if ((int)bytes != 1) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write ndims to output file.\n");
  }
  bytes = fwrite(&nvars,sizeof(int),1,out);
  if ((int)bytes != 1) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write nvars to output file.\n");
  }

  /* write dimensions */
  bytes = fwrite(dim,sizeof(int),ndims,out);
  if ((int)bytes != ndims) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write dimensions to output file.\n");
  }

  /* write grid */
  size = 0;
  for (d = 0; d < ndims; d++) size += dim[d];
  bytes = fwrite(x,sizeof(double),size,out);
  if ((int)bytes != size) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write grid to output file.\n");
  }

  /* write solution */
  size = 1;
  for (d = 0; d < ndims; d++) size *= dim[d]; size *= nvars;
  bytes = fwrite(u,sizeof(double),size,out);
  if ((int)bytes != size) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write solution to output file.\n");
  }

  fclose(out);
  return(0);
}
