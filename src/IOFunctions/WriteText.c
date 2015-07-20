/*! @file WriteText.c
    @author Debojyoti Ghosh
    @brief Write a vector field and its grid to a text file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

/*! Write a vector field and its associated grid to a text file. The data 
    is written to the text file in the following format:\n\n
    Each line of the file contains\n
    i0 i1 ... i{ndims-1} x0[i0] x1[i1] ... x{ndims-1}[i{ndims-1}] u0[p] u1[p] ... u{nvars-1}[p]\n\n
    where \n
    (i0,i1,...,i{ndims-1}) represents a ndims-dimensional grid index,\n
    x0, x1, ..., x{ndims-1} are the spatial dimensions (e.g. in 3D, ndims = 3, x0 = x, x1 = y, x2 = z)\n
    u0, u1, ..., u{nvars-1} are the components of the vector field u\n
    p = i0 + dim[0] * (i1 + dim[1] * (i2 + dim[2] * ( ... + dim[ndims-2] * i{ndims-1} ))) 
    (see #_ArrayIndexnD_)\n
    Thus, the number of lines is dim[0]*dim[1]*...*dim[ndims-1]
    \n\n
    This is a plain text file, so it can be visualized using, say MATLAB, or any software
    that can read data from a plain text file. For 1D (ndims = 1), the data can be easily
    plotted in Gnuplot (for example).
*/
int WriteText(
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
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int i, p;
    _ArrayIndex1D_(ndims,dim,index,0,p);
    for (i=0; i<ndims; i++) fprintf(out,"%4d ",index[i]);
    for (i=0; i<ndims; i++) {
      int j,offset = 0; for (j=0; j<i; j++) offset += dim[j];
      fprintf(out,"%+1.16E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+1.16E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return(0);
}
