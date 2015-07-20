/*! @file WriteTecplot3D.c
    @author Debojyoti Ghosh
    @brief Write a vector field and its grid to a 3D Tecplot file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

/*! Write a vector field and its associated grid to a 3D Tecplot file. This 
    file can then be visualized using Tecplot (http://www.tecplot.org) 
    (if available) \n
    \n
    Note: It's essentially a text file, and apart from the first two lines
    with Tecplot specific headers, the data is written out in the same 
    format as WriteText().
*/
int WriteTecplot3D(
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
  if (ndims !=3) {
    fprintf(stderr,"Error in WriteTecplot3D(): This functions is hardcoded for 3-dimensional ");
    fprintf(stderr,"problems only. Instead, ndims=%d.\n",ndims);
    return(1);
  }
  int i;
  int imax = dim[0];
  int jmax = dim[1];
  int kmax = dim[2];

  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* writing tecplot data file headers */
  fprintf(out,"VARIABLES=\"I\",\"J\",\"K\",\"X\",\"Y\",\"Z\",");
  char varname[3] = "00";
  for (i = 0; i < nvars; i++) {
    fprintf(out,"\"%s\",",varname);
    if (varname[1] == '9') { varname[0]++; varname[1] = '0'; }
    else                     varname[1]++;
  }
  fprintf(out,"\n");
  fprintf(out,"ZONE I=%d,J=%d,K=%d,F=POINT\n",imax,jmax,kmax);

  /* writing the data */
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
