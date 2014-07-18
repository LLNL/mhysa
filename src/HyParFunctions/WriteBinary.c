#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

int WriteBinary(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
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
