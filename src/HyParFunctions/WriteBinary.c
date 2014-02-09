#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

int WriteBinary(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  int size, d;
  printf("Writing solution file %s in binary format.\n",f);
  FILE *out;
  out = fopen(f,"wb");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* write ndims, nvars */
  fwrite(&ndims,sizeof(int),1,out);
  fwrite(&nvars,sizeof(int),1,out);

  /* write dimensions */
  fwrite(dim,sizeof(int),ndims,out);

  /* write grid */
  size = 0;
  for (d = 0; d < ndims; d++) size += dim[d];
  fwrite(x,sizeof(double),size,out);

  /* write solution */
  size = 1;
  for (d = 0; d < ndims; d++) size *= dim[d]; size *= nvars;
  fwrite(u,sizeof(double),size,out);

  fclose(out);
  return(0);
}
