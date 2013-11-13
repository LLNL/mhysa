#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

int WriteText(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  printf("Writing solution file %s in ASCII text format.\n",f);
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
      fprintf(out,"%+E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return(0);
}
