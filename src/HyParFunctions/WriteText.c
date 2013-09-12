#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

int WriteText(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  int ierr = 0;
  printf("Writing solution file %s in ASCII text format.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  int done = 0;
  ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int i;
    int p = ArrayIndex1D(ndims,dim,index,NULL,0);
    for (i=0; i<ndims; i++) fprintf(out,"%4d ",index[i]);
    for (i=0; i<ndims; i++) {
      int j,offset = 0; for (j=0; j<i; j++) offset += dim[j];
      fprintf(out,"%E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%E ",u[nvars*p+i]);
    fprintf(out,"\n");
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  fclose(out);
  return(0);
}
