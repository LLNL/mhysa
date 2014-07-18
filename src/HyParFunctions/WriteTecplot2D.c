#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>

int WriteTecplot2D(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  if (ndims !=2) {
    fprintf(stderr,"Error in WriteTecplot3D(): This functions is hardcoded for 2-dimensional ");
    fprintf(stderr,"problems only. Instead, ndims=%d.\n",ndims);
    return(1);
  }
  int i;
  int imax = dim[0];
  int jmax = dim[1];

  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* writing tecplot data file headers */
  fprintf(out,"VARIABLES=\"I\",\"J\",\"X\",\"Y\",");
  char varname[3] = "00";
  for (i = 0; i < nvars; i++) {
    fprintf(out,"\"%s\",",varname);
    if (varname[1] == '9') { varname[0]++; varname[1] = '0'; }
    else                     varname[1]++;
  }
  fprintf(out,"\n");
  fprintf(out,"ZONE I=%d,J=%d,F=POINT\n",imax,jmax);

  /* writing the data */
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
