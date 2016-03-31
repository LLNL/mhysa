/*
  This code extracts the midplane-slice of a specified 
  dimension in a nD solution (binary format) and writes 
  out a (n-1)D solution for that slice (in binary format).

  NOTE: make sure to make a subdirectory called "slices"
  before using it.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _ArraySetValue_(x,size,value)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);                       \
  }

#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

#define _ArrayIncrementIndex_(N,imax,i,done) \
  { \
    int arraycounter = 0; \
    while (arraycounter < (N)) { \
      if (i[arraycounter] == imax[arraycounter]-1) { \
        i[arraycounter] = 0; \
        arraycounter++; \
      } else { \
        i[arraycounter]++; \
        break; \
      } \
    } \
    if (arraycounter == (N)) done = 1; \
    else          done = 0; \
  }

void IncrementFilename(char *f)
{
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
}

int WriteBinary(int ndims, int nvars, int *dim, double *x, double *u, char *f)
{
  int size, d, index[ndims];
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

int ExtractSlice(char *f, int slicedim)
{
  FILE *in;
  in = fopen(f,"rb");
  if (!in) return(1);
  printf("Reading file %s.\n",f);
  int ndims, nvars;
  double *U,*x;

  /* read the file headers */
  fread(&ndims,sizeof(int),1,in);
  fread(&nvars,sizeof(int),1,in);

  /* read dimensions */
  int dims[ndims], d;
  fread(dims,sizeof(int),ndims,in);
  printf("Dimensions: ");
  for (d=0; d<ndims; d++) printf("%d ",dims[d]);
  printf("\n");
  printf("Nvars     : %d\n",nvars);
  
  /* allocate grid and solution arrays */
  int sizex = 0;      for (d=0; d<ndims; d++) sizex += dims[d];
  int sizeu = nvars;  for (d=0; d<ndims; d++) sizeu *= dims[d];
  x = (double*) calloc (sizex,sizeof(double));
  U = (double*) calloc (sizeu,sizeof(double));

  /* read grid and solution */
  fread(x,sizeof(double),sizex,in);
  fread(U,sizeof(double),sizeu,in);
  /* done reading */
  fclose(in);

  /* extract slice */
  int sdims[ndims-1], dd=0;
  sdims[0] = sdims[1] = 0;
  for (d=0; d<ndims; d++) {
    if (d != slicedim) {
      sdims[dd] = dims[d];
      dd++;
    }
  }
  printf("Slice dimensions: ");
  for (d=0; d<ndims-1; d++) printf("%d ",sdims[d]);
  printf("\n");
  int slice_sizex = 0;      for (d=0; d<ndims-1;d++) slice_sizex += sdims[d];
  int slice_sizeu = nvars;  for (d=0; d<ndims-1;d++) slice_sizeu *= sdims[d];
  double *SU, *Sx;
  Sx = (double*) calloc (slice_sizex,sizeof(double));
  SU = (double*) calloc (slice_sizeu,sizeof(double));

  int offset = 0, soffset = 0;
  for (d=0; d<ndims; d++) {
    if (d != slicedim) {
      int i;
      for (i = 0; i < dims[d]; i++) Sx[soffset+i] = x[offset+i];
      soffset += dims[d];
    }
    offset += dims[d];
  }

  int index[ndims], done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p;
    _ArrayIndex1D_(ndims,dims,index,0,p);
    if (index[slicedim] == dims[slicedim]/2) {
      int index_slice[ndims-1], pslice;
      int d, dd = 0;
      for (d=0; d<ndims; d++) {
        if (d != slicedim) {
          index_slice[dd] = index[d];
          dd++;
        }
      }
      _ArrayIndex1D_((ndims-1),sdims,index_slice,0,pslice);
      int v;
      for (v = 0; v < nvars; v++) SU[pslice*nvars+v] = U[p*nvars+v];
    }
    _ArrayIncrementIndex_(ndims,dims,index,done);
  }

  /* set filename */
  char slicefilename[100] = "slices/";
  strcat(slicefilename,f);
  
  /* write slice file */
  printf("Writing slice file %s.\n",slicefilename);
  WriteBinary(ndims-1,nvars,sdims,Sx,SU,slicefilename);
  
  /* clean up */
  free(U);
  free(x);
  free(SU);
  free(Sx);
  return(0);
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  double dt;
  int file_op_iter;
  char filename[50], op_file_format[50], overwrite[50];

  inputs = fopen("solver.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return(1);
  } else {
	  char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    fscanf(inputs,"%s",word);
   			if      (!strcmp(word, "dt"               ))  fscanf(inputs,"%lf",&dt           );
   			else if (!strcmp(word, "op_file_format"   ))  fscanf(inputs,"%s" ,op_file_format);
   			else if (!strcmp(word, "file_op_iter"     ))  fscanf(inputs,"%d" ,&file_op_iter  );
   			else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,overwrite      );
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  int slicedim;
  printf("Enter dimension along which to extract slice: ");
  scanf("%d",&slicedim);

  if (!strcmp(overwrite,"no")) {

    strcpy(filename,"op_00000.bin");
    while(1) {
      int flag = ExtractSlice(filename,slicedim);
      if (flag) { 
        printf("No more files found. Exiting.\n");
        break;
      } 
      IncrementFilename(filename);
    }
  } else if (!strcmp(overwrite,"yes")) {
    strcpy(filename,"op.bin");
    int flag = ExtractSlice(filename,slicedim);
    if (flag) printf("File not found. Exiting.\n");
  }

  return(0);
}
