/*
  This code converts the binary solution files for a 2D
  system to tecplot files
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

void WriteTecplot2D(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  if (ndims !=2) {
    fprintf(stderr,"Error in WriteTecplot3D(): This functions is hardcoded for 2-dimensional ");
    fprintf(stderr,"problems only. Instead, ndims=%d.\n",ndims);
    return;
  }
  int i;
  int imax = dim[0];
  int jmax = dim[1];

  printf("\tWriting tecplot solution file %s.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return;
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
  return;
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  double dt;
  int file_op_iter;
  char filename[50], op_file_format[50], tecfile[50];

  inputs = fopen("solver.inp","r");
  if (!in) {
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
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  strcpy(filename,"op_00000.bin");
  while(1) {
    in = fopen(filename,"rb");

    if (!in) {
      printf("No more files found. Exiting.\n");
      break;
    } else {
      printf("Reading file %s.\n",filename);
      int ndims, nvars, dims[2];
      double *U,*x;

      /* read the file headers */
      fread(&ndims,sizeof(int),1,in);
      fread(&nvars,sizeof(int),1,in);

      /* some checks */
      if (ndims != 2) {
        printf("Error: ndims in %s not equal to 2!\n",filename);
        return(0);
      }

      /* read dimensions */
      fread(dims,sizeof(int),ndims,in);
      printf("Dimensions: %d x %d\n",dims[0],dims[1]);
      printf("Nvars     : %d\n",nvars);

      /* allocate grid and solution arrays */
      int sizex = dims[0] + dims[1];
      int sizeu = dims[0] * dims[1] * nvars;
      x = (double*) calloc (sizex,sizeof(double));
      U = (double*) calloc (sizeu,sizeof(double));

      /* read grid and solution */
      fread(x,sizeof(double),sizex,in);
      fread(U,sizeof(double),sizeu,in);
      /* done reading */
      fclose(in);

      /* set filename */
      strcpy(tecfile,filename);
      tecfile[9]  = 'd';
      tecfile[10] = 'a';
      tecfile[11] = 't';

      /* write Tecplot file */
      int ind[2];
      printf("\tWriting file %s.\n",tecfile);
      WriteTecplot2D(2,nvars,dims,x,U,tecfile,&ind[0]);

      /* clean up */
      free(U);
      free(x);
    }

    IncrementFilename(filename);
  }

  return(0);
}
