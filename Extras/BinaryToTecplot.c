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
    fprintf(stderr,"Error in WriteTecplot2D(): This functions is hardcoded for 2-dimensional ");
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

int WriteTecplot3D(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
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

  printf("Writing tecplot solution file %s.\n",f);
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
      fprintf(out,"%+E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return(0);
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  double dt;
  int file_op_iter;
  char filename[50], op_file_format[50], tecfile[50], overwrite[50];

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

  if (!strcmp(overwrite,"no")) {
    strcpy(filename,"op_00000.bin");
    while(1) {
      in = fopen(filename,"rb");

      if (!in) {
        printf("No more files found. Exiting.\n");
        break;
      } else {
        printf("Reading file %s.\n",filename);
        int ndims, nvars;
        double *U,*x;

        /* read the file headers */
        fread(&ndims,sizeof(int),1,in);
        fread(&nvars,sizeof(int),1,in);

        /* some checks */
        if ((ndims != 2) && (ndims != 3)) {
          printf("Error: ndims in %s not equal to 2 or 3!\n",filename);
          return(0);
        }

        /* read dimensions */
        int dims[ndims];
        fread(dims,sizeof(int),ndims,in);
        if      (ndims == 2) printf("Dimensions: %d x %d\n",dims[0],dims[1]);
        else if (ndims == 3) printf("Dimensions: %d x %d x %d\n",dims[0],dims[1],dims[2]);
        printf("Nvars     : %d\n",nvars);

        /* allocate grid and solution arrays */
        int d;
        int sizex = 0;      for (d=0; d<ndims; d++) sizex += dims[d];
        int sizeu = nvars;  for (d=0; d<ndims; d++) sizeu *= dims[d];
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
        int ind[ndims];
        if      (ndims == 2) WriteTecplot2D(2,nvars,dims,x,U,tecfile,&ind[0]);
        else if (ndims == 3) WriteTecplot3D(3,nvars,dims,x,U,tecfile,&ind[0]);

        /* clean up */
        free(U);
        free(x);
      }

      IncrementFilename(filename);
    }
  } else if (!strcmp(overwrite,"yes")) {
    strcpy(filename,"op.bin");
    in = fopen(filename,"rb");
    if (!in) {
      printf("File no found. Exiting.\n");
    } else {
      printf("Reading file %s.\n",filename);
      int ndims, nvars;
      double *U,*x;

      /* read the file headers */
      fread(&ndims,sizeof(int),1,in);
      fread(&nvars,sizeof(int),1,in);

      /* some checks */
      if ((ndims != 2) && (ndims != 3)) {
        printf("Error: ndims in %s not equal to 2 or 3!\n",filename);
        return(0);
      }

      /* read dimensions */
      int dims[ndims];
      fread(dims,sizeof(int),ndims,in);
      if      (ndims == 2) printf("Dimensions: %d x %d\n",dims[0],dims[1]);
      else if (ndims == 3) printf("Dimensions: %d x %d x %d\n",dims[0],dims[1],dims[2]);
      printf("Nvars     : %d\n",nvars);

      /* allocate grid and solution arrays */
      int d;
      int sizex = 0;      for (d=0; d<ndims; d++) sizex += dims[d];
      int sizeu = nvars;  for (d=0; d<ndims; d++) sizeu *= dims[d];
      x = (double*) calloc (sizex,sizeof(double));
      U = (double*) calloc (sizeu,sizeof(double));

      /* read grid and solution */
      fread(x,sizeof(double),sizex,in);
      fread(U,sizeof(double),sizeu,in);
      /* done reading */
      fclose(in);

      /* set filename */
      strcpy(tecfile,filename);
      tecfile[3] = 'd';
      tecfile[4] = 'a';
      tecfile[5] = 't';

      /* write Tecplot file */
      int ind[ndims];
      if      (ndims == 2) WriteTecplot2D(2,nvars,dims,x,U,tecfile,&ind[0]);
      else if (ndims == 3) WriteTecplot3D(3,nvars,dims,x,U,tecfile,&ind[0]);

      /* clean up */
      free(U);
      free(x);
    }
  }

  return(0);
}
