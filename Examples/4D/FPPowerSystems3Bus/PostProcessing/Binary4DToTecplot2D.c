/*
  This code converts the 4D binary file for the 3-Bus
  power system simulations to two 2D tecplot files, one
  for each generator.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
  char varname[3] = "P";
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
  char filename[50], op_file_format[50], tecfile1[50], tecfile2[50];

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

  int count = 0;
  double samay = 0.0;
  strcpy(filename,"op_00000.bin");
  double energy0 = 1;
  while(1) {
    in = fopen(filename,"rb");

    if (!in) {
      printf("No more files found. Exiting.\n");
      break;
    } else {
      printf("Processing file %s.\n",filename);
      int ndims, nvars, dims[4];
      double *U,*x;

      /* read the file headers */
      fread(&ndims,sizeof(int),1,in);
      fread(&nvars,sizeof(int),1,in);

      /* some checks */
      if (ndims != 4) {
        printf("Error: ndims in %s not equal to 3!\n",filename);
        return(0);
      }
      if (nvars != 1) {
        printf("Error: nvars in %s not equal to 4!\n",filename);
        return(0);
      }

      /* read dimensions */
      fread(dims,sizeof(int),ndims,in);
      printf("Dimensions: %d x %d x %d x %d\n",dims[0],dims[1],dims[2],dims[3]);

      /* allocate grid and solution arrays */
      int sizex = dims[0] + dims[1] + dims[2] + dims[3];
      int sizeu = dims[0] * dims[1] * dims[2] * dims[3];
      x = (double*) calloc (sizex,sizeof(double));
      U = (double*) calloc (sizeu,sizeof(double));

      /* read grid and solution */
      fread(x,sizeof(double),sizex,in);
      fread(U,sizeof(double),sizeu,in);
      /* done reading */
      fclose(in);

      /* set filenames */
      strcpy(tecfile1,filename); tecfile1[0]='g'; tecfile1[1]='1';
      strcpy(tecfile2,filename); tecfile2[0]='g'; tecfile2[1]='2';
      tecfile1[9] = tecfile2[9] = 'd';
      tecfile1[10] = tecfile2[10] = 'a';
      tecfile1[11] = tecfile2[11] = 't';

      /* arrays for the 2D solution */
      double *x2d, *U2d;
      int    dim2d[2],NI,NJ,i,j,off1,off2,ind[2];

      /* Generator 1 */
      printf("\tExtracting data for Generator 1.\n");
      NI = dim2d[0] = dims[0]; NJ = dim2d[1] = dims[2];
      x2d = (double*) calloc (NI+NJ,sizeof(double));
      U2d = (double*) calloc (NI*NJ,sizeof(double));

      /* copy grid */
      for (i=0; i<NI; i++) x2d[i] = x[i];
      off1 = NI; off2 = dims[0]+dims[1];
      for (i=0; i<NJ; i++) x2d[i+off1] = x[i+off2];

      /* copy U - integrate along the other two dimensions */
      for (i=0; i<NI; i++) {
        for (j=0; j<NJ; j++) {
          int p = i + NI*j;
          U2d[p] = 0;
          int k,l,n;
          for (k=0; k<dims[1]; k++) {
            for (l=0; l<dims[3]; l++) {
              int q = i + dims[0] * (k + dims[1] * (j + dims[2] * l));
              double dx[4];
              int offset = 0, index[4];
              index[0] = i; index[1] = k; index[2] = j; index[3] = l;
              for (n=0; n<4; n++) {
                if      (index[n] == 0)         dx[n] = (x[index[n]+offset+1]-x[index[n]+offset]  );
                else if (index[n] == dims[n]-1) dx[n] = (x[index[n]+offset]  -x[index[n]+offset-1]);
                else                            dx[n] = (x[index[n]+offset+1]-x[index[n]+offset-1])/2.0;
                offset += dims[n];
              }
              U2d[p] += U[q] * dx[1]*dx[3];
            }
          }
        }
      }

      /* write to file */
      printf("\tWriting file %s.\n",tecfile1);
      WriteTecplot2D(2,1,dim2d,x2d,U2d,tecfile1,&ind[0]);

      /* Generator 2 */
      printf("\tExtracting data for Generator 2.\n");
      NI = dim2d[0] = dims[1]; NJ = dim2d[1] = dims[3];
      x2d = (double*) calloc (NI+NJ,sizeof(double));
      U2d = (double*) calloc (NI*NJ,sizeof(double));

      /* copy grid */
      for (i=0; i<NI; i++) x2d[i] = x[i+dims[0]];
      off1 = NI; off2 = dims[0]+dims[1]+dims[2];
      for (i=0; i<NJ; i++) x2d[i+off1] = x[i+off2];

      /* copy U - integrate along the other two dimensions */
      for (i=0; i<NI; i++) {
        for (j=0; j<NJ; j++) {
          int p = i + NI*j;
          U2d[p] = 0;
          int k,l,n;
          for (k=0; k<dims[0]; k++) {
            for (l=0; l<dims[2]; l++) {
              int q = k + dims[0] * (i + dims[1] * (l + dims[2] * j));
              double dx[4];
              int offset = 0, index[4];
              index[0] = k; index[1] = i; index[2] = l; index[3] = j;
              for (n=0; n<4; n++) {
                if      (index[n] == 0)         dx[n] = (x[index[n]+offset+1]-x[index[n]+offset]  );
                else if (index[n] == dims[n]-1) dx[n] = (x[index[n]+offset]  -x[index[n]+offset-1]);
                else                            dx[n] = (x[index[n]+offset+1]-x[index[n]+offset-1])/2.0;
                offset += dims[n];
              }
              U2d[p] += U[q] * dx[0]*dx[2];
            }
          }
        }
      }

      /* write to file */
      WriteTecplot2D(2,1,dim2d,x2d,U2d,tecfile2,&ind[0]);

      /* done */
      free(x2d);
      free(U2d);

      /* clean up */
      free(U);
      free(x);
    }

    count++;
    samay += (dt * (double)file_op_iter);
    IncrementFilename(filename);
  }

  return(0);
}
