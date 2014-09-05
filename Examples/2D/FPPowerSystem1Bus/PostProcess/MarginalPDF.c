/*

  This code extracts the marginal PDFs from the 2D 
  PDF output and writes them to text files.
  (Note: output file must be binary.)

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

int WriteText(int N,double *x,double *u,char *f)
{
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }
  int i;
  for (i=0; i<N; i++) fprintf(out, "%4d\t%1.16E\t%1.16E\n",i,x[i],u[i]);
  fclose(out);
  return(0);
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  char filename[50], op_file_format[50], tecfile1[50], tecfile2[50], op_overwrite[50];

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
   			if      (!strcmp(word, "op_file_format"   ))  fscanf(inputs,"%s" ,op_file_format);
   			else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,op_overwrite  );
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  if (!strcmp(op_overwrite,"no")) {
    strcpy(filename,"op_00000.bin");
    while(1) {
      in = fopen(filename,"rb");

      if (!in) {
        printf("No more files found. Exiting.\n");
        break;
      } else {
        printf("Processing file %s.\n",filename);
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
        if (nvars != 1) {
          printf("Error: nvars in %s not equal to 1!\n",filename);
          return(0);
        }

        /* read dimensions */
        fread(dims,sizeof(int),ndims,in);
        printf("Dimensions: %d x %d \n",dims[0],dims[1]);

        /* allocate grid and solution arrays */
        int sizex = dims[0] + dims[1];
        int sizeu = dims[0] * dims[1];
        x = (double*) calloc (sizex,sizeof(double));
        U = (double*) calloc (sizeu,sizeof(double));

        /* read grid and solution */
        fread(x,sizeof(double),sizex,in);
        fread(U,sizeof(double),sizeu,in);
        /* done reading */
        fclose(in);

        /* set filenames */
        strcpy(tecfile1,filename); tecfile1[0]='t'; tecfile1[1]='h';
        strcpy(tecfile2,filename); tecfile2[0]='o'; tecfile2[1]='m';
        tecfile1[9] = tecfile2[9] = 'd';
        tecfile1[10] = tecfile2[10] = 'a';
        tecfile1[11] = tecfile2[11] = 't';

        /* variables for the 1D solution */
        double *X, *Y, *U1d, integral;
        int    NI,NJ,i,j,ind;
        NI = dims[0]; NJ = dims[1];
        X = x;
        Y = x + NI;

        /* calculate integral of the 2D PDF */
        integral = 0.0;
        for (i=0; i<NI; i++) {
          for (j=0; j<NJ; j++) {
            int q = i + NI * j;
            double dx, dy;
            if      (i == 0)    dx = (X[i+1]-X[i]  );
            else if (i == NI-1) dx = (X[i]  -X[i-1]);
            else                dx = (X[i+1]-X[i-1])/2.0;
            if      (j == 0)    dy = (Y[j+1]-Y[j]  );
            else if (j == NJ-1) dy = (Y[j]  -Y[j-1]);
            else                dy = (Y[j+1]-Y[j-1])/2.0;
            integral += U[q] * dx * dy;
          }
        }
        printf("\tIntegral of 2D PDF: %1.16E\n",integral);

        /* Theta */
        printf("\tExtracting marginal PDF for theta.\n");
        U1d = (double*) calloc (NI,sizeof(double));

        /* copy U - integrate along the other dimension */
        for (i=0; i<NI; i++) {
          U1d[i] = 0;
          for (j=0; j<dims[1]; j++) {
            int q = i + NI * j;
            double dy;
            if      (j == 0)    dy = (Y[j+1]-Y[j]  );
            else if (j == NJ-1) dy = (Y[j]  -Y[j-1]);
            else                dy = (Y[j+1]-Y[j-1])/2.0;
            U1d[i] += U[q] * dy;
          }
        }

        /* calculating integral of marginal PDF */
        integral = 0;
        for (i=0; i<NI; i++) {
          double dx;
          if      (i == 0)    dx = (X[i+1]-X[i]  );
          else if (i == NI-1) dx = (X[i]  -X[i-1]);
          else                dx = (X[i+1]-X[i-1])/2.0;
          integral += U1d[i] * dx;
        }
        printf("\tIntegral of marginal PDF for theta: %1.16E\n",integral);

        /* write to file */
        printf("\tWriting file %s.\n",tecfile1);
        WriteText(NI,X,U1d,tecfile1);

        /* free arrays */
        free(U1d);

        /* Omega */
        printf("\tExtracting marginal PDF for omega.\n");
        U1d = (double*) calloc (NJ,sizeof(double));

        /* copy U - integrate along the other dimension */
        for (j=0; j<NJ; j++) {
          U1d[j] = 0;
          for (i=0; i<NI; i++) {
            int q = i + NI * j;
            double dx;
            if      (i == 0)    dx = (X[i+1]-X[i]  );
            else if (i == NI-1) dx = (X[i]  -X[i-1]);
            else                dx = (X[i+1]-X[i-1])/2.0;
            U1d[j] += U[q] * dx;
          }
        }

        /* calculating integral of marginal PDF */
        integral = 0;
        for (j=0; j<NJ; j++) {
          double dy;
          if      (j == 0)    dy = (Y[j+1]-Y[j]  );
          else if (j == NJ-1) dy = (Y[j]  -Y[j-1]);
          else                dy = (Y[j+1]-Y[j-1])/2.0;
          integral += U1d[j] * dy;
        }
        printf("\tIntegral of marginal PDF for omega: %1.16E\n",integral);

        /* write to file */
        printf("\tWriting file %s.\n",tecfile2);
        WriteText(NJ,Y,U1d,tecfile2);

        /* done */
        free(U1d);

        /* clean up */
        free(U);
        free(x);
      }

      IncrementFilename(filename);
    }
  } else {
    strcpy(filename,"op.bin");
    in = fopen(filename,"rb");

    if (!in) {
      printf("File op.bin not found. Exiting.\n");
    } else {
      printf("Processing file %s.\n",filename);
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
      if (nvars != 1) {
        printf("Error: nvars in %s not equal to 1!\n",filename);
        return(0);
      }

      /* read dimensions */
      fread(dims,sizeof(int),ndims,in);
      printf("Dimensions: %d x %d \n",dims[0],dims[1]);

      /* allocate grid and solution arrays */
      int sizex = dims[0] + dims[1];
      int sizeu = dims[0] * dims[1];
      x = (double*) calloc (sizex,sizeof(double));
      U = (double*) calloc (sizeu,sizeof(double));

      /* read grid and solution */
      fread(x,sizeof(double),sizex,in);
      fread(U,sizeof(double),sizeu,in);
      /* done reading */
      fclose(in);

      /* set filenames */
      strcpy(tecfile1,filename); tecfile1[0]='t'; tecfile1[1]='h';
      strcpy(tecfile2,filename); tecfile2[0]='o'; tecfile2[1]='m';
      tecfile1[3] = tecfile2[3] = 'd';
      tecfile1[4] = tecfile2[4] = 'a';
      tecfile1[5] = tecfile2[5] = 't';

      /* variables for the 1D solution */
      double *X, *Y, *U1d, integral;
      int    NI,NJ,i,j,ind;
      NI = dims[0]; NJ = dims[1];
      X = x;
      Y = x + NI;

      /* calculate integral of the 2D PDF */
      integral = 0.0;
      for (i=0; i<NI; i++) {
        for (j=0; j<NJ; j++) {
          int q = i + NI * j;
          double dx, dy;
          if      (i == 0)    dx = (X[i+1]-X[i]  );
          else if (i == NI-1) dx = (X[i]  -X[i-1]);
          else                dx = (X[i+1]-X[i-1])/2.0;
          if      (j == 0)    dy = (Y[j+1]-Y[j]  );
          else if (j == NJ-1) dy = (Y[j]  -Y[j-1]);
          else                dy = (Y[j+1]-Y[j-1])/2.0;
          integral += U[q] * dx * dy;
        }
      }
      printf("\tIntegral of 2D PDF: %1.16E\n",integral);

      /* Theta */
      printf("\tExtracting marginal PDF for theta.\n");
      U1d = (double*) calloc (NI,sizeof(double));

      /* copy U - integrate along the other dimension */
      for (i=0; i<NI; i++) {
        U1d[i] = 0;
        for (j=0; j<dims[1]; j++) {
          int q = i + NI * j;
          double dy;
          if      (j == 0)    dy = (Y[j+1]-Y[j]  );
          else if (j == NJ-1) dy = (Y[j]  -Y[j-1]);
          else                dy = (Y[j+1]-Y[j-1])/2.0;
          U1d[i] += U[q] * dy;
        }
      }

      /* calculating integral of marginal PDF */
      integral = 0;
      for (i=0; i<NI; i++) {
        double dx;
        if      (i == 0)    dx = (X[i+1]-X[i]  );
        else if (i == NI-1) dx = (X[i]  -X[i-1]);
        else                dx = (X[i+1]-X[i-1])/2.0;
        integral += U1d[i] * dx;
      }
      printf("\tIntegral of marginal PDF for theta: %1.16E\n",integral);

      /* write to file */
      printf("\tWriting file %s.\n",tecfile1);
      WriteText(NI,X,U1d,tecfile1);

      /* free arrays */
      free(U1d);

      /* Omega */
      printf("\tExtracting marginal PDF for omega.\n");
      U1d = (double*) calloc (NJ,sizeof(double));

      /* copy U - integrate along the other dimension */
      for (j=0; j<NJ; j++) {
        U1d[j] = 0;
        for (i=0; i<NI; i++) {
          int q = i + NI * j;
          double dx;
          if      (i == 0)    dx = (X[i+1]-X[i]  );
          else if (i == NI-1) dx = (X[i]  -X[i-1]);
          else                dx = (X[i+1]-X[i-1])/2.0;
          U1d[j] += U[q] * dx;
        }
      }

      /* calculating integral of marginal PDF */
      integral = 0;
      for (j=0; j<NJ; j++) {
        double dy;
        if      (j == 0)    dy = (Y[j+1]-Y[j]  );
        else if (j == NJ-1) dy = (Y[j]  -Y[j-1]);
        else                dy = (Y[j+1]-Y[j-1])/2.0;
        integral += U1d[j] * dy;
      }
      printf("\tIntegral of marginal PDF for omega: %1.16E\n",integral);

      /* write to file */
      printf("\tWriting file %s.\n",tecfile2);
      WriteText(NJ,Y,U1d,tecfile2);

      /* done */
      free(U1d);

      /* clean up */
      free(U);
      free(x);
    }
  }

  return(0);
}
