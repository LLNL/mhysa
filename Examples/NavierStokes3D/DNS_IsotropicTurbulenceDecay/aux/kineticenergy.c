/*

  This code calculates the kinetic energy as 
  a function of time. It reads in the solution
  files (assuming they are not overwritten and 
  are available as op_xxxxx.bin).

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

int main()
{
  FILE *out, *in, *inputs;
  double dt;
  int file_op_iter;
  char filename[50], op_file_format[50];

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
  out = fopen("energy.dat","w");
  strcpy(filename,"op_00000.bin");
  double energy0 = 1;
  while(1) {
    in = fopen(filename,"rb");

    if (!in) {
      printf("No more files found. Exiting.\n");
      break;
    } else {
      printf("File %s.\n",filename);
      int ndims, nvars, dims[3];
      double energy = 0.0, *U,*x,*y,*z;

      fread(&ndims,sizeof(int),1,in);
      fread(&nvars,sizeof(int),1,in);

      if (ndims != 3) {
        printf("Error: ndims in %s not equal to 3!\n",filename);
        return(0);
      }
      if (nvars != 5) {
        printf("Error: nvars in %s not equal to 4!\n",filename);
        return(0);
      }

      fread(dims,sizeof(int),ndims,in);
      printf("Dimensions: %d x %d x %d\n",dims[0],dims[1],dims[2]);
      int N = dims[0], N3 = N*N*N;
      if ((dims[1] != N) || (dims[2] != N)) {
        printf("Error: Unequal number of points in the 3 dimensions!\n");
        return(0);
      }

      double pi = 4*atan(1.0);
      double dx = 2*pi / (double)N;
      double dx3 = dx*dx*dx;

      x = (double*) calloc (N,sizeof(double));
      y = (double*) calloc (N,sizeof(double));
      z = (double*) calloc (N,sizeof(double));
      U = (double*) calloc (N3*nvars,sizeof(double));

      fread(x,sizeof(double),N,in);
      fread(y,sizeof(double),N,in);
      fread(z,sizeof(double),N,in);
      fread(U,sizeof(double),N3*nvars,in);

      int p;
      for (p = 0; p < N3; p++) {
        double rho, u, v, w;
        rho = U[nvars*p];
        u = U[nvars*p+1] / rho;
        v = U[nvars*p+2] / rho;
        w = U[nvars*p+3] / rho;
        energy += ((u*u + v*v + w*w) * dx3);
      }
      energy /= 2.0;

      /* non dimensionalize */
      if (!count) energy0 = energy;
      energy /= energy0;

      free(U);
      free(x);
      free(y);
      free(z);
      fprintf(out, "%1.16E\t%1.16E\n", samay, energy);
      fclose(in);
    }

    count++;
    samay += (dt * (double)file_op_iter);
    IncrementFilename(filename);
  }

  fclose(out);
  return(0);
}
