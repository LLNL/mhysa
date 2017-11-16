/*
 * This code extracts the minimum pressure and density at
 * the vortex core as a function of time for an unsteady 
 * simulation.
 * Output files have to be binary, and they can't be 
 * overwritten; i.e., it needs files named as 
 * op_xxxxx.bin
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
  double dt, samay;
  int file_op_iter, restart_iter=0, t, count;
  char filename[50], op_file_format[50];

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
   			else if (!strcmp(word, "restart_iter"     ))  fscanf(inputs,"%d" ,&restart_iter  );
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  samay = (double) restart_iter *dt;
  if (!restart_iter) {
    printf("Writing to core_pressure.dat (new file).\n");
    out = fopen("core_pressure.dat","w");
  } else {
    printf("Writing to core_pressure.dat (append).\n");
    out = fopen("core_pressure.dat","a");
  }
  strcpy(filename,"op_00000.bin");
  for (t=0; t<restart_iter; t++) {
    if ((t+1)%file_op_iter == 0) IncrementFilename(filename);
  }
  count = 0;
  while(1) {
    in = fopen(filename,"rb");

    if (!in) {
      printf("No more files found. Exiting.\n");
      break;
    } else {
      printf("File %s.\n",filename);
      int ndims, nvars, dims[2];
      double *U,*x,*y;

      fread(&ndims,sizeof(int),1,in);
      fread(&nvars,sizeof(int),1,in);

      if (ndims != 2) {
        printf("Error: ndims in %s not equal to 2!\n",filename);
        return(0);
      }
      if (nvars != 4) {
        printf("Error: nvars in %s not equal to 4!\n",filename);
        return(0);
      }

      fread(dims,sizeof(int),ndims,in);
      printf("Dimensions: %d x %d\n",dims[0],dims[1]);

      int NI = dims[0];
      int NJ = dims[1];
      int sizeu = NI*NJ*nvars;

      x = (double*) calloc (NI   ,sizeof(double));
      y = (double*) calloc (NJ   ,sizeof(double));
      U = (double*) calloc (sizeu,sizeof(double));

      fread(x,sizeof(double),NI   ,in);
      fread(y,sizeof(double),NJ   ,in);
      fread(U,sizeof(double),sizeu,in);

      int q;
      double min_p, min_rho;
      for (q = 0; q < NI*NJ; q++) {
        double rho, u, v, e, p;
        rho = U[nvars*q+0];
        u   = U[nvars*q+1] / rho;
        v   = U[nvars*q+2] / rho;
        e   = U[nvars*q+3];
        p = 0.4 * (e - 0.5*rho*(u*u+v*v));

        if (!q) { min_p = p; min_rho = rho; }
        
        if (rho < min_rho) min_rho = rho;
        if (p   < min_p  ) min_p   = p;
      }

      free(U);
      free(x);
      free(y);

      fprintf(out, "%1.16E\t%1.16E\t%1.16E\n", samay, min_p, min_rho);
      fclose(in);
    }

    count++;
    samay += (dt * (double)file_op_iter);
    IncrementFilename(filename);
  }

  fclose(out);
  return(0);
}
