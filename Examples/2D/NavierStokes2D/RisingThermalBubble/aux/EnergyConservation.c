/*
  Rising Thermal Bubble:-
  This code reads in the binary solution output
  files for the rising thermal bubble problem
  and computes the total (kinetic + potential)
  energy as a function of time.

  Note: Use this only if the solution was written
  out in unsteady form (i.e. op_overwrite is "no" 
  in solver.inp).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct _parameters_{
  double grav_x, grav_y, gamma;
} Parameters;

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

int CalculateEnergy(char *fname, double *kinetic_energy, double *potential_energy,
                    double *internal_energy, double *total_energy, void *p)
{
  Parameters *params = (Parameters*) p;
  FILE *in; in = fopen(fname,"rb");

  if (!in) return(-1);

  printf("Reading file %s.\n",fname);
  int ndims, nvars;
  double *U,*x;

  /* read the file headers */
  fread(&ndims,sizeof(int),1,in);
  fread(&nvars,sizeof(int),1,in);

  /* some checks */
  if (ndims != 2) {
    printf("Error: ndims in %s not equal to 2!\n",fname);
    return(1);
  }
  if (nvars != 4) {
    printf("Error: nvars in %s not equal to 4!\n",fname);
    return(1);
  }

  /* read dimensions */
  int dims[ndims];
  fread(dims,sizeof(int),ndims,in);
  printf("Dimensions: %d x %d\n",dims[0],dims[1]);
  printf("Nvars     : %d\n",nvars);

  /* allocate grid and solution arrays */
  x = (double*) calloc (dims[0]+dims[1]       ,sizeof(double));
  U = (double*) calloc (dims[0]*dims[1]*nvars ,sizeof(double));

  /* read grid and solution */
  fread(x,sizeof(double),dims[0]+dims[1]      ,in);
  fread(U,sizeof(double),dims[0]*dims[1]*nvars,in);
  /* done reading */
  fclose(in);

  int imax = dims[0];
  int jmax = dims[1];

  /* calculate primitive variables */
  double *X           = x;
  double *Y           = x+imax;
  double grav_x       = params->grav_x;
  double grav_y       = params->grav_y;
  double gamma        = params->gamma;

  int i, j;
  *kinetic_energy = 0;
  *potential_energy = 0;
  *internal_energy = 0;
  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      int p = i + imax*j;

      /* cell volume */
      double dx, dy;
      if (i == 0)             dx = X[1]-X[0];
      else if (i == imax-1)   dx = X[imax-1] - X[imax-2];
      else                    dx = 0.5 * (X[i+1]-X[i-1]);
      if (j == 0)             dy = Y[1] - Y[0];
      else if (j == jmax-1)   dy = Y[jmax-1] - Y[jmax-2];
      else                    dy = 0.5 * (Y[j+1]-Y[j-1]);
      double vol = dx * dy;

      double rho, uvel, vvel, E, P, theta;
      rho   = U[nvars*p+0];
      uvel  = U[nvars*p+1] / rho;
      vvel  = U[nvars*p+2] / rho;
      E     = U[nvars*p+3];
      P     = (gamma-1.0) * (E - 0.5*rho*(uvel*uvel+vvel*vvel));

      *kinetic_energy += (0.5*rho*(uvel*uvel+vvel*vvel)*vol);
      *potential_energy += (rho * (grav_x*X[i] + grav_y*Y[j]) * vol);
      *internal_energy += ( P/(gamma-1.0) * vol);
    }
  }
  *total_energy = *kinetic_energy + *potential_energy + *internal_energy;

  /* clean up */
  free(U);
  free(x);
}

int main()
{
  FILE *out, *in, *inputs;
  char filename[50], overwrite[50], op_file_format[50];
  double dt = 0;
  int file_op_iter = 0;

  printf("Reading solver.inp.\n");
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
        else if (!strcmp(word, "dt"               ))  fscanf(inputs,"%lf",&dt);
        else if (!strcmp(word, "file_op_iter"     ))  fscanf(inputs,"%d" ,&file_op_iter);
   			else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,overwrite);
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }
  if (!strcmp(overwrite,"yes")) {
    printf("Error: op_overwrite needs to be \"no\" in solver.inp.\n");
    return(0);
  }

  Parameters params;
  /* default values */
  params.grav_x   = 0.0;
  params.grav_y   = 9.8;
  params.gamma    = 1.4;
  /* read these parameters from file */
  printf("Reading physics.inp.\n");
  inputs = fopen("physics.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
	  char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    fscanf(inputs,"%s",word);
   			if      (!strcmp(word, "gamma"))    fscanf(inputs,"%lf",&params.gamma);
   			else if (!strcmp(word, "gravity")) {
          fscanf(inputs,"%lf",&params.grav_x);
          fscanf(inputs,"%lf",&params.grav_y);
        }
      }
    }
    fclose(inputs);
  }

  strcpy(filename,"op_00000.bin");
  out = fopen("energy.dat","w");
  int count = 0;
  while(1) {
    double kinetic_energy = 0,
           potential_energy = 0,
           internal_energy = 0,
           total_energy = 0,
           waqt = dt * file_op_iter * count;
    int err = CalculateEnergy(filename, &kinetic_energy, &potential_energy,
                              &internal_energy,&total_energy,&params);
    if (err == -1) {
      printf("No more files found. Exiting.\n");
      break;
    }

    fprintf(out,"%1.16e %1.16e %1.16e %1.16e %1.16e\n", waqt,
            kinetic_energy, potential_energy, internal_energy,
            total_energy);
    count++;
    IncrementFilename(filename);
  }
  fclose(out);
  return(0);
}
