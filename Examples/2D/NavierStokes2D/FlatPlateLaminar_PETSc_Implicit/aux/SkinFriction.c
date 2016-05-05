/* 
  This code calculates the skin friction
  as a function of Reynolds number for
  the flow over a flat plate.

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main()
{
  FILE *out, *in, *inputs;
  char filename[50], op_file_format[50];
  double Re_inf, M_inf;

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
   			if      (!strcmp(word, "op_file_format"))  fscanf(inputs,"%s" ,op_file_format);
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

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
   			if      (!strcmp(word, "Re"))   fscanf(inputs,"%lf",&Re_inf);
   			else if (!strcmp(word, "Minf")) fscanf(inputs,"%lf",&M_inf);
      }
    }
    fclose(inputs);
  }

  strcpy(filename,"op.bin");
  in  = fopen(filename,"rb");
  out = fopen("SkinFriction.dat","w");

  if (!in) {
    printf("File op.bin found. Exiting.\n");
    return(0);
  } else {
    printf("Reading file %s.\n",filename);
    int ndims, nvars, dims[3];
    double *U,*x,*y,*z;

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

    int i, j, q;
    j = 0;
    for (i = 1; i < NI-1; i++) {
      int q1,q2,q3;
      double rho1, u1, v1, e1, p1;
      double rho2, u2, v2, e2, p2;
      double rho3, u3, v3, e3, p3;

      /* calculate du_dy at the wall */

      q1   = i + NI*j;
      rho1 = U[nvars*q1+0];
      u1   = U[nvars*q1+1] / rho1;
      v1   = U[nvars*q1+2] / rho1;
      e1   = U[nvars*q1+3];
      p1   = 0.4 * (e1 - 0.5*rho1*(u1*u1+v1*v1));

      q2   = i + NI*(j+1);
      rho2 = U[nvars*q2+0];
      u2   = U[nvars*q2+1] / rho2;
      v2   = U[nvars*q2+2] / rho2;
      e2   = U[nvars*q2+3];
      p2   = 0.4 * (e2 - 0.5*rho2*(u2*u2+v2*v2));

      q3   = i + NI*(j+2);
      rho3 = U[nvars*q3+0];
      u3   = U[nvars*q3+1] / rho3;
      v3   = U[nvars*q3+2] / rho3;
      e3   = U[nvars*q3+3];
      p3   = 0.4 * (e3 - 0.5*rho3*(u3*u3+v3*v3));

      double u0 = -u1;

      double u_wall  = 0.5*u0 + 0.5*u1;
      double u_flow1 = 0.5*u1 + 0.5*u2;
      double u_flow2 = 0.5*u2 + 0.5*u3;
      double dy = y[j+1] - y[j];

      double du_dy = (-3.0*u_wall + 4.0*u_flow1 - u_flow2) / (2.0*dy);

      /* Calculate Re_x */
      double Re_x = Re_inf * x[i];

      /* Calculate exact solution */
      double exact_cf = 0.664 / sqrt(Re_x);

      /* Calculate the computed cf */
      double cf = (M_inf/Re_inf) * (2.0/(M_inf*M_inf)) * du_dy;

      if (x[i] > 0) fprintf(out,"%1.16E   %1.16E    %1.16E    %1.16E\n",Re_x,cf,exact_cf,du_dy);
    }

    free(U);
    free(x);
    free(y);

    fclose(in);
    fclose(out);
  }

  return(0);
}
