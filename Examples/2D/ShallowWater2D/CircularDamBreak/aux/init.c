/*
  Code to generate the initial solution for:
  Case: Small Peturbation 
  Model: ShallowWater2D

  Reference:
  Delis, Katsaounis, "Numerical solution of the two-dimensional
  shallow water equations by the application of relaxation methods",
  Applied Mathematical Modelling, 29 (2005), pp. 754--783
  http://dx.doi.org/10.1016/j.apm.2004.11.001
  Section 6.3

  Needs: solver.inp
  Writes out: initial.inp, topography.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(){
	int   NI = 101, NJ = 101, ndims = 2;
  char  ip_file_type[50]; 
  strcpy(ip_file_type,"ascii");

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        } else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d X %d\n",NI,NJ);

	int i,j;
  double length_x = 50.0, length_y = 50.0;
	double dx = length_x / ((double)NI-1);
	double dy = length_y / ((double)NJ-1);

	double *x, *y, *u0, *u1, *u2, *b;
	x  = (double*) calloc (NI   , sizeof(double));
	y  = (double*) calloc (NJ   , sizeof(double));
	u0 = (double*) calloc (NI*NJ, sizeof(double));
	u1 = (double*) calloc (NI*NJ, sizeof(double));
	u2 = (double*) calloc (NI*NJ, sizeof(double));
	b  = (double*) calloc (NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      double r = sqrt((x[i]-length_x/2)*(x[i]-length_x/2)+(y[j]-length_y/2)*(y[j]-length_y/2));
      int p = NI*j + i;
      double h, u, v;
      b[p] = 0.0;
      if      (r <= 11.0) h = 10.0;
      else                h = 1.0;
      u = v = 0.0;
      u0[p] = h;
      u1[p] = h*u;
      u2[p] = h*v;
	  }
	}

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
  	out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++) fprintf(out,"%1.16e ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++) fprintf(out,"%1.16e ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NI*j + i;
        fprintf(out,"%1.16e ",u0[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NI*j + i;
        fprintf(out,"%1.16e ",u1[p]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NI*j + i;
        fprintf(out,"%1.16e ",u2[p]);
      }
    }
    fprintf(out,"\n");
	  fclose(out);
    printf("Writing ASCII topography file topography.inp\n");
  	out = fopen("topography.inp","w");
    for (i = 0; i < NI; i++) fprintf(out,"%1.16e ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++) fprintf(out,"%1.16e ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++) {
	    for (i = 0; i < NI; i++) {
        int p = NI*j + i;
        fprintf(out,"%1.16e ",b[p]);
      }
    }
    fprintf(out,"\n");
	  fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    double *U = (double*) calloc (3*NI*NJ,sizeof(double));
    for (i=0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        int p = NI*j + i;
        U[3*p+0] = u0[p];
        U[3*p+1] = u1[p];
        U[3*p+2] = u2[p];
      }
    }
    fwrite(U,sizeof(double),3*NI*NJ,out);
    free(U);
    fclose(out);
    printf("Writing binary topography file topography.inp\n");
  	out = fopen("topography.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(b,sizeof(double),NI*NJ,out);
    fclose(out);
  }

	free(x);
	free(y);
	free(u0);
	free(u1);
	free(u2);
	free(b);

	return(0);
}
