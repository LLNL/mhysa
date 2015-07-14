/*
  Code to generate the initial solution for:
  Case: Limited Area Channel Flow
  Model: ShallowWater2D

  References:
  Zhu, Et. al., "Variational Data Assimilation with a Variable Resolution 
  Finite-Element Shallow-Water Equations Model", Monthly Weather Review,
  122, 1994, pp. 946--965
  http://dx.doi.org/10.1175/1520-0493(1994)122%3C0946:VDAWAV%3E2.0.CO;2
  Eqns. (4.1)-(4.3)

  Needs: solver.inp
  Writes out: initial.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

  /* some global parameters */
double L = 6000000.0,
       D = 4400000.0,
       g = 10.0,
       fhat = 1e-4,
       beta = 1.5e-11,
       H0 = 2000.0,
       H1 = 220.0,
       H2 = 133.0;

double hfunction(double x,double y) {
  double pi = 4.0*atan(1.0);
  double h = H0 + H1*tanh((9.0*(0.5*D-y))/(2.0*D))
             + H2 * (1.0/cosh(9.0*(0.5*D-y)/D)) * (1.0/cosh(9.0*(0.5*D-y)/D)) * sin(2*pi*x/L);
  return(h);
}

int main(){
	int   NI = 15, NJ = 12, ndims = 2;
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
	double dx = L / ((double)NI);
	double dy = D / ((double)NJ-1);

	double *x, *y, *u0, *u1, *u2;
	x  = (double*) calloc (NI   , sizeof(double));
	y  = (double*) calloc (NJ   , sizeof(double));
	u0 = (double*) calloc (NI*NJ, sizeof(double));
	u1 = (double*) calloc (NI*NJ, sizeof(double));
	u2 = (double*) calloc (NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NI*j + i;
      double h, u, v, hE, hW, hN, hS;
      h  = hfunction(x[i],y[j]);
      hE = hfunction(x[i]+dx,y[j]);
      hW = hfunction(x[i]-dx,y[j]);
      hN = hfunction(x[i],y[j]+dy);
      hS = hfunction(x[i],y[j]-dy);

      double f = fhat + beta * (y[j] - 0.5*D);
      double h_x, h_y;
      h_x = (hE - hW) / (2*dx);
      h_y = (hN - hS) / (2*dy);

      u = -(g/f) * h_y;
      v =  (g/f) * h_x;

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
  }

	free(x);
	free(y);
	free(u0);
	free(u1);
	free(u2);

	return(0);
}
