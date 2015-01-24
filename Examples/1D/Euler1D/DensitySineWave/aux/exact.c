/*
  Code to generate the initial and exact solution for:
  Case: Density Wave Advection
  Model: Euler1D

  Needs: solver.inp
  Writes out: initial.inp, exact.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main() {
  
	int     NI,ndims,niter;
  double  dt, pi = 4.0*atan(1.0), gamma = 1.4;
  FILE    *in, *out;

  /* default values */
  NI    = 40;
  ndims = 1;
  niter = 100;
  dt    = 0.01;

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if      (!strcmp(word, "ndims"))    fscanf(in,"%d" , &ndims);
        else if (!strcmp(word, "size"))     fscanf(in,"%d" ,    &NI);
        else if (!strcmp(word, "n_iter"))   fscanf(in,"%d" , &niter);
        else if (!strcmp(word, "dt"))       fscanf(in,"%lf",    &dt);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. Make sure the correct solver.inp is being used.\n");
    return(0);
  }
	printf("Grid: %d\n", NI);

	int i;
	double dx = 1.0 / ((double)NI);
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  double RHO,U,P;
  RHO   = 1.0;
  U     = 1.0;
  P     = 1.0/gamma;

	double *x, *rho,*rhou,*e;
	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
    double DRHO = 0.1*sin(2*pi*(x[i]-U*tf));
    rho[i]  = RHO + DRHO;
    rhou[i] = rho[i]*U;
    e[i]    = P/(gamma-1.0) + 0.5*rho[i]*U*U;
	}
	out = fopen("exact.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",rho[i]);						
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",rhou[i]);						
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",e[i]);						
  fprintf(out,"\n");
	fclose(out);

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
    double DRHO = 0.1*sin(2*pi*x[i]);
    rho[i]  = RHO + DRHO;
    rhou[i] = rho[i]*U;
    e[i]    = P/(gamma-1.0) + 0.5*rho[i]*U*U;
	}
	out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",rho[i]);						
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",rhou[i]);						
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",e[i]);						
  fprintf(out,"\n");
	fclose(out);

	free(x);
	free(rho);
	free(rhou);
	free(e);

	return(0);
}
