#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(){

  double eta = 0.0001;
  printf("Enter the perturbation constant (eta): ");
  scanf("%lf",&eta);
  
	int NI,ndims;
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims"))     fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) fscanf(in,"%d",&NI);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d\n",NI);

	int i;
	double dx = 1.0 / ((double)(NI));

	double *x, *rho,*rhou,*e;
	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
    double pi = 4.0 * atan(1.0);
    double phi = -sin(2*pi*x[i])/(2*pi);
    double RHO = exp(-phi);
    double U   = 0.0;
    double P   = exp(-phi) + eta * exp (-100*(x[i]-0.5)*(x[i]-0.5));
    rho[i]  = RHO;
    rhou[i] = RHO*U;
    e[i]    = P/0.4 + 0.5*RHO*U*U;
	}

  FILE *out;
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
