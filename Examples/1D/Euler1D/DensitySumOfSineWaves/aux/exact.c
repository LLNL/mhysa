#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double raiseto(double x,double a) 
{
  return(exp(a*log(x)));
}

int main() 
{
	int NI,ndims,niter;
  double dt, pi = 4.0*atan(1.0), gamma = 1.4;
  FILE *in, *out;

  double rho_inf = 1.0,
         u_inf   = 1.0,
         p_inf   = 1.0/gamma;

  /* default values */
  ndims = 1;
  NI    = 64;
  niter = 1000;
  dt    = 0.001;

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if      (!strcmp(word, "ndims"))  fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size"))   fscanf(in,"%d",&NI);
        else if (!strcmp(word, "n_iter")) fscanf(in,"%d",&niter);
        else if (!strcmp(word, "dt"))     fscanf(in,"%lf",&dt);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. Make sure the correct solver.inp file is present.\n");
    return(0);
  }
	printf("Grid: %d\n",NI);
  printf("Input maximum wavenumber (typically NI/2): ");
  int limit; scanf("%d",&limit);
  printf("Max wavenumber: %d\n",limit);

	int i,k;
	double dx = 1.0 / ((double)NI);
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  double factor = 0.01;

  srand(time(NULL));
  double *phi = (double*) calloc (limit+1,sizeof(double));
  for (k=1; k<=limit; k++) {
    phi[k] = -pi + 2*pi*(((double) rand()) / ((double) RAND_MAX));
//    phi[k] = pi/2.0;;
  }

	double *x, *rho,*rhou,*e;
	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = -0.5 + i*dx;
    double RHO,U,P,drho=0;
    for (k=1; k<=limit; k++) {
      double Ak = factor * raiseto(((double)k),-5.0/6.0);
      drho += (Ak * cos(2*pi*((double)k)*(x[i]-u_inf*tf)+phi[k]));
    }
    RHO = rho_inf + drho;
    U   = u_inf;
    P   = p_inf;
    rho[i]  = RHO;
    rhou[i] = RHO*U;
    e[i]    = P/0.4 + 0.5*RHO*U*U;
    printf("%d: %f\n",i,rho[i]);
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
	free(x);
	free(rho);
	free(rhou);
	free(e);

	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));
	for (i = 0; i < NI; i++){
		x[i] = -0.5 + i*dx;
    double RHO,U,P,drho=0;
    for (k=1; k<=limit; k++) {
      double Ak = factor * raiseto(((double)k),-5.0/6.0);
      drho += (Ak * cos(2*pi*((double)k)*(x[i])+phi[k]));
    }
    RHO = rho_inf + drho;
    U   = u_inf;
    P   = p_inf;
    rho[i]  = RHO;
    rhou[i] = RHO*U;
    e[i]    = P/0.4 + 0.5*RHO*U*U;
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

  free(phi);
	return(0);
}
