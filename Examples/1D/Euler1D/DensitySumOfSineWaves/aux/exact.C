#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

double raiseto(double x,double a) {
  return(exp(a*log(x)));
}

int main(){
  
	int NI,ndims,niter;
  double dt;
  std::ifstream in;
  std::cout << "Reading file \"solver.inp\"...\n";
  in.open("solver.inp");
  if (!in) {
    std::cout << "Error: Input file \"solver.inp\" not found. Default values will be used.\n";
  } else {
    char word[500];
    in >> word;
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        in >> word;
        if (!strcmp(word, "ndims"))     in >> ndims;
        else if (!strcmp(word, "size")) in >> NI;
        else if (!strcmp(word, "n_iter")) in >> niter;
        else if (!strcmp(word, "dt")) in >> dt;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 1) {
    std::cout << "ndims is not 1 in solver.inp. this code is to generate 1D exact conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << "\n";
  std::cout << "Input maximum wavenumber (typically NI/2): ";
  int limit; std::cin >> limit;

	int i,k;
	double dx = 1.0 / ((double)NI);
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

  double factor = 0.01;

  srand(time(NULL));
  double *phi = (double*) calloc (limit,sizeof(double));
  for (k=1; k<limit; k++) {
//    phi[k] = -pi + 2*pi*(((double) rand()) / ((double) RAND_MAX));
    phi[k] = pi/2.0;;
  }

	double *x, *rho,*rhou,*e;
  FILE *out;

	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));
	for (i = 0; i < NI; i++){
		x[i] = -0.5 + i*dx;
    double RHO,U,P,drho=0;
    for (k=1; k<limit; k++) {
      double Ak = factor * raiseto(((double)k),-5.0/6.0);
      drho += (Ak * cos(2*pi*((double)k)*(x[i]-tf)+phi[k]));
    }
    RHO = 1.0 + drho;
    U   = 1.0;
    P   = 1.0/1.4;
    rho[i]  = RHO;
    rhou[i] = RHO*U;
    e[i]    = P/0.4 + 0.5*RHO*U*U;
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
    for (k=1; k<limit; k++) {
      double Ak = factor * raiseto(((double)k),-5.0/6.0);
      drho += (Ak * cos(2*pi*((double)k)*(x[i])+phi[k]));
    }
    RHO = 1.0 + drho;
    U   = 1.0;
    P   = 1.0/1.4;
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
