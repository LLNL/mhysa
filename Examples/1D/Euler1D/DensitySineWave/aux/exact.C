#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

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

	int i;
	double dx = 1.0 / ((double)NI);
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);

	double *x, *rho,*rhou,*e;
	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
    double RHO,U,P;
    RHO = 1.0 + 0.1*sin(2*pi*(x[i]-tf));
    U   = 1.0;
    P   = 1.0/1.4;
    rho[i]  = RHO;
    rhou[i] = RHO*U;
    e[i]    = P/0.4 + 0.5*RHO*U*U;
	}

  FILE *out;
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

	return(0);
}
