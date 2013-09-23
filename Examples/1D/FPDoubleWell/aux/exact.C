#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

int main(){
  
	int NI,ndims;
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
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 1) {
    std::cout << "ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << "\n";

  double nu;
  std::cout << "Reading file \"physics.inp\"...\n";
  in.open("physics.inp");
  if (!in) {
    std::cout << "Error: Input file \"physics.inp\" not found. Default values will be used.\n";
    nu = 1.0;
  } else {
    char word[500];
    in >> word;
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        in >> word;
        if (!strcmp(word,"q"))     in >> nu;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
	std::cout << "q:\t\t\t" << nu << "\n";

	int i;
	double dx = 5.0 / ((double)(NI-1));

	double *x, *u;
	x = (double*) calloc (NI, sizeof(double));
	u = (double*) calloc (NI, sizeof(double));

  double N = 1.0/2622.75;
	for (i = 0; i < NI; i++){
		x[i] = -2.5 + i*dx;
    u[i] = N*exp(-(2.0*x[i]*x[i]*(x[i]*x[i]-2.0))/nu);
	}

  FILE *out;
	out = fopen("exact.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
  fprintf(out,"\n");
	for (i = 0; i < NI; i++)	fprintf(out,"%lf ",u[i]);						
  fprintf(out,"\n");
	fclose(out);

	free(x);
	free(u);

	return(0);
}
