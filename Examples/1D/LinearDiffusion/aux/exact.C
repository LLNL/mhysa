#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

int main(){
  
	int NI,ndims,niter;
  std::ifstream in;
  double nu, dt,final_time;

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
        else if (!strcmp(word, "dt")) in >> dt;
        else if (!strcmp(word, "n_iter")) in >> niter;
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

  final_time = (double)niter * dt;
	std::cout << "Final Time:\t\t\t" << final_time << "\n";


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
        if (!strcmp(word,"diffusion"))     in >> nu;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
	std::cout << "Diffusion Coeff:\t\t\t" << nu << "\n";

	int i;
	double dx = 1.0 / ((double)NI);

	double *x, *u;
	x = (double*) calloc (NI, sizeof(double));
	u = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
		u[i] = exp(-nu*4*pi*pi*final_time)*sin(2*pi*x[i]);
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
