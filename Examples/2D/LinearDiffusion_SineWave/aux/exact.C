#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

int main(){
  
	int NI,NJ,ndims,niter;
  std::ifstream in;
  double nu_x,nu_y,dt,final_time;

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
        else if (!strcmp(word, "size")) in >> NI >> NJ;
        else if (!strcmp(word, "dt")) in >> dt;
        else if (!strcmp(word, "n_iter")) in >> niter;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 2) {
    std::cout << "ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " X " << NJ << "\n";

  final_time = (double)niter * dt;
	std::cout << "Final Time:\t\t\t" << final_time << "\n";

  std::cout << "Reading file \"physics.inp\"...\n";
  in.open("physics.inp");
  if (!in) {
    std::cout << "Error: Input file \"physics.inp\" not found. Default values will be used.\n";
    nu_x = nu_y = 1.0;
  } else {
    char word[500];
    in >> word;
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        in >> word;
        if (!strcmp(word,"diffusion"))     in >> nu_x >> nu_y;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
	std::cout << "Diffusion Coeff:\t\t\t" << nu_x << ", " << nu_y << "\n";

	int i,j;
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);

	double *x, *y, *u;
	x = (double*) calloc (NI   , sizeof(double));
	y = (double*) calloc (NJ   , sizeof(double));
	u = (double*) calloc (NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NJ*i + j;
		  u[p] = exp(-pi*pi*(4*nu_x+4*nu_y)*final_time) * sin(2*pi*x[i]) * sin(2*pi*y[j]);
	  }
	}

  FILE *out;
	out = fopen("exact.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%lf ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%lf ",u[p]);
    }
  }
  fprintf(out,"\n");
	fclose(out);

	free(x);
	free(y);
	free(u);

	return(0);
}
