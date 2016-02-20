#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double GAMMA  = 1.4;

int main(){
  
	int NI,NJ,ndims;
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
        else if (!strcmp(word, "size")) in >> NI >> NJ;
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

	int i,j;
	double dx = 1.0 / ((double)(NI-1));
	double dy = 1.0 / ((double)(NJ-1));

	double *x, *y, *u0, *u1, *u2, *u3;
	x   = (double*) calloc (NI   , sizeof(double));
	y   = (double*) calloc (NJ   , sizeof(double));
	u0  = (double*) calloc (NI*NJ, sizeof(double));
	u1  = (double*) calloc (NI*NJ, sizeof(double));
	u2  = (double*) calloc (NI*NJ, sizeof(double));
	u3  = (double*) calloc (NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = -0.5 + i*dx;
	  	y[j] = -0.5 + j*dy;
      int p = NJ*i + j;

      double rho, u, v, P;
      if ((x[i] < 0) && (y[j] < 0)) {
        rho = 1.0;
        P   = 1.0;
        u   = -0.75;
        v   = 0.5;
      } else if ((x[i] >= 0) && (y[j] <  0)) {
        rho = 3.0;
        P   = 1.0;
        u   = -0.75;
        v   = -0.5;
      } else if ((x[i] <  0) && (y[j] >= 0)) {
        rho = 2.0;
        P   = 1.0;
        u   = 0.75;
        v   = 0.5;
      } else if ((x[i] >= 0) && (y[j] >= 0)) {
        rho = 1.0;
        P   = 1.0;
        u   = 0.75;
        v   = -0.5;
      }

      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
	  }
	}

  FILE *out;
	out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%lf ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%lf ",u0[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%lf ",u1[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%lf ",u2[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%lf ",u3[p]);
    }
  }
  fprintf(out,"\n");
	fclose(out);

	free(x);
	free(y);
	free(u0);
	free(u1);
	free(u2);
	free(u3);

	return(0);
}
