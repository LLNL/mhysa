#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi     = 4.0*atan(1.0);
const double GAMMA  = 1.4;

int main(){
  
	int NI,NJ,ndims;
  int n_iter;
  double tf, dt;
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
        else if (!strcmp(word, "n_iter")) in >> n_iter;
        else if (!strcmp(word, "dt")) in >> dt;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 2) {
    std::cout << "ndims is not 2 in solver.inp. this code is to generate 2D exact solution\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " X " << NJ << "\n";

	int i,j;
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);

  tf = (double)n_iter * dt;
  printf("Final time: %lf\n", tf);

  double rho_inf = 1.0;
  double drho = 0.1;
  double u_inf = 1.0;
  double v_inf = 1.0;
  double P_inf = 1.0/GAMMA;

	double *x, *y, *u0, *u1, *u2, *u3;
  FILE *out;

	x   = (double*) calloc (NI   , sizeof(double));
	y   = (double*) calloc (NJ   , sizeof(double));
	u0  = (double*) calloc (NI*NJ, sizeof(double));
	u1  = (double*) calloc (NI*NJ, sizeof(double));
	u2  = (double*) calloc (NI*NJ, sizeof(double));
	u3  = (double*) calloc (NI*NJ, sizeof(double));
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NJ*i + j;
      double rho, u, v, P;
      rho = rho_inf + drho * sin(2*pi*(x[i]-u_inf*tf)) * cos(2*pi*(y[j]-v_inf*tf));
      P   = P_inf;
      u   = u_inf;
      v   = v_inf;
      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
	  }
	}
	out = fopen("exact.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u0[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u1[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u2[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u3[p]);
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

	x   = (double*) calloc (NI   , sizeof(double));
	y   = (double*) calloc (NJ   , sizeof(double));
	u0  = (double*) calloc (NI*NJ, sizeof(double));
	u1  = (double*) calloc (NI*NJ, sizeof(double));
	u2  = (double*) calloc (NI*NJ, sizeof(double));
	u3  = (double*) calloc (NI*NJ, sizeof(double));
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NJ*i + j;
      double rho, u, v, P;
      rho = rho_inf + drho * sin(2*pi*x[i]) * cos(2*pi*y[j]);
      P   = P_inf;
      u   = u_inf;
      v   = v_inf;
      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
	  }
	}
	out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u0[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u1[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u2[p]);
    }
  }
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16E ",u3[p]);
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
