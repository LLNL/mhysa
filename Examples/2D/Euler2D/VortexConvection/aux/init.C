#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi     = 4.0*atan(1.0);
const double GAMMA  = 1.4;

double power(double x,double a)
{
  return(exp(a*log(x)));
}

int main(){
  
	int NI,NJ,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
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
        if (!strcmp(word, "ndims"))             in >> ndims;
        else if (!strcmp(word, "size"))         in >> NI >> NJ;
        else if (!strcmp(word, "ip_file_type")) in >> ip_file_type;
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
	double dx = 10.0 / ((double)NI);
	double dy = 10.0 / ((double)NJ);
  double u_inf = 0.5;
  double v_inf = 0.0;
  double b = u_inf;
  double x0 = 5.0, y0 = 5.0;

	double *x, *y, *u0, *u1, *u2, *u3;
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
      double rsq = (x[i]-x0)*(x[i]-x0) + (y[j]-y0)*(y[j]-y0);
      double rho, u, v, P;
      double du, dv;
      rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
      P   = power(rho,GAMMA);
      du  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * (y[j]-y0);
      dv  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * (x[i]-x0);
      u   = u_inf + du;
      v   = v_inf + dv;
      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
	  }
	}

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
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
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    double *U = (double*) calloc (4*NI*NJ,sizeof(double));
    for (i=0; i < NI; i++) {
      for (j = 0; j < NJ; j++) {
        int p = NJ*i + j;
        int q = NI*j + i;
        U[4*q+0] = u0[p];
        U[4*q+1] = u1[p];
        U[4*q+2] = u2[p];
        U[4*q+3] = u3[p];
      }
    }
    fwrite(U,sizeof(double),4*NI*NJ,out);
    free(U);
    fclose(out);
  }

	free(x);
	free(y);
	free(u0);
	free(u1);
	free(u2);
	free(u3);

	return(0);
}
