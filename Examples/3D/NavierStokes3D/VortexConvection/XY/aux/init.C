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
  
	int NI,NJ,NK,ndims;
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
        else if (!strcmp(word, "size")) in >> NI >> NJ >> NK;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 3) {
    std::cout << "ndims is not 3 in solver.inp. this code is to generate 3D initial conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " X " << NJ << " X " << NK << "\n";

	int i,j,k;
	double dx = 10.0 / ((double)NI);
	double dy = 10.0 / ((double)NJ);
	double dz = 10.0 / ((double)NK);
  double u_inf = 0.5;
  double v_inf = 0.0;
  double b = u_inf;
  double x0 = 5.0, y0 = 5.0;

	double *x, *y, *z, *u0, *u1, *u2, *u3, *u4;
	x   = (double*) calloc (NI      , sizeof(double));
	y   = (double*) calloc (NJ      , sizeof(double));
	z   = (double*) calloc (NK      , sizeof(double));
	u0  = (double*) calloc (NI*NJ*NK, sizeof(double));
	u1  = (double*) calloc (NI*NJ*NK, sizeof(double));
	u2  = (double*) calloc (NI*NJ*NK, sizeof(double));
	u3  = (double*) calloc (NI*NJ*NK, sizeof(double));
	u4  = (double*) calloc (NI*NJ*NK, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
    	for (k = 0; k < NK; k++){
	    	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;
        int p = NJ*NK*i + NK*j + k;
        double rsq = (x[i]-x0)*(x[i]-x0) + (y[j]-y0)*(y[j]-y0);
        double rho, u, v, w, P;
        double du, dv;
        rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
        P   = power(rho,GAMMA);
        du  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * (y[j]-y0);
        dv  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * (x[i]-x0);
        u   = u_inf + du;
        v   = v_inf + dv;
        w   = 0.0;
        u0[p] = rho;
        u1[p] = rho*u;
        u2[p] = rho*v;
        u3[p] = rho*w;
        u4[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}

  FILE *out;
	out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%lf ",y[j]);
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)  fprintf(out,"%lf ",z[k]);
  fprintf(out,"\n");
	for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%lf ",u0[p]);
      }
    }
  }
  fprintf(out,"\n");
	for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%lf ",u1[p]);
      }
    }
  }
  fprintf(out,"\n");
	for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%lf ",u2[p]);
      }
    }
  }
  fprintf(out,"\n");
	for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%lf ",u3[p]);
      }
    }
  }
  fprintf(out,"\n");
	for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%lf ",u4[p]);
      }
    }
  }
  fprintf(out,"\n");
	fclose(out);

	free(x);
	free(y);
	free(z);
	free(u0);
	free(u1);
	free(u2);
	free(u3);
	free(u4);

	return(0);
}
