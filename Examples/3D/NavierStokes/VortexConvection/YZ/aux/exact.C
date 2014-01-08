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
        else if (!strcmp(word, "size")) in >> NI >> NJ >> NK;
        else if (!strcmp(word, "n_iter")) in >> n_iter;
        else if (!strcmp(word, "dt")) in >> dt;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 3) {
    std::cout << "ndims is not 3 in solver.inp. this code is to generate 3D exact solution\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " X " << NJ << " X " << NK << "\n";

	int i,j,k;
	double dx = 10.0 / ((double)NI);
	double dy = 10.0 / ((double)NJ);
	double dz = 10.0 / ((double)NK);

  tf = (double)n_iter * dt; double tff = tf;
  while (tf > 20) tf -= 20; // Time period

  double v_inf = 0.5;
  double w_inf = 0.0;
  double b = v_inf;
  double y0 = 5.0+tf*v_inf, z0 = 5.0;
  if (y0 > 10) y0 -= 10; //periodic domain
  printf("Final time: %lf, Vortex center: %lf, %lf\n",tff,y0,z0);

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

        double ry, rz;
        ry = (y[j] - y0);
        rz = (z[k] - z0);
        if (ry < -5)      { ry += 10; }
        else if (ry > 5)  { ry -= 10; }

        double rsq = ry*ry + rz*rz;
        double rho, u, v, w, P;
        double dv, dw;
        rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
        P   = power(rho,GAMMA);
        dv  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * rz;
        dw  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * ry;
        u   = 0.0;
        v   = v_inf + dv;
        w   = w_inf + dw;
        u0[p] = rho;
        u1[p] = rho*u;
        u2[p] = rho*v;
        u3[p] = rho*w;
        u4[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}

  FILE *out;
	out = fopen("exact.inp","w");
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
	free(u0);
	free(u1);
	free(u2);
	free(u3);

	return(0);
}
