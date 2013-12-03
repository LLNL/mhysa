#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi     = 4.0*atan(1.0);
const double GAMMA  = 1.4;

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
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);
	double dz = 1.0 / ((double)NK);

  tf = (double)n_iter * dt;
  printf("Final time: %lf\n",tf);

  double u_inf = 1.0;
  double v_inf = 1.0;
  double w_inf = 1.0;
  double rho_inf = 1.0, drho = 0.1;
  double p_inf = 1.0/GAMMA;

	double *x, *y, *z, *u0, *u1, *u2, *u3, *u4;
  FILE *out;

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

        double rho, u, v, w, P;
        rho = rho_inf + drho * sin(2*pi*x[i]-2*pi*u_inf*tf) * sin(2*pi*y[j]-2*pi*v_inf*tf) * sin(2*pi*z[k]-2*pi*w_inf*tf);
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        u0[p] = rho;
        u1[p] = rho*u;
        u2[p] = rho*v;
        u3[p] = rho*w;
        u4[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
	out = fopen("exact.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)  fprintf(out,"%1.16E ",z[k]);
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u0[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u1[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u2[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u3[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u4[p]);
      }
    }
  }
  fprintf(out,"\n");
	fclose(out);

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

        double rho, u, v, w, P;
        rho = rho_inf + drho * sin(2*pi*x[i]) * sin(2*pi*y[j]) * sin(2*pi*z[k]);
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        u0[p] = rho;
        u1[p] = rho*u;
        u2[p] = rho*v;
        u3[p] = rho*w;
        u4[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
	out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)  fprintf(out,"%1.16E ",z[k]);
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u0[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u1[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u2[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u3[p]);
      }
    }
  }
  fprintf(out,"\n");
  for (k = 0; k < NK; k++)	{
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*NK*i + NK*j + k;
        fprintf(out,"%1.16E ",u4[p]);
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

	return(0);
}
