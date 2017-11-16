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
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
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
        else if (!strcmp(word, "ip_file_type")) in >> ip_file_type;
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

	double *x, *y, *z, *U;
  FILE *out;

 	x   = (double*) calloc (NI        , sizeof(double));
	y   = (double*) calloc (NJ        , sizeof(double));
	z   = (double*) calloc (NK        , sizeof(double));
	U   = (double*) calloc (5*NI*NJ*NK, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
  	  	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;
        int p = i + NK*j + NI*NJ*k;

        double rho, u, v, w, P;
        rho = rho_inf + drho * sin(2*pi*x[i]) * sin(2*pi*y[j]) * sin(2*pi*z[k]);
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        U[5*p+0] = rho;
        U[5*p+1] = rho*u;
        U[5*p+2] = rho*v;
        U[5*p+3] = rho*w;
        U[5*p+4] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
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
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+0]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+1]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+2]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+3]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+4]);
        }
      }
    }
    fprintf(out,"\n");
	  fclose(out);

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(z,sizeof(double),NK,out);
    fwrite(U,sizeof(double),5*NI*NJ*NK,out);
    fclose(out);

  }

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
  	  	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;
        int p = i + NK*j + NI*NJ*k;

        double rho, u, v, w, P;
        rho = rho_inf + drho * sin(2*pi*x[i]-2*pi*u_inf*tf) * sin(2*pi*y[j]-2*pi*v_inf*tf) * sin(2*pi*z[k]-2*pi*w_inf*tf);
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        U[5*p+0] = rho;
        U[5*p+1] = rho*u;
        U[5*p+2] = rho*v;
        U[5*p+3] = rho*w;
        U[5*p+4] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII exact solution file exact.inp\n");
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
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+0]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+1]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+2]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+3]);
        }
      }
    }
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)	{
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16E ",U[5*p+4]);
        }
      }
    }
    fprintf(out,"\n");
	  fclose(out);

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary exact solution file exact.inp\n");
  	out = fopen("exact.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(z,sizeof(double),NK,out);
    fwrite(U,sizeof(double),5*NI*NJ*NK,out);
    fclose(out);

  }

	free(x);
	free(y);
	free(z);
	free(U);

	return(0);
}
