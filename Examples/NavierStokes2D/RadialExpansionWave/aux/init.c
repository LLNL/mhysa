#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double power(double x,double a)
{
  return(exp(a*log(x)));
}

int main(){
  const double pi     = 4.0*atan(1.0);
  const double GAMMA  = 1.4;
	int NI=101,NJ=101,ndims=2;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
  FILE *in, *out;

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims"))     fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        } else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d X %d\n",NI,NJ);

	int i,j;
	double dx = 8.0 / ((double)(NI-1));
	double dy = 8.0 / ((double)(NJ-1));

	double *x, *y, *u0, *u1, *u2, *u3;
	x   = (double*) calloc (NI   , sizeof(double));
	y   = (double*) calloc (NJ   , sizeof(double));
	u0  = (double*) calloc (NI*NJ, sizeof(double));
	u1  = (double*) calloc (NI*NJ, sizeof(double));
	u2  = (double*) calloc (NI*NJ, sizeof(double));
	u3  = (double*) calloc (NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = -4.0 + i*dx;
	  	y[j] = -4.0 + j*dy;
      int p = NJ*i + j;
      double r = sqrt(x[i]*x[i] + y[j]*y[j]);

      double q;
      if      (r < 0.5)  q = 0.0;
      else if (r < 1.5)  q = (1.0/GAMMA) * (1.0 + tanh((r-1)/(0.25-(r-1)*(r-1))));
      else               q = 2.0/GAMMA;

      double rho, u, v, P, a;

      if (r == 0) u = v = 0;
      else {
        u = x[i]*q/r;
        v = y[j]*q/r;
      }

      a = 1.0 - (GAMMA-1.0)*q/2.0;

      rho = GAMMA*power(a,2/(GAMMA-1.0));
      P   = rho*a*a/GAMMA;

      u0[p] = rho;
      u1[p] = rho*u;
      u2[p] = rho*v;
      u3[p] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
	  }
	}

  if (!strcmp(ip_file_type,"ascii")) {
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
    printf("Error: Writing binary initial solution file not implemented. ");
    printf("Please choose ip_file_type in solver.inp as \"ascii\".\n");
  }

	free(x);
	free(y);
	free(u0);
	free(u1);
	free(u2);
	free(u3);

	return(0);
}
