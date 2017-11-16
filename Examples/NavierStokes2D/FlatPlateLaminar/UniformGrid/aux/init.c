#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double GAMMA  = 1.4;

int main()
{  
	int NI,NJ,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in= fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if      (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        } else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D solution\n");
    return(0);
  }
	printf("Grid:\t\t\t%d x %d\n",NI,NJ);

  double xmin = -0.25;
  double xmax =  1.0;
  double ymin =  0.0;
  double ymax =  0.25;

  double Lx = xmax - xmin;
  double Ly = ymax - ymin;

	int i,j;
	double dx = Lx / ((double)(NI-1));
	double dy = Ly / ((double)(NJ-1));

	double *x, *y, *U;
  FILE *out;

 	x   = (double*) calloc (NI     , sizeof(double));
	y   = (double*) calloc (NJ     , sizeof(double));
	U   = (double*) calloc (4*NI*NJ, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  x[i] = xmin + i*dx;
	    y[j] = ymin + j*dy;
      int p = i + NI*j;

      double rho, u, v, P;
      rho = 1.0;
      P   = 1.0/1.4;
      u   = 0.3;
      v   = 0.0;

      U[4*p+0] = rho;
      U[4*p+1] = rho*u;
      U[4*p+2] = rho*v;
      U[4*p+3] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v);
	  }
	}
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
  	out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = i + NI*j;
        fprintf(out,"%1.16E ",U[4*p+0]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = i + NI*j;
        fprintf(out,"%1.16E ",U[4*p+1]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = i + NI*j;
        fprintf(out,"%1.16E ",U[4*p+2]);
      }
    }
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = i + NI*j;
        fprintf(out,"%1.16E ",U[4*p+3]);
      }
    }
    fprintf(out,"\n");
	  fclose(out);

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(U,sizeof(double),4*NI*NJ,out);
    fclose(out);

  }

	free(x);
	free(y);
	free(U);

	return(0);
}
