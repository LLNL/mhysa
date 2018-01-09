#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{
  const double GAMMA  = 1.4;

	int NI,NJ,NK,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
        } else if (!strcmp(word, "size" )) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
          fscanf(in,"%d",&NK);
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in,"%s",ip_file_type);
        }
      }
    } else {
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
      return(0);
    }
  }
  fclose(in);
  if (ndims != 3) {
    printf("ndims is not 3 in solver.inp. this code is to generate 3D exact solution\n");
    return(0);
  }
	printf("Grid:\t\t\t%d x %d x %d.\n", NI, NJ, NK);

	int i,j,k;
	double dx = 1.0 / ((double)(NI-1));
	double dy = 1.0 / ((double)(NJ-1));
	double dz = (dx < dy ? dx : dy);

	double *x, *y, *z, *U;
	x  = (double*) calloc (NI         , sizeof(double));
	y  = (double*) calloc (NJ         , sizeof(double));
	z  = (double*) calloc (NK         , sizeof(double));
	U  = (double*) calloc (5*NI*NJ*NK , sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){

  	  	x[i] = -0.5 + i*dx;
  	  	y[j] = -0.5 + j*dy;
  	  	z[k] = k*dz;

        int p = i + NI*j + NI*NJ*k;
  
        double rho, u, v, w, P;
        if ((x[i] < 0) && (y[j] < 0)) {
          rho = 1.1;
          P   = 1.1;
          u   = 0.8939;
          v   = 0.8939;
        } else if ((x[i] >= 0) && (y[j] <  0)) {
          rho = 0.5065;
          P   = 0.35;
          u   = 0.0;
          v   = 0.8939;
        } else if ((x[i] <  0) && (y[j] >= 0)) {
          rho = 0.5065;
          P   = 0.35;
          u   = 0.8939;
          v   = 0.0;
        } else if ((x[i] >= 0) && (y[j] >= 0)) {
          rho = 1.1;
          P   = 1.1;
          u   = 0.0;
          v   = 0.0;
        }
        w = 0;

        U[5*p+0] = rho;
        U[5*p+1] = rho*u;
        U[5*p+2] = rho*v;
        U[5*p+3] = rho*w;
        U[5*p+4] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
  	  }
  	}
  }
  
  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {

    printf("Error: ascii output not implemented. Set ip_file_type to \"binary\" in solver.inp.\n");

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
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
