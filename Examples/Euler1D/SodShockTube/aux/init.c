#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(){
  
	int NI=101,ndims=1;
  FILE *in;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims"))     fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) fscanf(in,"%d",&NI);
        else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else { 
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
    }
  }
  fclose(in);
  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d\n",NI);

	int i;
	double dx = 1.0 / ((double)(NI-1));

	double *x, *rho,*rhou,*e;
	x    = (double*) calloc (NI, sizeof(double));
	rho  = (double*) calloc (NI, sizeof(double));
	rhou = (double*) calloc (NI, sizeof(double));
	e    = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
    double RHO,U,P;
    if (x[i] < 0.5) {
      RHO = 1.0;
      U   = 0.0;
      P   = 1.0;
    } else {
      RHO = 0.125;
      U   = 0;
      P   = 0.1;
    }
    rho[i]  = RHO;
    rhou[i] = RHO*U;
    e[i]    = P/0.4 + 0.5*RHO*U*U;
	}

  if (!strcmp(ip_file_type,"ascii")) {
    FILE *out;
	  out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",rho[i]);						
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",rhou[i]);						
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",e[i]);						
    fprintf(out,"\n");
	  fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Error: Writing binary initial solution file not implemented. ");
    printf("Please choose ip_file_type in solver.inp as \"ascii\".\n");
  }

	free(x);
	free(rho);
	free(rhou);
	free(e);

	return(0);
}
