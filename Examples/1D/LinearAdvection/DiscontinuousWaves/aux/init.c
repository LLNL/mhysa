#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double absolute(double x)
{
  return (x<0? -x : x);
}

int main(){
  
	int NI,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims"))             fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size"))         fscanf(in,"%d",&NI);
        else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d\n",NI);

	int i;
	double dx = 2.0 / ((double)NI);

	double *x, *u;
	x = (double*) calloc (NI, sizeof(double));
	u = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = -1.0 + i*dx;
    if      (x[i] < -0.8)  u[i] = 0.0;
    else if (x[i] < -0.6)  u[i] = exp(-log(2.0)*(x[i]+0.7)*(x[i]+0.7)/0.0009);
    else if (x[i] < -0.4)  u[i] = 0.0;
    else if (x[i] < -0.2)  u[i] = 1.0;
    else if (x[i] <  0  )  u[i] = 0.0;
    else if (x[i] <  0.2)  u[i] = 1.0 - absolute(10.0*(x[i]-0.1));
    else if (x[i] <  0.4)  u[i] = 0.0;
    else if (x[i] <  0.6)  u[i] = sqrt(1.0-100*(x[i]-0.5)*(x[i]-0.5));
    else                   u[i] = 0.0;
	}

  FILE *out;

  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
  	out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",u[i]);						
    fprintf(out,"\n");
    fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
    out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(u,sizeof(double),NI,out);
    fclose(out);
  }

	free(x);
	free(u);

	return(0);
}
