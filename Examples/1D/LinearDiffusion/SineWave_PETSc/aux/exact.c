#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(){
  
  double pi = 4.0*atan(1.0);
	int NI,ndims,niter;
  double nu, dt,final_time;
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
        if (!strcmp(word, "ndims"))       fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size"))   fscanf(in,"%d",&NI);
        else if (!strcmp(word, "dt"))     fscanf(in,"%lf",&dt);
        else if (!strcmp(word, "n_iter")) fscanf(in,"%d",&niter);
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

  final_time = (double)niter * dt;
	printf("Final Time:\t\t\t%lf\n",final_time);

  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    printf("Error: Input file \"physics.inp\" not found. Default values will be used.\n");
    nu = 1.0;
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word,"diffusion")) fscanf(in,"%lf",&nu);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
	printf("Diffusion Coeff:\t\t\t%lf\n",nu);

	int i;
	double dx = 1.0 / ((double)NI);

	double *x, *u;
	x = (double*) calloc (NI, sizeof(double));
	u = (double*) calloc (NI, sizeof(double));
  /* set grid */
	for (i = 0; i < NI; i++) x[i] = i*dx;

  /* Initial  solution */
	for (i = 0; i < NI; i++) u[i] = sin(2*pi*x[i]);
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

  /* Exact solution */
	for (i = 0; i < NI; i++) u[i] = exp(-nu*4*pi*pi*final_time)*sin(2*pi*x[i]);
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII exact solution file exact.inp\n");
	  out = fopen("exact.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",u[i]);						
    fprintf(out,"\n");
	  fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary exact solution file exact.inp\n");
    out = fopen("exact.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(u,sizeof(double),NI,out);
    fclose(out);
  }

	free(x);
	free(u);

	return(0);
}
