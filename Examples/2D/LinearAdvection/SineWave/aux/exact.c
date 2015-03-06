#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{  
  double pi = 4.0*atan(1.0);
	int NI,NJ,ndims,n_iter;
  double tf, dt;
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    fprintf(stderr,"Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
        } else if (!strcmp(word, "n_iter")) fscanf(in,"%d",&n_iter);
        else if (!strcmp(word, "dt")) fscanf(in,"%lf",&dt);
      }
    } else {
      fprintf(stderr,"Error: Illegal format in solver.inp. Crash and burn!\n");
      return(0);
    }
  }
  fclose(in);
  if (ndims != 2) {
    fprintf(stderr,"ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
	printf("Grid: %d, %d\n",NI,NJ);

  double ax, ay;
  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    fprintf(stderr,"Error: Input file \"physics.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "advection")) {
          fscanf(in,"%lf",&ax);
          fscanf(in,"%lf",&ay);
        } 
      }
    } else {
      fprintf(stderr,"Error: Illegal format in physics.inp. Crash and burn!\n");
      return(0);
    }
  }
  fclose(in);
	printf("Advection: %3.1f, %3.1f\n",ax,ay);
	
  int i,j;
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);

  tf = (double)n_iter * dt;
  printf("dt: %lf, n_iter: %d, Final time: %lf\n",dt,n_iter,tf);
	
  double *x, *y, *u;
	x = (double*) calloc (NI   , sizeof(double));
	y = (double*) calloc (NJ   , sizeof(double));
	u = (double*) calloc (NI*NJ, sizeof(double));

  FILE *out;
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NJ*i + j;
		  u[p] = sin(2*pi*(x[i]-ax*tf)) * cos(2*pi*(y[j]-ay*tf));
	  }
	}
	out = fopen("exact.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16e ",u[p]);
    }
  }
  fprintf(out,"\n");
	fclose(out);
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = i*dx;
	  	y[j] = j*dy;
      int p = NJ*i + j;
		  u[p] = sin(2*pi*x[i]) * cos(2*pi*y[j]);
	  }
	}
	out = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
  fprintf(out,"\n");
  for (j = 0; j < NJ; j++)	{
	  for (i = 0; i < NI; i++)	{
      int p = NJ*i + j;
      fprintf(out,"%1.16e ",u[p]);
    }
  }
  fprintf(out,"\n");
	fclose(out);


	free(x);
	free(y);
	free(u);

	return(0);
}
