/*
  Code to generate the initial solution for:
  Case: Dam Breaking over Rectangular Bump
  Model: ShallowWater1D

  Reference:
  Xing, Y., Shu, C.-W., "High order finite difference WENO 
  schemes with the exact conservation property for the shallow 
  water equations", Journal of Computational Physics, 208, 2005, 
  pp. 206-227. http://dx.doi.org/10.1016/j.jcp.2005.02.006
  Section 4.4

  Needs: solver.inp
  Writes out: initial.inp, topography.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double absolute(double x)
{
  return( x<0 ? -x : x);
}

int main() {
	int     NI,ndims;
  FILE    *in, *out;
  char    ip_file_type[50]; 

  /* default values */
  NI    = 500;
  ndims = 1;
  strcpy(ip_file_type,"ascii");

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if      (!strcmp(word, "ndims"))        fscanf(in,"%d" , &ndims);
        else if (!strcmp(word, "size"))         fscanf(in,"%d" ,    &NI);
        else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 1) {
    printf("ndims is not 1 in solver.inp. Make sure the correct solver.inp is being used.\n");
    return(0);
  }
	printf("Grid: %d\n", NI);

	int i;
  double length = 1500.0;
	double dx = length / ((double)NI);

	double *x,*h,*hu,*b;
	x  = (double*) calloc (NI, sizeof(double));
	h  = (double*) calloc (NI, sizeof(double));
	hu = (double*) calloc (NI, sizeof(double));
	b  = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
    if ( absolute(x[i]-750.0) <= 1500.0/8.0 ) b[i] = 8.0;
    else                                      b[i] = 0.0;
    if (x[i] <= 750.0)  h[i] = 20.0-b[i];
    else                h[i] = 15.0-b[i];
    hu[i] = 0;
	}

  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
	  out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",h[i]);						
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",hu[i]);						
    fprintf(out,"\n");
	  fclose(out);
    printf("Writing ASCII topography file topography.inp\n");
	  out = fopen("topography.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%1.16E ",b[i]);						
    fprintf(out,"\n");
	  fclose(out);
  } else {
    printf("Writing binary initial solution file initial.inp\n");
	  out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    double *U = (double*) calloc (2*NI,sizeof(double));
    for (i=0; i < NI; i++) {
      U[2*i+0] = h [i];
      U[2*i+1] = hu[i];
    }
    fwrite(U,sizeof(double),2*NI,out);
    free(U);
	  fclose(out);
    printf("Writing binary topography file topography.inp\n");
	  out = fopen("topography.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(b,sizeof(double),NI,out);
	  fclose(out);
  }

	free(x);
	free(h);
	free(hu);
	free(b);

	return(0);
}
