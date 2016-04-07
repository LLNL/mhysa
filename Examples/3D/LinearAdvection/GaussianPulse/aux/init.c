#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
	int NI=64,NJ=64,NK=64,ndims=3;
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
        if (!strcmp(word, "ndims"))     fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
          fscanf(in,"%d",&NK);
        } else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  if (ndims != 3) {
    printf("ndims is not 3 in solver.inp. this code is to generate 3D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d X %d X %d\n",NI,NJ,NK);

	int i,j,k;
	double dx = 12.0 / ((double)NI);
	double dy = 6.0 / ((double)NJ);
	double dz = 6.0 / ((double)NK);

	double *x, *y, *z, *u;
	x = (double*) calloc (NI      , sizeof(double));
	y = (double*) calloc (NJ      , sizeof(double));
	z = (double*) calloc (NK      , sizeof(double));
	u = (double*) calloc (NI*NJ*NK, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
      for (k = 0; k < NK; k++) {
	  	  x[i] = -6 + i*dx;
	  	  y[j] = -3 + j*dy;
	  	  z[k] = -3 + k*dz;
        int p = i + NI*j + NI*NJ*k;
		    u[p] = exp(-((x[i]*x[i]/2+y[j]*y[j]/2+z[k]*z[k]/2)));
      }
	  }
	}

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
	  out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16e ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%1.16e ",y[j]);
    fprintf(out,"\n");
    for (k = 0; k < NK; k++)  fprintf(out,"%1.16e ",z[k]);
    fprintf(out,"\n");
    for (k = 0; k < NK; k++) {
      for (j = 0; j < NJ; j++)	{
	      for (i = 0; i < NI; i++)	{
          int p = i + NK*j + NI*NJ*k;
          fprintf(out,"%1.16e ",u[p]);
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
    fwrite(u,sizeof(double),NI*NJ*NK,out);
    fclose(out);
  }

	free(x);
	free(y);
	free(z);
	free(u);

	return(0);
}
