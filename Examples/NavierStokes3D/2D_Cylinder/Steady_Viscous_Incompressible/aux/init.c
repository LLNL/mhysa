#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{
  
  const double  GAMMA  = 1.4;
	int           NI,NJ,NK,ndims;
  char          ip_file_type[50];
  FILE          *in;

  strcpy(ip_file_type,"ascii");

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
    printf("Error: ndims is not 3 in solver.inp.\n");
    return(0);
  }
  printf("Grid: %3d x %3d x %3d.\n",NI,NJ,NK);

  double Lx = 8.0;
  double Ly = 4.0;
  double sf_x1 = 1.15;
  double sf_x2 = 1.20;
  double sf_y = 1.4;

	int i,j,k;

  double  u   = 0.1,
          v   = 0.0,
          w   = 0.0, 
          rho = 1.0,
          P   = 1.0/GAMMA;

	double *x, *y, *z, *U;
  double dx, dy, dz;
  FILE *out;

 	x   = (double*) calloc (NI        , sizeof(double));
	y   = (double*) calloc (NJ        , sizeof(double));
	z   = (double*) calloc (NK        , sizeof(double));
	U   = (double*) calloc (5*NI*NJ*NK, sizeof(double));

	x[NI/8]   = -Lx/4;
	x[7*NI/8] = 3*Lx/4;
	dx = Lx / ((double)(6*NI/8)); 
	for (i = NI/8  ; i < 7*NI/8; i++) x[i] = x[NI/8] + dx * (i-NI/8);
	for (i = 7*NI/8; i < NI    ; i++)	x[i] = x[i-1] + sf_x2 * (x[i-1] - x[i-2]);
	for (i = NI/8-1; i >= 0    ; i--) x[i] = x[i+1] - sf_x1 * (x[i+2] - x[i+1]);

	y[NJ/8]   = -Ly/2;
	y[7*NJ/8] =  Ly/2;
	dy = Ly / ((double)(6*NJ/8)); 
	for (j = NJ/8  ; j < 7*NJ/8; j++) y[j] = y[NJ/8] + dy * (j-NJ/8);
	for (j = 7*NJ/8; j < NJ    ; j++)	y[j] = y[j-1] + sf_y * (y[j-1] - y[j-2]);
	for (j = NJ/8-1; j >= 0    ; j--)	y[j] = y[j+1] - sf_y * (y[j+2] - y[j+1]);

	dz = 0.5 * (dx +dy);
  for (k = 0; k < NK; k++) z[k] = (k-NK/2) * dz;

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
        int p = i + NI*j + NI*NJ*k;
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

	free(x);
	free(y);
	free(z);
	free(U);

	return(0);
}
