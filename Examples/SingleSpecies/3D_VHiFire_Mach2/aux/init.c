#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main()
{
  
  const double  gamma  = 1.4;
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

  /* bounding box around body, with some margin */
  double xbmin = -0.5;
  double xbmax =  1.5;
  double ybmin =  0.0;
  double ybmax =  0.15;
  double zbmin = -0.15;
  double zbmax =  0.15;

	int i,j,k;
	double *x, *y, *z;

 	x   = (double*) calloc (NI, sizeof(double));
	y   = (double*) calloc (NJ, sizeof(double));
	z   = (double*) calloc (NK, sizeof(double));

  double dx = (xbmax - xbmin) / ((double)(NI/2));
  x[NI/8] = xbmin;
  for (i = NI/8+1; i <= 5*NI/8; i++) {
    x[i] = xbmin + dx * (i - NI/8);
  }
  for (i = 5*NI/8+1; i < NI; i++) {
    x[i] = x[i-1] + 1.02*(x[i-1]-x[i-2]);
  }
  for (i = NI/8-1; i >= 0; i--) {
    x[i] = x[i+1] - 1.001*(x[i+2]-x[i+1]);
  }
  printf("xmin: %f, xmax: %f\n", x[0], x[NI-1]);

  double dy = (ybmax - ybmin) / ((double)(NJ/2));
  y[0] = ybmin;
  for (j = 1; j <= NJ/2; j++) {
    y[j] = ybmin + dy * j;
  }
  for (j = NJ/2+1; j < NJ; j++) {
    y[j] = y[j-1] + 1.08*(y[j-1]-y[j-2]);
  }
  printf("ymin: %f, ymax: %f\n", y[0], y[NJ-1]);

  double dz = (zbmax - zbmin) / ((double) (NK/2));
  z[NK/4] = zbmin;
  for (k = NK/4+1; k <= 3*NK/4; k++) {
    z[k] = zbmin + dz * (k - NK/4);
  }
  for (k = 3*NK/4+1; k < NK; k++) {
    z[k] = z[k-1] + 1.08*(z[k-1]-z[k-2]);
  }
  for (k = NK/4-1; k >= 0; k--) {
    z[k] = z[k+1] - 1.08*(z[k+2]-z[k+1]);
  }
  printf("zmin: %f, zmax: %f\n", z[0], z[NK-1]);


  double rho_inf  = 1.0;
  double p_inf    = 1.0/gamma;
  double u_inf    = 0.0;
  double v_inf    = 0.0;
  double w_inf    = 0.0;
  
	double *U = (double*) calloc (5*NI*NJ*NK, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){

        int p = i + NI*j + NI*NJ*k;

        double rho, u, v, w, P;
        rho = rho_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;
        P   = p_inf;

        U[5*p+0] = rho;
        U[5*p+1] = rho*u;
        U[5*p+2] = rho*v;
        U[5*p+3] = rho*w;
        U[5*p+4] = P/(gamma-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}

  FILE *out;
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
