#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
  
  const double pi     = 4.0*atan(1.0);
  const double GAMMA  = 1.4;

  const int nvars = 6;

	int NI, NJ, NK, ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  int n_iter;
  double tf, dt;

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
        } else if (!strcmp(word, "n_iter" )) {
          fscanf(in,"%d",&n_iter);
        } else if (!strcmp(word, "dt" )) {
          fscanf(in,"%lf",&dt);
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
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);
	double dz = 1.0 / ((double)NK);

  tf = (double)n_iter * dt;
  printf("Final time: %lf\n",tf);

  double u_inf = 1.0;
  double v_inf = 1.0;
  double w_inf = 1.0;
  double rho_inf = 1.0, drho = 0.1;
  double p_inf = 1.0/GAMMA;

  double mass_fraction_1 = 0.3,
         mass_fraction_2 = 1-mass_fraction_1;

	double *x, *y, *z, *U;
  FILE *out;

 	x   = (double*) calloc (NI        , sizeof(double));
	y   = (double*) calloc (NJ        , sizeof(double));
	z   = (double*) calloc (NK        , sizeof(double));
	U   = (double*) calloc (6*NI*NJ*NK, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
  	  	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;
        int p = i + NK*j + NI*NJ*k;

        double rho, rho1, rho2, u, v, w, P;
        rho = rho_inf + drho * sin(2*pi*x[i]) * sin(2*pi*y[j]) * sin(2*pi*z[k]);
        rho1= mass_fraction_1 * rho;
        rho2= mass_fraction_2 * rho;
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        U[nvars*p+0] = rho1;
        U[nvars*p+1] = rho2;
        U[nvars*p+2] = rho*u;
        U[nvars*p+3] = rho*v;
        U[nvars*p+4] = rho*w;
        U[nvars*p+5] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
  if (!strcmp(ip_file_type,"ascii")) {

    printf("Error: ascii output not available. Set ip_file_type to \"binary\" in solver.inp\n");

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(z,sizeof(double),NK,out);
    fwrite(U,sizeof(double),nvars*NI*NJ*NK,out);
    fclose(out);

  }

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
  	  	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;
        int p = i + NK*j + NI*NJ*k;

        double rho, rho1, rho2, u, v, w, P;
        rho = rho_inf + drho * sin(2*pi*x[i]-2*pi*u_inf*tf) * sin(2*pi*y[j]-2*pi*v_inf*tf) * sin(2*pi*z[k]-2*pi*w_inf*tf);
        rho1= mass_fraction_1 * rho;
        rho2= mass_fraction_2 * rho;
        P   = p_inf;
        u   = u_inf;
        v   = v_inf;
        w   = w_inf;

        U[nvars*p+0] = rho1;
        U[nvars*p+1] = rho2;
        U[nvars*p+2] = rho*u;
        U[nvars*p+3] = rho*v;
        U[nvars*p+4] = rho*w;
        U[nvars*p+5] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
  if (!strcmp(ip_file_type,"ascii")) {

    printf("Error: ascii output not available. Set ip_file_type to \"binary\" in solver.inp\n");

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary exact solution file exact.inp\n");
  	out = fopen("exact.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(z,sizeof(double),NK,out);
    fwrite(U,sizeof(double),nvars*NI*NJ*NK,out);
    fclose(out);

  }

	free(x);
	free(y);
	free(z);
	free(U);

	return(0);
}
