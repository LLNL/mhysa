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

  double Lx = 10.0; /* length of domain along x */
  double Ly = 10.0; /* length of domain along y */

  double xmin = -2.5; /* left x-boundary */
  double ymin = -5.0; /* left y-boundary */

  double xs = -2.0; /* location of initial shock */
  double Ms = 3.0;  /* shock Mach number */

	int i,j,k;
	double *x, *y, *z;
  double dx, dy, dz;

 	x   = (double*) calloc (NI, sizeof(double));
	y   = (double*) calloc (NJ, sizeof(double));
	z   = (double*) calloc (NK, sizeof(double));

  dx = Lx / ((double)(NI-1));
  for (i=0; i<NI; i++) x[i] = xmin + i*dx;
  dy = Ly / ((double)(NJ-1));
  for (j=0; j<NJ; j++) y[j] = ymin + j*dy;
	dz = 0.5 * (dx +dy);
  for (k = 0; k < NK; k++) z[k] = (k-NK/2) * dz;


  /* pre-shock conditions */
  double rho_inf  = 1.0;
  double p_inf    = 1.0;
  double u_inf    = 0.0; //- Ms * sqrt(gamma*p_inf/rho_inf);
  double v_inf    = 0.0;
  double w_inf    = 0.0;
  
  /* post-shock conditions */
  double rho_ps   = rho_inf * ((gamma+1.0)*Ms*Ms) / ((gamma-1.0)*Ms*Ms+2.0);
  double p_ps     = p_inf * (2.0*gamma*Ms*Ms-(gamma-1.0)) / (gamma+1.0);
  double M2       = sqrt (((gamma-1.0)*Ms*Ms+2.0) / (2.0*gamma*Ms*Ms-(gamma-1.0)));
  double u_ps     = - (M2 * sqrt(gamma*p_ps/rho_ps) - Ms * sqrt(gamma*p_inf/rho_inf));
  double v_ps     = 0.0;
  double w_ps     = 0.0;

  printf("Pre-shock conditions:-\n");
  printf("Density    = %1.16E\n", rho_inf);
  printf("Pressure   = %1.16E\n", p_inf  );
  printf("u-velocity = %1.16E\n", u_inf  );
  printf("v-velocity = %1.16E\n", v_inf  );
  printf("w-velocity = %1.16E\n", w_inf  );
  printf("Mach number= %1.16E\n", Ms     );
  printf("Post-shock conditions:-\n");
  printf("Density    = %1.16E\n", rho_ps);
  printf("Pressure   = %1.16E\n", p_ps  );
  printf("u-velocity = %1.16E\n", u_ps  );
  printf("v-velocity = %1.16E\n", v_ps  );
  printf("w-velocity = %1.16E\n", w_ps  );
  printf("Mach number= %1.16E\n", M2     );
  printf("Ratios:-\n");
  printf("p2/p1      = %1.16E\n",p_ps/p_inf);
  printf("rho2/rho1  = %1.16E\n",rho_ps/rho_inf);

	double *U = (double*) calloc (5*NI*NJ*NK, sizeof(double));
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
        int p = i + NI*j + NI*NJ*k;

        double rho, u, v, w, P;
        if  (x[i] > xs) {
          rho = rho_inf;
          u   = u_inf;
          v   = v_inf;
          w   = w_inf;
          P   = p_inf;
        } else {
          rho = rho_ps;
          u   = u_ps;
          v   = v_ps;
          w   = w_ps;
          P   = p_ps;
        }

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
