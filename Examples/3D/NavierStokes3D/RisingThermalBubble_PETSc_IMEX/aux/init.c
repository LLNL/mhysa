/*
  This code generates the initial solution for the rising 
  thermal bubble case for the 3D Navier-Stokes equations.

  Note: this code allocates the global domain, so be careful
  when using it for a large number of grid points.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double raiseto(double x, double a)
{
  return(exp(a*log(x)));
}

int main()
{
  double gamma    = 1.4;
  double R        = 287.058;
  double rho_ref  = 1.1612055171196529;
  double p_ref    = 100000.0;
  double grav_x   = 0.0;
  double grav_y   = 9.8;
  double grav_z   = 0.0;
  int    HB       = 0;
  
	int   NI,NJ,NK,ndims;
  char  ip_file_type[50]; strcpy(ip_file_type,"ascii");
  FILE *in;

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
          fscanf(in,"%d",&NK);
        } else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);

  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    printf("Error: Input file \"physics.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if      (!strcmp(word, "rho_ref")) fscanf(in,"%lf",&rho_ref);
        else if (!strcmp(word, "p_ref"  )) fscanf(in,"%lf",&p_ref  );
        else if (!strcmp(word, "gamma"  )) fscanf(in,"%lf",&gamma  );
        else if (!strcmp(word, "R"      )) fscanf(in,"%lf",&R      );
        else if (!strcmp(word, "HB"     )) fscanf(in,"%d" ,&HB     );
        else if (!strcmp(word, "gravity")) {
          fscanf(in,"%lf",&grav_x );
          fscanf(in,"%lf",&grav_y );
          fscanf(in,"%lf",&grav_z );
        }
      }
    } else printf("Error: Illegal format in physics.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 3) {
    printf("ndims is not 3 in solver.inp. this code is to generate 3D initial conditions\n");
    return(0);
  }
  if (HB != 2) {
    printf("Error: Specify \"HB\" as 2 in physics.inp.\n");
  }
  if ((grav_x != 0.0) || (grav_z != 0.0)) {
    printf("Error: gravity must be zero along x and z for this example.\n");
    return(0);
  }
	printf("Grid:\t\t\t%d X %d X %d\n",NI,NJ,NK);
  printf("Reference density and pressure: %lf, %lf.\n",rho_ref,p_ref);

	int i,j,k;
	double dx = 1000.0  / ((double)(NI-1));
	double dy = 1000.0  / ((double)(NJ-1));
	double dz = 1000.0  / ((double)(NK-1));

	double *x, *y, *z, *U;
	x   = (double*) calloc (NI        , sizeof(double));
	y   = (double*) calloc (NJ        , sizeof(double));
	z   = (double*) calloc (NK        , sizeof(double));
	U   = (double*) calloc (5*NI*NJ*NK, sizeof(double));

  /* Initial perturbation center */
  double xc = 500;
  double yc = 260;
  double zc = 500;
  double Cp = gamma * R / (gamma-1.0);

  /* initial perturbation parameters */
  double pi = 4.0*atan(1.0);
  double rc = 250.0;
  double T_ref = p_ref / (R * rho_ref);

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
      for (k = 0; k < NK; k++){
	  	  x[i] = i*dx;
	  	  y[j] = j*dy;
        z[k] = k*dz;
        int p = i + NI*j + NI*NJ*k;

        /* temperature peturbation */
        double r      = sqrt((x[i]-xc)*(x[i]-xc)+(y[j]-yc)*(y[j]-yc)+(z[k]-zc)*(z[k]-zc));
        double dtheta = (r>rc ? 0.0 : (0.5*(1.0+cos(pi*r/rc))) );
        double theta  = T_ref + dtheta;
        double Pexner = 1.0 - (grav_y*y[j])/(Cp*T_ref);

        double rho    = (p_ref/(R*theta)) * raiseto(Pexner, (1.0/(gamma-1.0)));
        double E      = rho * (R/(gamma-1.0)) * theta*Pexner;

        U[5*p+0] = rho;
        U[5*p+1] = 0.0;
        U[5*p+2] = 0.0;
        U[5*p+3] = 0.0;
        U[5*p+4] = E;
      }
	  }
	}

  FILE *out;
  if (!strcmp(ip_file_type,"ascii")) {

    printf("Error: sorry, ASCII format initial solution file not implemented. ");
    printf("Please choose \"ip_file_format\" in solver.inp as \"binary\".\n");
    return(0);

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
