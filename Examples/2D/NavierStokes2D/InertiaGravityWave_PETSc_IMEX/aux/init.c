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
  int    HB       = 0;
  double BV       = 0.0;
  
	int NI,NK,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf ("Error: Input file \"solver.inp\" not found.\n");
    return(1);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
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
        else if (!strcmp(word, "HB"     )) {
          fscanf(in,"%d" ,&HB     );
          if (HB == 3) fscanf(in, "%lf", &BV);
        } else if (!strcmp(word, "gravity")) {
          fscanf(in,"%lf",&grav_x );
          fscanf(in,"%lf",&grav_y );
        }
      }
    } else printf("Error: Illegal format in physics.inp. Crash and burn!\n");
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
  if (HB != 3) {
    printf("Error: Specify \"HB\" as 3 in physics.inp.\n");
  }
	if (grav_x != 0.0) {
    printf("Error: Gravity force along x must be zero for HB = 3.\n");
    return(0);
  }
  printf("Grid:\t\t\t%d X %d\n",NI,NK);
  printf("Reference density and pressure: %lf, %lf.\n",rho_ref,p_ref);

  /* Define the domain */
  double xmin, xmax, zmin, zmax;
  xmin    =  0.0;
  xmax    =  300000;
  zmin    =  0.0;
  zmax    =  10000.0 ;

  double Lx = xmax - xmin;
  double Lz = zmax - zmin;

	int i,k;
	double dx = Lx / ((double)NI-1);
	double dz = Lz / ((double)NK-1);

	double *x, *z, *U;
  FILE *out;

 	x   = (double*) calloc (NI     , sizeof(double));
	z   = (double*) calloc (NK     , sizeof(double));
	U   = (double*) calloc (4*NI*NK, sizeof(double));

  double inv_gamma_m1 = 1.0 / (gamma-1.0);
  double Cp = gamma * R * inv_gamma_m1;
  double Cv =         R * inv_gamma_m1;
  double T_ref = p_ref / (R*rho_ref);

  /* initial perturbation parameters */
  double pi = 4.0*atan(1.0);
  double tc = 0.01;
  double hc = 10000;
  double ac = 5000;
  double xc = 100000;
  double uc = 20.0;

	for (i = 0; i < NI; i++){
    for (k = 0; k < NK; k++){
    	x[i] = xmin + i*dx;
	   	z[k] = zmin + k*dz;
      int p = i + NI*k;

      /* temperature peturbation */
      double dtheta = tc * sin(pi*z[k]/hc) / (1.0 + ((x[i]-xc)/ac)*((x[i]-xc)/ac));

      double theta  = T_ref*exp(BV*BV*z[k]/grav_y) + dtheta;
      double Pexner = 1.0 + ((grav_y*grav_y)/(Cp*T_ref*BV*BV))*(exp(-BV*BV*z[k]/grav_y)-1.0);
      double rho    = (p_ref/(R*theta)) * raiseto(Pexner,inv_gamma_m1);
      double E      = Cv * theta * Pexner + 0.5 * (uc * uc);

      U[4*p+0] = rho;
      U[4*p+1] = rho*uc;
      U[4*p+2] = 0.0;
      U[4*p+3] = rho*E;
    }
	}

  if (!strcmp(ip_file_type,"ascii")) {
    printf("ASCII not supported. Use binary format\n");
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(z,sizeof(double),NK,out);
    fwrite(U,sizeof(double),4*NI*NK,out);
    fclose(out);

  }

	free(x);
	free(z);
	free(U);

	return(0);
}
