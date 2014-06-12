#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double raiseto(double x, double a)
{
  return(exp(a*log(x)));
}

int main(){
  
	int NI,NK,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");

  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) printf ("Error: Input file \"solver.inp\" not found. Default values will be used.\n");
  else {
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

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D exact solution\n");
    return(0);
  }
	printf("Grid:\t\t\t%d x %d\n", NI, NK);

  /* read in reference pressure and temperature */
  double P0    = 100000.0;
  double T0    = 300.0;
  double g     = 9.8;
  double R     = 287.058;
  double gamma = 1.4;
  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) printf ("Error: Input file \"physics.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if      (!strcmp(word, "Pref" )) fscanf(in,"%lf",&P0);
        else if (!strcmp(word, "Tref" )) fscanf(in,"%lf",&T0);
        else if (!strcmp(word, "gamma")) fscanf(in,"%lf",&gamma);
        else if (!strcmp(word, "g"    )) fscanf(in,"%lf",&g);
        else if (!strcmp(word, "R"    )) fscanf(in,"%lf",&R);
      }
    } else printf("Error: Illegal format in physics.inp. Crash and burn!\n");
  }
  fclose(in);
  double rho0 = P0 / (R*T0);
  double inv_gamma_m1 = 1.0 / (gamma-1.0);
  double Cp = gamma * inv_gamma_m1 * R;

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

  /* initial perturbation parameters */
  double pi = 4.0*atan(1.0);
  double tc = 0.01;
  double hc = 10000;
  double ac = 5000;
  double xc = 100000;
  double uc = 20.0;
  double BV = 0.01;

	double *x, *z, *U;
  FILE *out;

 	x   = (double*) calloc (NI     , sizeof(double));
	z   = (double*) calloc (NK     , sizeof(double));
	U   = (double*) calloc (4*NI*NK, sizeof(double));

	for (i = 0; i < NI; i++){
    for (k = 0; k < NK; k++){
    	x[i] = xmin + i*dx;
	   	z[k] = zmin + k*dz;
      int p = i + NI*k;

      /* temperature peturbation */
      double dtheta = tc * sin(pi*z[k]/hc) / (1.0 + ((x[i]-xc)/ac)*((x[i]-xc)/ac));

      double theta      = T0*exp(BV*BV*z[k]/g) + dtheta;
      double Pexner     = 1.0 + ((g*g)/(Cp*T0*BV*BV))*(exp(-BV*BV*z[k]/g)-1.0);
      double rho        = (P0/(R*theta)) * raiseto(Pexner,inv_gamma_m1);

      double theta_ref  = T0*exp(BV*BV*z[k]/g);
      double Pexner_ref = 1.0 + ((g*g)/(Cp*T0*BV*BV))*(exp(-BV*BV*z[k]/g)-1.0);
      double rho_ref    = (P0/(R*theta_ref)) * raiseto(Pexner,inv_gamma_m1);

      U[4*p+0] = rho - rho_ref;
      U[4*p+1] = rho * uc;
      U[4*p+2] = 0.0;
      U[4*p+3] = (rho*theta) - (rho_ref*theta_ref);
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
