#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double power(double x,double a)
{
  return(exp(a*log(x)));
}

int main()
{
  
  const double pi     = 4.0*atan(1.0);
  const double GAMMA  = 1.4;

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
	double dx = 10.0 / ((double)NI);
	double dy = 10.0 / ((double)NJ);
	double dz = (dx < dy ? dx : dy);

  double u_inf = 0.5;
  double v_inf = 0.0;
  double b = u_inf;
  double x0 = 5.0, y0 = 5.0;

	double *x, *y, *z, *U;
  FILE *out;

 	x   = (double*) calloc (NI        , sizeof(double));
	y   = (double*) calloc (NJ        , sizeof(double));
	z   = (double*) calloc (NK        , sizeof(double));
	U   = (double*) calloc (5*NI*NJ*NK, sizeof(double));

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){

  	  	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;

      }
	  }
	}

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
    	for (k = 0; k < NK; k++){

	    	x[i] = i*dx;
	    	y[j] = j*dy;
	    	z[k] = k*dz;

        int p = i + NI*j + NI*NJ*k;

        double rsq = (x[i]-x0)*(x[i]-x0) + (y[j]-y0)*(y[j]-y0);
        double rho, u, v, w, P;
        double du, dv;
        rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
        P   = power(rho,GAMMA);
        du  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * (y[j]-y0);
        dv  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * (x[i]-x0);
        u   = u_inf + du;
        v   = v_inf + dv;
        w   = 0.0;

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

  tf = (double)n_iter * dt; double tff = tf;
  while (tf > 20) tf -= 20; // Time period

  x0 = 5.0+tf*u_inf, y0 = 5.0;
  if (x0 > 10) x0 -= 10; //periodic domain
  printf("Final time: %lf, Vortex center: %lf, %lf\n",tff,x0,y0);

	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){

        int p = i + NI*j + NI*NJ*k;

        double rx, ry;
        rx = (x[i] - x0);
        ry = (y[j] - y0);
        if (rx < -5)      { rx += 10; }
        else if (rx > 5)  { rx -= 10; }

        double rsq = rx*rx + ry*ry;
        double rho, u, v, w, P;
        double du, dv;
        rho = power(1.0 - ((GAMMA-1.0)*b*b)/(8.0*GAMMA*pi*pi) * exp(1.0-rsq), 1.0/(GAMMA-1.0));
        P   = power(rho,GAMMA);
        du  = - b/(2.0*pi) * exp(0.5*(1.0-rsq)) * ry;
        dv  =   b/(2.0*pi) * exp(0.5*(1.0-rsq)) * rx;
        u   = u_inf + du;
        v   = v_inf + dv;
        w   = 0.0;

        U[5*p+0] = rho;
        U[5*p+1] = rho*u;
        U[5*p+2] = rho*v;
        U[5*p+3] = rho*w;
        U[5*p+4] = P/(GAMMA-1.0) + 0.5*rho*(u*u+v*v+w*w);
      }
	  }
	}
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII exact solution file exact.inp\n");
  	out = fopen("exact.inp","w");
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

    printf("Writing binary exact solution file exact.inp\n");
  	out = fopen("exact.inp","wb");
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
