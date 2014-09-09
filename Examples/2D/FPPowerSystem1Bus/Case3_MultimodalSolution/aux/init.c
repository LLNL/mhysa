#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main()
{
  int ierr = 0;
  FILE *out, *in;
  char ip_file_type[100];

  /* Initial condition / domain parameters */
  /* may also be read in                   */
  double xmin   = 0.0;
  double xmax   = 1.0;
  double ymin   = 0.99;
  double ymax   = 1.01;

  strcpy(ip_file_type,"ascii");

  /* position of initial dirac */
  double x0, y0;

  /* inputs to be read from solver.inp */
	int NI,NJ,ndims;

  /* reading grid inputs */
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in, "%s", word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in, "%s", word);
        if (!strcmp(word, "ndims")) {
          fscanf(in, "%d", &ndims);
        } else if (!strcmp(word, "size")) {
          fscanf(in, "%d", &NI);
          fscanf(in, "%d", &NJ);
        } else if (!strcmp(word, "ip_file_type")) {
          fscanf(in, "%s", ip_file_type);
        }
      }
    } else { 
      printf("Error: Illegal format in solver.inp. Crash and burn!\n");
      return(0);
    }
  }
  fclose(in);

  if (ndims != 2) {
    printf("ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%4d X %4d\n", NI, NJ);

  /* reading physics inputs */
  double Pm_avg, Pmax;
  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) {
    printf("Error: Input file \"physics.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in, "%s", word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in, "%s", word);
        if      (!strcmp(word, "xmin"  )) fscanf(in,"%lf",&xmin  );
        else if (!strcmp(word, "xmax"  )) fscanf(in,"%lf",&xmax  );
        else if (!strcmp(word, "ymin"  )) fscanf(in,"%lf",&ymin  );
        else if (!strcmp(word, "ymax"  )) fscanf(in,"%lf",&ymax  );
        else if (!strcmp(word, "Pm_avg")) fscanf(in,"%lf",&Pm_avg);
        else if (!strcmp(word, "Pmax"  )) fscanf(in,"%lf",&Pmax  );
      }
    } else { 
      printf("Error: Illegal format in physics.inp. Crash and burn!\n");
    }
  }
  fclose(in);

  /* calculating location of equilibrium */
  x0 = asin(Pm_avg/Pmax);
  y0 = 1.0;
  printf("Location of equilibrium: %1.16E %1.16E\n",x0,y0);

	int i,j;
	double dx = (xmax-xmin) / ((double)(NI));
	double dy = (ymax-ymin) / ((double)(NJ-1));

	double *x, *y, *u;
	x = (double*) calloc (NI   , sizeof(double));
	y = (double*) calloc (NJ   , sizeof(double));
	u = (double*) calloc (NI*NJ, sizeof(double));

  /* Generating grid and zero initial solution */
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
	  	x[i] = xmin + i*dx;
	  	y[j] = ymin + j*dy;
      int p = NJ*i + j;
      u[p] = 0;
	  }
	}

  /* Finding the nearest grid point to x0,y0 */
  printf("Finding nearest grid point to equilibrium.\n");
  double min_distance = sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin));
  int imin = -1, jmin = -1;
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
      double distance = sqrt((x[i]-x0)*(x[i]-x0) + (y[j]-y0)*(y[j]-y0));
      if (distance < min_distance) {
        min_distance = distance;
        imin = i;
        jmin = j;
      }
	  }
	}
  printf("Nearest grid point to equilibrium is (%d,%d) -> (%1.16E, %1.16E)\n",imin,jmin,x[imin],y[jmin]);
  if ((imin > -1) && (jmin > -1)) {
    x0 = x[imin];
    y0 = y[jmin];
    /* locate the dirac at the equilibrium */
    int p = NI*jmin + imin;
    u[p] = 1.0/(dx*dy);
    printf("Placing dirac at %d,%d.\n",imin,jmin);
  } else {
    printf("Error: location of dirac does not seem to be inside ");
    printf("the specified domain.\n");
    return(0);
  }

  /* calculating domain integral of initial PDF */
  double integral = 0;
	for (i = 1; i < NI-1; i++){
  	for (j = 1; j < NJ-1; j++){
      int p = NJ*i + j;
      integral += (u[p] * (dx*dy));
	  }
	}
  printf("Integral = %1.16E\n",integral);


  /* writing to file */

  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
	  out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%1.16E ",x[i]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)  fprintf(out,"%1.16E ",y[j]);
    fprintf(out,"\n");
    for (j = 0; j < NJ; j++)	{
	    for (i = 0; i < NI; i++)	{
        int p = NJ*i + j;
        fprintf(out,"%1.16E ",u[p]);
      }
    }
    fprintf(out,"\n");
	  fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
  	out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(y,sizeof(double),NJ,out);
    fwrite(u,sizeof(double),NI*NJ,out);
    fclose(out);
  } 

	free(x);
	free(y);
	free(u);

	return(0);
}
