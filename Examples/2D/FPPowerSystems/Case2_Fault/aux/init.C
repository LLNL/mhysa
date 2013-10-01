#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

int main(){
  int ierr = 0;
  FILE *out;
  std::ifstream in;

  /* Initial condition / domain parameters */
  /* may also be read in                   */
  double H   = 5.0;
  double O_s = 120*(4.0*atan(1.0));
  double D   = 5.0;
  double E   = 1.1358;
  double V   = 1.0;
  double g1  = 0.545;
  double g2  = 0.745;
  double Pm  = 0.9;
  double l   = 0.1;
  double q   = 0.0;
  double tf  = 0.1;
  double tcl = 0.2;
  double xmin =  0;
  double xmax =  pi;
  double ymin = 0.9;
  double ymax = 1.1;
  /* position of initial dirac */
  double x0, y0;

  /* inputs to be read from solver.inp */
	int NI,NJ,ndims;

  /* reading grid inputs */
  std::cout << "Reading file \"solver.inp\"...\n";
  in.open("solver.inp");
  if (!in) {
    std::cout << "Error: Input file \"solver.inp\" not found.\n";
    return(0);
  } else {
    char word[500];
    in >> word;
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        in >> word;
        if (!strcmp(word, "ndims"))     in >> ndims;
        else if (!strcmp(word, "size")) in >> NI >> NJ;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 2) {
    std::cout << "ndims is not 2 in solver.inp. this code is to generate 2D initial conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " X " << NJ << "\n";

  /* reading physics inputs */
  std::cout << "Reading file \"physics.inp\"...\n";
  in.open("physics.inp");
  if (!in) {
    std::cout << "Error: Input file \"physics.inp\" not found.\n";
    return(0);
  } else {
    char word[500];
    in >> word;
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        in >> word;
        if      (!strcmp(word, "inertia")) in >> H;
        else if (!strcmp(word, "omega_s")) in >> O_s;
        else if (!strcmp(word, "E"      )) in >> E;
        else if (!strcmp(word, "V"      )) in >> V;
        else if (!strcmp(word, "D"      )) in >> D;
        else if (!strcmp(word, "g1"     )) in >> g1;
        else if (!strcmp(word, "g2"     )) in >> g2;
        else if (!strcmp(word, "PM_min" )) in >> Pm;
        else if (!strcmp(word, "lambda" )) in >> l;
        else if (!strcmp(word, "q"      )) in >> q;
        else if (!strcmp(word, "tf"     )) in >> tf;
        else if (!strcmp(word, "tcl"    )) in >> tcl;
        else if (!strcmp(word, "xmin"   )) in >> xmin;
        else if (!strcmp(word, "xmax"   )) in >> xmax;
        else if (!strcmp(word, "ymin"   )) in >> ymin;
        else if (!strcmp(word, "ymax"   )) in >> ymax;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();

  /* calculating location of unstable equilibrium */
  double Pmax = E * V / g1;
  x0 = asin(Pm/Pmax);
  y0 = 1.0;

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
  double min_distance = sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin));
  int imin = -1, jmin = -1;
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
      double distance = (x[i]-x0)*(x[i]-x0)+(y[j]-y0)*(y[j]-y0);
      distance = sqrt(distance);
      if (distance < min_distance) {
        min_distance = distance;
        imin = i;
        jmin = j;
      }
	  }
	}
  if ((imin > -1) && (jmin > -1)) {
    /* corrected x0,y0 and Pm so that 
       the unstable equilibrium point
       lies on a grid point
    */
    x0 = x[imin];
    y0 = y[jmin];
    Pm = Pmax * sin(x0);
    /* locate the dirac here */
    int p = NJ*imin + jmin;
    u[p] = 1.0/(dx*dy);
    printf("Placing dirac at %d,%d.\n",imin,jmin);
  } else {
    printf("Error: location of dirac does not seem to be inside ");
    printf("the specified domain.\n");
    return(0);
  }

  /* rewrite physics.inp with the grid-corrected Pm */
  printf("Rewriting physics.inp\n");
  out = fopen("physics.inp","w");
  fprintf(out,"begin\n");
  fprintf(out,"\tinertia         %1.16E\n",H   );
  fprintf(out,"\tomega_s         %1.16E\n",O_s );
  fprintf(out,"\tE               %1.16E\n",E   );
  fprintf(out,"\tV               %1.16E\n",V   );
  fprintf(out,"\tD               %1.16E\n",D   );
  fprintf(out,"\tg1              %1.16E\n",g1  );
  fprintf(out,"\tg2              %1.16E\n",g2  );
  fprintf(out,"\tPM_min          %1.16E\n",Pm  );
  fprintf(out,"\tlambda          %1.16E\n",l   );
  fprintf(out,"\tq               %1.16E\n",q   );
  fprintf(out,"\ttf              %1.16E\n",tf  );
  fprintf(out,"\ttcl             %1.16E\n",tcl );
  fprintf(out,"\txmin            %1.16E\n",xmin);
  fprintf(out,"\txmax            %1.16E\n",xmax);
  fprintf(out,"\tymin            %1.16E\n",ymin);
  fprintf(out,"\tymax            %1.16E\n",ymax);
  fprintf(out,"end\n");
  fclose(out);
  
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

	free(x);
	free(y);
	free(u);

	return(0);
}
