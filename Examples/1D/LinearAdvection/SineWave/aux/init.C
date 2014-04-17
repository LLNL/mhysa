#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 4.0*atan(1.0);

int main(){
  
	int NI,ndims;
  char ip_file_type[50];
  strcpy(ip_file_type,"ascii");
  std::ifstream in;
  std::cout << "Reading file \"solver.inp\"...\n";
  in.open("solver.inp");
  if (!in) {
    std::cout << "Error: Input file \"solver.inp\" not found. Default values will be used.\n";
  } else {
    char word[500];
    in >> word;
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        in >> word;
        if (!strcmp(word, "ndims"))             in >> ndims;
        else if (!strcmp(word, "size"))         in >> NI;
        else if (!strcmp(word, "ip_file_type")) in >> ip_file_type;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 1) {
    std::cout << "ndims is not 1 in solver.inp. this code is to generate 1D initial conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << "\n";

	int i;
	double dx = 1.0 / ((double)NI);

	double *x, *u;
	x = (double*) calloc (NI, sizeof(double));
	u = (double*) calloc (NI, sizeof(double));

	for (i = 0; i < NI; i++){
		x[i] = i*dx;
		u[i] = sin(2*pi*x[i]);
	}

  FILE *out;

  if (!strcmp(ip_file_type,"ascii")) {
    printf("Writing ASCII initial solution file initial.inp\n");
  	out = fopen("initial.inp","w");
    for (i = 0; i < NI; i++)  fprintf(out,"%lf ",x[i]);
    fprintf(out,"\n");
	  for (i = 0; i < NI; i++)	fprintf(out,"%lf ",u[i]);						
    fprintf(out,"\n");
    fclose(out);
  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {
    printf("Writing binary initial solution file initial.inp\n");
    out = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),NI,out);
    fwrite(u,sizeof(double),NI,out);
    fclose(out);
  }

	free(x);
	free(u);

	return(0);
}
