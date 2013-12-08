#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

double abs(double x) {
  return (x < 0 ? -x : x);
}
void fourier_analysis(int,double*,int);

int main()
{
  int i,j,NI,NJ;
  std::ifstream in;
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
        if (!strcmp(word, "size")) in >> NI >> NJ;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  std::cout << "Grid size: " << NI << " x " << NJ << "\n";

  if (NI != NJ) {
    printf("Error: NI != NJ.\n");
    return(0);
  }
  int N = NI, N2 = N*N;

  std::ifstream solution;
  double *u;
  
  u = (double*) calloc (N2, sizeof(double));
  solution.open("op.dat");
  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      double temp; int ii, jj;
      solution >> ii >> jj >> temp >> temp >> u[j+N*i] >> temp >> temp >> temp;
      if ((ii != i) || (jj != j)) printf ("Error in reading file (%d,%d) != (%d,%d)!\n",ii,jj,i,j);
    }
  }
  solution.close();
  fourier_analysis(N,u,1);
  free(u);
 
  u = (double*) calloc (N2, sizeof(double));
  solution.open("op_exact.dat");
  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      double temp; int ii, jj;
      solution >> ii >> jj >> temp >> temp >> u[j+N*i] >> temp >> temp >> temp;
      if ((ii != i) || (jj != j)) printf ("Error in reading file (%d,%d) != (%d,%d)!\n",ii,jj,i,j);
    }
  }
  solution.close();
  fourier_analysis(N,u,0);
  free(u);
 
  return(0);
}

void fourier_analysis(int N, double *u, int flag)
{
	double PI = 4.0*atan(1.0);
	int i,j;
	fftw_complex *in, *out;
	fftw_plan p;

	in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N*N));
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N*N));
	p = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i < N; i++){
	  for (j = 0; j < N; j++){
		  in[i*N+j][0] = u[i*N+j];
		  in[i][1] = 0;
    }
	}

	fftw_execute(p);

	for (i = 0; i < N*N; i++){
		out[i][0] /= (N*N);
		out[i][1] /= (N*N);
	}

  std::ofstream out_file;
	if (flag == 0){
		out_file.open("spectrum_exact.dat");
	}else{
		out_file.open("spectrum_final.dat");
	}
	for (i = 1; i < N/2; i++){
		for (j = 1; j < i; j++){
			int isq = i*i, jsq = j*j;
			double kk = sqrt(isq + jsq);
			double Eng = 0;  
			Eng += (out[(j+N*i)][0]*out[(j+N*i)][0] + out[(j+N*i)][1]*out[(j+N*i)][1]);
			Eng += (out[(i+N*j)][0]*out[(i+N*j)][0] + out[(i+N*j)][1]*out[(i+N*j)][1]);
		  out_file << kk << "\t" << Eng << "\n";
		}
    for (j = i; j < i+1; j++) {
			int isq = i*i, jsq = j*j;
			double kk = sqrt(isq + jsq);
			double Eng = 0;  
			Eng += (out[(j+N*i)][0]*out[(j+N*i)][0] + out[(j+N*i)][1]*out[(j+N*i)][1]);
			Eng += (out[(i+N*j)][0]*out[(i+N*j)][0] + out[(i+N*j)][1]*out[(i+N*j)][1]);
		  out_file << kk << "\t" << Eng << "\n";
    }
	}
	out_file.close();
                                                                                                                                                        
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}
