#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

double abs(double x) {
  return (x < 0 ? -x : x);
}
void fourier_analysis(int,int,double*,int);

int main()
{
  int i,NI;
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
        if (!strcmp(word, "size")) in >> NI;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  std::cout << "Grid size: " << NI << "\n";

  std::ifstream solution;
  double *u;
  
  u = (double*) calloc (NI, sizeof(double));
  solution.open("op_exact.dat");
  for (i = 0; i < NI; i++) {
    double temp; int ii;
    solution >> ii >> temp >> u[i] >> temp >> temp;
    if (ii != i) printf ("Error in reading file (%d,%d)!\n",ii,i);
  }
  solution.close();
  fourier_analysis(NI,0,u,0);
  free(u);
  
  u = (double*) calloc (NI, sizeof(double));
  solution.open("op.dat");
  for (i = 0; i < NI; i++) {
    double temp; int ii;
    solution >> ii >> temp >> u[i] >> temp >> temp;
    if (ii != i) printf ("Error in reading file (%d,%d)!\n",ii,i);
  }
  solution.close();
  fourier_analysis(NI,0,u,1);
  free(u);
}

void fourier_analysis(int N, int G, double *u, int flag){
	double PI = 4.0*atan(1.0);
	int i;
	fftw_complex *in, *out;
	fftw_plan p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i < N; i++){
		in[i][0] = u[i+G];
		in[i][1] = 0;
	}

	fftw_execute(p);
	for (i = 0; i < N; i++){
		out[i][0] /= (N);
		out[i][1] /= (N);
	}

	std::ofstream out_file;
	if (flag == 0){
		out_file.open("spectrum_exact.dat");
	}else{
		out_file.open("spectrum_final.dat");
	}
	for (i = 1; i < N; i++){
		double term = 0;
		double phase = 0;
		double xx, yy;
		xx = abs(out[i][0]);
		yy = abs(out[i][1]);

		double small = 1e-5;
		if (xx < small){
			if (yy < small)	phase = 0;
			else		phase = PI/2.0;
		}else{
			phase = atan2(yy,xx);
		}
		
		term += (out[i][0]*out[i][0]     +  out[i][1]*out[i][1]);
		term += (out[N-i][0]*out[N-i][0] +  out[N-i][1]*out[N-i][1]);
		out_file << i << "\t" << term << "\t" << phase << "\n";
	}
	out_file.close();

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}
