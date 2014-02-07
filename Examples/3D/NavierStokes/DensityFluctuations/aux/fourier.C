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
  int i,j,k,NI,NJ,NK;
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
        if (!strcmp(word, "size")) in >> NI >> NJ >> NK;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  std::cout << "Grid size: " << NI << " x " << NJ << " x " << NK << "\n";

  int N = NI, N3 = N*N*N;

  std::ifstream solution;
  double *u;
  
  u = (double*) calloc (N3, sizeof(double));
  solution.open("op.dat");
  if (solution) {
    for (k = 0; k < N; k++) {
      for (j = 0; j < N; j++) {
        for (i = 0; i < N; i++) {
          double temp; int ii, jj, kk;
          solution >> ii >> jj >> kk >> temp >> temp >> temp >> u[k+N*(j+N*i)] >> temp >> temp >> temp >> temp;
//          if ((ii != i) || (jj != j) || (kk != k)) printf ("Error in reading file (%d,%d,%d) != (%d,%d,%d)!\n",ii,jj,kk,i,j,k);
        }
      }
    }
    solution.close();
    fourier_analysis(N,u,1);
  } else std::cout << "File not found: op.dat\n";
  free(u);
 
  u = (double*) calloc (N3, sizeof(double));
  solution.open("op_exact.dat");
  if (solution) {
    for (k = 0; k < N; k++) {
      for (j = 0; j < N; j++) {
        for (i = 0; i < N; i++) {
          double temp; int ii, jj, kk;
          solution >> ii >> jj >> kk >> temp >> temp >> temp >> u[k+N*(j+N*i)] >> temp >> temp >> temp >> temp;
          if ((ii != i) || (jj != j) || (kk != k)) printf ("Error in reading file (%d,%d,%d) != (%d,%d,%d)!\n",ii,jj,kk,i,j,k);
        }
      }
    }
    solution.close();
    fourier_analysis(N,u,0);
  } else std::cout << "File not found: op_exact.dat\n";
  free(u);
 
  return(0);
}

void fourier_analysis(int N, double *u, int flag)
{
	double PI = 4.0*atan(1.0);
	int i,j,k;
	fftw_complex *in, *out;
	fftw_plan p;
  int N3 = N*N*N;

	in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3);
	p = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i < N; i++){
	  for (j = 0; j < N; j++){
	    for (k = 0; k < N; k++){
		    in[(i*N+j)*N+k][0] = u[(i*N+j)*N+k];
		    in[(i*N+j)*N+k][1] = 0.0;
      }
    }
	}

	fftw_execute(p);

	for (i = 0; i < N3; i++){
		out[i][0] /= N3;
		out[i][1] /= N3;
	}

  std::ofstream out_file;
	if (flag == 0){
		out_file.open("spectrum_exact.dat");
	}else{
		out_file.open("spectrum_final.dat");
	}
  int kkmax = 3*(N*N/4) + 1;
  double *Eng = (double*) calloc (kkmax, sizeof(double));
  double *frq = (double*) calloc (kkmax, sizeof(double));
  double *cnt = (double*) calloc (kkmax, sizeof(double));
  for (i = 1; i < N/2; i++) {
    for (j = 1; j < N/2; j++) {
      for (k = 1; k < N/2; k++) {
        int f2 = i*i + j*j + k*k;
        Eng[f2] += (out[(i*N+j)*N+k][0]*out[(i*N+j)*N+k][0] + out[(i*N+j)*N+k][1]*out[(i*N+j)*N+k][1]);
        frq[f2] = sqrt((double)f2);
        cnt[f2] += 1.0;
      }
    }
  }
  for (i=0; i < kkmax; i++) {
    if (Eng[i] > 0) out_file << frq[i] << "\t" << Eng[i]/cnt[i] << "\n";
  }
  free(Eng);
  free(frq);
  free(cnt);
	out_file.close();
                                                                                                                                                        
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}
