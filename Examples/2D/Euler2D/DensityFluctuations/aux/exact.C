#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

double raiseto(double x, double a){
	return exp(a*log(x));
}

double magnitude(double a, double b){
	return sqrt(a*a + b*b);
}

int main(){
	const double PI = 4 * atan(1.0);
	std::ofstream out;

	int NI,NJ,ndims,niter;
  double dt;
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
        if (!strcmp(word, "ndims"))     in >> ndims;
        else if (!strcmp(word, "size")) in >> NI >> NJ;
        else if (!strcmp(word, "n_iter")) in >> niter;
        else if (!strcmp(word, "dt")) in >> dt;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 2) {
    std::cout << "ndims is not 2 in solver.inp. this code is to generate 2D exact conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " x " << NJ << "\n";
  std::cout << "Input maximum wavenumber (typically N/2): ";
  int limit; std::cin >> limit;

  if (NI != NJ) { printf("Error: NI != NJ.\n"); return(0); }
	int N = NI;
	int i,j;
	double dx = 1.0 / ((double)NI);
	double dy = 1.0 / ((double)NJ);
  double tf = ((double)niter) * dt;
  printf("Final Time: %lf\n",tf);
  double factor = 0.01;

	double kk = sqrt(2 * limit * limit);
	int kkmax = (int) kk;

	fftw_complex *uhat = (fftw_complex*) fftw_malloc(N*N * sizeof(fftw_complex));	
	fftw_complex *u = (fftw_complex*) fftw_malloc(N*N * sizeof(fftw_complex));	

	fftw_plan inv_trans_u, trans_u;
	inv_trans_u = fftw_plan_dft_2d(N, N, uhat, u, 1, FFTW_MEASURE);
	trans_u     = fftw_plan_dft_2d(N, N, u, uhat,-1, FFTW_MEASURE);

  if (limit > N/2) limit = N/2;

	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			uhat[(j+N*i)][0] = 0;
			uhat[(j+N*i)][1] = 0;
		}
	}
	for (i = 1; i <= limit; i++){
		for (j = 1; j <= limit; j++){
			double kk = sqrt(i*i + j*j);
      double Ak = factor * raiseto(kk,-5.0/6.0);
			uhat[(j+N*i)][0] = Ak/sqrt(2.0);
			uhat[(j+N*i)][1] = Ak/sqrt(2.0);
		}
	}
	for (i = 1+N/2; i < N; i++){
		for (j = 0; j <= N/2; j++){
			uhat[(j+N*i)][0] = uhat[(j+N*(N-i))][0];
			uhat[(j+N*i)][1] = -uhat[(j+N*(N-i))][1];
		}
	}
	for (i = 0; i < N; i++){
		for (j = N/2+1; j < N; j++){
			uhat[(j+N*i)][0] = uhat[((N-j)+N*i)][0];
			uhat[(j+N*i)][1] = -uhat[((N-j)+N*i)][1];
		}
	}

#if 0
	out.open("spectrum_test1.dat");
	for (i = 1; i < N/2; i++){
		for (j = 1; j < i; j++){
			int isq = i*i, jsq = j*j;
			double kk = sqrt(isq + jsq);
			double Eng = 0;  
			Eng += (uhat[(j+N*i)][0]*uhat[(j+N*i)][0] + uhat[(j+N*i)][1]*uhat[(j+N*i)][1]);
			Eng += (uhat[(i+N*j)][0]*uhat[(i+N*j)][0] + uhat[(i+N*j)][1]*uhat[(i+N*j)][1]);
		  out << kk << "\t" << Eng << "\n";
		}
    for (j = i; j < i+1; j++) {
			int isq = i*i, jsq = j*j;
			double kk = sqrt(isq + jsq);
			double Eng = 0;  
			Eng += (uhat[(j+N*i)][0]*uhat[(j+N*i)][0] + uhat[(j+N*i)][1]*uhat[(j+N*i)][1]);
			Eng += (uhat[(i+N*j)][0]*uhat[(i+N*j)][0] + uhat[(i+N*j)][1]*uhat[(i+N*j)][1]);
		  out << kk << "\t" << Eng << "\n";
    }
	}
	out.close();
#endif

	fftw_execute(inv_trans_u);

#if 0
  for (i=0; i<N*N; i++)  u[i][1] = 0.0;

	fftw_execute(trans_u);
  for (i=0; i<N*N; i++) {
    uhat[i][0] /= (N*N);
    uhat[i][1] /= (N*N);
  }

	out.open("spectrum_test2.dat");
	for (i = 1; i < N/2; i++){
		for (j = 1; j < i; j++){
			int isq = i*i, jsq = j*j;
			double kk = sqrt(isq + jsq);
			double Eng = 0;  
			Eng += (uhat[(j+N*i)][0]*uhat[(j+N*i)][0] + uhat[(j+N*i)][1]*uhat[(j+N*i)][1]);
			Eng += (uhat[(i+N*j)][0]*uhat[(i+N*j)][0] + uhat[(i+N*j)][1]*uhat[(i+N*j)][1]);
		  out << kk << "\t" << Eng << "\n";
		}
    for (j = i; j < i+1; j++) {
			int isq = i*i, jsq = j*j;
			double kk = sqrt(isq + jsq);
			double Eng = 0;  
			Eng += (uhat[(j+N*i)][0]*uhat[(j+N*i)][0] + uhat[(j+N*i)][1]*uhat[(j+N*i)][1]);
			Eng += (uhat[(i+N*j)][0]*uhat[(i+N*j)][0] + uhat[(i+N*j)][1]*uhat[(i+N*j)][1]);
		  out << kk << "\t" << Eng << "\n";
    }
	}
	out.close();
#endif
	fftw_destroy_plan(inv_trans_u);
	fftw_destroy_plan(trans_u);

  double *x,*y;
	x    = (double*) calloc (NI, sizeof(double));
	y    = (double*) calloc (NJ, sizeof(double));

	double *rho, *rhou, *rhov, *e;
	rho  = (double*) calloc (N*N, sizeof(double));
	rhou = (double*) calloc (N*N, sizeof(double));
	rhov = (double*) calloc (N*N, sizeof(double));
	e    = (double*) calloc (N*N, sizeof(double));

	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
		  x[i] = -0.5 + i*dx;
		  y[j] = -0.5 + j*dy;
      int p = i*NJ + j;
      double RHO,U,V,P;
      RHO = 1.0 + u[j+N*i][0];
      U   = 1.0;
      V   = 1.0;
      P   = 1.0/1.4;
      rho[p]  = RHO;
      rhou[p] = RHO*U;
      rhov[p] = RHO*V;
      e[p]    = P/0.4 + 0.5*RHO*(U*U+V*V);
		}
	}

	fftw_free(uhat);
	fftw_free(u);

  FILE *op;
	op = fopen("initial.inp","w");
  for (i = 0; i < NI; i++)  fprintf(op,"%1.16E ",x[i]);
  fprintf(op,"\n");
  for (j = 0; j < NJ; j++)  fprintf(op,"%1.16E ",y[j]);
  fprintf(op,"\n");
	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NJ; j++) {
      int p = i*NJ + j;
      fprintf(op,"%1.16E ",rho[p]);
    }
  }
  fprintf(op,"\n");
	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NJ; j++) {
      int p = i*NJ + j;
      fprintf(op,"%1.16E ",rhou[p]);
    }
  }
  fprintf(op,"\n");
	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NJ; j++) {
      int p = i*NJ + j;
      fprintf(op,"%1.16E ",rhov[p]);
    }
  }
  fprintf(op,"\n");
	for (i = 0; i < NI; i++) {
	  for (j = 0; j < NJ; j++) {
      int p = i*NJ + j;
      fprintf(op,"%1.16E ",e[p]);
    }
  }
  fprintf(op,"\n");
	fclose(op);
	free(x);
	free(y);
	free(rho);
	free(rhou);
	free(rhov);
	free(e);

	return(0);
}
