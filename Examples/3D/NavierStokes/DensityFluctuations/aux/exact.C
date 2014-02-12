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

	int NI,NJ,NK,ndims;
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
        if (!strcmp(word, "ndims"))     in >> ndims;
        else if (!strcmp(word, "size")) in >> NI >> NJ >> NK;
        else if (!strcmp(word, "ip_file_type")) in >> ip_file_type;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  if (ndims != 3) {
    std::cout << "ndims is not 3 in solver.inp. this code is to generate 3D exact conditions\n";
    return(0);
  }
	std::cout << "Grid:\t\t\t" << NI << " x " << NJ << " x " << NK << "\n";
  std::cout << "Input maximum wavenumber (typically N/2): ";
  int limit; std::cin >> limit;

  if ((NI != NJ) || (NI != NK) || (NJ != NK)) { printf("Error: NI,NJ,NK not equal. Bye!\n"); return(0); }
	int N = NI, N3 = N*N*N;
	int i,j,k;
	double dx = 1.0 / ((double)N);
  double factor = 0.00001;

	fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
	fftw_complex *u    = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));
	fftw_plan inv_trans_u = fftw_plan_dft_3d(N, N, N, uhat, u, 1, FFTW_MEASURE);

  if (limit > N/2) limit = N/2;

	for (i = 0; i < N3; i++){
		uhat[i][0] = 0;
		uhat[i][1] = 0;
	}
	for (i = 1; i <= limit; i++){
		for (j = 1; j <= limit; j++){
		  for (k = 1; k <= limit; k++){
			  double kk = sqrt(i*i + j*j + k*k);
        double Ak = raiseto(kk,-5.0/6.0);
			  uhat[k+N*(j+N*i)][0] = Ak/sqrt(2.0);
			  uhat[k+N*(j+N*i)][1] = Ak/sqrt(2.0);
      }
		}
	}
	for (i = 1+N/2; i < N; i++){
		for (j = 0; j <= N/2; j++){
		  for (k = 0; k <= N/2; k++){
			  uhat[k+N*(j+N*i)][0] =  uhat[k+N*(j+N*(N-i))][0];
			  uhat[k+N*(j+N*i)][1] = -uhat[k+N*(j+N*(N-i))][1];
      }
		}
	}
	for (i = 0; i < N; i++){
		for (j = N/2+1; j < N; j++){
      for (k = 0; k <= N/2; k++) {
			  uhat[k+N*(j+N*i)][0] =  uhat[k+N*((N-j)+N*i)][0];
			  uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*i)][1];
      }
		}
	}
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
      for (k = N/2+1; k < N; k++) {
			  uhat[k+N*(j+N*i)][0] =  uhat[(N-k)+N*(j+N*i)][0];
			  uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*i)][1];
      }
		}
	}

	fftw_execute(inv_trans_u);
	fftw_destroy_plan(inv_trans_u);

  double *x,*y,*z;
	x    = (double*) calloc (N, sizeof(double));
	y    = (double*) calloc (N, sizeof(double));
	z    = (double*) calloc (N, sizeof(double));

	double *U;
	U  = (double*) calloc (5*N3, sizeof(double));

	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
  		  x[i] = -0.5 + i*dx;
	  	  y[j] = -0.5 + j*dx;
	  	  z[k] = -0.5 + k*dx;
        double RHO, uvel, vvel, wvel, P;
        RHO = 1.0 + factor * u[k+N*(j+N*i)][0];
        if (RHO < 0) std::cout << "Warning: negative density calculated!\n";
        uvel   = 1.0;
        vvel   = 1.0;
        wvel   = 1.0;
        P   = 1.0/1.4;
        int p = i + N*j + N*N*k;
        U[5*p+0] = RHO;
        U[5*p+1] = RHO*uvel;
        U[5*p+2] = RHO*vvel;
        U[5*p+3] = RHO*wvel;
        U[5*p+4] = P/0.4 + 0.5*RHO*(uvel*uvel+vvel*vvel+wvel*wvel);
      }
		}
	}

	fftw_free(uhat);
	fftw_free(u);

  FILE *op;
  if (!strcmp(ip_file_type,"ascii")) {

    printf("Writing ASCII initial solution file initial.inp\n");
  	op = fopen("initial.inp","w");
    for (i = 0; i < N; i++)  fprintf(op,"%1.16E ",x[i]);
    fprintf(op,"\n");
    for (j = 0; j < N; j++)  fprintf(op,"%1.16E ",y[j]);
    fprintf(op,"\n");
    for (k = 0; k < N; k++)  fprintf(op,"%1.16E ",z[k]);
    fprintf(op,"\n");
	  for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
	      for (k = 0; k < N; k++) {
          int p = i + N*j + N*N+k;
          fprintf(op,"%1.16E ",U[5*p+0]);
        }
      }
    }
    fprintf(op,"\n");
  	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
	      for (k = 0; k < N; k++) {
          int p = i + N*j + N*N+k;
          fprintf(op,"%1.16E ",U[5*p+1]);
        }
      }
    }
    fprintf(op,"\n");
  	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
	      for (k = 0; k < N; k++) {
          int p = i + N*j + N*N+k;
          fprintf(op,"%1.16E ",U[5*p+2]);
        }
      }
    }
    fprintf(op,"\n");
	  for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
	      for (k = 0; k < N; k++) {
          int p = i + N*j + N*N+k;
          fprintf(op,"%1.16E ",U[5*p+3]);
        }
      }
    }
    fprintf(op,"\n");
  	for (i = 0; i < N; i++) {
	    for (j = 0; j < N; j++) {
	      for (k = 0; k < N; k++) {
          int p = i + N*j + N*N+k;
          fprintf(op,"%1.16E ",U[5*p+4]);
        }
      }
    }
    fprintf(op,"\n");
	  fclose(op);

  } else if ((!strcmp(ip_file_type,"binary")) || (!strcmp(ip_file_type,"bin"))) {

    printf("Writing binary initial solution file initial.inp\n");
  	op = fopen("initial.inp","wb");
    fwrite(x,sizeof(double),N,op);
    fwrite(y,sizeof(double),N,op);
    fwrite(z,sizeof(double),N,op);
    fwrite(U,sizeof(double),5*N3,op);
    fclose(op);

  }

	free(x);
	free(y);
	free(z);
  free(U);

	return(0);
}
