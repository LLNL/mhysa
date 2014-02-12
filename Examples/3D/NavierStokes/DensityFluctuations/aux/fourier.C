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
  char op_file_format[50];
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
        else if (!strcmp(word, "op_file_format")) in >> op_file_format;
      }
    }else{ 
      std::cout << "Error: Illegal format in solver.inp. Crash and burn!\n";
    }
  }
  in.close();
  std::cout << "Grid size: " << NI << " x " << NJ << " x " << NK << "\n";

  int N = NI, N3 = N*N*N;

  if (!strcmp(op_file_format,"text")) {
    std::ifstream solution;
    double *u;
  
    u = (double*) calloc (N3, sizeof(double));
    printf("Reading ASCII solution file op.dat\n");
    solution.open("op.dat");
    if (solution) {
      for (k = 0; k < N; k++) {
        for (j = 0; j < N; j++) {
          for (i = 0; i < N; i++) {
            double temp; int ii, jj, kk;
            solution >> ii >> jj >> kk >> temp >> temp >> temp >> u[k+N*(j+N*i)] >> temp >> temp >> temp >> temp;
          }
        }
      }
      solution.close();
      fourier_analysis(N,u,1);
    } else std::cout << "File not found: op.dat\n";
    free(u);
 
    u = (double*) calloc (N3, sizeof(double));
    printf("Reading ASCII solution file op_exact.dat\n");
    solution.open("op_exact.dat");
    if (solution) {
      for (k = 0; k < N; k++) {
        for (j = 0; j < N; j++) {
          for (i = 0; i < N; i++) {
            double temp; int ii, jj, kk;
            solution >> ii >> jj >> kk >> temp >> temp >> temp >> u[k+N*(j+N*i)] >> temp >> temp >> temp >> temp;
          }
        }
      }
      solution.close();
      fourier_analysis(N,u,0);
    } else std::cout << "File not found: op_exact.dat\n";
    free(u);

  } else if ((!strcmp(op_file_format,"binary")) || (!strcmp(op_file_format,"bin"))) {

    double *x = (double*) calloc (N,sizeof(double));
    double *y = (double*) calloc (N,sizeof(double));
    double *z = (double*) calloc (N,sizeof(double));
    double *U = (double*) calloc (5*N3,sizeof(double));

    int ndims, nvars, dims[3];

    double *u;
    FILE *inp;

    u = (double*) calloc (N3, sizeof(double));
    printf("Reading binary solution file op.bin\n");
    inp = fopen("op.bin","rb");
    fread(&ndims,sizeof(int),1,inp);
    fread(&nvars,sizeof(int),1,inp);
    fread(dims,sizeof(int),3,inp);
    fread(x,sizeof(double),N,inp);
    fread(y,sizeof(double),N,inp);
    fread(z,sizeof(double),N,inp);
    fread(U,sizeof(double),5*N3,inp);
    fclose(inp);
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          u[k+N*(j+N*i)] = U[5*(i+N*(j+N*k))];
        }
      }
    }
    fourier_analysis(N,u,1);
    free(u);

    u = (double*) calloc (N3, sizeof(double));
    printf("Reading binary solution file op_exact.bin\n");
    inp = fopen("op_exact.bin","rb");
    fread(&ndims,sizeof(int),1,inp);
    fread(&nvars,sizeof(int),1,inp);
    fread(dims,sizeof(int),3,inp);
    fread(x,sizeof(double),N,inp);
    fread(y,sizeof(double),N,inp);
    fread(z,sizeof(double),N,inp);
    fread(U,sizeof(double),5*N3,inp);
    fclose(inp);
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          u[k+N*(j+N*i)] = U[5*(i+N*(j+N*k))];
        }
      }
    }
    fourier_analysis(N,u,0);
    free(u);

    free(x);
    free(y);
    free(z);
    free(U);
  }
 
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
