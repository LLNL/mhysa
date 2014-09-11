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
  /* initial spectrum parameters */
	double kp = 16.0;
	double u0 = 0.2 * sqrt(2.0);

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
        if (!strcmp(word, "ndims"))             in >> ndims;
        else if (!strcmp(word, "size"))         in >> NI >> NJ >> NK;
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
  std::cout << "Input wavenumber of peak: ";
  std::cin >> kp;

  if ((NI != NJ) || (NI != NK) || (NJ != NK)) { 
    printf("Error: NI,NJ,NK not equal. Bye!\n"); 
    return(0); 
  }
	int N = NI, N3 = N*N*N;
  if (limit > N/2) limit = N/2;
	int i,j,k;
	double dx = 2*PI / ((double)N);

  /* Calculating the velocity components through a Fourier transform */

	double kk = sqrt(3 * (N/2)*(N/2));
	int kkmax = (int) kk;

	fftw_complex *uhat = (fftw_complex*) fftw_malloc(N*N*N * sizeof(fftw_complex));	
	fftw_complex *vhat = (fftw_complex*) fftw_malloc(N*N*N * sizeof(fftw_complex));	
	fftw_complex *what = (fftw_complex*) fftw_malloc(N*N*N * sizeof(fftw_complex));	
	fftw_complex *u = (fftw_complex*) fftw_malloc(N*N*N * sizeof(fftw_complex));	
	fftw_complex *v = (fftw_complex*) fftw_malloc(N*N*N * sizeof(fftw_complex));	
	fftw_complex *w = (fftw_complex*) fftw_malloc(N*N*N * sizeof(fftw_complex));	

	fftw_plan transform_u, transform_v, transform_w;
	fftw_plan inv_trans_u, inv_trans_v, inv_trans_w;

	transform_u = fftw_plan_dft_3d(N, N, N, u, uhat, -1, FFTW_MEASURE);
	transform_v = fftw_plan_dft_3d(N, N, N, v, vhat, -1, FFTW_MEASURE);
	transform_w = fftw_plan_dft_3d(N, N, N, w, what, -1, FFTW_MEASURE);
	inv_trans_u = fftw_plan_dft_3d(N, N, N, uhat, u, 1, FFTW_MEASURE);
	inv_trans_v = fftw_plan_dft_3d(N, N, N, vhat, v, 1, FFTW_MEASURE);
	inv_trans_w = fftw_plan_dft_3d(N, N, N, what, w, 1, FFTW_MEASURE);

  /* Specifying the velocities in Fourier space */
	for (i = 1; i <= N/2; i++){
		for (j = 0; j <= N/2; j++){
			for (k = 0; k <= N/2; k++){
				double kk = sqrt(i*i + j*j + k*k);
				double th1 = 2*PI * ((double)rand())/((double)RAND_MAX);
				double th2 = 2*PI * ((double)rand())/((double)RAND_MAX);
				double phi1 = 2*PI * ((double)rand())/((double)RAND_MAX);

				double E = 16.0 * sqrt(2.0/PI) * (u0*u0/kp) * raiseto(kk/kp, 4.0) * exp(-2.0*(kk/kp)*(kk/kp));
				double alfa_real = sqrt(E/(4*PI*kk*kk))*cos(th1)*cos(phi1);
				double alfa_imag = sqrt(E/(4*PI*kk*kk))*sin(th1)*cos(phi1);
				double beta_real = sqrt(E/(4*PI*kk*kk))*cos(th2)*sin(phi1);
				double beta_imag = sqrt(E/(4*PI*kk*kk))*sin(th2)*sin(phi1);

				uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
				vhat[k+N*(j+N*i)][0] = (beta_real*j*k+alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
				what[k+N*(j+N*i)][0] = (beta_real*sqrt(i*i+j*j))/kk;
				uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));
				vhat[k+N*(j+N*i)][1] = (beta_imag*j*k+alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));
				what[k+N*(j+N*i)][1] = (beta_imag*sqrt(i*i+j*j))/kk;
			}
		}
	}
	for (i = 0; i < 1; i++){
		for (k = 0; k <= N/2; k++){
			for (j = 1; j <= N/2; j++){
				double kk = sqrt(i*i + j*j + k*k);
				double th1 = 2*PI * ((double)rand())/((double)RAND_MAX);
				double th2 = 2*PI * ((double)rand())/((double)RAND_MAX);
				double phi1 = 2*PI * ((double)rand())/((double)RAND_MAX);

				double E = 16.0 * sqrt(2.0/PI) * (u0*u0/kp) * raiseto(kk/kp, 4.0) * exp(-2.0*(kk/kp)*(kk/kp));
				double alfa_real = sqrt(E/(4*PI*kk*kk))*cos(th1)*cos(phi1);
				double alfa_imag = sqrt(E/(4*PI*kk*kk))*sin(th1)*cos(phi1);
				double beta_real = sqrt(E/(4*PI*kk*kk))*cos(th2)*sin(phi1);
				double beta_imag = sqrt(E/(4*PI*kk*kk))*sin(th2)*sin(phi1);

				uhat[k+N*(j+N*i)][0] = (alfa_real*kk*j+beta_real*i*k)/(kk*sqrt(i*i+j*j));
				vhat[k+N*(j+N*i)][0] = (beta_real*j*k+alfa_real*kk*i)/(kk*sqrt(i*i+j*j));
				what[k+N*(j+N*i)][0] = (beta_real*sqrt(i*i+j*j))/kk;
				uhat[k+N*(j+N*i)][1] = (alfa_imag*kk*j+beta_imag*i*k)/(kk*sqrt(i*i+j*j));
				vhat[k+N*(j+N*i)][1] = (beta_imag*j*k+alfa_imag*kk*i)/(kk*sqrt(i*i+j*j));
				what[k+N*(j+N*i)][1] = (beta_imag*sqrt(i*i+j*j))/kk;
			}
		}
	}
	for (i = 0; i < 1; i++){
		for (j = 0; j < 1; j++){
			for (k = 0; k <= N/2; k++){
				uhat[k+N*(j+N*i)][0] = 0;
				vhat[k+N*(j+N*i)][0] = 0;
				what[k+N*(j+N*i)][0] = 0;
				uhat[k+N*(j+N*i)][1] = 0;
				vhat[k+N*(j+N*i)][1] = 0;
				what[k+N*(j+N*i)][1] = 0;
			}
		}
	}
	for (i = 1+N/2; i < N; i++){
		for (j = 0; j <= N/2; j++){
			for (k = 0; k <= N/2; k++){
				uhat[k+N*(j+N*i)][0] = uhat[k+N*(j+N*(N-i))][0];
				vhat[k+N*(j+N*i)][0] = vhat[k+N*(j+N*(N-i))][0];
				what[k+N*(j+N*i)][0] = what[k+N*(j+N*(N-i))][0];
				uhat[k+N*(j+N*i)][1] = -uhat[k+N*(j+N*(N-i))][1];
				vhat[k+N*(j+N*i)][1] = -vhat[k+N*(j+N*(N-i))][1];
				what[k+N*(j+N*i)][1] = -what[k+N*(j+N*(N-i))][1];
			}
		}
	}
	for (i = 0; i < N; i++){
		for (j = N/2+1; j < N; j++){
			for (k = 0; k <= N/2; k++){
				uhat[k+N*(j+N*i)][0] = uhat[k+N*((N-j)+N*i)][0];
				vhat[k+N*(j+N*i)][0] = vhat[k+N*((N-j)+N*i)][0];
				what[k+N*(j+N*i)][0] = what[k+N*((N-j)+N*i)][0];
				uhat[k+N*(j+N*i)][1] = -uhat[k+N*((N-j)+N*i)][1];
				vhat[k+N*(j+N*i)][1] = -vhat[k+N*((N-j)+N*i)][1];
				what[k+N*(j+N*i)][1] = -what[k+N*((N-j)+N*i)][1];
			}
		}
	}
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = N/2+1; k < N; k++){
				uhat[k+N*(j+N*i)][0] = uhat[(N-k)+N*(j+N*i)][0];
				vhat[k+N*(j+N*i)][0] = vhat[(N-k)+N*(j+N*i)][0];
				what[k+N*(j+N*i)][0] = what[(N-k)+N*(j+N*i)][0];
				uhat[k+N*(j+N*i)][1] = -uhat[(N-k)+N*(j+N*i)][1];
				vhat[k+N*(j+N*i)][1] = -vhat[(N-k)+N*(j+N*i)][1];
				what[k+N*(j+N*i)][1] = -what[(N-k)+N*(j+N*i)][1];
			}
		}
	}

	double *freq = (double*) calloc(kkmax+1, sizeof(double));
	double *Eng = (double*) calloc(kkmax+1, sizeof(double));
	double total_energy = 0.0;
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = 0; k < N; k++){
				int isq, jsq, ksq;
				if (i > N/2)	isq = (i-N) * (i-N);
				else		isq = i*i;
				if (j > N/2)	jsq = (j-N) * (j-N);
				else		jsq = j*j;
				if (k > N/2)	ksq = (k-N) * (k-N);
				else		ksq = k*k;
				double kk = sqrt(isq + jsq + ksq);
				freq[(int)kk] = kk;
				Eng[(int)kk] = Eng[(int)kk] 
					+ ((uhat[k+N*(j+N*i)][0]*uhat[k+N*(j+N*i)][0] + uhat[k+N*(j+N*i)][1]*uhat[k+N*(j+N*i)][1])
					+  (vhat[k+N*(j+N*i)][0]*vhat[k+N*(j+N*i)][0] + vhat[k+N*(j+N*i)][1]*vhat[k+N*(j+N*i)][1])
					+  (what[k+N*(j+N*i)][0]*what[k+N*(j+N*i)][0] + what[k+N*(j+N*i)][1]*what[k+N*(j+N*i)][1]));
				total_energy = total_energy 
					+ ((uhat[k+N*(j+N*i)][0]*uhat[k+N*(j+N*i)][0] + uhat[k+N*(j+N*i)][1]*uhat[k+N*(j+N*i)][1])
					+  (vhat[k+N*(j+N*i)][0]*vhat[k+N*(j+N*i)][0] + vhat[k+N*(j+N*i)][1]*vhat[k+N*(j+N*i)][1])
					+  (what[k+N*(j+N*i)][0]*what[k+N*(j+N*i)][0] + what[k+N*(j+N*i)][1]*what[k+N*(j+N*i)][1]));
			}
		}
	}
	std::cout << "Total Energy: " << total_energy << "\n";
	out.open("spectrum_initial.dat");
	for (i = 0; i < kkmax; i++){
		out << freq[i] << "\t" << Eng[i]/total_energy << "\n";
	}
	out.close();
	free(freq);
	free(Eng);

  /* Inverse Fourier transform */
	fftw_execute(inv_trans_u);
	fftw_execute(inv_trans_v);
	fftw_execute(inv_trans_w);
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = 0; k < N; k++){
				u[k+N*(j+N*i)][1] = 0.0;
				v[k+N*(j+N*i)][1] = 0.0;
				w[k+N*(j+N*i)][1] = 0.0;
			}
		}
	}

#if 0
  /* check */
	fftw_execute(transform_u);
	fftw_execute(transform_v);
	fftw_execute(transform_w);
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = 0; k < N; k++){
				uhat[k+N*(j+N*i)][0] /= (N*N*N);
				vhat[k+N*(j+N*i)][0] /= (N*N*N);
				what[k+N*(j+N*i)][0] /= (N*N*N);
				uhat[k+N*(j+N*i)][1] /= (N*N*N);
				vhat[k+N*(j+N*i)][1] /= (N*N*N);
				what[k+N*(j+N*i)][1] /= (N*N*N);
			}
		}
	}
	freq = (double*) calloc(kkmax+1, sizeof(double));
	Eng = (double*) calloc(kkmax+1, sizeof(double));
	double total_energy2 = 0.0;
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = 0; k < N; k++){
				int isq, jsq, ksq;
				if (i > N/2)	isq = (i-N) * (i-N);
				else		isq = i*i;
				if (j > N/2)	jsq = (j-N) * (j-N);
				else		jsq = j*j;
				if (k > N/2)	ksq = (k-N) * (k-N);
				else		ksq = k*k;
				double kk = sqrt(isq + jsq + ksq);
				freq[(int)kk] = kk;
				Eng[(int)kk] = Eng[(int)kk] 
					+ ((uhat[k+N*(j+N*i)][0]*uhat[k+N*(j+N*i)][0] + uhat[k+N*(j+N*i)][1]*uhat[k+N*(j+N*i)][1])
					+  (vhat[k+N*(j+N*i)][0]*vhat[k+N*(j+N*i)][0] + vhat[k+N*(j+N*i)][1]*vhat[k+N*(j+N*i)][1])
					+  (what[k+N*(j+N*i)][0]*what[k+N*(j+N*i)][0] + what[k+N*(j+N*i)][1]*what[k+N*(j+N*i)][1]));
				total_energy2 = total_energy2 
					+ ((uhat[k+N*(j+N*i)][0]*uhat[k+N*(j+N*i)][0] + uhat[k+N*(j+N*i)][1]*uhat[k+N*(j+N*i)][1])
					+  (vhat[k+N*(j+N*i)][0]*vhat[k+N*(j+N*i)][0] + vhat[k+N*(j+N*i)][1]*vhat[k+N*(j+N*i)][1])
					+  (what[k+N*(j+N*i)][0]*what[k+N*(j+N*i)][0] + what[k+N*(j+N*i)][1]*what[k+N*(j+N*i)][1]));
			}
		}
	}
	std::cout << "Total Energy (after one round of transformations): " << total_energy2 << "\n";
	out.open("spectrum_initial_test2.dat");
	for (i = 0; i < kkmax; i++){
		out << freq[i] << "\t" << Eng[i] << "\n";
	}
	out.close();
	free(freq);
	free(Eng);
#endif

	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(what);

	fftw_destroy_plan(transform_u);
	fftw_destroy_plan(transform_v);
	fftw_destroy_plan(transform_w);
	fftw_destroy_plan(inv_trans_u);
	fftw_destroy_plan(inv_trans_v);
	fftw_destroy_plan(inv_trans_w);

	double rms_velocity = 0;
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = 0; k < N; k++){
        double uu, vv, ww;
				uu = u[k+N*(j+N*i)][0];
				vv = v[k+N*(j+N*i)][0];
				ww = w[k+N*(j+N*i)][0];
				rms_velocity += (uu*uu + vv*vv + ww*ww);
			}
		}
	}
	rms_velocity = sqrt(rms_velocity / (3*N*N*N));
	std::cout << "RMS velocity (component-wise): " << rms_velocity << "\n";

  /* grid and solution in conserved variable form */
  double *x,*y,*z,*U;
	x    = (double*) calloc (N, sizeof(double));
	y    = (double*) calloc (N, sizeof(double));
	z    = (double*) calloc (N, sizeof(double));
	U  = (double*) calloc (5*N3, sizeof(double));
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
      for (k = 0; k < N; k++){
  		  x[i] = i*dx;
	  	  y[j] = j*dx;
	  	  z[k] = k*dx;
        double RHO, uvel, vvel, wvel, P;
        RHO = 1.0;
				uvel = u[k+N*(j+N*i)][0];
				vvel = v[k+N*(j+N*i)][0];
				wvel = w[k+N*(j+N*i)][0];
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

	fftw_free(u);
	fftw_free(v);
	fftw_free(w);

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
