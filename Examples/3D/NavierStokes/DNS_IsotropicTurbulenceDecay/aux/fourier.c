/* 

  Code to calculate the energy spectrum of the solution
  at a given instance.

  Create a symbolic link named "op.bin" to the desired
  solution file "op_xxxxx.bin" if not overwriting solu-
  tion files.

  It will write a file "spectrum_final.dat" which has
  two columns - k and E(k)

  It needs the FFTW3 library installed. To compile it, 
  make sure the fftw3.h and libfftw3 are available in
  the include and linking paths.

*/


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

void fourier_analysis(int,double*,double*,double*);

int main()
{
  int i,j,k,NI,NJ,NK;
  char op_file_format[50];
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "size")) {
          fscanf(in,"%d",&NI);
          fscanf(in,"%d",&NJ);
          fscanf(in,"%d",&NK);
        } else if (!strcmp(word, "op_file_format")) fscanf(in,"%s",op_file_format);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  printf("Grid size: %d x %d x %d.\n",NI,NJ,NK);
  int N = NI, N3 = N*N*N;

  if ((!strcmp(op_file_format,"binary")) || (!strcmp(op_file_format,"bin"))) {

    FILE *inp;
    double *x = (double*) calloc (N,sizeof(double));
    double *y = (double*) calloc (N,sizeof(double));
    double *z = (double*) calloc (N,sizeof(double));
    double *U = (double*) calloc (5*N3,sizeof(double));
    int ndims, nvars, dims[3];

    double *u,*v,*w;
    u = (double*) calloc (N3, sizeof(double));
    v = (double*) calloc (N3, sizeof(double));
    w = (double*) calloc (N3, sizeof(double));

    printf("Reading binary solution file op.bin\n");
    inp = fopen("op.bin","rb");

    if (!inp) {
      printf("Output file op.bin not found!\n");
      return(0);
    }

    fread(&ndims,sizeof(int),1,inp);
    fread(&nvars,sizeof(int),1,inp);
    fread(dims,sizeof(int),3,inp);
    if ((dims[0] != N) || (dims[1] != N) || (dims[2] != N)) {
      printf("Error: incorrect dimensions read from solution file. N=%d, dims=%d,%d,%d.\n",
              N,dims[0],dims[1],dims[2]);
      return(0);
    }
    fread(x,sizeof(double),N,inp);
    fread(y,sizeof(double),N,inp);
    fread(z,sizeof(double),N,inp);
    fread(U,sizeof(double),5*N3,inp);
    fclose(inp);
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        for (k=0; k<N; k++) {
          double rho = U[5*(i+N*(j+N*k))];
          u[k+N*(j+N*i)] = U[5*(i+N*(j+N*k))+1]/rho;
          v[k+N*(j+N*i)] = U[5*(i+N*(j+N*k))+2]/rho;
          w[k+N*(j+N*i)] = U[5*(i+N*(j+N*k))+3]/rho;
        }
      }
    }
    fourier_analysis(N,u,v,w);
    free(u);
    free(v);
    free(w);

    free(x);
    free(y);
    free(z);
    free(U);

  } else {
    printf("Error: Unsupported output file type. Use binary output only!\n");
  }
 
  return(0);
}

void fourier_analysis(int N, double *uu, double *vv, double *ww)
{
	int i,j,k;
  int N3 = N*N*N;

	fftw_complex *u    = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));	
	fftw_complex *v    = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));	
	fftw_complex *w    = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));	
	fftw_complex *uhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));	
	fftw_complex *vhat = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));	
	fftw_complex *what = (fftw_complex*) fftw_malloc(N3 * sizeof(fftw_complex));	

	fftw_plan transform_u, transform_v, transform_w;
	transform_u = fftw_plan_dft_3d(N, N, N, u, uhat, -1, FFTW_MEASURE);
	transform_v = fftw_plan_dft_3d(N, N, N, v, vhat, -1, FFTW_MEASURE);
	transform_w = fftw_plan_dft_3d(N, N, N, w, what, -1, FFTW_MEASURE);

	double rms_velocity = 0;
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			for (k = 0; k < N; k++){
				u[k+N*(j+N*i)][0] = uu[k+N*(j+N*i)];
				v[k+N*(j+N*i)][0] = vv[k+N*(j+N*i)];
				w[k+N*(j+N*i)][0] = ww[k+N*(j+N*i)];
				u[k+N*(j+N*i)][1] = 0;
				v[k+N*(j+N*i)][1] = 0;
				w[k+N*(j+N*i)][1] = 0;

				rms_velocity +=  (uu[k+N*(j+N*i)]*uu[k+N*(j+N*i)] 
						            + vv[k+N*(j+N*i)]*vv[k+N*(j+N*i)]
						            + ww[k+N*(j+N*i)]*ww[k+N*(j+N*i)]);
			}
		}
	}
	rms_velocity = sqrt(rms_velocity / (3*N3));
	printf("RMS velocity (component-wise): %1.16E\n",rms_velocity);

	fftw_execute(transform_u);
	fftw_execute(transform_v);
	fftw_execute(transform_w);

	fftw_free(u);
	fftw_free(v);
	fftw_free(w);

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

	double kk = sqrt(3 * (N/2)*(N/2));
	int kkmax = (int) kk;
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
					+ 0.5 * ( (uhat[k+N*(j+N*i)][0]*uhat[k+N*(j+N*i)][0] + uhat[k+N*(j+N*i)][1]*uhat[k+N*(j+N*i)][1])
					        + (vhat[k+N*(j+N*i)][0]*vhat[k+N*(j+N*i)][0] + vhat[k+N*(j+N*i)][1]*vhat[k+N*(j+N*i)][1])
					        + (what[k+N*(j+N*i)][0]*what[k+N*(j+N*i)][0] + what[k+N*(j+N*i)][1]*what[k+N*(j+N*i)][1]) );
				total_energy = total_energy 
					+ 0.5 * ( (uhat[k+N*(j+N*i)][0]*uhat[k+N*(j+N*i)][0] + uhat[k+N*(j+N*i)][1]*uhat[k+N*(j+N*i)][1])
					        + (vhat[k+N*(j+N*i)][0]*vhat[k+N*(j+N*i)][0] + vhat[k+N*(j+N*i)][1]*vhat[k+N*(j+N*i)][1])
					        + (what[k+N*(j+N*i)][0]*what[k+N*(j+N*i)][0] + what[k+N*(j+N*i)][1]*what[k+N*(j+N*i)][1]) );
			}
		}
	}
	printf("Total Energy: %1.16E\n",total_energy);
	FILE *out;
	out = fopen("spectrum_final.dat","w");
	for (i = 1; i < kkmax; i++) fprintf(out,"%1.16E\t%1.16E\n",freq[i],Eng[i]/total_energy);
	fclose(out);
	free(freq);
	free(Eng);

	fftw_destroy_plan(transform_u);
	fftw_destroy_plan(transform_v);
	fftw_destroy_plan(transform_w);

	fftw_free(uhat);
	fftw_free(vhat);
	fftw_free(what);
}
