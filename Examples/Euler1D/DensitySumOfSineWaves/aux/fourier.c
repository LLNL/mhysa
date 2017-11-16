#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

double absolute(double x) {
  return (x < 0 ? -x : x);
}
void fourier_analysis(int,int,double*,int);

int main()
{
  int i,NI;
  FILE *in;
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")) {
      while (strcmp(word, "end")) {
        fscanf(in,"%s",word);
        if (!strcmp(word, "size")) fscanf(in,"%d",&NI);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  printf("Grid size: %d\n",NI);

  FILE   *solution;
  double *u;
  
  u = (double*) calloc (NI, sizeof(double));
  solution = fopen("op_exact.dat","r");
  if (solution) {
    for (i = 0; i < NI; i++) {
      double temp; int ii;
      fscanf(solution,"%d",&ii);
      fscanf(solution,"%lf",&temp);
      fscanf(solution,"%lf",&u[i]);
      fscanf(solution,"%lf",&temp);
      fscanf(solution,"%lf",&temp);
      if (ii != i) printf ("Exact solution: Error in reading file (%d,%d)!\n",ii,i);
    }
    fclose(solution);
    fourier_analysis(NI,0,u,0);
  } else printf("Error: exact solution file not found.\n");
  free(u);
  
  u = (double*) calloc (NI, sizeof(double));
  solution = fopen("op.dat","r");
  if (solution) {
    for (i = 0; i < NI; i++) {
      double temp; int ii;
      fscanf(solution,"%d",&ii);
      fscanf(solution,"%lf",&temp);
      fscanf(solution,"%lf",&u[i]);
      fscanf(solution,"%lf",&temp);
      fscanf(solution,"%lf",&temp);
      if (ii != i) printf ("Computed Solution: Error in reading file (%d,%d)!\n",ii,i);
    }
    fclose(solution);
    fourier_analysis(NI,0,u,1);
  } else printf("Error: Solution file not found.\n");
  free(u);

  return(0);
}

void fourier_analysis(int N, int G, double *u, int flag)
{
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

	FILE *out_file;
	if (!flag) out_file = fopen("spectrum_exact.dat","w");
	else       out_file = fopen("spectrum_final.dat","w");
	for (i = 1; i < N; i++){
		double term = 0;
		double phase = 0;
		double xx, yy;
		xx = absolute(out[i][0]);
		yy = absolute(out[i][1]);

		double small = 1e-5;
		if (xx < small){
			if (yy < small)	phase = 0;
			else		phase = PI/2.0;
		}else{
			phase = atan2(yy,xx);
		}
		
		term += (out[i][0]*out[i][0]     +  out[i][1]*out[i][1]);
		term += (out[N-i][0]*out[N-i][0] +  out[N-i][1]*out[N-i][1]);
		fprintf(out_file,"%5d  %1.16E  %1.16E\n",i,term,phase);
	}
	fclose(out_file);

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
}
