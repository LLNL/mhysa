#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <mpivars.h>
#include <hypar.h>

double FPPowerSystem3BusComputeCFL        (void*,void*,double,double);
double FPPowerSystem3BusComputeDiffNumber (void*,void*,double,double);
int    FPPowerSystem3BusAdvection         (double*,double*,int,void*,double);
int    FPPowerSystem3BusDiffusion         (double*,double*,int,int,void*,double);
int    FPPowerSystem3BusUpwind            (double*,double*,double*,double*,
                                           double*,double*,int,void*,double);
int    FPPowerSystem3BusPostStep          (double*,void*,void*,double);
int    FPPowerSystem3BusPrintStep         (void*,void*,double);

int LUDecomp(double *A, double *rhs, int N)
{
	int i, j, k;
  for (i = 0; i < N; i++) {
    if (A[i*N+i] == 0){
      fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): Zero encountered on main diagonal!\n");
			return(1);
		}
		for (j = i+1; j < N; j++){
			double factor = A[j*N+i] / A[i*N+i];
			A[j*N+i] = 0;
			for (k = i+1; k < N; k++) A[j*N+k] -= factor * A[i*N+k];
			rhs[j] -= factor * rhs[i];
		}
	}
	for (i = N-1; i >=0; i--){
		double sum = 0;
		for (j = i+1; j < N; j++) sum += A[i*N+j] * rhs[j];
		rhs[i] = (rhs[i] - sum) / A[i*N+i];
	}
	return(0);
}

int MatInverse(double *A, double *B, int N)
{
	int i, j, k;
	double *r = (double*) calloc(N,sizeof(double));
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (j == i)	r[j] = 1.0;
			else		    r[j] = 0.0;
		}
		int ierr = LUDecomp(A,r,N); if (ierr) return(ierr);
		for (j = 0; j < N; j++)	B[j*N+i] = r[j];
	}
  free(r);
	return(0);
}

int FPPowerSystem3BusInitialize(void *s,void *m)
{
  HyPar               *solver  = (HyPar*)             s;
  MPIVariables        *mpi     = (MPIVariables*)      m; 
  FPPowerSystem3Bus   *physics = (FPPowerSystem3Bus*) solver->physics;
  int                 ferr, N;
  _DECLARE_IERR_;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  physics->N    = N = 6;
  physics->a    = (double*) calloc ((N/2)*(N/2),sizeof(double));
  physics->b    = (double*) calloc ((N/2)*(N/2),sizeof(double));
  physics->Ainv = (double*) calloc (N*N        ,sizeof(double));

  /* default values of model parameters */
  physics->PM1     = 0;
  physics->PM2     = 0;
  physics->H1      = 0;
  physics->H2      = 0;
  physics->omegaB  = 0;
  physics->D1      = 0;
  physics->D2      = 0;
  physics->E1      = 0;
  physics->E2      = 0;
  physics->Xd1     = 0;
  physics->Xd2     = 0;

  physics->a[0*N/2+0] = 0;
  physics->a[0*N/2+1] = 0;
  physics->a[0*N/2+2] = 0;
  physics->a[1*N/2+0] = 0;
  physics->a[1*N/2+1] = 0;
  physics->a[1*N/2+2] = 0;
  physics->a[2*N/2+0] = 0;
  physics->a[2*N/2+1] = 0;
  physics->a[2*N/2+2] = 0;

  physics->b[0*N/2+0] = 0;
  physics->b[0*N/2+1] = 0;
  physics->b[0*N/2+2] = 0;
  physics->b[1*N/2+0] = 0;
  physics->b[1*N/2+1] = 0;
  physics->b[1*N/2+2] = 0;
  physics->b[2*N/2+0] = 0;
  physics->b[2*N/2+1] = 0;
  physics->b[2*N/2+2] = 0;


  /* reading physical model specific inputs - all processes */
  FILE *in;
  if (!mpi->rank) printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) {
    if (!mpi->rank) fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
    char word[_MAX_STRING_SIZE_];
    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      }
	  } else {
    	if (!mpi->rank) fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
  }
  fclose(in);

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = FPPowerSystem3BusComputeCFL;
  solver->ComputeDiffNumber  = FPPowerSystem3BusComputeDiffNumber;
  solver->FFunction          = FPPowerSystem3BusAdvection;
  solver->HFunction          = FPPowerSystem3BusDiffusion;
  solver->Upwind             = FPPowerSystem3BusUpwind;
  solver->PostStep           = FPPowerSystem3BusPostStep;
  solver->PrintStep          = FPPowerSystem3BusPrintStep;

  /* Calculate and print the PDF integral of the initial solution */
  IERR FPPowerSystem3BusPostStep(solver->u,solver,mpi,0.0);        CHECKERR(ierr);
  if (!mpi->rank) IERR FPPowerSystem3BusPrintStep(solver,mpi,0.0); CHECKERR(ierr);

  /* Some initial calculations of the physical parameters */
  double A[N*N];
  A[0*N+0] =  physics->a[0*N/2+0];
  A[0*N+1] = -physics->b[0*N/2+0];
  A[1*N+0] =  physics->b[0*N/2+0];
  A[1*N+1] =  physics->a[0*N/2+0];
  A[0*N+2] =  physics->a[0*N/2+1];
  A[0*N+3] = -physics->b[0*N/2+1];
  A[1*N+2] =  physics->b[0*N/2+1];
  A[1*N+3] =  physics->a[0*N/2+1];
  A[0*N+4] =  physics->a[0*N/2+2];
  A[0*N+5] = -physics->b[0*N/2+2];
  A[1*N+4] =  physics->b[0*N/2+2];
  A[1*N+5] =  physics->a[0*N/2+2];
  A[2*N+0] =  physics->a[1*N/2+0];
  A[2*N+1] = -physics->b[1*N/2+0];
  A[3*N+0] =  physics->b[1*N/2+0];
  A[3*N+1] =  physics->a[1*N/2+0];
  A[2*N+2] =  physics->a[1*N/2+1];
  A[2*N+3] = -physics->b[1*N/2+1];
  A[3*N+2] =  physics->b[1*N/2+1];
  A[3*N+3] =  physics->a[1*N/2+1];
  A[2*N+4] =  physics->a[1*N/2+2];
  A[2*N+5] = -physics->b[1*N/2+2];
  A[3*N+4] =  physics->b[1*N/2+2];
  A[3*N+5] =  physics->a[1*N/2+2];
  A[4*N+0] =  physics->a[2*N/2+0];
  A[4*N+1] = -physics->b[2*N/2+0];
  A[5*N+0] =  physics->b[2*N/2+0];
  A[5*N+1] =  physics->a[2*N/2+0];
  A[4*N+2] =  physics->a[2*N/2+1];
  A[4*N+3] = -physics->b[2*N/2+1];
  A[5*N+2] =  physics->b[2*N/2+1];
  A[5*N+3] =  physics->a[2*N/2+1];
  A[4*N+4] =  physics->a[2*N/2+2];
  A[4*N+5] = -physics->b[2*N/2+2];
  A[5*N+4] =  physics->b[2*N/2+2];
  A[5*N+5] =  physics->a[2*N/2+2];
  
  int inverr = MatInverse(A,physics->Ainv,N);
  if (inverr) {
    if (!mpi->rank) fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): Unable to invert matrix!\n");
    return(inverr);
  }

  return(0);
}

