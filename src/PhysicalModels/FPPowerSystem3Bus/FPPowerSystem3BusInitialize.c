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
	int i, j;
	double *r  = (double*) calloc(N  ,sizeof(double));
	double *AA = (double*) calloc(N*N,sizeof(double));
	for (i = 0; i < N; i++){
		for (j = 0; j < N*N; j++) AA[j] = A[j];
		for (j = 0; j < N; j++){
			if (j == i)	r[j] = 1.0;
			else		    r[j] = 0.0;
		}
		int ierr = LUDecomp(AA,r,N); if (ierr) return(ierr);
		for (j = 0; j < N; j++)	B[j*N+i] = r[j];
	}
  free(r);
  free(AA);
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

  double PI = 4.0*atan(1.0);

  physics->N    = N = 6;
  physics->G    = (double*) calloc ((N/2)*(N/2),sizeof(double));
  physics->B    = (double*) calloc ((N/2)*(N/2),sizeof(double));
  physics->Ainv = (double*) calloc (N*N        ,sizeof(double));

  /* default values of model parameters */
  physics->PM1     = 3.109260511864138;
  physics->PM2     = 1.0;
  physics->H1      = 1000.640;
  physics->H2      = 3.64;
  physics->omegaB  = 2*PI*60;
  physics->D1      = 450;
  physics->D2      = 1.0;
  physics->E1      = 1.044225672060674;
  physics->E2      = 1.034543707656856;
  physics->Xd1     = 0.02;
  physics->Xd2     = 0.2;
  physics->alpha   = 4.386890097147679;
  physics->beta    = -1.096722524286920;

  physics->lambda[0][0] = 0.1;
  physics->lambda[0][1] = 0.0;
  physics->lambda[1][0] = 0.0;
  physics->lambda[1][1] = 0.1;

  physics->sigma[0][0] = 1.0;
  physics->sigma[0][1] = 0.0;
  physics->sigma[1][0] = 0.0;
  physics->sigma[1][1] = 1.0;

  physics->G[0*N/2+0] =  7.631257631257632;
  physics->G[0*N/2+1] = -3.815628815628816;
  physics->G[0*N/2+2] = -3.815628815628816;
  physics->G[1*N/2+0] = -3.815628815628816;
  physics->G[1*N/2+1] =  6.839334669523348;
  physics->G[1*N/2+2] = -3.023705853894533;
  physics->G[2*N/2+0] = -3.815628815628816;
  physics->G[2*N/2+1] = -3.023705853894533;
  physics->G[2*N/2+2] =  6.839334669523348;

  physics->B[0*N/2+0] = -38.053788156288157;
  physics->B[0*N/2+1] =  19.078144078144078;
  physics->B[0*N/2+2] =  19.078144078144078;
  physics->B[1*N/2+0] =  19.078144078144078;
  physics->B[1*N/2+1] = -34.081673347616743;
  physics->B[1*N/2+2] =  15.118529269472663;
  physics->B[2*N/2+0] =  19.078144078144078;
  physics->B[2*N/2+1] =  15.118529269472663;
  physics->B[2*N/2+2] = -34.081673347616743;


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
        if      (!strcmp(word,"PM1"    ))  {ferr=fscanf(in,"%lf",&physics->PM1   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H1"     ))  {ferr=fscanf(in,"%lf",&physics->H1    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"D1"     ))  {ferr=fscanf(in,"%lf",&physics->D1    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E1"     ))  {ferr=fscanf(in,"%lf",&physics->E1    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Xd1"    ))  {ferr=fscanf(in,"%lf",&physics->Xd1   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"PM2"    ))  {ferr=fscanf(in,"%lf",&physics->PM2   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H2"     ))  {ferr=fscanf(in,"%lf",&physics->H2    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"D2"     ))  {ferr=fscanf(in,"%lf",&physics->D2    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E2"     ))  {ferr=fscanf(in,"%lf",&physics->E2    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Xd2"    ))  {ferr=fscanf(in,"%lf",&physics->Xd2   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"omegaB" ))  {ferr=fscanf(in,"%lf",&physics->omegaB) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"alpha"  ))  {ferr=fscanf(in,"%lf",&physics->alpha ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"beta"   ))  {ferr=fscanf(in,"%lf",&physics->beta  ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"sigma"  ))  {
          ferr=fscanf(in,"%lf",&physics->sigma[0][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->sigma[0][1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->sigma[1][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->sigma[1][1]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"lambda"))  {
          ferr=fscanf(in,"%lf",&physics->lambda[0][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->lambda[0][1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->lambda[1][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->lambda[1][1]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"G"))  {
          ferr=fscanf(in,"%lf",&physics->G[0*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[0*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[0*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[1*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[1*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[1*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[2*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[2*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[2*N/2+2]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"B"))  {
          ferr=fscanf(in,"%lf",&physics->B[0*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[0*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[0*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[1*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[1*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[1*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[2*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[2*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[2*N/2+2]) ;if(ferr!=1)return(1);
        }
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
  A[0*N+0] =  physics->G[0*N/2+0];
  A[0*N+1] = -physics->B[0*N/2+0] + 1.0/physics->Xd1;
  A[1*N+0] =  physics->B[0*N/2+0] - 1.0/physics->Xd1;
  A[1*N+1] =  physics->G[0*N/2+0];
  A[0*N+2] =  physics->G[0*N/2+1];
  A[0*N+3] = -physics->B[0*N/2+1];
  A[1*N+2] =  physics->B[0*N/2+1];
  A[1*N+3] =  physics->G[0*N/2+1];
  A[0*N+4] =  physics->G[0*N/2+2];
  A[0*N+5] = -physics->B[0*N/2+2];
  A[1*N+4] =  physics->B[0*N/2+2];
  A[1*N+5] =  physics->G[0*N/2+2];
  A[2*N+0] =  physics->G[1*N/2+0];
  A[2*N+1] = -physics->B[1*N/2+0];
  A[3*N+0] =  physics->B[1*N/2+0];
  A[3*N+1] =  physics->G[1*N/2+0];
  A[2*N+2] =  physics->G[1*N/2+1];
  A[2*N+3] = -physics->B[1*N/2+1] + 1.0/physics->Xd2;
  A[3*N+2] =  physics->B[1*N/2+1] - 1.0/physics->Xd2;
  A[3*N+3] =  physics->G[1*N/2+1];
  A[2*N+4] =  physics->G[1*N/2+2];
  A[2*N+5] = -physics->B[1*N/2+2];
  A[3*N+4] =  physics->B[1*N/2+2];
  A[3*N+5] =  physics->G[1*N/2+2];
  A[4*N+0] =  physics->G[2*N/2+0];
  A[4*N+1] = -physics->B[2*N/2+0];
  A[5*N+0] =  physics->B[2*N/2+0];
  A[5*N+1] =  physics->G[2*N/2+0];
  A[4*N+2] =  physics->G[2*N/2+1];
  A[4*N+3] = -physics->B[2*N/2+1];
  A[5*N+2] =  physics->B[2*N/2+1];
  A[5*N+3] =  physics->G[2*N/2+1];
  A[4*N+4] =  physics->G[2*N/2+2] + physics->alpha;
  A[4*N+5] = -physics->B[2*N/2+2] - physics->beta;
  A[5*N+4] =  physics->B[2*N/2+2] + physics->beta;
  A[5*N+5] =  physics->G[2*N/2+2] + physics->alpha;
  
  int inverr = MatInverse(A,physics->Ainv,N);
  if (inverr) {
    if (!mpi->rank) fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): Unable to invert matrix!\n");
    return(inverr);
  }

#if 0
  /* Verifying Ainv is correct */
  int M = N, i, j, k ;
  double *Ainv = physics->Ainv;
  double eye[M*M];
  for (j=0; j<M; j++) {
    for (k=0; k<M; k++) {
      eye[j*M+k] = 0.0;
      for (i=0; i<M; i++) eye[j*M+k] += (A[j*M+i] * Ainv[i*M+k]);
    }
  }
  printf("\n");
  printf("Verifying Ainv is correct. A*Ainv = \n");
  for (j=0; j<M; j++) {
    for (k=0; k<M; k++) printf("%+1.16E  ",eye[j*M+k]);
    printf("\n");
  }
  printf("\n");
#endif
  return(0);
}

