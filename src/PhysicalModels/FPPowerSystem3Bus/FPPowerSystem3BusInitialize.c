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
int    FPPowerSystem3BusCalculateAInv     (double*,double*,void*,double*);

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
  physics->Gf   = (double*) calloc ((N/2)*(N/2),sizeof(double));
  physics->Bf   = (double*) calloc ((N/2)*(N/2),sizeof(double));
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

  physics->Gf[0*N/2+0] =  7.631257631257632;
  physics->Gf[0*N/2+1] = -3.815628815628816;
  physics->Gf[0*N/2+2] = -3.815628815628816;
  physics->Gf[1*N/2+0] = -3.815628815628816;
  physics->Gf[1*N/2+1] =  6.839334669523348;
  physics->Gf[1*N/2+2] = -3.023705853894533;
  physics->Gf[2*N/2+0] = -3.815628815628816;
  physics->Gf[2*N/2+1] = -3.023705853894533;
  physics->Gf[2*N/2+2] =  1006.839334669523348;

  physics->Bf[0*N/2+0] = -38.053788156288157;
  physics->Bf[0*N/2+1] =  19.078144078144078;
  physics->Bf[0*N/2+2] =  19.078144078144078;
  physics->Bf[1*N/2+0] =  19.078144078144078;
  physics->Bf[1*N/2+1] = -34.081673347616743;
  physics->Bf[1*N/2+2] =  15.118529269472663;
  physics->Bf[2*N/2+0] =  19.078144078144078;
  physics->Bf[2*N/2+1] =  15.118529269472663;
  physics->Bf[2*N/2+2] = -34.081673347616743;

  physics->tf  = 0.1;
  physics->tcl = 0.3;

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
        else if (!strcmp(word,"tf"     ))  {ferr=fscanf(in,"%lf",&physics->tf    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"tcl"    ))  {ferr=fscanf(in,"%lf",&physics->tcl   ) ;if(ferr!=1)return(1);}
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
        } else if (!strcmp(word,"Gf"))  {
          ferr=fscanf(in,"%lf",&physics->Gf[0*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[0*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[0*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[1*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[1*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[1*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[2*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[2*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Gf[2*N/2+2]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"Bf"))  {
          ferr=fscanf(in,"%lf",&physics->Bf[0*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[0*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[0*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[1*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[1*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[1*N/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[2*N/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[2*N/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->Bf[2*N/2+2]) ;if(ferr!=1)return(1);
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

  /* calculate the inverse of the impedance matrix */
  IERR FPPowerSystem3BusCalculateAInv(physics->G,physics->B,physics,physics->Ainv); CHECKERR(ierr);

  return(0);
}

