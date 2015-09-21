/*! @file FPPowerSystem3BusInitialize.c
    @author Debojyoti Ghosh
    @brief Function to initialize the 3-bus power system model
*/

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

/*! Initialize the 3-bus power system model:
    Sets the default parameters, read in and set physics-related parameters,
    and set the physics-related function pointers in #HyPar.
*/
int FPPowerSystem3BusInitialize(
                                void *s, /*!< Solver object of type #HyPar */
                                void *m  /*!< MPI object of type #MPIVariables */
                               )
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

  double pi = 4.0*atan(1.0);
  /* default values of model parameters */
  physics->N            = 2;
  physics->Pm1_avg      = 0.8;
  physics->Pm2_avg      = 1.6;
  physics->Pmref_avg    = 0.79330781761651;
  physics->H1           = 3.20;
  physics->H2           = 6.40;
  physics->Href         = 13.60;
  physics->E1           = 1.01556070860155;
  physics->E2           = 1.0491099265981;
  physics->Eref         = 1.05623172878954;
  physics->omegaB       = 2*pi*60.0;
  physics->sigma[0][0]  = 0.0125;
  physics->sigma[0][1]  = 0.0;
  physics->sigma[1][0]  = 0.0;
  physics->sigma[1][1]  = 0.0125;
  physics->lambda[0][0] = 10.0/physics->omegaB;
  physics->lambda[0][1] = 0.0;
  physics->lambda[1][0] = 0.0;
  physics->lambda[1][1] = 10.0/physics->omegaB;
  physics->gamma        = 0.25;

  physics->G   = (double*) calloc ((N+1)*(N+1),sizeof(double));
  physics->B   = (double*) calloc ((N+1)*(N+1),sizeof(double));
  
  physics->G[0*(N+1)+0] = 0.276805493111691;
  physics->G[0*(N+1)+1] = 0.213024867595501;
  physics->G[0*(N+1)+2] = 0.209205876527443;
  physics->G[1*(N+1)+0] = 0.213024867595501;
  physics->G[1*(N+1)+1] = 0.419642083051144;
  physics->G[1*(N+1)+2] = 0.286592141665043;
  physics->G[2*(N+1)+0] = 0.209205876527443;
  physics->G[2*(N+1)+1] = 0.286592141665044;
  physics->G[2*(N+1)+2] = 0.844559256324453;
  
  physics->B[0*(N+1)+0] = -2.36794416971567;
  physics->B[0*(N+1)+1] =  1.08817493992579;
  physics->B[0*(N+1)+2] =  1.22601259339234;
  physics->B[1*(N+1)+0] =  1.08817493992579;
  physics->B[1*(N+1)+1] = -2.72352378723346;
  physics->B[1*(N+1)+2] =  1.51348094527252;
  physics->B[2*(N+1)+0] =  1.22601259339234;
  physics->B[2*(N+1)+1] =  1.51348094527252;
  physics->B[2*(N+1)+2] = -2.98729895217208;

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
        if      (!strcmp(word,"Pm1_avg"   ))  {ferr=fscanf(in,"%lf",&physics->Pm1_avg   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Pm2_avg"   ))  {ferr=fscanf(in,"%lf",&physics->Pm2_avg   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Pmref_avg" ))  {ferr=fscanf(in,"%lf",&physics->Pmref_avg ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H1"        ))  {ferr=fscanf(in,"%lf",&physics->H1        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H2"        ))  {ferr=fscanf(in,"%lf",&physics->H2        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Href"      ))  {ferr=fscanf(in,"%lf",&physics->Href      ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E1"        ))  {ferr=fscanf(in,"%lf",&physics->E1        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E2"        ))  {ferr=fscanf(in,"%lf",&physics->E2        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Eref"      ))  {ferr=fscanf(in,"%lf",&physics->Eref      ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"omegaB"    ))  {ferr=fscanf(in,"%lf",&physics->omegaB    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"gamma"     ))  {ferr=fscanf(in,"%lf",&physics->gamma     ) ;if(ferr!=1)return(1);}
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
          ferr=fscanf(in,"%lf",&physics->G[0*(N+1)+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[0*(N+1)+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[0*(N+1)+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[1*(N+1)+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[1*(N+1)+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[1*(N+1)+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[2*(N+1)+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[2*(N+1)+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->G[2*(N+1)+2]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"B"))  {
          ferr=fscanf(in,"%lf",&physics->B[0*(N+1)+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[0*(N+1)+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[0*(N+1)+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[1*(N+1)+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[1*(N+1)+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[1*(N+1)+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[2*(N+1)+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[2*(N+1)+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&physics->B[2*(N+1)+2]) ;if(ferr!=1)return(1);
        }
      }
	  } else {
    	if (!mpi->rank) fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
  }
  fclose(in);

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in FPPowerSystem3BusInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = FPPowerSystem3BusComputeCFL;
  solver->ComputeDiffNumber  = FPPowerSystem3BusComputeDiffNumber;
  solver->FFunction          = FPPowerSystem3BusAdvection;
  solver->HFunction          = FPPowerSystem3BusDiffusion;
  solver->Upwind             = FPPowerSystem3BusUpwind;

  return(0);
}
