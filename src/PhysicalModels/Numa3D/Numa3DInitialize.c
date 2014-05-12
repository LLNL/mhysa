#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <boundaryconditions.h>
#include <physicalmodels/numa3d.h>
#include <mpivars.h>
#include <hypar.h>

double Numa3DComputeCFL        (void*,void*,double,double);
int    Numa3DFlux              (double*,double*,int,void*,double);
int    Numa3DSource            (double*,double*,void*,double);
int    Numa3DRusanov           (double*,double*,double*,double*,double*,double*,int,void*,double);

static int Numa3DStandardAtmosphere(void*,double*,int);

int Numa3DInitialize(void *s,void *m)
{
  HyPar           *solver  = (HyPar*)         s;
  MPIVariables    *mpi     = (MPIVariables*)  m; 
  Numa3D          *physics = (Numa3D*)        solver->physics;
  int             ferr     = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in Numa3DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in Numa3DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma  = 1.4; 
  physics->R      = 287.058;        /* J kg^{-1} K^{-1} */
  physics->Omega  = 7.2921150E-05;  /* rad s^{-1}       */
  physics->g      = 9.8;            /* m s^{-2}         */

  physics->Pref   = 101327.0;       /* N m^{-2}         */
  physics->Tref   = 288.15;         /* Kelvin           */

  /* allocate rho0(z), P0(z), T0(z) */
  physics->rho0 = (double*) calloc (solver->dim_local[_ZDIR_]+2*solver->ghosts,sizeof(double));
  physics->P0   = (double*) calloc (solver->dim_local[_ZDIR_]+2*solver->ghosts,sizeof(double));
  physics->T0   = (double*) calloc (solver->dim_local[_ZDIR_]+2*solver->ghosts,sizeof(double));

  /* reading physical model specific inputs - all processes */
  FILE *in;
  if (!mpi->rank) printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) {
    if (!mpi->rank) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
  } else {
    char word[_MAX_STRING_SIZE_];
    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
        if (!strcmp(word, "gamma")) { 
          ferr = fscanf(in,"%lf",&physics->gamma); if (ferr != 1) return(1);
        } else if (!strcmp(word,"R")) {
          ferr = fscanf(in,"%lf",&physics->R); if (ferr != 1) return(1);
        } else if (!strcmp(word,"g")) {
          ferr = fscanf(in,"%lf",&physics->g); if (ferr != 1) return(1);
        } else if (!strcmp(word,"Omega")) {
          ferr = fscanf(in,"%lf",&physics->Omega); if (ferr != 1) return(1);
        } else if (!strcmp(word,"Pref")) {
          ferr = fscanf(in,"%lf",&physics->Pref); if (ferr != 1) return(1);
        } else if (!strcmp(word,"Tref")) {
          ferr = fscanf(in,"%lf",&physics->Tref); if (ferr != 1) return(1);
        } else if (strcmp(word,"end")) {
          char useless[_MAX_STRING_SIZE_];
          ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
          printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",word,useless);
          printf("recognized or extraneous. Ignoring.\n");
        }
      }
	  } else {
    	fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
  }
  fclose(in);

  /* calculate the mean hydrostatic atmosphere as a function of altitude */
  double *zcoord = solver->x + (solver->dim_local[0]+2*solver->ghosts)
                             + (solver->dim_local[1]+2*solver->ghosts);
  IERR Numa3DStandardAtmosphere(physics,zcoord,(solver->dim_local[2]+2*solver->ghosts)); 
  CHECKERR(ierr);

  /* initializing physical model-specific functions */
  solver->ComputeCFL  = Numa3DComputeCFL;
  solver->FFunction   = Numa3DFlux;
  solver->SFunction   = Numa3DSource;
  solver->Upwind      = Numa3DRusanov;

  /* set the value of gamma in all the boundary objects */
  int n;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  for (n = 0; n < solver->nBoundaryZones; n++)  boundary[n].gamma = physics->gamma;

  return(0);
}

int Numa3DStandardAtmosphere(void *p,double *z,int N)
{
  Numa3D *physics = (Numa3D*) p;

  double R      = physics->R;
  double gamma  = physics->gamma;
  double g      = physics->g;

  /* reference quantities at zero altitude */
  double rho_ref, P_ref, T_ref; 
  P_ref   = physics->Pref;
  T_ref   = physics->Tref;
  rho_ref = P_ref/(R*T_ref);

  double *rho = physics->rho0;
  double *P   = physics->P0;
  double *T   = physics->T0;

  int i;
  for (i=0; i<N; i++) {
    double zcoord = z[i];
    double pi = 1.0 - ((gamma-1.0)/gamma)*(g/(R*T_ref))*zcoord;
    double term = gamma/(gamma-1.0);
    
    rho[i] = rho_ref * raiseto(pi,term);
    P[i]   = P_ref   * raiseto(pi,term);
    T[i]   = T_ref;
  }
  return(0);
}
