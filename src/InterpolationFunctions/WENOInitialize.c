/*! @file WENOInitialize.c
    @brief Initializes the WENO-type schemes
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int WENOFifthOrderCalculateWeights    (double*,double*,double*,int,void*,void*);
int WENOFifthOrderCalculateWeightsChar(double*,double*,double*,int,void*,void*);
int WENOFifthOrderInitializeWeights   (int,void*,void*);

/*!
  This function initializes the WENO-type methods.
  + Sets the parameters to default values.
  + Reads in the parameters from optional input file "weno.inp", if available.
  + Allocates memory for and initializes the nonlinear weights used by WENO-type
    schemes.
*/
int WENOInitialize(
                    void *s,      /*!< Solver object of type #HyPar */
                    void *m,      /*!< MPI object of type #MPIVariables */
                    char *scheme, /*!< Name of scheme */
                    char *type    /*!< Type of interpolation */
                  )
{
  HyPar           *solver = (HyPar*) s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;

  int nvars = solver->nvars;
  int ndims = solver->ndims;

  /* default parameters */
  weno->mapped      = 0;
  weno->borges      = 0;
  weno->yc          = 0;
  weno->no_limiting = 0;
  weno->eps         = 1e-6;
  weno->p           = 2.0;

  weno->rc          = 0.3;
  weno->xi          = 0.001;
  weno->tol         = 1e-16;

  if (!mpi->rank) {
    FILE *in;
    int ferr;
    in = fopen("weno.inp","r");
    if (!in) printf("Warning: File weno.inp not found. Using default parameters for WENO5/CRWENO5/HCWENO5 scheme.\n");
    else {
      printf("Reading WENO parameters from weno.inp.\n");
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if      (!strcmp(word,"mapped"     )) { ferr = fscanf(in,"%d" ,&weno->mapped     ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"borges"     )) { ferr = fscanf(in,"%d" ,&weno->borges     ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"yc"         )) { ferr = fscanf(in,"%d" ,&weno->yc         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"no_limiting")) { ferr = fscanf(in,"%d" ,&weno->no_limiting); if (ferr != 1) return(1); }
          else if (!strcmp(word,"epsilon"    )) { ferr = fscanf(in,"%lf",&weno->eps        ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"p"          )) { ferr = fscanf(in,"%lf",&weno->p          ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"rc"         )) { ferr = fscanf(in,"%lf",&weno->rc         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"xi"         )) { ferr = fscanf(in,"%lf",&weno->xi         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"tol"        )) { ferr = fscanf(in,"%lf",&weno->tol        ); if (ferr != 1) return(1); }
          else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"weno.inp\" with value %s not ",word,useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
	    } else {
    	  fprintf(stderr,"Error: Illegal format in file \"weno.inp\".\n");
        return(1);
	    }
      fclose(in);
    }
  }

  int     integer_data[4];
  double  real_data[5];
  if (!mpi->rank) {
    integer_data[0] = weno->mapped;
    integer_data[1] = weno->borges;
    integer_data[2] = weno->yc;
    integer_data[3] = weno->no_limiting;
    real_data[0]    = weno->eps;
    real_data[1]    = weno->p;
    real_data[2]    = weno->rc;
    real_data[3]    = weno->xi;
    real_data[4]    = weno->tol;
  }
  MPIBroadcast_integer(integer_data,4,0,&mpi->world);
  MPIBroadcast_double (real_data   ,5,0,&mpi->world);
  
  weno->mapped      = integer_data[0];
  weno->borges      = integer_data[1];
  weno->yc          = integer_data[2];
  weno->no_limiting = integer_data[3];
  weno->eps         = real_data   [0];
  weno->p           = real_data   [1];
  weno->rc          = real_data   [2];
  weno->xi          = real_data   [3];
  weno->tol         = real_data   [4];

  /* WENO weight calculation is hard-coded for p=2, so return error if p != 2 in
   * user input file, so that there's no confusion */
  if (weno->p != 2.0) {
    if (!mpi->rank) printf("Warning from WENOInitialize(): \"p\" parameter is 2.0. Any other value will be ignored!\n");
  }

  weno->offset = (int*) calloc (ndims,sizeof(int));
  int dir,d;
  for (dir=0; dir<ndims; dir++) {
    weno->offset[dir] = 0;
    for (d=0; d<dir; d++) {
      int size = nvars, i;
      for (i=0; i<ndims; i++) 
        size *= ( i==d ? solver->dim_local[i]+1 : solver->dim_local[i] );
      weno->offset[dir] += size;
    }
  }

  int total_size = 0;
  for (d=0; d<ndims; d++) {
    int size = nvars, i;
    for (i=0; i<ndims; i++) 
      size *= ( i==d ? solver->dim_local[i]+1 : solver->dim_local[i] );
    total_size += size;
  }
  weno->size = total_size;

  weno->w1 = (double*) calloc (4*total_size,sizeof(double));
  weno->w2 = (double*) calloc (4*total_size,sizeof(double));
  weno->w3 = (double*) calloc (4*total_size,sizeof(double));

  if ((!strcmp(type,_CHARACTERISTIC_)) && (nvars > 1)) 
    solver->SetInterpLimiterVar = WENOFifthOrderCalculateWeightsChar;
  else solver->SetInterpLimiterVar = WENOFifthOrderCalculateWeights;

  /* initialize WENO weights to their optimal values */
  for (d=0; d<ndims; d++) WENOFifthOrderInitializeWeights(d,solver,mpi);

  return(0);
}
