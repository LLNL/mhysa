#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int WENOInitialize(void *s,void *m, char *scheme)
{
  HyPar           *solver = (HyPar*) s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;

  /* default parameters */
  weno->mapped      = 0;
  weno->borges      = 0;
  weno->yc          = 0;
  weno->no_limiting = 0;
  weno->eps         = 1e-6;
  weno->p           = 2.0;

  weno->rc          = 0.3;
  weno->xi          = 0.001;

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
    }
  }

  int     integer_data[4];
  double  real_data[4];
  if (!mpi->rank) {
    integer_data[0] = weno->mapped;
    integer_data[1] = weno->borges;
    integer_data[2] = weno->yc;
    integer_data[3] = weno->no_limiting;
    real_data[0]    = weno->eps;
    real_data[1]    = weno->p;
    real_data[2]    = weno->rc;
    real_data[3]    = weno->xi;
  }
  MPIBroadcast_integer(integer_data,4,0,&mpi->world);
  MPIBroadcast_double (real_data   ,4,0,&mpi->world);
  
  weno->mapped      = integer_data[0];
  weno->borges      = integer_data[1];
  weno->yc          = integer_data[2];
  weno->no_limiting = integer_data[3];
  weno->eps         = real_data   [0];
  weno->p           = real_data   [1];
  weno->rc          = real_data   [2];
  weno->xi          = real_data   [3];

  if (   (!strcmp(scheme,_FIFTH_ORDER_CRWENO_))
      || (!strcmp(scheme,_FIFTH_ORDER_HCWENO_)) ) {
    int size = 1, d;
    for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+1);
    size *= solver->nvars;
    if (!strcmp(solver->interp_type,_CHARACTERISTIC_)) size *= solver->nvars;

    weno->A = (double*) calloc (size, sizeof(double));
    weno->B = (double*) calloc (size, sizeof(double));
    weno->C = (double*) calloc (size, sizeof(double));
    weno->R = (double*) calloc (size, sizeof(double));

    weno->sendbuf = (double*) calloc (size, sizeof(double));
    weno->recvbuf = (double*) calloc (size, sizeof(double));

  } else {

    weno->A = weno->B = weno->C = weno->R = NULL;
    weno->sendbuf = weno->recvbuf = NULL;

  }

  return(0);
}
