#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ExactSolutionSerial    (void*, void*, double*, int*);
static int ExactSolutionParallel  (void*, void*, double*, int*);

int ExactSolution(void *s, void *m, double *uex, int *flag)
{
  HyPar  *solver = (HyPar*) s;
  if      (!strcmp(solver->input_mode,"serial"))    return(ExactSolutionSerial    (s,m,uex,flag));
  else if (!strcmp(solver->input_mode,"parallel"))  return(ExactSolutionParallel  (s,m,uex,flag));
  else {
    fprintf(stderr,"Error: Illegal value (%s) for input_mode (may be \"serial\" or \"parallel\"\n",
            solver->input_mode);
    return(1);
  }
}

int ExactSolutionSerial(void *s, void *m, double *uex, int *exact_flag)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i, d, ferr;
  double        *ug = NULL, *xg = NULL; 
  _DECLARE_IERR_;

  *exact_flag = 0;
  /* Only root process reads in exact solution file */
  if (!mpi->rank) {

    if (!strcmp(solver->ip_file_type,"ascii")) {

      FILE *in; in = fopen("exact.inp","r");
      if (!in) *exact_flag = 0;
      else {
        *exact_flag = 1;
        /* Reading exact solution */
        printf("Reading exact solution from ASCII file \"exact.inp\" (Serial mode).\n");
        int size,offset;
        /* allocate global solution array */
        size  = solver->npoints_global*solver->nvars;
        ug    = (double*) calloc(size,sizeof(double));
        size  = 0; for (d=0; d<solver->ndims; d++) size += solver->dim_global[d];
        xg    = (double*) calloc(size,sizeof(double));

        /* read grid (not necessary but to keep initial and exact solutions file formats same */
        offset = 0;
        for (d = 0; d < solver->ndims; d++) {
          for (i = 0; i < solver->dim_global[d]; i++) ferr = fscanf(in,"%lf",&xg[i+offset]);
          offset += solver->dim_global[d];
        }

        /* read solution */
        for (i = 0; i < solver->nvars; i++) {
          int *index = solver->index;
          int done = 0; _ArraySetValue_(index,solver->ndims,0);
          while (!done) {
            int p; _ArrayIndex1D_(solver->ndims,solver->dim_global,index,0,p);
            ferr = fscanf(in,"%lf",&ug[p*solver->nvars+i]);
            if (ferr != 1) return(1);
            _ArrayIncrementIndex_(solver->ndims,solver->dim_global,index,done);
          }
        }

        fclose(in);
      }

    } else if ((!strcmp(solver->ip_file_type,"bin")) || (!strcmp(solver->ip_file_type,"binary"))) {

      FILE *in; in = fopen("exact.inp","rb");
      if (!in) *exact_flag = 0;
      else {
        *exact_flag = 1;
        printf("Reading exact solution from binary file \"exact.inp\" (Serial mode).\n");
        size_t bytes;
        int size,offset;
        /* allocate global solution array */
        size = solver->npoints_global*solver->nvars;
        ug = (double*) calloc(size,sizeof(double));
        size = 0; for (d=0; d<solver->ndims; d++) size += solver->dim_global[d];
        xg      = (double*) calloc(size,sizeof(double));

        /* read grid */
        size = 0;
        for (d = 0; d < solver->ndims; d++) size += solver->dim_global[d];
        bytes = fread(xg, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ExactSolution(): Unable to read grid. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        /* read solution */
        size = 1;
        for (d = 0; d < solver->ndims; d++) size *= solver->dim_global[d]; size *= solver->nvars;
        bytes = fread(ug, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ExactSolution(): Unable to read solution. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        fclose(in);
      }

    }

  }

  /* Broadcast exact_flag to all processes */
  IERR MPIBroadcast_integer(exact_flag,1,0,&mpi->world); CHECKERR(ierr);

  if (*exact_flag) {
    /* partition global exact solution array to all processes */
    IERR MPIPartitionArraynD(solver->ndims,mpi,(mpi->rank?NULL:ug),uex,
                             solver->dim_global,solver->dim_local,
                             solver->ghosts,solver->nvars); CHECKERR(ierr);
    /* free global arrays */
    if (!mpi->rank) {
      free(ug);
      free(xg);
    }
  }

  return(0);
}

int ExactSolutionParallel(void *s, void *m, double *uex, int *exact_flag)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  *exact_flag = 0;
  if (!strcmp(solver->ip_file_type,"ascii")) {

    fprintf(stderr,"Error in ExactSolutionParallel(): Reading ASCII exact solution file in parallel not yet supported. Choose the binary option.\n");
    return(1);

  } else if ((!strcmp(solver->ip_file_type,"binary")) || (!strcmp(solver->ip_file_type,"bin"))) {

    FILE *in; in = fopen("exact_par.inp","rb");
    if (!in) *exact_flag = 0;
    else {
      *exact_flag = 1;
      int myrank = mpi->rank, ndims  = solver->ndims, ghosts = solver->ghosts;
      int rank[ndims+1], size[ndims], flag = 0, nvars, ferr, total, n;

      double *x, *u;
      total = 0; for (n = 0; n < ndims; n++) total += solver->dim_local[n];
      x = (double*) calloc (total, sizeof(double));
      total = solver->nvars; for (n = 0; n < ndims; n++) total *= solver->dim_local[n];
      u = (double*) calloc (total, sizeof(double));

      if (!myrank) printf("Reading exact solution from binary file \"exact_par.inp\" (Parallel mode).\n");
      int done = 0;
      while((!feof(in)) && (!done)) {
        ferr = fread(rank,sizeof(int),ndims+1,in); 
        if (ferr != (ndims+1)) {
          fprintf(stderr,"Error (1) in reading binary file exact_par.inp in parallel mode on rank %d.\n", myrank);
          return(1);
        }
        ferr = fread(size,sizeof(int),ndims,in);
        if (ferr != ndims) {
          fprintf(stderr,"Error (2) in reading binary file exact_par.inp in parallel mode on rank %d.\n", myrank);
          return(1);
        }
        ferr = fread(&nvars,sizeof(int),1,in);
        if (ferr != 1) {
          fprintf(stderr,"Error (3) in reading binary file exact_par.inp in parallel mode on rank %d.\n", myrank);
          return(1);
        }
        if (rank[ndims] == myrank) {
          /* checks */
          int check = 1;
          for (n = 0; n < ndims; n++) {
            if (rank[n] != mpi->ip[n])            check = 0;
            if (size[n] != solver->dim_local[n])  check = 0;
          }
          if (nvars != solver->nvars) check = 0;
          if (!check) {
            fprintf(stderr,"Error in ExactSolutionParallel(): Inconsistent data read on rank %d.\n",myrank);
            return(1);
          }
          /* read grid */
          total = 0; for (n = 0; n < ndims; n++) total += size[n];
          ferr = fread(x,sizeof(double),total,in);
          if (ferr != total) {
            fprintf(stderr,"Error (4) in reading binary file exact_par.inp in parallel mode on rank %d.\n", myrank);
            return(1);
          }
          /* read exact solution */
          total = nvars; for (n = 0; n < ndims; n++) total *= size[n];
          ferr = fread(u,sizeof(double),total,in);
          if (ferr != total) {
            fprintf(stderr,"Error (5) in reading binary file exact_par.inp in parallel mode on rank %d.\n", myrank);
            return(1);
          }
          flag = 1;
          done = 1;
        } else {
          int n, offset1 = 0, offset2 = nvars;
          for (n = 0; n < ndims; n++) {
            offset1 += size[n];
            offset2 *= size[n];
          }
          ferr = fseek(in,sizeof(double)*(offset1+offset2),SEEK_CUR);
          if (ferr) {
            fprintf(stderr,"Error (6) in reading binary file exact_par.inp in parallel mode on rank %d.\n", myrank);
            return(1);
          }
        }
      }
      fclose(in);

      int index[ndims];
      IERR ArrayCopynD(ndims,u,uex,solver->dim_local,0,ghosts,index,solver->nvars); CHECKERR(ierr);

      free(x);
      free(u);
    }

    return(0);

  } else {

    fprintf(stderr,"Error in ExactSolutionParallel(): Illegal value (%s) for ip_file type. May be \"ascii\", \"binary\" or \"bin\".\n",
            solver->ip_file_type);
    return(1);

  }

}
