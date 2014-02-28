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
        int size;
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

    int     ndims = solver->ndims, d, n, *recv_int;
    double  *recv_double;

    /* allocate arrays to receive data in */
    int sizex = 0;             for (d=0; d<ndims; d++) sizex += solver->dim_local[d];
    int sizeu = solver->nvars; for (d=0; d<ndims; d++) sizeu *= solver->dim_local[d];
    int total_size = sizex + sizeu;
    recv_int    = (int*)    calloc (2*(ndims+1),sizeof(int  ));
    recv_double = (double*) calloc (total_size,sizeof(double));

    /* check if exact solution file exists */
    if (!mpi->rank) {
      FILE *in;
      in = fopen("exact_par.inp","rb");
      if (in) { *exact_flag = 1; fclose(in); }
      else    *exact_flag = 0;
    }
    /* Broadcast exact_flag to all processes */
    IERR MPIBroadcast_integer(exact_flag,1,0,&mpi->world); CHECKERR(ierr);

    if (*exact_flag) {

      if (!mpi->rank) {

        /* rank 0 opens the file, reads in the data and sends it to corresponding process */
        FILE  *in;
        int   ferr;
        int *flag = (int*) calloc (mpi->nproc, sizeof(int));
        for (n=0; n<mpi->nproc; n++) flag[n] = 0;
        in = fopen("exact_par.inp","rb");

        int count = 0;
        printf("Reading exact solution from binary file exact_par.inp (parallel mode).\n");
        while ((!feof(in)) && (count < mpi->nproc)) {
          int rank[ndims+1],size[ndims],nvars;
          ferr = fread(rank,sizeof(int),ndims+1,in);
          if (ferr != (ndims+1)) {
            fprintf(stderr,"Error (1) in reading binary file exact_par.inp in parallel mode on rank %d.\n", mpi->rank);
            return(1);
          }
          ferr = fread(size,sizeof(int),ndims,in);
          if (ferr != ndims) {
            fprintf(stderr,"Error (2) in reading binary file exact_par.inp in parallel mode on rank %d.\n", mpi->rank);
            return(1);
          }
          ferr = fread(&nvars,sizeof(int),1,in);
          if (ferr != 1) {
            fprintf(stderr,"Error (3) in reading binary file exact_par.inp in parallel mode on rank %d.\n", mpi->rank);
            return(1);
          }
          int sizex = 0;    for (d=0; d<ndims; d++) sizex += size[d];
          int sizeu = nvars;for (d=0; d<ndims; d++) sizeu *= size[d];
          int total_size = sizex + sizeu;
          if (!rank[ndims]) {
            /* self */
            ferr = fread(recv_double,sizeof(double),total_size,in);
            if (ferr != total_size) {
              fprintf(stderr,"Error (4) in reading binary file exact_par.inp in parallel mode on rank %d.\n", mpi->rank);
              return(1);
            }
            for (d=0; d<ndims+1; d++) recv_int[d]         = rank[d];
            for (d=0; d<ndims  ; d++) recv_int[d+ndims+1] = size[d];
            recv_int[2*ndims+1] = nvars;
          } else {
#ifndef serial
            /* allocate data array to send */
            int    *send_int     = (int*)     calloc (2*(ndims+1),sizeof(int   ));
            double *send_double  = (double*)  calloc (total_size ,sizeof(double));
            ferr = fread(send_double,sizeof(double),total_size,in);
            if (ferr != total_size) {
              fprintf(stderr,"Error (4) in reading binary file exact_par.inp in parallel mode on rank %d.\n", mpi->rank);
              return(1);
            }
            for (d=0; d<ndims+1; d++) send_int[d]         = rank[d];
            for (d=0; d<ndims  ; d++) send_int[d+ndims+1] = size[d];
            send_int[2*ndims+1] = nvars;
            /* send to the corresponding rank */
            MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
            MPI_Isend(send_int   ,2*(ndims+1),MPI_INT   ,rank[ndims],1247,mpi->world,&req[0]);
            MPI_Isend(send_double,total_size ,MPI_DOUBLE,rank[ndims],1248,mpi->world,&req[1]);
            MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);
            free(send_int);
            free(send_double);
#else
            fprintf(stderr,"Error in ExactSolutionParallel(): Code should not have reached here!\n");
#endif
          }
          flag[rank[ndims]] = 1;
          count++;
        }

        /* check if data for all the processes have been read */
        int check = 1;
        for (n = 0; n < mpi->nproc; n++) {
          if (!flag[n]) {
            fprintf(stderr,"Error in ExactSolutionParallel(): Data for rank %d not found.\n",n);
            check = 0;
          }
        }
        if (!check) return(1);
        free(flag);

        fclose(in);

      } else {

        /* receive the stuff from rank 0 */
        MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
        MPI_Irecv(recv_int   ,2*(ndims+1),MPI_INT   ,0,1247,mpi->world,&req[0]);
        MPI_Irecv(recv_double,total_size ,MPI_DOUBLE,0,1248,mpi->world,&req[1]);
        MPI_Waitall(2,&req[0],MPI_STATUS_IGNORE);

      }

      /* checks */
      int check = 1;
      for (n = 0; n < ndims; n++) {
        if (recv_int[n]         != mpi->ip[n]          )  check = 0;
        if (recv_int[n+ndims+1] != solver->dim_local[n])  check = 0;
      }
      if (recv_int[2*ndims+1] != solver->nvars) check = 0;
      if (!check) {
        fprintf(stderr,"Error in ExactSolutionParallel(): Inconsistent data read on rank %d.\n",mpi->rank);
        return(1);
      }

      int index[ndims];
      IERR ArrayCopynD(ndims,recv_double+sizex,uex,solver->dim_local,0,solver->ghosts,index,solver->nvars); CHECKERR(ierr);

      free(recv_int);
      free(recv_double);

    }
    return(0);

  } else {

    fprintf(stderr,"Error in ExactSolutionParallel(): Illegal value (%s) for ip_file type. May be \"ascii\", \"binary\" or \"bin\".\n",
            solver->ip_file_type);
    return(1);

  }
}
