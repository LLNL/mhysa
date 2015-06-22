/*! @file ExactSolution.c
    @author Debojyoti Ghosh
    @brief Read in exact solution, if available.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ExactSolutionSerial    (void*, void*, double*, int*);
#ifndef serial
static int ExactSolutionParallel  (void*, void*, double*, int*);
static int ExactSolutionMPI_IO    (void*, void*, double*, int*);
#endif

/*! Read in the exact solution, if available: wrapper function that calls 
    the appropriate function depending on input mode (#HyPar::input_mode).
*/
int ExactSolution(
                    void    *s,     /*!< Solver object of type #HyPar */
                    void    *m,     /*!< MPI object of type #MPIVariables */
                    double  *uex,   /*!< Array to hold the exact solution, if available */
                    int     *flag   /*!< Flag to indicate if exact solution was available */
                 )
{
  HyPar  *solver = (HyPar*) s;
  if      (!strcmp(solver->input_mode,"serial"))    return(ExactSolutionSerial    (s,m,uex,flag));
#ifndef serial
  else if (!strcmp(solver->input_mode,"parallel"))  return(ExactSolutionParallel  (s,m,uex,flag));
  else if (!strcmp(solver->input_mode,"mpi-io"  ))  return(ExactSolutionMPI_IO    (s,m,uex,flag));
#endif
  else {
    fprintf(stderr,"Error: Illegal value (%s) for input_mode.\n",solver->input_mode);
    return(1);
  }
}

/*! Read the exact solution in a serial fashion: For multi-processor simulation, only rank 0
    reads in the entire solution from the file, and then distributes the relevant portions
    to each of the processors. This involves memory allocation for the global domain on rank
    0. Thus, do not use for large domains. This approach is not very scalable either, if running
    with a very large number of processors (> 1000).
    \n
    + Supports both binary and ASCII formats for the exact solution file.

    For serial runs, this is the only input mode (of course!).
    \sa InitialSolutionSerial(), OutputSolutionSerial()
*/
int ExactSolutionSerial(
                          void    *s,           /*!< Solver object of type #HyPar */
                          void    *m,           /*!< MPI object of type #MPIVariables */
                          double  *uex,         /*!< Array to hold the exact solution, if available */
                          int     *exact_flag   /*!< Flag to indicate if exact solution was available */
                       )
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

#ifndef serial

/*! Read in the exact solution in a parallel fashion: The number of MPI ranks participating in file I/O
    is specified as an input. All the MPI ranks are divided into that many I/O groups, with one rank in 
    each group as the "leader" that does the file reading and writing. For reading in the solution,
    the leader of an I/O group reads its own file and distributes the solution to the processors in 
    its group. The number of I/O group is typically specified as the number of I/O nodes available 
    on the HPC platform, given the number of compute nodes the code is running on. This is a good 
    balance between all the processors serially reading from the same file, and having as many 
    files (with the local solution) as the number of processors. This approach has been observed to 
    be very scalable (up to ~ 100,000 - 1,000,000 processors).
    \n
    + The file read by each I/O leader must contain the partitioned solution for each processor in
      its group. Use Extras/ParallelInput.c to generate these files from the standard exact/initial
      solution file containing the entire solution.
    + Supports only binary format for exact solution file.

    \sa InitialSolutionParallel(), OutputSolutionParallel()
*/
int ExactSolutionParallel(
                          void    *s,           /*!< Solver object of type #HyPar */
                          void    *m,           /*!< MPI object of type #MPIVariables */
                          double  *uex,         /*!< Array to hold the exact solution, if available */
                          int     *exact_flag   /*!< Flag to indicate if exact solution was available */
                         )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           proc,d;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars; 
  int ghosts = solver->ghosts;
  int *dim_local = solver->dim_local;

  *exact_flag = 1;

  /* check for existence of exact solution file */
  if (mpi->IOParticipant) {
    FILE *in;
    char filename[_MAX_STRING_SIZE_];
    char filename_root[_MAX_STRING_SIZE_] = "exact_par.inp";
    MPIGetFilename(filename_root,&mpi->IOWorld,filename);
    in = fopen(filename,"rb");
    if (!in)  *exact_flag = 0;
    else {
      *exact_flag = 1;
      fclose(in);
    }
  }
  IERR MPIMin_integer(exact_flag,exact_flag,1,&mpi->world);

  if (*exact_flag) {

    if (!mpi->rank) printf("Reading exact solution from binary file exact_par.inp.xxx (parallel mode).\n");

    /* calculate size of the local grid on this rank */
    int sizex = 0;     for (d=0; d<ndims; d++) sizex += dim_local[d];
    int sizeu = nvars; for (d=0; d<ndims; d++) sizeu *= dim_local[d];

    /* allocate buffer arrays to read in grid and solution */
    double *buffer = (double*) calloc (sizex+sizeu, sizeof(double));

    if (mpi->IOParticipant) {

      /* if this rank is responsible for file I/O */
      double *read_buffer = NULL;
      int     read_size_x, read_size_u, read_total_size;
      int     is[ndims], ie[ndims];

      /* open the file */
      FILE *in;
      int  bytes;
      char filename[_MAX_STRING_SIZE_];
      char filename_root[_MAX_STRING_SIZE_] = "exact_par.inp";
      MPIGetFilename(filename_root,&mpi->IOWorld,filename);

      in = fopen(filename,"rb");
      if (!in) {
        fprintf(stderr,"Error in ExactSolutionParallel(): File %s could not be opened.\n",filename);
        return(1);
      }

      /* Read own data */
      bytes = fread(buffer,sizeof(double),(sizex+sizeu),in);
      if (bytes != (sizex+sizeu)) {
        fprintf(stderr,"Error in ExactSolutionParallel(): File %s contains insufficient data.\n",filename);
        return(1);
      }

      /* read and send the data for the other processors in this IO rank's group */
      for (proc=mpi->GroupStartRank+1; proc<mpi->GroupEndRank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(ndims,proc,mpi,solver->dim_global,is,ie);
        /* calculate the size of its local data and allocate read buffer */
        read_size_x = 0;      for (d=0; d<ndims; d++) read_size_x += (ie[d]-is[d]);
        read_size_u = nvars;  for (d=0; d<ndims; d++) read_size_u *= (ie[d]-is[d]);
        read_total_size = read_size_x + read_size_u;
        read_buffer = (double*) calloc (read_total_size, sizeof(double));
        /* read the data */
        bytes = fread(read_buffer,sizeof(double),read_total_size,in);
        if (bytes != read_total_size) {
          fprintf(stderr,"Error in ExactSolutionParallel(): File %s contains insufficient data.\n",filename);
          return(1);
        }
        /* send the data */
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(read_buffer,read_total_size,MPI_DOUBLE,proc,1100,mpi->world,&req);
        MPI_Wait(&req,MPI_STATUS_IGNORE);
        free(read_buffer);
      }

      /* close the file */
      fclose(in);

    } else {

      /* all other processes, just receive the data from
       * the rank responsible for file I/O */
      MPI_Request req = MPI_REQUEST_NULL;
      MPI_Irecv(buffer,(sizex+sizeu),MPI_DOUBLE,mpi->IORank,1100,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

    }

    /* copy the solution */
    int index[ndims];
    IERR ArrayCopynD(ndims,(buffer+sizex),uex,dim_local,0,ghosts,index,nvars); 
    CHECKERR(ierr);

    /* free buffers */
    free(buffer);
  }

  return(0);

}

/*! Read in the exact solution in a parallel fashion using MPI-IO: Similar to ExactSolutionParallel(),
    except that the I/O leaders read from the file using the MPI I/O routines. These are constantly 
    being developed to be scalable on the latest and greatest HPC platforms.
    \n
    + Read the documentation for ExactSolutionParallel() first.
    + The file read by each I/O leader must contain the partitioned solution for each processor in
      its group. Use Extras/MPIInput.c to generate these files from the standard exact/initial
      solution file containing the entire solution.
    + Supports only binary format for exact solution file.

    \sa InitialSolutionMPI_IO
*/
int ExactSolutionMPI_IO(
                          void    *s,           /*!< Solver object of type #HyPar */
                          void    *m,           /*!< MPI object of type #MPIVariables */
                          double  *uex,         /*!< Array to hold the exact solution, if available */
                          int     *exact_flag   /*!< Flag to indicate if exact solution was available */
                       )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           proc,d;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars; 
  int ghosts = solver->ghosts;
  int *dim_local = solver->dim_local;

  *exact_flag = 0;

  /* check for existence of exact solution file */
  if (!mpi->rank) {
    FILE *in;
    in = fopen("exact_mpi.inp","rb");
    if (!in)  *exact_flag = 0;
    else {
      *exact_flag = 1;
      fclose(in);
    }
  }
  IERR MPIBroadcast_integer(exact_flag,1,0,&mpi->world);

  if (*exact_flag) {

    if (!mpi->rank) printf("Reading exact solution from binary file exact_mpi.inp (MPI-IO mode).\n");

    /* calculate size of the local grid on this rank */
    int sizex = 0;     for (d=0; d<ndims; d++) sizex += dim_local[d];
    int sizeu = nvars; for (d=0; d<ndims; d++) sizeu *= dim_local[d];

    /* allocate buffer arrays to read in grid and solution */
    double *buffer = (double*) calloc (sizex+sizeu, sizeof(double));

    if (mpi->IOParticipant) {

      /* if this rank is responsible for file I/O */
      double *read_buffer = NULL;
      int     read_size_x, read_size_u, read_total_size;
      int     is[ndims], ie[ndims], size;

      /* calculate offset */
      long long offset = 0;
      for (proc=0; proc < mpi->rank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(ndims,proc,mpi,solver->dim_global,is,ie);
        /* calculate the size of its local grid */
        size = 0; for (d=0; d<ndims; d++) size += (ie[d]-is[d]);
        offset += size;
        /* calculate the size of the local solution */
        size = nvars; for (d=0; d<ndims; d++) size *= (ie[d]-is[d]);
        offset += size;
      }
    
      /* open the file */
      MPI_Status  status;
      MPI_File    in;
      int         error;
      error = MPI_File_open(mpi->IOWorld,"exact_mpi.inp",MPI_MODE_RDONLY,MPI_INFO_NULL,&in);
      if (error != MPI_SUCCESS) {
        fprintf(stderr,"Error in ExactSolutionMPI_IO(): Unable to open exact_mpi.inp.\n");
        return(1);
      }

      /* set offset */
      MPI_Offset FileOffset = (MPI_Offset) (offset * sizeof(double));
      MPI_File_seek(in,FileOffset,MPI_SEEK_SET);

      /* Read own data */
      MPI_File_read(in,buffer,(sizex+sizeu)*sizeof(double),MPI_BYTE,&status);

      /* read and send the data for the other processors in this IO rank's group */
      for (proc=mpi->GroupStartRank+1; proc<mpi->GroupEndRank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(ndims,proc,mpi,solver->dim_global,is,ie);
        /* calculate the size of its local data and allocate read buffer */
        read_size_x = 0;      for (d=0; d<ndims; d++) read_size_x += (ie[d]-is[d]);
        read_size_u = nvars;  for (d=0; d<ndims; d++) read_size_u *= (ie[d]-is[d]);
        read_total_size = read_size_x + read_size_u;
        read_buffer = (double*) calloc (read_total_size, sizeof(double));
        /* read the data */
        MPI_File_read(in,read_buffer,read_total_size*sizeof(double),MPI_BYTE,&status);
        /* send the data */
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(read_buffer,read_total_size,MPI_DOUBLE,proc,1100,mpi->world,&req);
        MPI_Wait(&req,MPI_STATUS_IGNORE);
        free(read_buffer);
      }

      /* close the file */
      MPI_File_close(&in);

    } else {

      /* all other processes, just receive the data from
       * the rank responsible for file I/O */
      MPI_Request req = MPI_REQUEST_NULL;
      MPI_Irecv(buffer,(sizex+sizeu),MPI_DOUBLE,mpi->IORank,1100,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

    }

    /* copy the solution */
    int index[ndims];
    IERR ArrayCopynD(ndims,(buffer+sizex),uex,dim_local,0,ghosts,index,nvars); 
    CHECKERR(ierr);

    /* free buffers */
    free(buffer);
  }

  return(0);

}

#endif
