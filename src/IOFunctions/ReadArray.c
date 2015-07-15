/*! @file ReadArray.c
    @author Debojyoti Ghosh
    @brief Read in array data from file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ReadArraySerial    (int,int,int*,int*,int,void*,void*,double*,char*,int*);
#ifndef serial
static int ReadArrayParallel  (int,int,int*,int*,int,void*,void*,double*,char*,int*);
static int ReadArrayMPI_IO    (int,int,int*,int*,int,void*,void*,double*,char*,int*);
#endif

/*! Read in an array data from file, if available: wrapper function that calls 
    the appropriate function depending on input mode (#HyPar::input_mode).\n\n
    This function can be used to read in a variable array that is required by
    a physical model, but is not a part of the solution (for example, the 
    topography in #ShallowWater1D). The array to be read will be stored with
    the same memory layout as the solution; however, global and local dimensions, 
    number of ghost points, number of spatial dimensions, and number of variables
    per grid point may be different.
    The mode and typeof input is same as 
    that specified for reading in the initial (and exact, if available) solution(s) 
    (#HyPar::input_mode, #HyPar::ip_file_type).
*/
int ReadArray(
              int     ndims,        /*!< Number of spatial dimensions */
              int     nvars,        /*!< Number of variables per grid point */
              int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
              int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
              int     ghosts,       /*!< Number of ghost points */
              void    *s,           /*!< Solver object of type #HyPar */
              void    *m,           /*!< MPI object of type #MPIVariables */
              double  *uex,         /*!< Array to hold the variable */
              char    *fname_root,  /*!< Filename root (If #HyPar::input_mode is "serial", the filename is <fname_root>.inp; 
                                         if #Hypar::input_mode is "parallel", the filename is <fname_root>_par.inp.xxx; 
                                         and if #HyPar::input_mode is "mpi-io", the filename is <fname_root>_mpi.inp.xxx; 
                                         where "xxx" is a 3 digit index denoting the MPI rank which will read the file.) */
              int     *read_flag    /*!< Flag to indicate if the file was read */
             )
{
  HyPar  *solver = (HyPar*) s;
  if      (!strcmp(solver->input_mode,"serial"))    return(ReadArraySerial  (ndims,nvars,dim_global,dim_local,ghosts,s,m,uex,fname_root,read_flag));
#ifndef serial
  else if (!strcmp(solver->input_mode,"parallel"))  return(ReadArrayParallel(ndims,nvars,dim_global,dim_local,ghosts,s,m,uex,fname_root,read_flag));
  else if (!strcmp(solver->input_mode,"mpi-io"  ))  return(ReadArrayMPI_IO  (ndims,nvars,dim_global,dim_local,ghosts,s,m,uex,fname_root,read_flag));
#endif
  else {
    fprintf(stderr,"Error: Illegal value (%s) for input_mode.\n",solver->input_mode);
    return(1);
  }
}

/*! Read an array in a serial fashion: For multi-processor simulation, only rank 0
    reads in the entire solution from the file, and then distributes the relevant portions
    to each of the processors. This involves memory allocation for the global domain on rank
    0. Thus, do not use for large domains. This approach is not very scalable either, if running
    with a very large number of processors (> 1000).
    \n
    + Supports both binary and ASCII formats.
    + See ReadArray() on what this is used for.

    For serial runs, this is the only input mode (of course!).
    \sa InitialSolutionSerial(), OutputSolutionSerial()
*/
int ReadArraySerial(
                      int     ndims,        /*!< Number of spatial dimensions */
                      int     nvars,        /*!< Number of variables per grid point */
                      int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     ghosts,       /*!< Number of ghost points */
                      void    *s,           /*!< Solver object of type #HyPar */
                      void    *m,           /*!< MPI object of type #MPIVariables */
                      double  *uex,         /*!< Array to hold the variable */
                      char    *fname_root,  /*!< Filename root (If #HyPar::input_mode is "serial", the filename is <fname_root>.inp; 
                                                 if #Hypar::input_mode is "parallel", the filename is <fname_root>_par.inp.xxx; 
                                                 and if #HyPar::input_mode is "mpi-io", the filename is <fname_root>_mpi.inp.xxx; 
                                                 where "xxx" is a 3 digit index denoting the MPI rank which will read the file.) */
                      int     *read_flag    /*!< Flag to indicate if the file was read */
                    )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i, d, ferr, index[ndims];
  double        *ug = NULL, *xg = NULL; 
  _DECLARE_IERR_;

  *read_flag = 0;
  /* Only root process reads from the file */
  if (!mpi->rank) {

    if (!strcmp(solver->ip_file_type,"ascii")) {
      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,fname_root);
      strcat(filename,".inp");
      FILE *in; in = fopen(filename,"r");
      if (!in) *read_flag = 0;
      else {
        *read_flag = 1;
        /* Reading from file */
        printf("Reading array from ASCII file %s (Serial mode).\n",filename);
        int size,offset;
        /* allocate global solution array */
        size  = 1; for (d=0; d<ndims; d++) size *= dim_global[d]; size *= nvars;
        ug    = (double*) calloc(size,sizeof(double));
        size  = 0; for (d=0; d<ndims; d++) size += dim_global[d];
        xg    = (double*) calloc(size,sizeof(double));

        /* read grid (not necessary but to keep the format same as initial and exact solutions files) */
        offset = 0;
        for (d = 0; d < ndims; d++) {
          for (i = 0; i < dim_global[d]; i++) ferr = fscanf(in,"%lf",&xg[i+offset]);
          offset += dim_global[d];
        }

        /* read solution */
        for (i = 0; i < nvars; i++) {
          int done = 0; _ArraySetValue_(index,ndims,0);
          while (!done) {
            int p; _ArrayIndex1D_(ndims,dim_global,index,0,p);
            ferr = fscanf(in,"%lf",&ug[p*nvars+i]);
            if (ferr != 1) return(1);
            _ArrayIncrementIndex_(ndims,dim_global,index,done);
          }
        }

        fclose(in);
      }
    } else if ((!strcmp(solver->ip_file_type,"bin")) || (!strcmp(solver->ip_file_type,"binary"))) {

      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,fname_root);
      strcat(filename,".inp");
      FILE *in; in = fopen(filename,"rb");
      if (!in) *read_flag = 0;
      else {
        *read_flag = 1;
        printf("Reading array from binary file %s (Serial mode).\n",filename);
        size_t bytes;
        int size;
        /* allocate global solution array */
        size  = 1; for (d=0; d<ndims; d++) size *= dim_global[d]; size *= nvars;
        ug = (double*) calloc(size,sizeof(double));
        size = 0; for (d=0; d<ndims; d++) size += dim_global[d];
        xg      = (double*) calloc(size,sizeof(double));

        /* read grid */
        size = 0;
        for (d = 0; d < ndims; d++) size += dim_global[d];
        bytes = fread(xg, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read grid. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        /* read solution */
        size = 1;
        for (d = 0; d < ndims; d++) size *= dim_global[d]; size *= nvars;
        bytes = fread(ug, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read solution. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        fclose(in);
      }

    }

  }

  /* Broadcast read_flag to all processes */
  IERR MPIBroadcast_integer(read_flag,1,0,&mpi->world); CHECKERR(ierr);

  if (*read_flag) {
    /* partition global array to all processes */
    IERR MPIPartitionArraynD(ndims,mpi,(mpi->rank?NULL:ug),uex,dim_global,
                             dim_local,ghosts,nvars); CHECKERR(ierr);
    /* free global arrays */
    if (!mpi->rank) {
      free(ug);
      free(xg);
    }
  }

  return(0);
}

#ifndef serial

/*! Read in an array in a parallel fashion: The number of MPI ranks participating in file I/O
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
      its group. Use Extras/ParallelInput.c to generate these files from the standard file containing 
      the entire solution.
    + Supports only binary format.
    + See ReadArray() on what this is used for.

    \sa InitialSolutionParallel(), OutputSolutionParallel()
*/
int ReadArrayParallel(
                      int     ndims,        /*!< Number of spatial dimensions */
                      int     nvars,        /*!< Number of variables per grid point */
                      int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     ghosts,       /*!< Number of ghost points */
                      void    *s,           /*!< Solver object of type #HyPar */
                      void    *m,           /*!< MPI object of type #MPIVariables */
                      double  *uex,         /*!< Array to hold the variable */
                      char    *fname_root,  /*!< Filename root (If #HyPar::input_mode is "serial", the filename is <fname_root>.inp; 
                                                 if #Hypar::input_mode is "parallel", the filename is <fname_root>_par.inp.xxx; 
                                                 and if #HyPar::input_mode is "mpi-io", the filename is <fname_root>_mpi.inp.xxx; 
                                                 where "xxx" is a 3 digit index denoting the MPI rank which will read the file.) */
                      int     *read_flag    /*!< Flag to indicate if file was read */
                     )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           proc,d;
  _DECLARE_IERR_;

  *read_flag = 1;
  char filename_root[_MAX_STRING_SIZE_];
  strcpy(filename_root,fname_root);
  strcat(filename_root,"_par.inp");

  /* check for existence of the file */
  if (mpi->IOParticipant) {
    FILE *in;
    char filename[_MAX_STRING_SIZE_];
    MPIGetFilename(filename_root,&mpi->IOWorld,filename);
    in = fopen(filename,"rb");
    if (!in)  *read_flag = 0;
    else {
      *read_flag = 1;
      fclose(in);
    }
  }
  IERR MPIMin_integer(read_flag,read_flag,1,&mpi->world);

  if (*read_flag) {

    if (!mpi->rank) printf("Reading from binary file %s.xxx (parallel mode).\n",filename_root);

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
      MPIGetFilename(filename_root,&mpi->IOWorld,filename);

      in = fopen(filename,"rb");
      if (!in) {
        fprintf(stderr,"Error in ReadArrayParallel(): File %s could not be opened.\n",filename);
        return(1);
      }

      /* Read own data */
      bytes = fread(buffer,sizeof(double),(sizex+sizeu),in);
      if (bytes != (sizex+sizeu)) {
        fprintf(stderr,"Error in ReadArrayParallel(): File %s contains insufficient data.\n",filename);
        return(1);
      }

      /* read and send the data for the other processors in this IO rank's group */
      for (proc=mpi->GroupStartRank+1; proc<mpi->GroupEndRank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie);
        /* calculate the size of its local data and allocate read buffer */
        read_size_x = 0;      for (d=0; d<ndims; d++) read_size_x += (ie[d]-is[d]);
        read_size_u = nvars;  for (d=0; d<ndims; d++) read_size_u *= (ie[d]-is[d]);
        read_total_size = read_size_x + read_size_u;
        read_buffer = (double*) calloc (read_total_size, sizeof(double));
        /* read the data */
        bytes = fread(read_buffer,sizeof(double),read_total_size,in);
        if (bytes != read_total_size) {
          fprintf(stderr,"Error in ReadArrayParallel(): File %s contains insufficient data.\n",filename);
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

/*! Read in an array in a parallel fashion using MPI-IO: Similar to ReadArrayParallel(),
    except that the I/O leaders read from the file using the MPI I/O routines. These are constantly 
    being developed to be scalable on the latest and greatest HPC platforms.
    \n
    + Read the documentation for ReadArrayParallel() first.
    + The file read by each I/O leader must contain the partitioned solution for each processor in
      its group. Use Extras/MPIInput.c to generate these files from the standard file containing 
      the entire data.
    + Supports only binary format.
    + See ReadArray() on what this is used for.

    \sa InitialSolutionMPI_IO
*/
int ReadArrayMPI_IO(
                      int     ndims,        /*!< Number of spatial dimensions */
                      int     nvars,        /*!< Number of variables per grid point */
                      int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     ghosts,       /*!< Number of ghost points */
                      void    *s,           /*!< Solver object of type #HyPar */
                      void    *m,           /*!< MPI object of type #MPIVariables */
                      double  *uex,         /*!< Array to hold the variable */
                      char    *fname_root,  /*!< Filename root (If #HyPar::input_mode is "serial", the filename is <fname_root>.inp; 
                                                 if #Hypar::input_mode is "parallel", the filename is <fname_root>_par.inp.xxx; 
                                                 and if #HyPar::input_mode is "mpi-io", the filename is <fname_root>_mpi.inp.xxx; 
                                                 where "xxx" is a 3 digit index denoting the MPI rank which will read the file.) */
                      int     *read_flag    /*!< Flag to indicate if file was read */
                   )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           proc,d;
  _DECLARE_IERR_;

  *read_flag = 0;
  char filename[_MAX_STRING_SIZE_];
  strcpy(filename,fname_root);
  strcat(filename,"_mpi.inp");

  /* check for existence of file */
  if (!mpi->rank) {
    FILE *in;
    in = fopen(filename,"rb");
    if (!in)  *read_flag = 0;
    else {
      *read_flag = 1;
      fclose(in);
    }
  }
  IERR MPIBroadcast_integer(read_flag,1,0,&mpi->world);

  if (*read_flag) {

    if (!mpi->rank) printf("Reading from binary file %s (MPI-IO mode).\n",filename);

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
        IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie);
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
      error = MPI_File_open(mpi->IOWorld,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&in);
      if (error != MPI_SUCCESS) {
        fprintf(stderr,"Error in ReadArrayMPI_IO(): Unable to open %s.\n",filename);
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
        IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie);
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
