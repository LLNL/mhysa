/*! @file WriteArray.c
    @author Debojyoti Ghosh
    @brief Write an array to file

    Contains functions to write out an array (which is not
    a part of the solution) to a file in the same format 
    and mode as the solution output file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <timeintegration.h>
#include <hypar.h>

/* Function declarations */
static int WriteArraySerial   (int,int,int*,int*,int,double*,double*,void*,void*,char*);
static int WriteArrayParallel (int,int,int*,int*,int,double*,double*,void*,void*,char*);

/*! Write out an array data to file: wrapper function that calls 
    the appropriate function depending on output mode (#HyPar::output_mode).\n\n
    This function can be used to write out a variable array of a physical 
    quantity that is not a part of the solution (for example, the 
    topography in #ShallowWater1D). The array is written out in the same 
    format (#HyPar::op_file_format) and using the same mode (#HyPar::output_mode)
    as the solution. The global and local dimensions, number of ghost points, 
    number of spatial dimensions, and number of variables per grid point may be 
    different from the solution. 
    \sa OutputSolution()
*/
int WriteArray(
                int     ndims,        /*!< Number of spatial dimensions */
                int     nvars,        /*!< Number of variables per grid point */
                int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                int     ghosts,       /*!< Number of ghost points */
                double  *x,           /*!< Array of spatial coordinates */
                double  *u,           /*!< Array to write */
                void    *s,           /*!< Solver object of type #HyPar */
                void    *m,           /*!< MPI object of type #MPIVariables */
                char*   fname_root    /*!< Filename root (extension is added automatically). For unsteady output,
                                           a numerical index is added that is the same as for the solution output files. */
              )
{
  HyPar         *solver = (HyPar*)       s;
  MPIVariables  *mpi    = (MPIVariables*)m;
  _DECLARE_IERR_;
  
  /* if WriteOutput() is NULL, then return */
  if (!solver->WriteOutput) return(0);

  if (!strcmp(solver->output_mode,"serial")) {
    IERR WriteArraySerial(ndims,nvars,dim_global,dim_local,ghosts,x,u,
                          solver,mpi,fname_root); CHECKERR(ierr);
  } else {
    IERR WriteArrayParallel(ndims,nvars,dim_global,dim_local,ghosts,x,u,
                            solver,mpi,fname_root); CHECKERR(ierr);
  }
  
  return(0);
}

/*!
  Function to write out an array to a file in serial mode.
  It will allocate the global domain on rank 0, so do not
  use for big problems for which the entire global domain 
  will not fit on one node. This approach is also not very
  scalable.\n\n
  + See WriteArray() for what this is used for.
  + Supports both binary, and ascii-type formats.
  \sa ReadArraySerial(), OutputSolutionSerial()
*/
int WriteArraySerial(
                      int     ndims,        /*!< Number of spatial dimensions */
                      int     nvars,        /*!< Number of variables per grid point */
                      int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     ghosts,       /*!< Number of ghost points */
                      double  *x,           /*!< Array of spatial coordinates */
                      double  *u,           /*!< Array to write */
                      void    *s,           /*!< Solver object of type #HyPar */
                      void    *m,           /*!< MPI object of type #MPIVariables */
                      char*   fname_root    /*!< Filename root (extension is added automatically). For unsteady output,
                                                 a numerical index is added that is the same as for the solution output files. */
                    )
{
  HyPar         *solver = (HyPar*)       s;
  MPIVariables  *mpi    = (MPIVariables*)m;
  int           d;
  _DECLARE_IERR_;

  /* root process: allocate global output arrays */
  double *ug, *xg;
  if (!mpi->rank) {
    int size_global;

    size_global = 1;
    for (d=0; d<ndims; d++) size_global *= dim_global[d];
    ug = (double*) calloc (size_global*nvars,sizeof(double));
    _ArraySetValue_(ug,size_global*nvars,0.0);

    size_global = 0;
    for (d=0; d<ndims; d++) size_global += dim_global[d];
    xg = (double*) calloc (size_global,sizeof(double));
    _ArraySetValue_(xg,size_global,0.0); CHECKERR(ierr);

  } else {

    /* null pointers on non-root processes */
    ug = xg = NULL;

  }

  /* Assemble the local output arrays into the global output arrays */
  IERR MPIGatherArraynD(ndims,mpi,ug,u,dim_global,dim_local,
                        ghosts,nvars);  CHECKERR(ierr);
  int offset_global, offset_local;
  offset_global = offset_local = 0;
  for (d=0; d<ndims; d++) {
    IERR MPIGatherArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),
                            &x[offset_local+ghosts],
                            mpi->is[d],mpi->ie[d],dim_local[d],0); CHECKERR(ierr);
    offset_global += dim_global[d];
    offset_local  += dim_local [d] + 2*ghosts;
  }

  if (!mpi->rank) {
    /* write output file to disk */
    char filename[_MAX_STRING_SIZE_] = "";
    strcat(filename,fname_root);
    if (!strcmp(solver->op_overwrite,"no")) {
      strcat(filename,"_");
      strcat(filename,solver->filename_index);
    }
    strcat(filename,solver->solnfilename_extn);
    printf("Writing solution file %s.\n",filename);
    IERR solver->WriteOutput(ndims,nvars,dim_global,xg,ug,filename,
                             solver->index); CHECKERR(ierr);

    /* Clean up output arrays */
    free(xg);
    free(ug);
  }

  return(0);
}

/*!
  Function to write an array to file in parallel. See InitialSolutionParallel()
  to understand the logic of parallel I/O. This function uses the same I/O groups
  and ranks to write out the files. This approach to file output has been observed
  to be very scalable. \n
  + Supports only binary output
  + See WriteArray() to understand what this is used for.
  \sa InitialSolutionParallel, OutputSolutionParallel()
*/
int WriteArrayParallel(
                        int     ndims,        /*!< Number of spatial dimensions */
                        int     nvars,        /*!< Number of variables per grid point */
                        int     *dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                        int     *dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                        int     ghosts,       /*!< Number of ghost points */
                        double  *x,           /*!< Array of spatial coordinates */
                        double  *u,           /*!< Array to write */
                        void    *s,           /*!< Solver object of type #HyPar */
                        void    *m,           /*!< MPI object of type #MPIVariables */
                        char*   fname_root    /*!< Filename root (extension is added automatically). For unsteady output,
                                                   a numerical index is added that is the same as for the solution output files. */
                      )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           proc,d;
  _DECLARE_IERR_;

  static int count = 0;

  char filename_root[_MAX_STRING_SIZE_];
  strcpy(filename_root,fname_root);
  strcat(filename_root,solver->solnfilename_extn);
  if (!mpi->rank) printf("Writing solution file %s.xxxx (parallel mode).\n",filename_root);

  /* calculate size of the local grid on this rank */
  int sizex = 0;     for (d=0; d<ndims; d++) sizex += dim_local[d];
  int sizeu = nvars; for (d=0; d<ndims; d++) sizeu *= dim_local[d];

  /* allocate buffer arrays to write grid and solution */
  double *buffer = (double*) calloc (sizex+sizeu, sizeof(double));

  /* copy the grid to buffer */
  int offset1 = 0, offset2 = 0;
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_((x+offset1+ghosts),(buffer+offset2),dim_local[d]);
    offset1 += (dim_local[d]+2*ghosts);
    offset2 +=  dim_local[d];
  }

  /* copy the solution */
  int index[ndims];
  IERR ArrayCopynD(ndims,u,(buffer+sizex),dim_local,ghosts,0,index,nvars); CHECKERR(ierr);

  if (mpi->IOParticipant) {

    /* if this rank is responsible for file I/O */
    double *write_buffer = NULL;
    int     write_size_x, write_size_u, write_total_size;
    int     is[ndims], ie[ndims];

    /* open the file */
    FILE *out;
    int  bytes;
    char filename[_MAX_STRING_SIZE_];
    MPIGetFilename(filename_root,&mpi->IOWorld,filename);

    if (!strcmp(solver->op_overwrite,"no")) {
      if ((!count) && (!solver->restart_iter)) {
        /* open a new file, since this function is being called the first time 
           and this is not a restart run*/
        out = fopen(filename,"wb");
        if (!out) {
          fprintf(stderr,"Error in WriteArrayParallel(): File %s could not be opened for writing.\n",filename);
          return(1);
        }
      } else {
        /* append to existing file */
        out = fopen(filename,"ab");
        if (!out) {
          fprintf(stderr,"Error in WriteArrayParallel(): File %s could not be opened for appending.\n",filename);
          return(1);
        }
      }
    } else {
      /* write a new file / overwrite existing file */
      out = fopen(filename,"wb");
      if (!out) {
        fprintf(stderr,"Error in WriteArrayParallel(): File %s could not be opened for writing.\n",filename);
        return(1);
      }
    }
    count++;

    /* Write own data and free buffer */
    bytes = fwrite(buffer,sizeof(double),(sizex+sizeu),out);
    if (bytes != (sizex+sizeu)) {
      fprintf(stderr,"Error in WriteArrayParallel(): Failed to write data to file %s.\n",filename);
      return(1);
    }
    free(buffer);

    /* receive and write the data for the other processors in this IO rank's group */
    for (proc=mpi->GroupStartRank+1; proc<mpi->GroupEndRank; proc++) {
      /* get the local domain limits for process proc */
      IERR MPILocalDomainLimits(ndims,proc,mpi,dim_global,is,ie);
      /* calculate the size of its local data and allocate write buffer */
      write_size_x = 0;      for (d=0; d<ndims; d++) write_size_x += (ie[d]-is[d]);
      write_size_u = nvars;  for (d=0; d<ndims; d++) write_size_u *= (ie[d]-is[d]);
      write_total_size = write_size_x + write_size_u;
      write_buffer = (double*) calloc (write_total_size, sizeof(double));
      /* receive the data */
      MPI_Request req = MPI_REQUEST_NULL;
      MPI_Irecv(write_buffer,write_total_size,MPI_DOUBLE,proc,1449,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
      /* write the data */
      bytes = fwrite(write_buffer,sizeof(double),write_total_size,out);
      if (bytes != write_total_size) {
        fprintf(stderr,"Error in WriteArrayParallel(): Failed to write data to file %s.\n",filename);
        return(1);
      }
      free(write_buffer);
    }

    /* close the file */
    fclose(out);

  } else {

    /* all other processes, just send the data to the rank responsible for file I/O */
    MPI_Request req = MPI_REQUEST_NULL;
    MPI_Isend(buffer,(sizex+sizeu),MPI_DOUBLE,mpi->IORank,1449,mpi->world,&req);
    MPI_Wait(&req,MPI_STATUS_IGNORE);
    free(buffer);

  }

  return(0);
}
