#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

int VolumeIntegral(double*,double*,void*,void*);

static int InitialSolutionSerial    (void*, void*);
#ifndef serial
static int InitialSolutionParallel  (void*, void*);
static int InitialSolutionMPI_IO    (void*, void*);
#endif

int InitialSolution(void *s, void *m)
{
  HyPar  *solver = (HyPar*) s;
  
  if      (!strcmp(solver->input_mode,"serial"  ))  return(InitialSolutionSerial  (s,m));
#ifndef serial
  else if (!strcmp(solver->input_mode,"parallel"))  return(InitialSolutionParallel(s,m));
  else if (!strcmp(solver->input_mode,"mpi-io"  ))  return(InitialSolutionMPI_IO  (s,m));
#endif
  else {
    fprintf(stderr,"Error: Illegal value (%s) for input_mode.\n",
            solver->input_mode);
    return(1);
  }
}

int InitialSolutionSerial(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i,d, ferr;
  int           offset_global, offset_local, total_size;
  _DECLARE_IERR_;

  /* Only root process reads in initial solution file */
  double *ug,*xg,*dxinvg; /* global solution vector and grid arrays */
  if (!mpi->rank) {
    int size,offset;
    /* allocate global solution array */
    size = solver->npoints_global*solver->nvars;
    ug = (double*) calloc(size,sizeof(double));
    size = 0;
    for (d=0; d<solver->ndims; d++) size += solver->dim_global[d];
    xg      = (double*) calloc(size,sizeof(double));
    dxinvg  = (double*) calloc(size,sizeof(double));

    if (!strcmp(solver->ip_file_type,"ascii")) {

      /* Reading grid and initial solution */
      printf("Reading grid and initial conditions from ASCII file \"initial.inp\" (Serial mode).\n");
      FILE *in; in = fopen("initial.inp","r");
      if (!in) {
        fprintf(stderr,"Error: initial solution file \"initial.inp\" not found.\n");
        return(1);
      }

      /* read grid and calculate dxinv*/
      offset = 0;
      for (d = 0; d < solver->ndims; d++) {
        for (i = 0; i < solver->dim_global[d]; i++) ferr = fscanf(in,"%lf",&xg[i+offset]);
        for (i = 0; i < solver->dim_global[d]; i++) {
          if      (i == 0)                        dxinvg[i+offset] = 1.0/(xg[i+1+offset]-xg[i  +offset]);
          else if (i == solver->dim_global[d]-1)  dxinvg[i+offset] = 1.0/(xg[i  +offset]-xg[i-1+offset]);
          else                                    dxinvg[i+offset] = 2.0/(xg[i+1+offset]-xg[i-1+offset]);
        }
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

    } else if ((!strcmp(solver->ip_file_type,"bin")) || (!strcmp(solver->ip_file_type,"binary"))) {

      size_t bytes;
      /* Reading grid and initial solution */
      printf("Reading grid and initial conditions from binary file \"initial.inp\" (Serial mode).\n");
      FILE *in; in = fopen("initial.inp","rb");
      if (!in) {
        fprintf(stderr,"Error: initial solution file \"initial.inp\" not found.\n");
        return(1);
      }

      /* read grid */
      total_size = 0;
      for (d = 0; d < solver->ndims; d++) total_size += solver->dim_global[d];
      bytes = fread(xg, sizeof(double), total_size, in);
      if ((int)bytes != total_size) {
        fprintf(stderr,"Error in InitialSolution(): Unable to read grid. Expected %d, Read %d.\n",
                total_size, (int)bytes);
      }

      /* read solution */
      total_size = 1;
      for (d = 0; d < solver->ndims; d++) total_size *= solver->dim_global[d]; total_size *= solver->nvars;
      bytes = fread(ug, sizeof(double), total_size, in);
      if ((int)bytes != total_size) {
        fprintf(stderr,"Error in InitialSolution(): Unable to read solution. Expected %d, Read %d.\n",
                total_size, (int)bytes);
      }

      /* calculate dxinv*/
      offset = 0;
      for (d = 0; d < solver->ndims; d++) {
        for (i = 0; i < solver->dim_global[d]; i++) {
          if      (i == 0)                        dxinvg[i+offset] = 1.0/(xg[i+1+offset]-xg[i  +offset]);
          else if (i == solver->dim_global[d]-1)  dxinvg[i+offset] = 1.0/(xg[i  +offset]-xg[i-1+offset]);
          else                                    dxinvg[i+offset] = 2.0/(xg[i+1+offset]-xg[i-1+offset]);
        }
        offset += solver->dim_global[d];
      }

      fclose(in);
    }

  } else {
    ug      = NULL;
    xg      = NULL;
    dxinvg  = NULL;
  }

  /* partition initial solution across the processes */
  IERR MPIPartitionArraynD(solver->ndims,mpi,(mpi->rank?NULL:ug),solver->u,
                             solver->dim_global,solver->dim_local,
                             solver->ghosts,solver->nvars); CHECKERR(ierr);

  /* partition x vector across the processes */
  offset_global = offset_local = 0;
  for (d=0; d<solver->ndims; d++) {
    IERR MPIPartitionArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),
                                    &solver->x[offset_local+solver->ghosts],
                                    mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    IERR MPIPartitionArray1D(mpi,(mpi->rank?NULL:&dxinvg[offset_global]),
                                    &solver->dxinv[offset_local+solver->ghosts],
                                    mpi->is[d],mpi->ie[d],solver->dim_local[d],0); CHECKERR(ierr);
    offset_global += solver->dim_global[d];
    offset_local  += solver->dim_local [d] + 2*solver->ghosts;
  }

  /* exchange MPI-boundary values of x between processors */
  offset_local = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->x[offset_local],solver->dim_local[d],
                                   solver->ghosts,d,solver->ndims); CHECKERR(ierr);
    IERR MPIExchangeBoundaries1D(mpi,&solver->dxinv[offset_local],solver->dim_local[d],
                                   solver->ghosts,d,solver->ndims); CHECKERR(ierr);
    offset_local  += solver->dim_local [d] + 2*solver->ghosts;
  }

  /* fill in ghost values of x and dxinv at physical boundaries by extrapolation */
  offset_local = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *x     = &solver->x    [offset_local];
    double *dxinv = &solver->dxinv[offset_local];
    int    ghosts = solver->ghosts;
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) {
        int delta = ghosts - i;
        dxinv[i] = dxinv[ghosts];
        x[i] = x[ghosts] + ((double) delta) * (x[ghosts]-x[ghosts+1]);
      }
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
        int delta = i - (dim[d]+ghosts-1);
        dxinv[i] = dxinv[dim[d]+ghosts-1];
        x[i] =  x[dim[d]+ghosts-1] 
              + ((double) delta) * (x[dim[d]+ghosts-1]-x[dim[d]+ghosts-2]);
      }
    }
    offset_local  += dim[d] + 2*ghosts;
  }

  if (!mpi->rank) {
    free(ug);
    free(xg);
    free(dxinvg);
  }

  /* calculate volume integral of the initial solution */
  IERR VolumeIntegral(solver->VolumeIntegralInitial,solver->u,solver,mpi); CHECKERR(ierr);
  if (!mpi->rank) {
    printf("Volume integral of the initial solution:\n");
    for (d=0; d<solver->nvars; d++) printf("%4d:\t%1.16E\n",d,solver->VolumeIntegralInitial[d]);
  }
  /* Set initial total boundary flux integral to zero */
  _ArraySetValue_(solver->TotalBoundaryIntegral,solver->nvars,0);

  return(0);
}

#ifndef serial

int InitialSolutionParallel(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i,proc,d;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int ghosts = solver->ghosts;
  int *dim_local = solver->dim_local;

  if (!strcmp(solver->ip_file_type,"ascii")) {

    fprintf(stderr,"Error in InitialSolutionParallel(): Input file type must be binary.\n");
    return(1);

  } else if ((!strcmp(solver->ip_file_type,"binary")) || (!strcmp(solver->ip_file_type,"bin"))) {

    if (!mpi->rank) printf("Reading initial solution from binary file initial_par.inp.xxxx (parallel mode).\n");

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

      /* open the file */
      FILE *in;
      int  bytes;
      char filename[_MAX_STRING_SIZE_];
      MPIGetFilename("initial_par.inp",&mpi->IOWorld,filename);

      in = fopen(filename,"rb");
      if (!in) {
        fprintf(stderr,"Error in InitialSolutionParallel(): File %s could not be opened.\n",filename);
        return(1);
      }

      /* Read own data */
      bytes = fread(buffer,sizeof(double),(sizex+sizeu),in);
      if (bytes != (sizex+sizeu)) {
        fprintf(stderr,"Error in InitialSolutionParallel(): File %s contains insufficient data.\n",filename);
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
          fprintf(stderr,"Error in InitialSolutionParallel(): File %s contains insufficient data.\n",filename);
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

    /* copy the grid */
    int offset1 = 0, offset2 = 0;
    for (d = 0; d < ndims; d++) {
      _ArrayCopy1D_((buffer+offset2),(solver->x+offset1+ghosts),dim_local[d]);
      offset1 += (solver->dim_local[d]+2*ghosts);
      offset2 +=  solver->dim_local[d];
    }

    /* copy the solution */
    int index[ndims];
    IERR ArrayCopynD(ndims,(buffer+sizex),solver->u,solver->dim_local,0,ghosts,index,solver->nvars); CHECKERR(ierr);

    /* free buffers */
    free(buffer);

  } else {

    fprintf(stderr,"Error in InitialSolutionParallel(): Illegal value (%s) for ip_file type. May be \"ascii\", \"binary\" or \"bin\".\n",
            solver->ip_file_type);
    return(1);

  }

  int offset;

  /* exchange MPI-boundary values of x between processors */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->x[offset],solver->dim_local[d],
                                   ghosts,d,solver->ndims); CHECKERR(ierr);
    offset  += (solver->dim_local [d] + 2*ghosts);
  }
  /* fill in ghost values of x at physical boundaries by extrapolation */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *x     = &solver->x    [offset];
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) {
        int delta = ghosts - i;
        x[i] = x[ghosts] + ((double) delta) * (x[ghosts]-x[ghosts+1]);
      }
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
        int delta = i - (dim[d]+ghosts-1);
        x[i] =  x[dim[d]+ghosts-1] 
              + ((double) delta) * (x[dim[d]+ghosts-1]-x[dim[d]+ghosts-2]);
      }
    }
    offset  += (dim[d] + 2*ghosts);
  }

  /* calculate dxinv */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    for (i = 0; i < solver->dim_local[d]; i++) 
      solver->dxinv[i+offset+ghosts] = 2.0/(solver->x[i+1+offset+ghosts]-solver->x[i-1+offset+ghosts]);
    offset += (solver->dim_local[d] + 2*ghosts);
  }

  /* exchange MPI-boundary values of dxinv between processors */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->dxinv[offset],solver->dim_local[d],
                                   ghosts,d,solver->ndims); CHECKERR(ierr);
    offset  += (solver->dim_local[d] + 2*ghosts);
  }

  /* fill in ghost values of dxinv at physical boundaries by extrapolation */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *dxinv = &solver->dxinv[offset];
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) dxinv[i] = dxinv[ghosts];
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) dxinv[i] = dxinv[dim[d]+ghosts-1];
    }
    offset  += (dim[d] + 2*ghosts);
  }

  /* calculate volume integral of the initial solution */
  IERR VolumeIntegral(solver->VolumeIntegralInitial,solver->u,solver,mpi); CHECKERR(ierr);
  if (!mpi->rank) {
    printf("Volume integral of the initial solution:\n");
    for (d=0; d<solver->nvars; d++) printf("%2d:  %1.16E\n",d,solver->VolumeIntegralInitial[d]);
  }
  /* Set initial total boundary flux integral to zero */
  _ArraySetValue_(solver->TotalBoundaryIntegral,solver->nvars,0);

  return(0);
}

int InitialSolutionMPI_IO(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i,proc,d;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int nvars = solver->nvars;
  int ghosts = solver->ghosts;
  int *dim_local = solver->dim_local;

  if (!strcmp(solver->ip_file_type,"ascii")) {

    fprintf(stderr,"Error in InitialSolutionMPI_IO(): Input file type must be binary.\n");
    return(1);

  } else if ((!strcmp(solver->ip_file_type,"binary")) || (!strcmp(solver->ip_file_type,"bin"))) {

    if (!mpi->rank) printf("Reading initial solution from binary file initial_mpi.inp (MPI-IO mode).\n");

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
      error = MPI_File_open(mpi->IOWorld,"initial_mpi.inp",MPI_MODE_RDONLY,MPI_INFO_NULL,&in);
      if (error != MPI_SUCCESS) {
        fprintf(stderr,"Error in InitialSolutionMPI_IO(): File initial_mpi.inp could not be opened.\n");
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

    /* copy the grid */
    int offset1 = 0, offset2 = 0;
    for (d = 0; d < ndims; d++) {
      _ArrayCopy1D_((buffer+offset2),(solver->x+offset1+ghosts),dim_local[d]);
      offset1 += (solver->dim_local[d]+2*ghosts);
      offset2 +=  solver->dim_local[d];
    }

    /* copy the solution */
    int index[ndims];
    IERR ArrayCopynD(ndims,(buffer+sizex),solver->u,solver->dim_local,0,ghosts,index,solver->nvars); CHECKERR(ierr);

    /* free buffers */
    free(buffer);

  } else {

    fprintf(stderr,"Error in InitialSolutionParallel(): Illegal value (%s) for ip_file type. May be \"ascii\", \"binary\" or \"bin\".\n",
            solver->ip_file_type);
    return(1);

  }

  int offset;

  /* exchange MPI-boundary values of x between processors */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->x[offset],solver->dim_local[d],
                                   ghosts,d,solver->ndims); CHECKERR(ierr);
    offset  += (solver->dim_local [d] + 2*ghosts);
  }
  /* fill in ghost values of x at physical boundaries by extrapolation */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *x     = &solver->x    [offset];
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) {
        int delta = ghosts - i;
        x[i] = x[ghosts] + ((double) delta) * (x[ghosts]-x[ghosts+1]);
      }
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
        int delta = i - (dim[d]+ghosts-1);
        x[i] =  x[dim[d]+ghosts-1] 
              + ((double) delta) * (x[dim[d]+ghosts-1]-x[dim[d]+ghosts-2]);
      }
    }
    offset  += (dim[d] + 2*ghosts);
  }

  /* calculate dxinv */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    for (i = 0; i < solver->dim_local[d]; i++) 
      solver->dxinv[i+offset+ghosts] = 2.0/(solver->x[i+1+offset+ghosts]-solver->x[i-1+offset+ghosts]);
    offset += (solver->dim_local[d] + 2*ghosts);
  }

  /* exchange MPI-boundary values of dxinv between processors */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    IERR MPIExchangeBoundaries1D(mpi,&solver->dxinv[offset],solver->dim_local[d],
                                   ghosts,d,solver->ndims); CHECKERR(ierr);
    offset  += (solver->dim_local[d] + 2*ghosts);
  }

  /* fill in ghost values of dxinv at physical boundaries by extrapolation */
  offset = 0;
  for (d = 0; d < solver->ndims; d++) {
    double *dxinv = &solver->dxinv[offset];
    int    *dim   = solver->dim_local;
    if (mpi->ip[d] == 0) {
      /* fill left boundary along this dimension */
      for (i = 0; i < ghosts; i++) dxinv[i] = dxinv[ghosts];
    }
    if (mpi->ip[d] == mpi->iproc[d]-1) {
      /* fill right boundary along this dimension */
      for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) dxinv[i] = dxinv[dim[d]+ghosts-1];
    }
    offset  += (dim[d] + 2*ghosts);
  }

  /* calculate volume integral of the initial solution */
  IERR VolumeIntegral(solver->VolumeIntegralInitial,solver->u,solver,mpi); CHECKERR(ierr);
  if (!mpi->rank) {
    printf("Volume integral of the initial solution:\n");
    for (d=0; d<solver->nvars; d++) printf("%2d:  %1.16E\n",d,solver->VolumeIntegralInitial[d]);
  }
  /* Set initial total boundary flux integral to zero */
  _ArraySetValue_(solver->TotalBoundaryIntegral,solver->nvars,0);

  return(0);
}

#endif
