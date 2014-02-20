#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int InitialSolutionSerial    (void*, void*);
static int InitialSolutionParallel  (void*, void*);

int InitialSolution(void *s, void *m)
{
  HyPar  *solver = (HyPar*) s;
  if      (!strcmp(solver->input_mode,"serial"))    return(InitialSolutionSerial    (s,m));
  else if (!strcmp(solver->input_mode,"parallel"))  return(InitialSolutionParallel  (s,m));
  else {
    fprintf(stderr,"Error: Illegal value (%s) for input_mode (may be \"serial\" or \"parallel\"\n",
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
  return(0);
}

int InitialSolutionParallel(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i,d,offset;
  int           ghosts = solver->ghosts;
  _DECLARE_IERR_;

  if (!strcmp(solver->ip_file_type,"ascii")) {

    fprintf(stderr,"Error in InitialSolutionParallel(): Reading ASCII initial solution file in parallel not yet supported. Choose the binary option.\n");
    return(1);

  } else if ((!strcmp(solver->ip_file_type,"binary")) || (!strcmp(solver->ip_file_type,"bin"))) {

    int myrank = mpi->rank;
    int ndims  = solver->ndims;
    int rank[ndims+1], size[ndims], nvars;
    int ferr, total, n;
    FILE *in;

    double *x, *u;
    total = 0; for (n = 0; n < ndims; n++) total += solver->dim_local[n];
    x = (double*) calloc (total, sizeof(double));
    total = solver->nvars; for (n = 0; n < ndims; n++) total *= solver->dim_local[n];
    u = (double*) calloc (total, sizeof(double));

    in = fopen("initial_par.inp","rb");
    if (!in) {
      fprintf(stderr,"Error in InitialSolutionParallel(): Could not open file initial_par.inp for reading.\n");
      return(1);
    }
    if (!myrank) printf("Reading grid and initial conditions from binary file \"initial_par.inp\" (Parallel mode).\n");
    int done = 0;
    while((!feof(in) && (!done))) {
      ferr = fread(rank,sizeof(int),ndims+1,in); 
      if (ferr != (ndims+1)) {
        fprintf(stderr,"Error (1) in reading binary file initial_par.inp in parallel mode on rank %d.\n", myrank);
        return(1);
      }
      ferr = fread(size,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error (2) in reading binary file initial_par.inp in parallel mode on rank %d.\n", myrank);
        return(1);
      }
      ferr = fread(&nvars,sizeof(int),1,in);
      if (ferr != 1) {
        fprintf(stderr,"Error (3) in reading binary file initial_par.inp in parallel mode on rank %d.\n", myrank);
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
          fprintf(stderr,"Error in InitialSolutionParallel(): Inconsistent data read on rank %d.\n",myrank);
          return(1);
        }
        /* read grid */
        total = 0; for (n = 0; n < ndims; n++) total += size[n];
        ferr = fread(x,sizeof(double),total,in);
        if (ferr != total) {
          fprintf(stderr,"Error (4) in reading binary file initial_par.inp in parallel mode on rank %d.\n", myrank);
          return(1);
        }
        /* read initial solution */
        total = nvars; for (n = 0; n < ndims; n++) total *= size[n];
        ferr = fread(u,sizeof(double),total,in);
        if (ferr != total) {
          fprintf(stderr,"Error (5) in reading binary file initial_par.inp in parallel mode on rank %d.\n", myrank);
          return(1);
        }
        done = 1;
      } else {
        int n, offset1 = 0, offset2 = nvars;
        for (n = 0; n < ndims; n++) {
          offset1 += size[n];
          offset2 *= size[n];
        }
        ferr = fseek(in,sizeof(double)*(offset1+offset2),SEEK_CUR);
        if (ferr) {
          fprintf(stderr,"Error (6) in reading binary file initial_par.inp in parallel mode on rank %d.\n", myrank);
          return(1);
        }
      }
    }
    fclose(in);

    if (!done) {
      fprintf(stderr,"Error in InitialSolutionParallel(): data for rank %d not found!\n",myrank);
      return(1);
    }

    int offset1 = 0, offset2 = 0;
    for (n = 0; n < ndims; n++) {
      int i;
      for (i = 0; i < solver->dim_local[n]; i++) solver->x[i+offset1+ghosts] = x[i+offset2];
      offset1 += (solver->dim_local[n]+2*ghosts);
      offset2 += solver->dim_local[n];
    }

    int index[ndims];
    IERR ArrayCopynD(ndims,u,solver->u,solver->dim_local,0,ghosts,index,solver->nvars); CHECKERR(ierr);

    free(x);
    free(u);

  } else {

    fprintf(stderr,"Error in InitialSolutionParallel(): Illegal value (%s) for ip_file type. May be \"ascii\", \"binary\" or \"bin\".\n",
            solver->ip_file_type);
    return(1);

  }

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
  
  return(0);
}
