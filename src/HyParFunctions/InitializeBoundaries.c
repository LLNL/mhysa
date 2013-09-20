#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>

static int CalculateLocalExtent(void*,void*);

int InitializeBoundaries(void *s,void *m)
{
  HyPar           *solver   = (HyPar*)        s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  DomainBoundary  *boundary = NULL;
  int             n,ierr    = 0;

  /* root process reads boundary condition file */
  if (!mpi->rank) {
    printf("Reading boundary conditions from \"boundary.inp\".\n");
    FILE *in;
    in = fopen("boundary.inp","r");
    if (!in) {
      fprintf(stderr,"Error: boundary condition file \"boundary.inp\" not found.\n");
      return(1);
    }

    /* read number of boundary conditions and allocate */
    ierr = fscanf(in,"%d",&solver->nBoundaryZones); if (ierr != 1) return(1);
    boundary = (DomainBoundary*) calloc (solver->nBoundaryZones,sizeof(DomainBoundary));

    /* read each boundary condition */
    for (n = 0; n < solver->nBoundaryZones; n++) {
      int d;
      boundary[n].xmin = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
      boundary[n].xmax = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */

      ierr = fscanf(in,"%s",boundary[n].bctype); if (ierr != 1) return(1);
      ierr = fscanf(in,"%d",&boundary[n].var  ); if (ierr != 1) return(1);
      ierr = fscanf(in,"%d",&boundary[n].dim  ); if (ierr != 1) return(1);
      ierr = fscanf(in,"%d",&boundary[n].face ); if (ierr != 1) return(1);
      for (d=0; d < solver->ndims; d++) {
        ierr = fscanf(in,"%lf %lf", &boundary[n].xmin[d], &boundary[n].xmax[d]);
        if (ierr != 2) return(1);
      }
      /* some checks */
      if (boundary[n].dim >= solver->ndims) {
        fprintf(stderr,"Error in reading boundary condition %d: dim %d is invalid (ndims = %d).\n",
                n,boundary[n].dim,solver->ndims);
        return(1);
      }
      if (boundary[n].var >= solver->nvars) {
        fprintf(stderr,"Error in reading boundary condition %d: var %d is invalid (nvars = %d).\n",
                n,boundary[n].var,solver->nvars);
        return(1);
      }
      printf("  Boundary %10s:  Variable %2d, along dimension %2d and face %+1d\n",
                boundary[n].bctype,boundary[n].var,boundary[n].dim,boundary[n].face);
    }
    fclose(in);
    printf("%d boundary condition(s) read.\n",solver->nBoundaryZones);
  }

  /* tell other processes how many BCs are there and let them allocate */
  ierr = MPIBroadcast_integer(&solver->nBoundaryZones,1,0); CHECKERR(ierr);
  if (mpi->rank) {
    boundary = (DomainBoundary*) calloc (solver->nBoundaryZones,sizeof(DomainBoundary));
    for (n = 0; n < solver->nBoundaryZones; n++) {
      boundary[n].xmin = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
      boundary[n].xmax = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
    }
  }
  /* communicate BC data to other processes */
  for (n = 0; n < solver->nBoundaryZones; n++) {
    MPIBroadcast_character(boundary[n].bctype,_MAX_STRING_SIZE_,0);
    MPIBroadcast_integer  (&boundary[n].var  ,1                ,0);
    MPIBroadcast_integer  (&boundary[n].dim  ,1                ,0);
    MPIBroadcast_integer  (&boundary[n].face ,1                ,0);
    MPIBroadcast_double   (boundary[n].xmin  ,solver->ndims    ,0);
    MPIBroadcast_double   (boundary[n].xmax  ,solver->ndims    ,0);
  }

  solver->boundary = boundary;

  /* each process calculates its own part of these boundaries */
  ierr = CalculateLocalExtent(solver,mpi); CHECKERR(ierr);

  /* initialize function pointers for each boundary condition */
  for (n = 0; n < solver->nBoundaryZones; n++) {
    ierr = BCInitialize(&boundary[n]); CHECKERR(ierr);
  }

  return(0);
}


int CalculateLocalExtent(void *s,void *m)
{
  HyPar           *solver   = (HyPar*)        s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  DomainBoundary  *boundary = solver->boundary;

  int n;
  for (n = 0; n < solver->nBoundaryZones; n++) {
    /* allocate */
    boundary[n].is   = (int*) calloc (solver->ndims,sizeof(int)); /* deallocated in BCCleanup.c */
    boundary[n].ie   = (int*) calloc (solver->ndims,sizeof(int)); /* deallocated in BCCleanup.c */

    int d,dim = boundary[n].dim;

    if (boundary[n].face == 1) {

      if (mpi->ip[dim] == 0) {
        boundary[n].on_this_proc = 1;
        int offset = 0;
        for (d=0; d<solver->ndims; d++) {
          if (d == dim) {
            boundary[n].is[d] = -solver->ghosts;
            boundary[n].ie[d] = 0;
          } else {
            int is, ie;
            FindInterval(boundary[n].xmin[d],boundary[n].xmax[d],
                         &solver->x[offset+solver->ghosts],
                         solver->dim_local[d],&is,&ie);
            boundary[n].is[d] = is;
            boundary[n].ie[d] = ie;
            if ((ie-is) <= 0) boundary[n].on_this_proc = 0;
          }
          offset += solver->dim_local[d] + 2*solver->ghosts;
        }
      } else  boundary[n].on_this_proc = 0;

    } else if (boundary[n].face == -1) {

      if (mpi->ip[dim] == mpi->iproc[dim]-1) {
        boundary[n].on_this_proc = 1;
        int offset = 0;
        for (d=0; d<solver->ndims; d++) {
          if (d == dim) {
            boundary[n].is[d] = solver->dim_local[dim];
            boundary[n].ie[d] = solver->dim_local[dim] + solver->ghosts;
          } else {
            int is, ie;
            FindInterval(boundary[n].xmin[d],boundary[n].xmax[d],
                         &solver->x[offset+solver->ghosts],
                         solver->dim_local[d],&is,&ie);
            boundary[n].is[d] = is;
            boundary[n].ie[d] = ie;
            if ((ie-is) <= 0) boundary[n].on_this_proc = 0;
          }
          offset += solver->dim_local[d] + 2*solver->ghosts;
        }
      } else  boundary[n].on_this_proc = 0;
    }

  }
  
  return(0);
}
