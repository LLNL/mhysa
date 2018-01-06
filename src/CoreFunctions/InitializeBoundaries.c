/*! @file InitializeBoundaries.c
    @author Debojyoti Ghosh
    @brief Initialize the boundary implementation
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <boundaryconditions.h>
#include <hypar.h>

static int CalculateLocalExtent(void*,void*);

/*! This function initializes the variables and functions related to implementing
    the boundary conditions.
    + Rank 0 reads in the boundary conditions file and broadcasts the information
      to all processors.
    + Depending on the type of boundary, additional information is read in. For
      example, for Dirichlet boundary, the Dirichlet value is read in.
    + Allocate and initialize arrays and variables related to implementing the 
      boundary conditions.
    + Each rank finds out if the subdomain it owns abuts any of the boundaries
      specified.

    Note that boundary conditions are implemented as boundary objects of the
    type #DomainBoundary.
*/
int InitializeBoundaries(
                          void *s,  /*!< Solver object of type #HyPar */
                          void *m   /*!< MPI object of type #MPIVariables */
                        )
{
  HyPar           *solver   = (HyPar*)        s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  DomainBoundary  *boundary = NULL;
  int             n, ferr;
  _DECLARE_IERR_;

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
    ferr = fscanf(in,"%d",&solver->nBoundaryZones); if (ferr != 1) return(1);
    boundary = (DomainBoundary*) calloc (solver->nBoundaryZones,sizeof(DomainBoundary));
    for (n = 0; n < solver->nBoundaryZones; n++) {
      boundary[n].DirichletValue = boundary[n].SpongeValue 
                                 = boundary[n].FlowDensity
                                 = boundary[n].FlowVelocity 
                                 = boundary[n].UnsteadyDirichletData
                                 = NULL;
      boundary[n].UnsteadyDirichletSize = NULL;
    }

    /* read each boundary condition */
    for (n = 0; n < solver->nBoundaryZones; n++) {
      int d, v;
      boundary[n].xmin = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
      boundary[n].xmax = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */

      ferr = fscanf(in,"%s",boundary[n].bctype); if (ferr != 1) return(1);
      ferr = fscanf(in,"%d",&boundary[n].dim  ); if (ferr != 1) return(1);
      ferr = fscanf(in,"%d",&boundary[n].face ); if (ferr != 1) return(1);
      for (d=0; d < solver->ndims; d++) {
        ferr = fscanf(in,"%lf %lf", &boundary[n].xmin[d], &boundary[n].xmax[d]);
        if (ferr != 2) return(1);
      }

      /* read in boundary type-specific additional data if required */

      if (!strcmp(boundary[n].bctype,_DIRICHLET_)) {
        boundary[n].DirichletValue = (double*) calloc (solver->nvars,sizeof(double)); 
                                     /* deallocated in BCCleanup.c */
        /* read the Dirichlet value for each variable on this boundary */
        for (v = 0; v < solver->nvars; v++) ferr = fscanf(in,"%lf",&boundary[n].DirichletValue[v]);
      }

      if (!strcmp(boundary[n].bctype,_SPONGE_)) {
        boundary[n].SpongeValue = (double*) calloc (solver->nvars,sizeof(double)); 
                                     /* deallocated in BCCleanup.c */
        /* read the sponge value for each variable on this boundary */
        for (v = 0; v < solver->nvars; v++) ferr = fscanf(in,"%lf",&boundary[n].SpongeValue[v]);
      }

      if (    (!strcmp(boundary[n].bctype,_SLIP_WALL_)) 
          ||  (!strcmp(boundary[n].bctype,_NOSLIP_WALL_)) ) {
        boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        /* read the wall velocity */
        for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowVelocity[v]);
      }

      if (!strcmp(boundary[n].bctype,_SUBSONIC_INFLOW_)) {
        boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        boundary[n].FlowDensity  = (double*) calloc (solver->nspecies,sizeof(double));
        /* read in the inflow density and velocity */
        for (v = 0; v < solver->nspecies; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowDensity[v]);
        for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowVelocity[v]);
      }

      if (!strcmp(boundary[n].bctype,_SUBSONIC_OUTFLOW_)) {
        /* read in the outflow pressure */
        ferr = fscanf(in,"%lf",&boundary[n].FlowPressure);
      }

      if (!strcmp(boundary[n].bctype,_SUBSONIC_AMBIVALENT_)) {
        boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        boundary[n].FlowDensity  = (double*) calloc (solver->nspecies,sizeof(double));
        /* read in the inflow density, velocity, and pressure */
        for (v = 0; v < solver->nspecies; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowDensity[v]);
        for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowVelocity[v]);
        ferr = fscanf(in,"%lf",&boundary[n].FlowPressure);
      }

      if (!strcmp(boundary[n].bctype,_SUPERSONIC_INFLOW_)) {
        boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        boundary[n].FlowDensity  = (double*) calloc (solver->nspecies,sizeof(double));
        /* read in the inflow density, velocity and pressure */
        for (v = 0; v < solver->nspecies; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowDensity[v]);
        for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowVelocity[v]);
        ferr = fscanf(in,"%lf",&boundary[n].FlowPressure);
      }

      if (!strcmp(boundary[n].bctype,_TURBULENT_SUPERSONIC_INFLOW_)) {
        boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
        boundary[n].FlowDensity  = (double*) calloc (solver->nspecies,sizeof(double));
        /* read in the inflow density, velocity and pressure */
        for (v = 0; v < solver->nspecies; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowDensity[v]);
        for (v = 0; v < solver->ndims; v++) ferr = fscanf(in,"%lf",&boundary[n].FlowVelocity[v]);
        ferr = fscanf(in,"%lf",&boundary[n].FlowPressure);
        ferr = fscanf(in,"%s" , boundary[n].UnsteadyDirichletFilename);
      }

      /* if boundary is periodic, let the MPI and HyPar know */
      if (!strcmp(boundary[n].bctype,_PERIODIC_)) {
        solver->isPeriodic[boundary[n].dim] = 1;
      }
      /* 
        The MPI function to exchange internal (MPI) boundary information will handle
        periodic boundaries ONLY IF number of process along that dimension is more 
        than 1.
      */
      if ((!strcmp(boundary[n].bctype,_PERIODIC_)) && (mpi->iproc[boundary[n].dim] > 1)) {
        mpi->bcperiodic[boundary[n].dim] = 1;
      }

      /* some checks */
      if (boundary[n].dim >= solver->ndims) {
        fprintf(stderr,"Error in reading boundary condition %d: dim %d is invalid (ndims = %d).\n",
                n,boundary[n].dim,solver->ndims);
        return(1);
      }
      printf("  Boundary %30s:  Along dimension %2d and face %+1d\n",
                boundary[n].bctype,boundary[n].dim,boundary[n].face);
    }

    fclose(in);
    printf("%d boundary condition(s) read.\n",solver->nBoundaryZones);
  }

  /* tell other processes how many BCs are there and let them allocate */
  IERR MPIBroadcast_integer(&solver->nBoundaryZones,1,0,&mpi->world); CHECKERR(ierr);
  if (mpi->rank) {
    boundary = (DomainBoundary*) calloc (solver->nBoundaryZones,sizeof(DomainBoundary));
    for (n = 0; n < solver->nBoundaryZones; n++) {
      boundary[n].xmin = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
      boundary[n].xmax = (double*) calloc (solver->ndims,sizeof(double)); /* deallocated in BCCleanup.c */
      boundary[n].DirichletValue = boundary[n].SpongeValue 
                                 = boundary[n].FlowVelocity 
                                 = boundary[n].UnsteadyDirichletData
                                 = NULL;
      boundary[n].UnsteadyDirichletSize = NULL;
    }
  }

  /* communicate BC data to other processes */
  for (n = 0; n < solver->nBoundaryZones; n++) {
    IERR MPIBroadcast_character(boundary[n].bctype,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
    IERR MPIBroadcast_integer  (&boundary[n].dim  ,1                ,0,&mpi->world); CHECKERR(ierr);
    IERR MPIBroadcast_integer  (&boundary[n].face ,1                ,0,&mpi->world); CHECKERR(ierr);
    IERR MPIBroadcast_double   (boundary[n].xmin  ,solver->ndims    ,0,&mpi->world); CHECKERR(ierr);
    IERR MPIBroadcast_double   (boundary[n].xmax  ,solver->ndims    ,0,&mpi->world); CHECKERR(ierr);
  }
  IERR MPIBroadcast_integer(solver->isPeriodic,solver->ndims,0,&mpi->world);CHECKERR(ierr);

  /* broadcast periodic boundary info for MPI to all processes */
  IERR MPIBroadcast_integer(mpi->bcperiodic,solver->ndims,0,&mpi->world);CHECKERR(ierr);

  /* On other processes, if necessary, allocate and receive boundary-type-specific data */
  for (n = 0; n < solver->nBoundaryZones; n++) {
    if (!strcmp(boundary[n].bctype,_DIRICHLET_)) {
      if (mpi->rank)  boundary[n].DirichletValue = (double*) calloc (solver->nvars,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].DirichletValue,solver->nvars,0,&mpi->world); CHECKERR(ierr);
    }

    if (!strcmp(boundary[n].bctype,_SPONGE_)) {
      if (mpi->rank)  boundary[n].SpongeValue = (double*) calloc (solver->nvars,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].SpongeValue,solver->nvars,0,&mpi->world); CHECKERR(ierr);
    }

    if (    (!strcmp(boundary[n].bctype,_SLIP_WALL_)) 
        ||  (!strcmp(boundary[n].bctype,_NOSLIP_WALL_)) ) {
      if (mpi->rank) boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
    }

    if (!strcmp(boundary[n].bctype,_SUBSONIC_INFLOW_)) {
      if (mpi->rank) boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
      if (mpi->rank) boundary[n].FlowDensity  = (double*) calloc (solver->nspecies,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].FlowDensity,solver->nspecies,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(boundary[n].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
    }

    if (!strcmp(boundary[n].bctype,_SUBSONIC_OUTFLOW_)) {
      IERR MPIBroadcast_double(&boundary[n].FlowPressure,1,0,&mpi->world); CHECKERR(ierr);
    }

    if (!strcmp(boundary[n].bctype,_SUBSONIC_AMBIVALENT_)) {
      if (mpi->rank) boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
      if (mpi->rank) boundary[n].FlowDensity = (double*) calloc (solver->nspecies,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].FlowDensity,solver->nspecies,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(boundary[n].FlowVelocity,solver->ndims,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(&boundary[n].FlowPressure,1,0,&mpi->world); CHECKERR(ierr);
    }

    if (!strcmp(boundary[n].bctype,_SUPERSONIC_INFLOW_)) {
      if (mpi->rank) boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
      if (mpi->rank) boundary[n].FlowDensity = (double*) calloc (solver->nspecies,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].FlowDensity ,solver->nspecies,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(boundary[n].FlowVelocity ,solver->ndims,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(&boundary[n].FlowPressure,1            ,0,&mpi->world); CHECKERR(ierr);
    }

    if (!strcmp(boundary[n].bctype,_TURBULENT_SUPERSONIC_INFLOW_)) {
      if (mpi->rank) boundary[n].FlowVelocity = (double*) calloc (solver->ndims,sizeof(double));
      if (mpi->rank) boundary[n].FlowDensity = (double*) calloc (solver->nspecies,sizeof(double));
      IERR MPIBroadcast_double(boundary[n].FlowDensity ,solver->nspecies,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(boundary[n].FlowVelocity ,solver->ndims,0,&mpi->world); CHECKERR(ierr);
      IERR MPIBroadcast_double(&boundary[n].FlowPressure,1            ,0,&mpi->world); CHECKERR(ierr);
      /* allocate arrays and read in unsteady boundary data */
      IERR BCReadTurbulentInflowData(&boundary[n],mpi,solver->ndims,solver->nvars,solver->dim_local); CHECKERR(ierr);
    }

  }

  solver->boundary = boundary;

  /* each process calculates its own part of these boundaries */
  IERR CalculateLocalExtent(solver,mpi); CHECKERR(ierr);

  /* initialize function pointers for each boundary condition */
  for (n = 0; n < solver->nBoundaryZones; n++) {
    IERR BCInitialize(&boundary[n]); CHECKERR(ierr);
  }

  /* set pointer to physics object in each boundary zone */
  for (n = 0; n < solver->nBoundaryZones; n++) {
    boundary[n].physics = &(solver->physics);
  }

  return(0);
}

/*! For each of the boundary conditions, compute its extent on the
    local sub-domain of each rank (if at all this subdomain abuts
    that boundary), and accordingly set the bounding grid indices.
*/
int CalculateLocalExtent(
                          void *s, /*!< Solver object of type #HyPar */
                          void *m  /*!< MPI object of type #MPIVariables */
                        )
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

    if (!strcmp(boundary[n].bctype,_SPONGE_)) {
      /* Sponge boundary condition */
      boundary[n].on_this_proc = 1;
      int offset = 0;
      for (d=0; d<solver->ndims; d++) {
        int is, ie;
        FindInterval(boundary[n].xmin[d],boundary[n].xmax[d],
                     &solver->x[offset+solver->ghosts],
                     solver->dim_local[d],&is,&ie);
        boundary[n].is[d] = is;
        boundary[n].ie[d] = ie;
        if ((ie-is) <= 0) boundary[n].on_this_proc = 0;
        offset += solver->dim_local[d] + 2*solver->ghosts;
      }
    } else {
      /* other boundary conditions */
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

  }
  
  return(0);
}
