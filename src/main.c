/*! @file main.c
 *  @brief Main driver.
 * The main driver function that calls the initialization, solving, and cleaning up functions.
 *  @author Debojyoti Ghosh
*/

/*! @mainpage

  @author Debojyoti Ghosh [\b Email: (first name) (dot) (last name) (at) gmail (dot) com, \b Website: http://debog.github.io/]

  MHYSA: Multispecies Hypersonic Flow Simulator
  ---------------------------------------------

  MHYSA is a conservative finite-difference code to solve the 1D, 2D, and 3D multispecies
  Navier-Stokes equations for hypersonic flow applications. It supports Cartesian grids 
  only and has the option of including immersed boundaries.

  MHYSA is based on HyPar (http://hypar.github.io/), a finite-difference algorithm to solve 
  hyperbolic-parabolic partial differential equations (with source terms). See the HyPar
  website and documentation for more details.

  Features
  --------
  + Solves <B>1D Euler</B> and <B>2D and 3D Navier-Stokes equations</B>.
  + Solves the PDEs over <B>Cartesian</B> grids.
  + Written entirely in C and uses the MPICH library. It also uses OpenMP threads 
    but this is a work-in-progress.
  + Can be <B>compiled with PETSc</B> (http://www.mcs.anl.gov/petsc/), if available, where 
    it can use PETSc's time integration module TS (http://www.mcs.anl.gov/petsc/petsc-current/src/ts/).
  + For 3-dimensional simulations, the <B>immersed boundaries</B> can be used to
    solve over non-Cartesian geometries.

  MHYSA is based on HyPar, which has been developed to be scalable, and apart from the usual 
  functionalities to solve a system of PDEs on distributed architectures, it provides scalable 
  file I/O functions. It has been tested on several platforms, including DOE Leadership-class
  supercomputers, with up to ~0.5 million MPI ranks.

  Download
  --------

  Please contact the author to obtain a copy.

  Documentation
  -------------
  To generate a local copy of this documentation, run "doxygen Doxyfile" in $(root_dir). The folder $(root_dir)/doc
  should contain the generated documentation in HTML format.

  Compiling
  ---------

  To compile MHYSA, follow these steps in the root directory:
  
        autoreconf -i
        [CFLAGS="..."] ./configure [options]
        make
        make install

  CFLAGS should include all the compiler flags.

  \b Note: Default installation target is its own directory, and thus "make install" should not require
           administrative privileges. The binary will be placed in \a bin/ subdirectory.

  The configure options can include options such as BLAS/LAPACK location, MPI directory, etc. Type "./configure --help"
  to see a full list. The options specific to MHYSA are:
  + --with-mpi-dir: Specify path where mpicc is installed.
  + --enable-omp: Enable OpenMP threads.
  + --enable-scalapack: Enable ScaLAPACK (this will make available a tridiagonal solver using ScaLAPACK).
  + --with-blas-dir: Specify path where BLAS is installed (relevant only if --enable-scalapack is specified).
  + --with-lapack-dir: Specify path where LAPACK is installed (relevant only if --enable-scalapack is specified).
  + --with-scalapack-dir: Specify path where ScaLAPACK is installed (relevant only if --enable-scalapack is specified).
  + --with-fortran-lib: Specify path where FORTRAN libraries are installed (for ScaLAPACK) (relevant only if --enable-scalapack 
    is specified).

  \b Compiling \b with \b PETSc:
  Install PETSc and make sure the environment variables \b PETSC_DIR and \b PETSC_ARCH are defined. Please PETSc's
  installation instructions for this. Once these environment variables are present, MHYSA will use them to compile
  itself with PETSc functionalities.

  Notes
  -----
  + This package has been tested using the GNU and IBM C compilers. The configuration script is designed to look for these 
    compilers only.
  + Feel free to contact me about anything regarding this (doubts/difficulties/suggestions).
  + Feel free to use and modify the code in any way.

  Running
  -------
  + It's best to start with some examples. See the section on examples.
  + To run more cases, see the section in input files for a complete description of input files required.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#ifdef with_petsc
#include <petscinterface.h>
#endif
#include <mpivars.h>
#include <hypar.h>

static const char help[] = "MHYSA - Multispecies Hypersonic Flow Simulator";

/*!
 * \brief Main driver
 *
 * The main driver function that calls the initialization, solving, and cleaning up functions.
*/
int main(int argc,char **argv)
{
  MPIVariables    mpi;
  HyPar           solver;
  int             ierr = 0, d;
  struct timeval  main_start, solve_start;
  struct timeval  main_end  , solve_end  ;
#ifdef with_petsc
  PetscBool       use_petscts;
#endif

#ifdef serial
  mpi.rank  = 0;
  mpi.nproc = 1;
  mpi.world = 0;
  mpi.comm  = NULL;
  printf("MHYSA - Serial Version\n");
#else
  MPI_Init(&argc,&argv);
  MPI_Comm_dup (MPI_COMM_WORLD,&mpi.world);
  MPI_Comm_rank(mpi.world,&mpi.rank );
  MPI_Comm_size(mpi.world,&mpi.nproc);
  if (!mpi.rank) printf("MHYSA - Parallel (MPI) version with %d processes\n",mpi.nproc);
#endif

#ifdef with_petsc
  PetscInitialize(&argc,&argv,(char*)0,help);
  if (!mpi.rank) printf("Compiled with PETSc time integration.\n");

  use_petscts = PETSC_FALSE; /* default value */
  ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-use-petscts" ,&use_petscts ,PETSC_NULL); CHKERRQ(ierr);
  solver.use_petscTS  = use_petscts;
#endif

  gettimeofday(&main_start,NULL);

  /* Read Inputs */
  ierr = ReadInputs(&solver,&mpi);
  if (ierr) {
    printf("Error: ReadInputs() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize and allocate arrays */
  ierr = Initialize(&solver,&mpi);
  if (ierr) {
    printf("Error: Initialize() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* read and set grid & initial solution */
  ierr = InitialSolution(&solver,&mpi);
  if (ierr) {
    printf("Error: InitialSolution() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize solvers */
  ierr = InitializeSolvers(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeSolvers() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize physics */
  ierr = InitializePhysics(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializePhysics() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize domain boundaries */
  ierr = InitializeBoundaries(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeBoundaries() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize immersed boundaries */
  ierr = InitializeImmersedBoundaries(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeImmersedBoundaries() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initializations complete */
  
  /* Run the solver */
#ifndef serial
  MPI_Barrier(mpi.world);
#endif
  gettimeofday(&solve_start,NULL);
#ifdef with_petsc
  if (solver.use_petscTS == PETSC_TRUE) {
    /* Use PETSc time-integration */
    ierr = SolvePETSc(&solver,&mpi);
    if (ierr) {
      printf("Error: SolvePETSc() returned with status %d on process %d.\n",ierr,mpi.rank);
      return(ierr);
    }
  } else {
    /* Use native time-integration */
    ierr = Solve(&solver,&mpi);
    if (ierr) {
      printf("Error: Solve() returned with status %d on process %d.\n",ierr,mpi.rank);
      return(ierr);
    }
  }
#else 
  /* Use native time-integration */
  ierr = Solve(&solver,&mpi);
  if (ierr) {
    printf("Error: Solve() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
#endif
  gettimeofday(&solve_end,NULL);
#ifndef serial
  MPI_Barrier(mpi.world);
#endif
  gettimeofday(&main_end,NULL);

  /* calculate solver and total runtimes */
  long long walltime;
  walltime = (  (main_end.tv_sec * 1000000   + main_end.tv_usec  ) 
              - (main_start.tv_sec * 1000000 + main_start.tv_usec));
  double main_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&main_runtime,&main_runtime,1,&mpi.world); if(ierr) return(ierr);
  walltime = (  (solve_end.tv_sec * 1000000   + solve_end.tv_usec  ) 
              - (solve_start.tv_sec * 1000000 + solve_start.tv_usec));
  double solver_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&solver_runtime,&solver_runtime,1,&mpi.world); if(ierr) return(ierr);

  if (!mpi.rank) {
    FILE *out; 
    /* write out solution errors and wall times to file */
    out = fopen("errors.dat","w");
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",solver.dim_global[d]);
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",mpi.iproc[d]);
    fprintf(out,"%1.16E  ",solver.dt);
    fprintf(out,"%1.16E %1.16E %1.16E   ",solver.error[0],solver.error[1],solver.error[2]);
    fprintf(out,"%1.16E %1.16E\n",solver_runtime,main_runtime);
    fclose(out);
    /* write out conservation errors to file */
    out = fopen("conservation.dat","w");
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",solver.dim_global[d]);
    for (d=0; d<solver.ndims; d++) fprintf(out,"%4d ",mpi.iproc[d]);
    fprintf(out,"%1.16E  ",solver.dt);
    for (d=0; d<solver.nvars; d++) fprintf(out,"%1.16E ",solver.ConservationError[d]);
    fprintf(out,"\n");
    fclose(out);
    /* write out function call counts to file */
    out = fopen("function_counts.dat","w");
    fprintf(out,"%d\n",solver.n_iter);
    fprintf(out,"%d\n",solver.count_hyp);
    fprintf(out,"%d\n",solver.count_par);
    fprintf(out,"%d\n",solver.count_sou);
#ifdef with_petsc
    fprintf(out,"%d\n",solver.count_RHSFunction);
    fprintf(out,"%d\n",solver.count_IFunction);
    fprintf(out,"%d\n",solver.count_IJacobian);
    fprintf(out,"%d\n",solver.count_IJacFunction);
#endif
    fclose(out);
    /* print solution errors, conservation errors, and wall times to screen */
    printf("Computed errors:\n");
    printf("  L1         Error           : %1.16E\n",solver.error[0]);
    printf("  L2         Error           : %1.16E\n",solver.error[1]);
    printf("  Linfinity  Error           : %1.16E\n",solver.error[2]);
    printf("Conservation Errors:\n");
    for (d=0; d<solver.nvars; d++) printf("\t%1.16E\n",solver.ConservationError[d]);
    printf("Solver runtime (in seconds): %1.16E\n",solver_runtime);
    printf("Total  runtime (in seconds): %1.16E\n",main_runtime);
  }

  /* Cleaning up */
  ierr = Cleanup(&solver,&mpi);
  if (ierr) {
    printf("Error: CleanUp() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  if (!mpi.rank) printf("Finished.\n");

#ifdef with_petsc
  PetscFinalize();
#endif

#ifndef serial
  MPI_Comm_free(&mpi.world);
  MPI_Finalize();
#endif
  return(0);
}
