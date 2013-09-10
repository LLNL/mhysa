#include <stdio.h>
#include <sys/time.h>
#include <mpivars.h>
#include <hypar.h>

int main(int argc,char **argv)
{
  MPIVariables    mpi;
  HyPar           solver;
  int             ierr = 0;
  struct timeval  main_start, solve_start;
  struct timeval  main_end  , solve_end  ;

#ifdef serial
  mpi.rank  = 0;
  mpi.nproc = 0;
  printf("HyPar - Serial Version\n");
#else
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi.rank );
  MPI_Comm_size(MPI_COMM_WORLD,&mpi.nproc);
  if (!mpi.rank) printf("HyPar - Parallel (MPI) version with %d processes\n",mpi.nproc);
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
//  ierr = InitialSolution(&solver,&mpi);
  if (ierr) {
    printf("Error: InitialSolution() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize domain boundaries */
//  ierr = InitializeBoundaries(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeBoundaries() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initialize solvers */
//  ierr = InitializeSolvers(&solver,&mpi);
  if (ierr) {
    printf("Error: InitializeSolvers() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  /* Initializations complete */

  /* Run the solver */
  gettimeofday(&solve_start,NULL);
//  ierr = Solve(&solver,&mpi);
  gettimeofday(&solve_end,NULL);
  if (ierr) {
    printf("Error: Solve() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }

  /* Write output */
//  ierr = OutputSolution(&solver,&mpi);
  if (ierr) {
    printf("Error: OutputSolution() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }

  /* Cleaning up */
//  ierr = CleanUp(&solver);
  if (ierr) {
    printf("Error: CleanUp() returned with status %d on process %d.\n",ierr,mpi.rank);
    return(ierr);
  }
  if (!mpi.rank) printf("Finished.\n");

  gettimeofday(&main_end,NULL);

  /* calculate solver and total runtimes */
  long long walltime;
  walltime = (  (main_end.tv_sec * 1000000   + main_end.tv_usec  ) 
              - (main_start.tv_sec * 1000000 + main_start.tv_usec));
  double main_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&main_runtime,&main_runtime,1); if(ierr) return(ierr);
  walltime = (  (solve_end.tv_sec * 1000000   + solve_end.tv_usec  ) 
              - (solve_start.tv_sec * 1000000 + solve_start.tv_usec));
  double solver_runtime = (double) walltime / 1000000.0;
  ierr = MPIMax_double(&solver_runtime,&solver_runtime,1); if(ierr) return(ierr);
  if (!mpi.rank) printf("Solver runtime (in seconds): %E\n",solver_runtime);
  if (!mpi.rank) printf("Total  runtime (in seconds): %E\n",main_runtime);

#ifndef serial
  MPI_Finalize();
#endif
  return(0);
}
