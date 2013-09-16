#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/* Function declarations */
int WriteText               (int,int,int*,double*,double*,char*,int*);
int ApplyBoundaryConditions (void*,void*,double*);

int InitializeSolvers(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;

  if (!mpi->rank) printf("Initializing solvers.\n");

  solver->ApplyBoundaryConditions = ApplyBoundaryConditions;

  /* Time integration */
  if (solver->time_scheme == _FORWARD_EULER_) solver->TimeIntegrate = TimeForwardEuler;
  else {
    fprintf(stderr,"Error: %d is a not a supported time-integration scheme.\n",
            solver->time_scheme);
    return(1);
  }

  /* Solution output function */
  if (!strcmp(solver->op_file_format,"text")) {
    solver->WriteOutput = WriteText;
    if (!strcmp(solver->op_overwrite,"no")) strcpy(solver->op_filename,"op_00000.dat");
    else                                    strcpy(solver->op_filename,"op.dat");
  } else if (!strcmp(solver->op_file_format,"none")) {
    solver->WriteOutput = NULL;
  } else {
    fprintf(stderr,"Error: %s is not a supported file format.\n",solver->op_file_format);
    return(1);
  }

  return(0);
}
