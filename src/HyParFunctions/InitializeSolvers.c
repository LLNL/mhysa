#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpivars.h>
#include <hypar.h>

/* Output functions */
int WriteText  (int,int,int*,double*,double*,char*,int*);

int InitializeSolvers(void *s, void *m)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
//  int           ierr    = 0;

  if (!mpi->rank) printf("Initializing solvers.\n");

  /* Solution output function */
  if (!strcmp(solver->op_file_format,"text")) {
    solver->WriteOutput = WriteText;
    if (!strcmp(solver->op_overwrite,"no")) strcpy(solver->op_filename,"op_00000.dat");
    else                                    strcpy(solver->op_filename,"op.dat");
  } else if (!strcmp(solver->op_file_format,"none")) {
    solver->WriteOutput = NULL;
  } else {
    fprintf(stderr,"Error: %s is not an unsupported file format.\n",solver->op_file_format);
    return(1);
  }

  return(0);
}
