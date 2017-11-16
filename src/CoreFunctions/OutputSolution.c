/*! @file OutputSolution.c
    @author Debojyoti Ghosh
    @brief Write out the solution to file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <io.h>
#include <timeintegration.h>
#include <hypar.h>

/* Function declarations */
void IncrementFilenameIndex       (char*,int);

/*! Write out the solution to file */
int OutputSolution(void *s, void *m)
{
  HyPar         *solver = (HyPar*)       s;
  MPIVariables  *mpi    = (MPIVariables*)m;
  _DECLARE_IERR_;
  
  /* if WriteOutput() is NULL, then return */
  if (!solver->WriteOutput) return(0);

  /* time integration module may have auxiliary arrays to write out, so get them */
  int     NSolutions = 0;
  double  *uaux = NULL;
  IERR TimeGetAuxSolutions(&NSolutions,NULL,solver,-1); CHECKERR(ierr);
  if (NSolutions > 10) NSolutions = 10;

  int  n;
  char fname_root[3]     = "op";
  char aux_fname_root[4] = "ts0";

  for (n=0; n<NSolutions; n++) {
    IERR TimeGetAuxSolutions(&NSolutions,&uaux,solver,n); CHECKERR(ierr);
    IERR WriteArray(solver->ndims,solver->nvars,solver->dim_global,solver->dim_local,
                    solver->ghosts,solver->x,uaux,solver,mpi,aux_fname_root); CHECKERR(ierr);
    aux_fname_root[2]++;
  }
  IERR WriteArray(solver->ndims,solver->nvars,solver->dim_global,solver->dim_local,
                  solver->ghosts,solver->x,solver->u,solver,mpi,fname_root); CHECKERR(ierr);

  /* increment the index string, if required */
  if ((!strcmp(solver->output_mode,"serial")) && (!strcmp(solver->op_overwrite,"no")))
      IncrementFilenameIndex(solver->filename_index,solver->index_length);
  
  return(0);
}
