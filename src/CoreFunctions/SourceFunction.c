/*! @file SourceFunction.c
    @author Debojyoti Ghosh
    @brief Evaluate the source term
*/

#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>
#include <hypar.h>

/*! Evaluate the source term \f${\bf S}\left({\bf u}\right)\f$ in the governing equation,
    if the physical model specifies one. In addition, if the simulation requires a sponge
    boundary treatment, the sponge BC function is called.
*/
int SourceFunction(
                    double  *source,  /*!< the computed source term */
                    double  *u,       /*!< solution */
                    void    *s,       /*!< solver object of type #HyPar */
                    void    *m,       /*!< MPI object of type #MPIVariables */
                    double  t         /*!< Current simulation time */
                  )
{
  HyPar           *solver   = (HyPar*)        s;
  MPIVariables    *mpi      = (MPIVariables*) m;
  int             n,d;
  _DECLARE_IERR_;

  /* extract boundary information to check for and implement sponge BC */
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  int             nb        = solver->nBoundaryZones;

  int     nvars   = solver->nvars;
  int     ghosts  = solver->ghosts;
  int     ndims   = solver->ndims;
  int     *dim    = solver->dim_local;

  /* initialize to zero */
  int size = 1;
  for (d=0; d<ndims; d++) size *= (dim[d] + 2*ghosts);
  _ArraySetValue_(source,size*nvars,0.0);

  /* call the source function of the physics model, if available */
  if (solver->SFunction) {
    IERR solver->SFunction(source,u,solver,mpi,t); CHECKERR(ierr);
    solver->count_sou++;
  }

  /* Apart from other source terms, implement sponge BC as a source */
  for (n = 0; n < nb; n++) {
    if (!strcmp(boundary[n].bctype,_SPONGE_)) {
      IERR BCSpongeSource(&boundary[n],ndims,nvars,ghosts,dim,solver->x,u,source); 
      CHECKERR(ierr);
    }
  }


  return(0);
}
