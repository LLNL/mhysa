/*! @file TimeGLMGEE.c
    @brief General Linear Methods with Global Error Estimators
    @author Debojyoti Ghosh
*/

#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Advance the ODE given by
  \f{equation}{
    \frac{d{\bf u}}{dt} = {\bf F} \left({\bf u}\right)
  \f}
  by one time step of size #TimeIntegration::dt using the \f$s\f$-stage General Linear Method with
  Global Error Estimation (GLM-GEE), given by
  \f{align}{
    {\bf U}^{\left(i\right)} &= c_{i0}{\bf u}_n + \sum_{j=1}^{r-1} c_{ij} \tilde{\bf u}_n^j 
                               + \Delta t \sum_{j=1}^{i-1} a_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), i=1,\cdots,s\\
    {\bf u}_{n+1} &= d_{00} {\bf u}_n + \sum_{j=1}^{r-1} d_{0j} \tilde{\bf u}_n^j 
                      + \Delta t \sum_{j=1}^s b_{0j} {\bf F}\left({\bf U}^{\left(j\right)}\right), \\
    \tilde{\bf u}_{n+1}^i &= d_{i0} {\bf u}_n + \sum_{j=1}^{r-1} d_{ij} \tilde{\bf u}_n^j 
                      + \Delta t \sum_{j=1}^s b_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), i=1,\cdots,r-1
  \f}
  where the superscripts in parentheses represent stages, the subscripts represent the time level, the 
  superscripts without parentheses represent auxiliary solutions, \f$\Delta t\f$ is the
  time step size #TimeIntegration::dt, \f$\tilde{\bf u}^i, i=1,\cdots,r-1\f$ are the auxiliary solutions,
  and \f${\bf F}\left({\bf u}\right)\f$ is computed by #TimeIntegration::RHSFunction. The coefficients 
  defining this methods are:
  + \f$a_{ij}, i=1,\cdots,s, j=1,\cdots,s\f$ (#GLMGEEParameters::A)
  + \f$b_{ij}, i=0,\cdots,r-1, j=1,\cdots,s\f$ (#GLMGEEParameters::B)
  + \f$c_{ij}, i=1,\cdots,s, j=0,\cdots,r-1\f$ (#GLMGEEParameters::C)
  + \f$d_{ij}, i=0,\cdots,r-1, j=0,\cdots,r-1\f$ (#GLMGEEParameters::D)
  
  where \f$s\f$ is the number of stages (#GLMGEEParameters::nstages) and \f$r\f$ is the number of auxiliary solutions
  propagated with the solution (#GLMGEEParameters::r).

  Note: In the code #TimeIntegration::Udot is equivalent to \f${\bf F}\left({\bf u}\right)\f$.

  References:
  + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
int TimeGLMGEE(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration   *TS     = (TimeIntegration*) ts;
  HyPar             *solver = (HyPar*)           TS->solver;
  MPIVariables      *mpi    = (MPIVariables*)    TS->mpi;
  GLMGEEParameters  *params = (GLMGEEParameters*)solver->msti;
  int               i, j;
  _DECLARE_IERR_;

  int    s      = params->nstages;
  int    r      = params->r;
  int    ndims  = solver->ndims;
  int    nvars  = solver->nvars;
  int    size   = nvars * solver->npoints_local_wghosts;
  double dt     = TS->dt;
  double *A=params->A, *B=params->B, *C=params->C, *D=params->D, *c=params->c,
         **U = TS->U, **Udot = TS->Udot, **Uaux = &TS->U[r];

  /* Calculate stage values */
  for (j=0; j<s; j++) {

    double stagetime = TS->waqt + c[j]*dt;

    _ArrayScaleCopy1D_(solver->u,C[j*r+0],U[0],size);
    for (i=1;i<r;i++) _ArrayAXPY_(Uaux[i-1],C[j*r+i]   ,U[0],size);
    for (i=0;i<j;i++) _ArrayAXPY_(Udot[i]  ,dt*A[j*s+i],U[0],size); 

    if (solver->PreStage) { IERR solver->PreStage(j,U,solver,mpi,stagetime); CHECKERR(ierr); }
    IERR TS->RHSFunction(Udot[j],U[0],solver,mpi,stagetime);
    if (solver->PostStage) { IERR solver->PostStage(U[j],solver,mpi,stagetime); CHECKERR(ierr); }

    _ArraySetValue_(TS->BoundaryFlux[j],2*ndims*nvars,0.0);
    _ArrayCopy1D_(solver->StageBoundaryIntegral,TS->BoundaryFlux[j],2*ndims*nvars);
  }

  /* Step completion */
  for (j=0; j<r; j++) {
    _ArrayScaleCopy1D_(solver->u,D[j*r+0],U[j],size);
    for (i=1; i<r; i++) _ArrayAXPY_(Uaux[i-1],D[j*r+i]   ,U[j],size);
    for (i=0; i<s; i++) _ArrayAXPY_(Udot[i]  ,dt*B[j*s+i],U[j],size);
  }

  for (i=0; i<s; i++) _ArrayAXPY_(TS->BoundaryFlux[i],dt*B[0*s+i],solver->StepBoundaryIntegral,
                                  2*ndims*nvars);

  _ArrayCopy1D_(U[0],solver->u,size);
  for (i=1; i<r; i++) _ArrayCopy1D_(U[i],Uaux[i-1],size);

  return(0);
}

