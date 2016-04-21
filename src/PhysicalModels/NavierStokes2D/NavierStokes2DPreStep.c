/*! @file NavierStokes2DPreStep.c
    @brief Pre-step function for 2D Navier Stokes equations
    @author Debojyoti Ghosh
*/
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*! Pre-step function for the 2D Navier Stokes equations: This function
    is called at the beginning of each time step.
    + The solution at the beginning of the step is copied into #NavierStokes2D::solution
      for the linearized flux partitioning.
    + The linearized "fast" Jacobian (representing the acoustic modes) #NavierStokes2D::fast_jac
      is computed as:
      \f{equation}{
        A_f\left({\bf u}\right) = R\left({\bf u}\right)\Lambda_f\left({\bf u}\right)L\left({\bf u}\right)
      \f}
      where \f$R\f$ and \f$L\f$ are the matrices of right and left eigenvectors, and,
      \f{equation}{
        \Lambda_f = diag\left[0,0,u+c,u-c \right]
      \f}
      with \f$c\f$ being the speed of sound.
    \n\n

    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      To appear.
*/
int NavierStokes2DPreStep(
                            double  *u,   /*!< Solution vector */
                            void    *s,   /*!< Solver object of type #HyPar */
                            void    *m,   /*!< MPI object of type #MPIVariables */
                            double  waqt  /*!< Current simulation time */
                         )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes2D    *param  = (NavierStokes2D*) solver->physics;
  int               *dim    = solver->dim_local;
  int               ghosts  = solver->ghosts, dir, p;
  double            *A;
  static const int  ndims   = _MODEL_NDIMS_;
  static const int  JacSize = _MODEL_NVARS_*_MODEL_NVARS_;
  static int        index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];
  static double     D[_MODEL_NVARS_*_MODEL_NVARS_],L[_MODEL_NVARS_*_MODEL_NVARS_],
                    R[_MODEL_NVARS_*_MODEL_NVARS_],DL[_MODEL_NVARS_*_MODEL_NVARS_];

  /* set bounds for array index to include ghost points */
  _ArrayAddCopy1D_(dim,(2*ghosts),bounds,ndims);
  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);
  /* copy the solution to act as a reference for linearization */
  _ArrayCopy1D_(u,param->solution,(solver->npoints_local_wghosts*_MODEL_NVARS_));

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

    dir = _XDIR_;
    A = (param->fast_jac + 2*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes2DEigenvalues_      ((u+_MODEL_NVARS_*p),D,param,dir); 
    _NavierStokes2DLeftEigenvectors_ ((u+_MODEL_NVARS_*p),L,param,dir);
    _NavierStokes2DRightEigenvectors_((u+_MODEL_NVARS_*p),R,param,dir);
    /* remove the entropy modes (corresponding to eigenvalues u) */
    D[2*_MODEL_NVARS_+2] = D[3*_MODEL_NVARS_+3] = 0.0;
    /* assemble the Jacobian */
    MatMult4(_MODEL_NVARS_,DL,D,L );
    MatMult4(_MODEL_NVARS_,A ,R,DL);

    dir = _YDIR_;
    A = (param->fast_jac + 2*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes2DEigenvalues_      ((u+_MODEL_NVARS_*p),D,param,dir)
    _NavierStokes2DLeftEigenvectors_ ((u+_MODEL_NVARS_*p),L,param,dir);
    _NavierStokes2DRightEigenvectors_((u+_MODEL_NVARS_*p),R,param,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[1*_MODEL_NVARS_+1] = D[3*_MODEL_NVARS_+3] = 0.0;
    /* assemble the Jacobian */
    MatMult4(_MODEL_NVARS_,DL,D,L );
    MatMult4(_MODEL_NVARS_,A ,R,DL);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
