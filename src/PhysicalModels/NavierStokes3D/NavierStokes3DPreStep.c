/*! @file NavierStokes3DPreStep.c
    @brief Pre-step function for 3D Navier Stokes equations
    @author Debojyoti Ghosh
*/
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*! Pre-step function for the 3D Navier Stokes equations: This function
    is called at the beginning of each time step.
    + The solution at the beginning of the step is copied into #NavierStokes3D::solution
      for the linearized flux partitioning.
    + The linearized "fast" Jacobian (representing the acoustic modes) #NavierStokes3D::fast_jac
      is computed as:
      \f{equation}{
        A_f\left({\bf u}\right) = R\left({\bf u}\right)\Lambda_f\left({\bf u}\right)L\left({\bf u}\right)
      \f}
      where \f$R\f$ and \f$L\f$ are the matrices of right and left eigenvectors, and,
      \f{equation}{
        \Lambda_f = diag\left[0,0,0,u+c,u-c \right]
      \f}
      with \f$c\f$ being the speed of sound.
    \n\n

    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      To appear.
*/
int NavierStokes3DPreStep(
                            double  *u,   /*!< Solution vector */
                            void    *s,   /*!< Solver object of type #HyPar */
                            void    *m,   /*!< MPI object of type #MPIVariables */
                            double  waqt  /*!< Current simulation time */
                         )
{
  HyPar             *solver = (HyPar*)   s;
  NavierStokes3D    *param  = (NavierStokes3D*) solver->physics;
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
    int q = _MODEL_NVARS_*p;

    dir = _XDIR_;
    A = (param->fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+q),D,param,dir); 
    _NavierStokes3DLeftEigenvectors_ ((u+q),L,param,dir);
    _NavierStokes3DRightEigenvectors_((u+q),R,param,dir);
    /* remove the entropy modes (corresponding to eigenvalues u) */
    D[0] = D[12] = D[18] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    dir = _YDIR_;
    A = (param->fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+q),D,param,dir)
    _NavierStokes3DLeftEigenvectors_ ((u+q),L,param,dir);
    _NavierStokes3DRightEigenvectors_((u+q),R,param,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[0] = D[6] = D[18] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    dir = _ZDIR_;
    A = (param->fast_jac + _MODEL_NDIMS_*JacSize*p + dir*JacSize);
    /* get the eigenvalues, and left & right eigenvectors */
    _NavierStokes3DEigenvalues_      ((u+q),D,param,dir)
    _NavierStokes3DLeftEigenvectors_ ((u+q),L,param,dir);
    _NavierStokes3DRightEigenvectors_((u+q),R,param,dir);
    /* remove the entropy modes (corresponding to eigenvalues v) */
    D[0] = D[6] = D[12] = 0.0;
    /* assemble the Jacobian */
    MatMult5(_MODEL_NVARS_,DL,D,L );
    MatMult5(_MODEL_NVARS_,A ,R,DL);

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
