/*! @file NavierStokes3DSource.c
    @author Debojyoti Ghosh
    @brief Compute the gravitational source term for the 3D Navier Stokes system
*/
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

static int NavierStokes3DSourceFunction (double*,double*,double*,void*,void*,double,int);
static int NavierStokes3DSourceUpwind   (double*,double*,double*,double*,int,void*,double);

/*! Computes the gravitational source term using a well-balanced formulation: the source term
    is rewritten as follows:
    \f{equation}{
      \left[\begin{array}{c} 0 \\ -\rho {\bf g}\cdot{\bf \hat{i}} \\ -\rho {\bf g}\cdot{\bf \hat{j}} \\ -\rho {\bf g}\cdot{\bf \hat{k}}  \\ -\rho u {\bf g}\cdot{\bf \hat{i}} - \rho v {\bf g}\cdot{\bf \hat{j}} - \rho w {\bf g}\cdot{\bf \hat{k}} \end{array}\right] 
      = 
      \left[\begin{array}{ccccc} 0 & p_0 \varrho^{-1} & 0 & 0 &  p_0 u \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial x}\left[\begin{array}{c} 0 \\ \varphi \\ 0 \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{ccccc} 0 & 0 & p_0 \varrho^{-1} & 0 &  p_0 v \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ \varphi \\ 0 \\ \varphi \end{array}\right]
      +
      \left[\begin{array}{ccccc} 0 & 0 & 0 &  p_0 \varrho^{-1} & p_0 w \varrho^{-1} \end{array}\right] \cdot \frac{\partial}{\partial y}\left[\begin{array}{c} 0 \\ 0 \\ 0 \\ \varphi \\ \varphi \end{array}\right]
    \f}
    where \f$\varphi = \varphi\left(x,y\right)\f$ and \f$\varrho = \varrho\left(x,y\right)\f$ are computed in 
    NavierStokes3DGravityField(). The derivatives are computed in an exactly identical manner as the hyperbolic
    flux (i.e., using a conservative finite-difference formulation) (see HyperbolicFunction()).
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, http://dx.doi.org/10.2514/1.J054580.
*/
int NavierStokes3DSource(
                          double  *source,  /*!< Array to hold the computed source */
                          double  *u,       /*!< Solution vector array */
                          void    *s,       /*!< Solver object of type #HyPar */
                          void    *m,       /*!< MPI object of type #MPIVariables */
                          double  t         /*!< Current simulation time */
                        )
{
  HyPar           *solver = (HyPar* )         s;
  MPIVariables    *mpi    = (MPIVariables*)   m;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  if ((param->grav_x == 0.0) && (param->grav_y == 0.0) && (param->grav_z == 0.0))  
    return(0); /* no gravitational forces */

  int     v, done, p, p1, p2, dir;
  double  *SourceI = solver->fluxI; /* interace source term       */
  double  *SourceC = solver->fluxC; /* cell-centered source term  */
  double  *SourceL = solver->fL;
  double  *SourceR = solver->fR;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  double  *x      = solver->x;
  double  *dxinv  = solver->dxinv;
  double  RT      =  param->p0 / param->rho0;
  static int index[_MODEL_NDIMS_],index1[_MODEL_NDIMS_],index2[_MODEL_NDIMS_],dim_interface[_MODEL_NDIMS_];
  static double grav[_MODEL_NDIMS_];

  grav[_XDIR_] = param->grav_x; 
  grav[_YDIR_] = param->grav_y; 
  grav[_ZDIR_] = param->grav_z; 
  for (dir = 0; dir < _MODEL_NDIMS_; dir++) {
    if (grav[dir] != 0.0) {
      /* set interface dimensions */
      _ArrayCopy1D_(dim,dim_interface,_MODEL_NDIMS_); dim_interface[dir]++;
      /* calculate the split source function exp(-phi/RT) */
      IERR NavierStokes3DSourceFunction(SourceC,u,x,solver,mpi,t,dir); CHECKERR(ierr);
      /* calculate the left and right interface source terms */
      IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,dir,solver,mpi,0); CHECKERR(ierr);
      IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,dir,solver,mpi,0); CHECKERR(ierr);
      /* calculate the final interface source term */
      IERR NavierStokes3DSourceUpwind(SourceI,SourceL,SourceR,u,dir,solver,t);
      /* calculate the final cell-centered source term */
      done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
      while (!done) {
        _ArrayCopy1D_(index,index1,_MODEL_NDIMS_);
        _ArrayCopy1D_(index,index2,_MODEL_NDIMS_); index2[dir]++;
        _ArrayIndex1D_(_MODEL_NDIMS_,dim          ,index ,ghosts,p );
        _ArrayIndex1D_(_MODEL_NDIMS_,dim_interface,index1,0     ,p1);
        _ArrayIndex1D_(_MODEL_NDIMS_,dim_interface,index2,0     ,p2);

        double dx_inverse; _GetCoordinate_(dir,index[dir],dim,ghosts,dxinv,dx_inverse);
        double rho, vel[_MODEL_NDIMS_], e, P; _NavierStokes3DGetFlowVar_((u+_MODEL_NVARS_*p),rho,vel[0],vel[1],vel[2],e,P,param);
        double term[_MODEL_NVARS_] = {0.0, rho*RT*(dir==_XDIR_), rho*RT*(dir==_YDIR_), rho*RT*(dir==_ZDIR_), rho*RT*vel[dir]};
        for (v=0; v<_MODEL_NVARS_; v++) {
          source[_MODEL_NVARS_*p+v] += (  (term[v]*param->grav_field_f[p]) 
                                        * (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse );
        }
        vel[0] = P; /* useless statement to avoid compiler warnings */
        _ArrayIncrementIndex_(_MODEL_NDIMS_,dim,index,done);
      }
    }
  }

  return(0);
}

/*! Compute the source function in the well-balanced treatment of the source terms. The source
    function is:
    \f{equation}{
      dir = x \rightarrow \left[\begin{array}{c}0 \\ \varphi\left(x,y,z\right) \\ 0 \\ 0 \\ \varphi\left(x,y,z\right) \end{array}\right],
      \ dir = y \rightarrow \left[\begin{array}{c}0 \\ 0 \\ \varphi\left(x,y,z\right) \\ 0 \\ \varphi\left(x,y,z\right) \end{array}\right],
      \ dir = z \rightarrow \left[\begin{array}{c}0 \\ 0 \\ 0 \\ \varphi\left(x,y,z\right) \\ \varphi\left(x,y,z\right) \end{array}\right]
    \f}
    where \f$\varphi\f$ (#NavierStokes3D::grav_field_g) is computed in NavierStokes3D::GravityField().
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, AIAA Journal, http://dx.doi.org/10.2514/1.J054580.
*/
int NavierStokes3DSourceFunction(
                                  double  *f, /*!< Array to hold the computed source function */
                                  double  *u, /*!< Solution vector array */
                                  double  *x, /*!< Array of spatial coordinates (grid) */
                                  void    *s, /*!< Solver object of type #HyPar */
                                  void    *m, /*!< MPI object of type #MPIVariables */
                                  double  t,  /*!< Current simulation time */
                                  int     dir /*!< Spatial dimension (x, y, or z) */
                                )
{
  HyPar           *solver = (HyPar* )         s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->physics;

  int ghosts  = solver->ghosts;
  int *dim    = solver->dim_local;
  static int index[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_], offset[_MODEL_NDIMS_];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_);
  int i; for (i=0; i<_MODEL_NDIMS_; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,_MODEL_NDIMS_,-ghosts);

  int done = 0; _ArraySetValue_(index,_MODEL_NDIMS_,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,index,offset,ghosts,p);
    (f+_MODEL_NVARS_*p)[0] = 0.0;
    (f+_MODEL_NVARS_*p)[1] = param->grav_field_g[p] * (dir == _XDIR_);
    (f+_MODEL_NVARS_*p)[2] = param->grav_field_g[p] * (dir == _YDIR_);
    (f+_MODEL_NVARS_*p)[3] = param->grav_field_g[p] * (dir == _ZDIR_);
    (f+_MODEL_NVARS_*p)[4] = param->grav_field_g[p];
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,index,done);
  }

  return(0);
}

/*! Compute the "upwind" source function value at the interface: the upwinding is just the
    arithmetic average of the left and right biased interpolated values of the source function
    at the interface.
    \n\n
    References:
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, Submitted
*/
int NavierStokes3DSourceUpwind(
                                double  *fI,  /*!< Array to hold the computed "upwind" interface source function */
                                double  *fL,  /*!< Interface source function value computed using left-biased interpolation */ 
                                double  *fR,  /*!< Interface source function value computed using right-biased interpolation */
                                double  *u,   /*!< Solution vector array */
                                int     dir,  /*!< Spatial dimension (x,y, or z) */
                                void    *s,   /*!< Solver object of type #HyPar */
                                double  t     /*!< Current simulation time */
                              )
{
  HyPar *solver = (HyPar*) s;
  int   done,k, *dim  = solver->dim_local;
  _DECLARE_IERR_;


  static int index_outer[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_], 
             bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[dir] += 1;

  done = 0; _ArraySetValue_(index_outer,_MODEL_NDIMS_,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p;  _ArrayIndex1D_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      for (k = 0; k < _MODEL_NVARS_; k++) 
        (fI+_MODEL_NVARS_*p)[k] = 0.5 * ((fL+_MODEL_NVARS_*p)[k] + (fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}
