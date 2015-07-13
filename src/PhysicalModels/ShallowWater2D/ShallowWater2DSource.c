/*! @file ShallowWater2DSource.c
    @author Debojyoti Ghosh
    @brief Contains the functions to compute the source terms for the 2D shallow water equations.
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater2d.h>
#include <mpivars.h>
#include <hypar.h>

static int ShallowWater2DSourceFunction1  (double*,double*,double*,void*,void*,double,int);
static int ShallowWater2DSourceFunction2  (double*,double*,double*,void*,void*,double,int);
int        ApplyBoundaryConditions        (void*,void*,double*,double*,int,double);

/*! Compute the source terms for the 2D shallow water equations. The source term
    is computed according to the balanced formulation introduced in the reference below.
    The source term is reformulated and "discretized" in a similar fashion as the hyperbolic
    flux to ensure that the hydrostatic balance is maintained to machine precision.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater2DSource(
                          double  *source, /*!< Computed source terms (array size & layout same as u) */
                          double  *u,      /*!< Solution (conserved variables) */
                          void    *s,      /*!< Solver object of type #HyPar */
                          void    *m,      /*!< MPI object of type #MPIVariables */
                          double  t        /*!< Current solution time */
                        )
{
  HyPar          *solver = (HyPar* ) s;
  MPIVariables   *mpi = (MPIVariables*) m;
  ShallowWater2D *param  = (ShallowWater2D*) solver->physics;

  int     v, done, p, p1, p2;
  double  *SourceI = solver->fluxI; /* interace source term       */
  double  *SourceC = solver->fluxC; /* cell-centered source term  */
  double  *SourceL = solver->fL;
  double  *SourceR = solver->fR;

  int     ndims   = solver->ndims;
  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  double  *x      = solver->x;
  double  *dxinv  = solver->dxinv;
  int     index[ndims],index1[ndims],index2[ndims],dim_interface[ndims];

  /* Along X */
  
  /* set interface dimensions */
  _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_XDIR_]++;

  /* calculate the first source function */
  IERR ShallowWater2DSourceFunction1(SourceC,u,x,solver,mpi,t,_XDIR_); CHECKERR(ierr);
  /* calculate the left and right interface source terms */
  IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
  /* calculate the final interface source term */
  IERR param->SourceUpwind(SourceI,SourceL,SourceR,u,_XDIR_,solver,t);
  /* calculate the final cell-centered source term */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_XDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
    double dx_inverse;   _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dx_inverse);
    for (v=0; v<_MODEL_NVARS_; v++)
      source[_MODEL_NVARS_*p+v] += (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse;
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  /* calculate the second source function */
  IERR ShallowWater2DSourceFunction2(SourceC,u,x,solver,mpi,t,_XDIR_); CHECKERR(ierr);
  /* calculate the left and right interface source terms */
  IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_XDIR_,solver,mpi,0); CHECKERR(ierr);
  /* calculate the final interface source term */
  IERR param->SourceUpwind(SourceI,SourceL,SourceR,u,_XDIR_,solver,t);
  /* calculate the final cell-centered source term */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_XDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
    double dx_inverse;  _GetCoordinate_(_XDIR_,index[_XDIR_],dim,ghosts,dxinv,dx_inverse);
    double h, uvel, vvel; _ShallowWater2DGetFlowVar_((u+_MODEL_NVARS_*p),h,uvel,vvel);
    double term[_MODEL_NVARS_] = { 0.0, -param->g * (h + param->b[p]), 0.0 };
    for (v=0; v<_MODEL_NVARS_; v++)
      source[_MODEL_NVARS_*p+v] += term[v]*(SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dx_inverse;
    uvel = vvel = h; /* useless statement to avoid compiler warning */
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  /* Along Y */

  /* set interface dimensions */
  _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[_YDIR_]++;

  /* calculate the first source function */
  IERR ShallowWater2DSourceFunction1(SourceC,u,x,solver,mpi,t,_YDIR_); CHECKERR(ierr);
  /* calculate the left and right interface source terms */
  IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_YDIR_,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_YDIR_,solver,mpi,0); CHECKERR(ierr);
  /* calculate the final interface source term */
  IERR param->SourceUpwind(SourceI,SourceL,SourceR,u,_YDIR_,solver,t);
  /* calculate the final cell-centered source term */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_YDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
    double dy_inverse;   _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,dxinv,dy_inverse);
    for (v=0; v<_MODEL_NVARS_; v++)
      source[_MODEL_NVARS_*p+v] += (SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dy_inverse;
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  /* calculate the second source function */
  IERR ShallowWater2DSourceFunction2(SourceC,u,x,solver,mpi,t,_YDIR_); CHECKERR(ierr);
  /* calculate the left and right interface source terms */
  IERR solver->InterpolateInterfacesHyp(SourceL,SourceC,u,x, 1,_YDIR_,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(SourceR,SourceC,u,x,-1,_YDIR_,solver,mpi,0); CHECKERR(ierr);
  /* calculate the final interface source term */
  IERR param->SourceUpwind(SourceI,SourceL,SourceR,u,_YDIR_,solver,t);
  /* calculate the final cell-centered source term */
  done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index,index1,ndims);
    _ArrayCopy1D_(index,index2,ndims); index2[_YDIR_]++;
    _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p );
    _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
    _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
    double dy_inverse;  _GetCoordinate_(_YDIR_,index[_YDIR_],dim,ghosts,dxinv,dy_inverse);
    double h, uvel, vvel; _ShallowWater2DGetFlowVar_((u+_MODEL_NVARS_*p),h,uvel,vvel);
    double term[_MODEL_NVARS_] = { 0.0, 0.0, -param->g * (h + param->b[p]) };
    for (v=0; v<_MODEL_NVARS_; v++)
      source[_MODEL_NVARS_*p+v] += term[v]*(SourceI[_MODEL_NVARS_*p2+v]-SourceI[_MODEL_NVARS_*p1+v])*dy_inverse;
    uvel = vvel = h; /* useless statement to avoid compiler warning */
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(0);
}

/*! Compute the first source function that is then "discretized" in a way similar to 
    the hyperbolic flux function for the balanced formulation introduced in the reference below.
    The source term is reformulated and "discretized" in a similar fashion as the hyperbolic
    flux to ensure that the hydrostatic balance is maintained to machine precision.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
   \n
   This function computes:
   \f{equation}{
     {\bf S}_2 = 
        \left\{
          \begin{array}{cc} 
            \left[ \begin{array}{c} 0 \\ \frac{1}{2}gb^2 \\ 0 \end{array}\right] & {\rm dir} = x \\
            \left[ \begin{array}{c} 0 \\ 0 \\ \frac{1}{2}gb^2 \end{array}\right] & {\rm dir} = y 
          \end{array} 
        \right.
   \f}
*/
int ShallowWater2DSourceFunction1(
                                  double  *f, /*!< Computed source function (array size and layout same as u) */
                                  double  *u, /*!< Solution (conserved variables) */
                                  double  *x, /*!< Spatial coordinates */
                                  void    *s, /*!< Solver object of type #HyPar */
                                  void    *m, /*!< MPI object of type #MPIVariables */
                                  double  t,  /*!< Current solution time */
                                  int     dir /*!< Spatial dimension */
                                 )
{
  HyPar          *solver = (HyPar* ) s;
  ShallowWater2D *param  = (ShallowWater2D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    (f+_MODEL_NVARS_*p)[0] = 0.0;
    (f+_MODEL_NVARS_*p)[1] = (dir == _XDIR_ ? 0.5 * param->g * param->b[p] * param->b[p] : 0.0 );
    (f+_MODEL_NVARS_*p)[2] = (dir == _YDIR_ ? 0.5 * param->g * param->b[p] * param->b[p] : 0.0 );
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}

/*! Compute the second source function that is then "discretized" in a way similar to 
    the hyperbolic flux function for the balanced formulation introduced in the reference below.
    The source term is reformulated and "discretized" in a similar fashion as the hyperbolic
    flux to ensure that the hydrostatic balance is maintained to machine precision.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
   \n
   This function computes:
   \f{equation}{
     {\bf S}_2 = 
        \left\{
          \begin{array}{cc} 
            \left[ \begin{array}{c} 0 \\ b \\ 0 \end{array}\right] & {\rm dir} = x \\
            \left[ \begin{array}{c} 0 \\ 0 \\ b \end{array}\right] & {\rm dir} = y 
          \end{array} 
        \right.
   \f}
*/
int ShallowWater2DSourceFunction2(
                                  double  *f, /*!< Computed source function (array size and layout same as u) */
                                  double  *u, /*!< Solution (conserved variables) */
                                  double  *x, /*!< Spatial coordinates */
                                  void    *s, /*!< Solver object of type #HyPar */
                                  void    *m, /*!< MPI object of type #MPIVariables */
                                  double  t,  /*!< Current solution time */
                                  int     dir /*!< Spatial dimension */
                                 )
{
  HyPar          *solver = (HyPar* ) s;
  ShallowWater2D *param  = (ShallowWater2D*) solver->physics;

  int     ghosts  = solver->ghosts;
  int     *dim    = solver->dim_local;
  int     ndims   = solver->ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    (f+_MODEL_NVARS_*p)[0] = 0.0;
    (f+_MODEL_NVARS_*p)[1] = (dir == _XDIR_ ? param->b[p] : 0.0 );
    (f+_MODEL_NVARS_*p)[2] = (dir == _YDIR_ ? param->b[p] : 0.0 );
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
