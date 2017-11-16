/*! @file HyperbolicFunction.c
    @author Debojyoti Ghosh
    @brief Compute the hyperbolic term of the governing equations
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ReconstructHyperbolic (double*,double*,double*,double*,int,void*,void*,double,int,
                                  int(*)(double*,double*,double*,double*,double*,double*,int,void*,double));
static int DefaultUpwinding      (double*,double*,double*,double*,double*,double*,int,void*,double);

/*! This function computes the hyperbolic term of the governing equations, expressed as follows:- 
    \f{equation}{
      {\bf F} \left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {\partial {\bf f}_d\left({\bf u}\right)} {\partial x_d}
    \f}
    using a conservative finite-difference discretization is space:
    \f{equation}{
      \hat{\bf F} \left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {1} {\Delta x_d} \left[ \hat{\bf f}_{d,j+1/2} - \hat{\bf f}_{d,j-1/2} \right]
    \f}
    where \f$d\f$ denotes the spatial dimension, \f$D\f$ denotes the total number of spatial dimensions, the hat denotes 
    the discretized quantity, and \f$j\f$ is the grid coordinate along dimension \f$d\f$.
    The approximation to the flux function \f${\bf f}_d\f$ at the interfaces \f$j\pm1/2\f$, denoted by \f$\hat{\bf f}_{d,j\pm 1/2}\f$,
    are computed using the function ReconstructHyperbolic().
*/
int HyperbolicFunction(
                        double  *hyp, /*!< Array to hold the computed hyperbolic term (shares the same layout as u */
                        double  *u,   /*!< Solution array */
                        void    *s,   /*!< Solver object of type #HyPar */
                        void    *m,   /*!< MPI object of type #MPIVariables */
                        double  t,    /*!< Current simulation time */
                        int     LimFlag,  /*!< Flag to indicate if the nonlinear coefficients for solution-dependent 
                                               interpolation method should be recomputed (see ReconstructHyperbolic() for
                                               an explanation on why this is needed) */
                        /*! Function pointer to the flux function for the hyperbolic term */
                        int(*FluxFunction)(double*,double*,int,void*,double), 
                        /*! Function pointer to the upwinding function for the hyperbolic term */
                        int(*UpwindFunction)(double*,double*,double*,double*,double*,double*,int,void*,double) 
                      )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           d, v, i, done;
  double        *FluxI  = solver->fluxI; /* interface flux     */
  double        *FluxC  = solver->fluxC; /* cell centered flux */
  _DECLARE_IERR_;

  int     ndims  = solver->ndims;
  int     nvars  = solver->nvars;
  int     ghosts = solver->ghosts;
  int     *dim   = solver->dim_local;
  int     size   = solver->npoints_local_wghosts;
  double  *x     = solver->x;
  double  *dxinv = solver->dxinv;
  int     index[ndims], index1[ndims], index2[ndims], dim_interface[ndims];

  LimFlag = (LimFlag && solver->flag_nonlinearinterp && solver->SetInterpLimiterVar);

  _ArraySetValue_(hyp,size*nvars,0.0);
  _ArraySetValue_(solver->StageBoundaryIntegral,2*ndims*nvars,0.0);
  if (!FluxFunction) return(0); /* zero hyperbolic term */
  solver->count_hyp++;

  int offset = 0;
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(dim,dim_interface,ndims); dim_interface[d]++;
    int size_cellcenter = 1; for (i = 0; i < ndims; i++) size_cellcenter *= (dim[i] + 2*ghosts);
    int size_interface = 1; for (i = 0; i < ndims; i++) size_interface *= dim_interface[i];

    /* evaluate cell-centered flux */
    IERR FluxFunction(FluxC,u,d,solver,t); CHECKERR(ierr);
    /* compute interface fluxes */
    IERR ReconstructHyperbolic(FluxI,FluxC,u,x+offset,d,solver,mpi,t,LimFlag,UpwindFunction); 
    CHECKERR(ierr);

    /* calculate the first derivative */
    done = 0; _ArraySetValue_(index,ndims,0);
    int p, p1, p2;
    while (!done) {
      _ArrayCopy1D_(index,index1,ndims);
      _ArrayCopy1D_(index,index2,ndims); index2[d]++;
      _ArrayIndex1D_(ndims,dim          ,index ,ghosts,p);
      _ArrayIndex1D_(ndims,dim_interface,index1,0     ,p1);
      _ArrayIndex1D_(ndims,dim_interface,index2,0     ,p2);
      for (v=0; v<nvars; v++) hyp[nvars*p+v] += dxinv[offset+ghosts+index[d]] 
                                              * (FluxI[nvars*p2+v]-FluxI[nvars*p1+v]);
      /* boundary flux integral */
      if (index[d] == 0) 
        for (v=0; v<nvars; v++) solver->StageBoundaryIntegral[(2*d+0)*nvars+v] -= FluxI[nvars*p1+v];
      if (index[d] == dim[d]-1) 
        for (v=0; v<nvars; v++) solver->StageBoundaryIntegral[(2*d+1)*nvars+v] += FluxI[nvars*p2+v];

      _ArrayIncrementIndex_(ndims,dim,index,done);
    }

    offset += dim[d] + 2*ghosts;
  }

  if (solver->flag_ib) _ArrayBlockMultiply_(hyp,solver->iblank,size,nvars);
  return(0);
}

/*! This function computes the numerical flux \f$\hat{\bf f}_{j+1/2}\f$ at the interface from the cell-centered 
    flux function \f${\bf f}_j\f$. This happens in two steps:-

    \b Interpolation: High-order accurate approximations to the flux at the interface \f$j+1/2\f$ are computed from
    the cell-centered flux with left- and right-biased interpolation methods. This is done by the
    #HyPar::InterpolateInterfacesHyp function. This can be expressed as follows:
    \f{align}{
      \hat{\bf f}^L_{j+1/2} &= \mathcal{I}\left({\bf f}_j,+1\right), \\
      \hat{\bf f}^R_{j+1/2} &= \mathcal{I}\left({\bf f}_j,-1\right),
    \f}
    where the \f$\pm 1\f$ indicates the interpolation bias, and \f$\mathcal{I}\f$ is the interpolation operator 
    pointed to by #HyPar::InterpolateInterfacesHyp (see \b src/InterpolationFunctions for all the available operators).
    The interface values of the solution are similarly computed:
    \f{align}{
      \hat{\bf u}^L_{j+1/2} &= \mathcal{I}\left({\bf u}_j,+1\right), \\
      \hat{\bf u}^R_{j+1/2} &= \mathcal{I}\left({\bf u}_j,-1\right),
    \f}
    The specific choice of \f$\mathcal{I}\f$ is set based on #HyPar::spatial_scheme_hyp.

    \b Upwinding: The final flux at the interface is computed as
    \f{equation}{
      \hat{\bf f}_{j+1/2} = \mathcal{U}\left( \hat{\bf f}^L_{j+1/2}, \hat{\bf f}^R_{j+1/2}, \hat{\bf u}^L_{j+1/2}, \hat{\bf u}^R_{j+1/2} \right),
    \f}
    where \f$\mathcal{U}\f$ denotes the upwinding function UpwindFunction() passed as an argument (if NULL, DefaultUpwinding() is used). The
    upwinding function is specified by the physical model.

    \b Note:
    Solution-dependent, nonlinear interpolation methods (such as WENO, CRWENO) are implemented in a way that separates the calculation of the 
    nonlinear interpolation weights (based on, say, the smoothness of the flux function), and the actual evaluation of the interpolant, into 
    different functions. This allows the flexibility to choose if and when the nonlinear coefficients are evaluated (or previously computed 
    values are reused). Some possible scenarios are:
    + For explicit time integration, they are computed every time the hyperbolic flux term is being computed.
    + For implicit time integration, consistency or linearization may require that they be computed and "frozen" 
      at the beginning of a time step or stage. 

    The argument \b LimFlag controls this behavior:
    + LimFlag = 1 means recompute the nonlinear coefficients.
    + LimFlag = 0 means reuse the the previously computed coefficients.
*/
int ReconstructHyperbolic(
                            double  *fluxI,     /*!< Array to hold the computed interface fluxes. This array does not
                                                     have ghost points. The dimensions are the same as those of u without
                                                     ghost points in all dimensions, except along dir, where it is one more */
                            double  *fluxC,     /*!< Array of the flux function computed at the cell centers 
                                                     (same layout as u) */
                            double  *u,         /*!< Solution array */
                            double  *x,         /*!< Array of spatial coordinates */
                            int     dir,        /*!< Spatial dimension along which to reconstruct the interface fluxes */
                            void    *s,         /*!< Solver object of type #HyPar */
                            void    *m,         /*!< MPI object of type #MPIVariables */
                            double  t,          /*!< Current solution time */
                            int     LimFlag,    /*!< Flag to indicate if the nonlinear coefficients for solution-dependent
                                                     interpolation method should be recomputed */
                            /*! Function pointer to the upwinding function for the interface flux computation. If NULL, 
                                DefaultUpwinding() will be used. */
                            int(*UpwindFunction)(double*,double*,double*,double*,double*,double*,int,void*,double)
                          )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  double *uC     = NULL;
  double *uL     = solver->uL;
  double *uR     = solver->uR;
  double *fluxL  = solver->fL;
  double *fluxR  = solver->fR;

  /* 
    precalculate the non-linear interpolation coefficients if required 
    else reuse the weights previously calculated
  */
  if (LimFlag) IERR solver->SetInterpLimiterVar(fluxC,u,x,dir,solver,mpi);

  /* if defined, calculate the modified u-function to be used for upwinding
     e.g.: used in well-balanced schemes for Euler/Navier-Stokes with gravity
     otherwise, just copy u to uC */
  if (solver->UFunction) {
    uC = solver->uC;
    IERR solver->UFunction(uC,u,dir,solver,mpi,t); CHECKERR(ierr);
  } else uC = u;

  /* Interpolation -> to calculate left and right-biased interface flux and state variable*/
  IERR solver->InterpolateInterfacesHyp(uL   ,uC   ,u,x, 1,dir,solver,mpi,1); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(uR   ,uC   ,u,x,-1,dir,solver,mpi,1); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(fluxL,fluxC,u,x, 1,dir,solver,mpi,0); CHECKERR(ierr);
  IERR solver->InterpolateInterfacesHyp(fluxR,fluxC,u,x,-1,dir,solver,mpi,0); CHECKERR(ierr);

  /* Upwind -> to calculate the final interface flux */
  if (UpwindFunction) { IERR UpwindFunction   (fluxI,fluxL,fluxR,uL  ,uR  ,u   ,dir,solver,t); CHECKERR(ierr); }
  else                { IERR DefaultUpwinding (fluxI,fluxL,fluxR,NULL,NULL,NULL,dir,solver,t); CHECKERR(ierr); }

  return(0);
}

/*! If no upwinding scheme is specified, this function defines the "upwind" flux as the 
    arithmetic mean of the left- and right-biased fluxes. */
int DefaultUpwinding(
                      double  *fI,  /*!< Computed upwind interface flux */
                      double  *fL,  /*!< Left-biased reconstructed interface flux */
                      double  *fR,  /*!< Right-biased reconstructed interface flux */
                      double  *uL,  /*!< Left-biased reconstructed interface solution */
                      double  *uR,  /*!< Right-biased reconstructed interface solution */
                      double  *u,   /*!< Cell-centered solution */
                      int     dir,  /*!< Spatial dimension */
                      void    *s,   /*!< Solver object of type #HyPar */
                      double  t     /*!< Current solution time */
                    )
{
  HyPar *solver = (HyPar*)    s;
  int   done;

  int *dim  = solver->dim_local;
  int ndims = solver->ndims;
  int nvars = solver->nvars;

  int bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  done = 0; int index_outer[ndims], index_inter[ndims];
  _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      int v; for (v=0; v<nvars; v++) fI[nvars*p+v] = 0.5 * (fL[nvars*p+v]+fR[nvars*p+v]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
