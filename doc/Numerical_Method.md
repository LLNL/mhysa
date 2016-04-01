Numerical Method
================

HyPar solves the following partial differential equation (PDE) using a conservative finite-difference
algorithm on a Cartesian grid.
\f{equation}{
  \frac {\partial {\bf u}} {\partial t} = {\bf F}_{\rm hyp}\left({\bf u}\right) + {\bf F}_{\rm par}\left({\bf u}\right) + {\bf F}_{\rm sou}\left({\bf u}\right)
\f}
where \f${\bf F}_{\rm hyp}\f$ is the hyperbolic term, \f${\bf F}_{\rm par}\f$ is the parabolic term, and 
\f${\bf F}_{\rm sou}\f$ is the source term. Each of these is discretized in space as described below, to 
obtain the following semi-discrete ordinary differential equation (ODE) in time:
\f{equation}{
  \frac {d {\bf u}} {d t} = \hat{\bf F}_{\rm hyp}\left({\bf u}\right) + \hat{\bf F}_{\rm par}\left({\bf u}\right) + \hat{\bf F}_{\rm sou}\left({\bf u}\right)
\f}
where \f$\hat{\left(\cdot\right)}\f$ represents the spatially discretized terms. The governing PDE can be
of any space dimension. The semi-discrete ODE is integrated in time by
+ Using native time integrators (see timeintegration.h, Solve()).
+ PETSc (https://www.mcs.anl.gov/petsc/) - specifically the TS module (see petscinterface.h, SolvePETSc()).

Hyperbolic term
---------------

The hyperbolic term is of the following form:
\f{equation}{
  {\bf F}_{\rm hyp}\left({\bf u}\right) = -\sum_{d=0}^{D-1} \frac{\partial {\bf f}_d\left({\bf u}\right)}{\partial x_d}
\f}
and is discretized as:
\f{equation}{
  {\bf F}_{\rm hyp}\left({\bf u}\right) \approx - \sum_{d=0}^{D-1} \frac{1}{\Delta x_d} \left[ \hat{\bf f}_{d,j+1/2} - \hat{\bf f}_{d,j-1/2} \right]
\f}
where \f$d\f$ is the spatial dimension, \f$D\f$ is the total number of spatial dimensions, \f$j\f$ denotes the grid index along \f$d\f$. This
is implemented in HyperbolicFunction().

The numerical approximation \f$\hat{\bf f}_{d,j+1/2}\f$ of the primitive of the flux \f${\bf f}_d\left({\bf u}\right)\f$ at the grid interfaces \f$j+1/2\f$ is expressed
as
\f{equation}{
  \hat{\bf f}_{d,j+1/2} = \mathcal{U}\left(\hat{\bf f}^L_{d,j+1/2},\hat{\bf f}^R_{d,j+1/2},\hat{\bf u}^L_{d,j+1/2},\hat{\bf u}^R_{d,j+1/2}\right)
\f}
where \f$\mathcal{U}\f$ is an upwinding function. #HyPar::Upwind points to the physical model-specific upwinding function that implements
\f$\mathcal{U}\f$, and is set by the initialization function of a specific physical model (for example Euler1DInitialize()). The physical model
is specified by setting #HyPar::model (read from \a "solver.inp" in ReadInputs()).

\f$\hat{\bf f}^{L,R}_{d,j+1/2}\f$ are the left- and right-biased numerically interpolated values of the primitive of the flux \f${\bf f}_d\left({\bf u}\right)\f$
at the grid interfaces and are computed using #HyPar::InterpolateInterfacesHyp. They are initialized in InitializeSolvers() based on the value of 
#HyPar::spatial_scheme_hyp (read from \a "solver.inp" in ReadInputs()).

#HyPar::HyperbolicFunction points to HyperbolicFunction().

Parabolic term
--------------

The parabolic term can take two different forms as described below:-

\b No \b cross-derivatives: In this form, the parabolic term is of the following form:
\f{equation}{
  {\bf F}_{\rm par}\left({\bf u}\right) = \sum_{d=0}^{D-1} \frac {\partial^2 {\bf g}_d\left({\bf u}\right)} {\partial x_d^2}
\f}
where \f$d\f$ is the spatial dimension, and \f$D\f$ is the total number of spatial dimensions. If the parabolic function is 
in this form, then the physical model must specify #HyPar::GFunction which must point to the function that computes 
\f${\bf g}_d\left({\bf u}\right)\f$. In this case, the spatial discretization is carried out in one of two ways:-

+ Conservative discretization, implemented in ParabolicFunctionCons1Stage() and invoked by setting #HyPar::spatial_type_par 
  to #_CONS_1STAGE_: The parabolic term is discretized as
\f{equation}{
  {\bf F}_{\rm par}\left({\bf u}\right) \approx \sum_{d=0}^{D-1} \frac {1}{\Delta x_d^2} \left[ \hat{\bf g}_{d,j+1/2} - \hat{\bf g}_{d,j-1/2} \right]
\f}
  where \f$j\f$ denotes the grid index along \f$d\f$. \f$\hat{\bf g}_{d,j+1/2}\f$ is the numerical approximation to the second primitive
  of \f${\bf g}_d\left({\bf u}\right)\f$ and is computed using #HyPar::InterpolateInterfacesPar.

+ Non-conservative, 1-Stage discretization, implemented in ParabolicFunctionNC1Stage() and invoked by setting #HyPar::spatial_type_par
  to #_NC_1STAGE_: The parabolic term is discretized as
\f{equation}{
  {\bf F}_{\rm par}\left({\bf u}\right) \approx \sum_{d=0}^{D-1} \frac {1}{\Delta x_d^2} \left[ \mathcal{L}_d\left({\bf g}_d\right) \right]
\f}
  where \f$\mathcal{L}\f$ represents the finite-difference Laplacian operator, and is computed using #HyPar::SecondDerivativePar.

\b With \b cross-derivatives: In this form, the parabolic term is of the following form:
\f{equation}{
  {\bf F}_{\rm par}\left({\bf u}\right) = \sum_{d1=0}^{D-1}\sum_{d2=0}^{D-1} \frac {\partial^2 h_{d1,d2}\left({\bf u}\right)} {\partial x_{d1} \partial x_{d2}}
\f}
where \f$d1,d2\f$ are spatial dimensions, \f$D\f$ is the total number of spatial dimensions. If the parabolic function is
in this form, then the physical model must specify #HyPar::HFunction which must point to the function that computes
\f${\bf h}_{d1,d2}\left({\bf u}\right)\f$. In this case, the spatial discretization is carried out in one of two ways:-

+ Non-conservative 2-stage discretization, implemented in ParabolicFunctionNC2Stage() and invoked by setting #HyPar::spatial_type_par
  to #_NC_2STAGE_: The parabolic term is discretized as
\f{equation}{
  {\bf F}_{\rm par}\left({\bf u}\right) \approx \sum_{d1=0}^{D-1}\sum_{d2=0}^{D-1} \frac {1}{\Delta x_{d1} \Delta x_{d2}} \left[ \mathcal{D}_{d1}\left(\mathcal{D}_{d2}\left({\bf g}_d\right)\right) \right]
\f}
  where \f$\mathcal{D}_d\f$ denotes the finite-difference first derivative operator along spatial dimension \f$d\f$, and is computed using #HyPar::FirstDerivativePar.

+ Non-conservative, "1.5-Stage" discretization, implemented in ParabolicFunctionNC1_5Stage() and invoked by setting #HyPar::spatial_type_par
  to #_NC_1_5STAGE_: The parabolic term is discretized as
\f{equation}{
  {\bf F}_{\rm par}\left({\bf u}\right) \approx \sum_{d1=0}^{D-1}\sum_{d2=0,d2 \ne d1}^{D-1} \frac {1}{\Delta x_{d1} \Delta x_{d2}} \left[ \mathcal{D}_{d1}\left(\mathcal{D}_{d2}\left({\bf g}_d\right)\right) \right]
                  + \sum_{d=0}^{D-1} \frac {1}{\Delta x_d^2} \left[ \mathcal{L}_d\left({\bf g}_d\right) \right]
\f}
  where \f$\mathcal{D}_d\f$ denotes the finite-difference first derivative operator along spatial dimension \f$d\f$ (computed using #HyPar::FirstDerivativePar).
  and \f$\mathcal{L}\f$ represents the finite-difference Laplacian operator (computed using #HyPar::SecondDerivativePar).

The function pointers #HyPar::InterpolateInterfacesPar, #HyPar::SecondDerivativePar, and #HyPar::FirstDerivativePar are set in InitializeSolvers() based on the value of 
#HyPar::spatial_scheme_par (read from \a "solver.inp" in ReadInputs()).

Depending on which of the above forms are used, #HyPar::ParabolicFunction points to either of ParabolicFunctionCons1Stage(), ParabolicFunctionNC1Stage(),
ParabolicFunctionNC2Stage(), or ParabolicFunctionNC1_5Stage() (set in InitializeSolvers()).

Source term
-----------

#HyPar::SourceFunction points to SourceFunction(). There is no discretization involved in general, and this function
just calls the physical model-specific source function to which #HyPar::SFunction points.
  
