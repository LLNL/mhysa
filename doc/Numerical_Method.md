Numerical Method
================

MHYSA solves the following partial differential equation (PDE) using a conservative finite-difference
algorithm on a Cartesian grid.
\f{equation}{
  \frac {\partial {\bf u}} {\partial t} = {\bf F}_{\rm hyp}\left({\bf u}\right) + {\bf F}_{\rm par}\left({\bf u}\right) + {\bf F}_{\rm sou}\left({\bf u}\right)
\f}
where \f${\bf F}_{\rm hyp}\f$ is the hyperbolic term, \f${\bf F}_{\rm par}\f$ is the parabolic term, and 
\f${\bf F}_{\rm sou}\f$ is the source term. Each of these is discretized in space as described below (in
the section "Spatial Discretization"), to 
obtain the following semi-discrete ordinary differential equation (ODE) in time:
\f{equation}{
  \frac {d {\bf u}} {d t} = \hat{\bf F}_{\rm hyp}\left({\bf u}\right) + \hat{\bf F}_{\rm par}\left({\bf u}\right) + \hat{\bf F}_{\rm sou}\left({\bf u}\right)
\f}
where \f$\hat{\left(\cdot\right)}\f$ represents the spatially discretized terms. The governing PDE can be
of any space dimension. 

\section spatial_discretization Spatial Discretization

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
#HyPar::spatial_scheme_hyp (read from \a "solver.inp" in ReadInputs()). See interpolation.h for all the spatial discretization schemes implemented.

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
#HyPar::spatial_scheme_par (read from \a "solver.inp" in ReadInputs()). See interpolation.h, firstderivative.h, and secondderivative.h for the spatial
disretization methods implemented.

Depending on which of the above forms are used, #HyPar::ParabolicFunction points to either of ParabolicFunctionCons1Stage(), ParabolicFunctionNC1Stage(),
ParabolicFunctionNC2Stage(), or ParabolicFunctionNC1_5Stage() (set in InitializeSolvers()).

Source term
-----------

#HyPar::SourceFunction points to SourceFunction(). There is no discretization involved in general, and this function
just calls the physical model-specific source function to which #HyPar::SFunction points.


\section time_integration Time Integration

Native Time Integrators
-----------------------

The semi-discrete ODE is integrated in time using explicit multi-stage time integration methods. The ODE can be written as:
\f{equation}{
  \frac {d {\bf u}} {d t} = {\bf F}\left({\bf u}\right)
\f}
where
\f{equation}{
  {\bf F}\left({\bf u}\right) = \hat{\bf F}_{\rm hyp}\left({\bf u}\right) + \hat{\bf F}_{\rm par}\left({\bf u}\right) + \hat{\bf F}_{\rm sou}\left({\bf u}\right)
\f}
The following explicit time integration methods are implemented in MHYSA (see timeintegration.h):
+ Forward Euler - TimeForwardEuler(), #_FORWARD_EULER_
+ Explicit Runge-Kutta - TimeRK(), #_RK_
+ Explicit General Linear Methods with Global Error Estimation - TimeGLMGEE(), #_GLM_GEE_

\sa Solve()


PETSc Time Integrators
----------------------

If compiled with PETSc (https://www.mcs.anl.gov/petsc/), MHYSA can use all the time integration methods and features implemented in the \b TS module of
PETSc. See the following for relevant documentation of PETSc time integrators:
+ http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/index.html
+ http://www.mcs.anl.gov/petsc/petsc-current/src/ts/examples/tutorials/index.html
+ http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf (Chapter 6: Scalable ODE and DAE Solvers)

In addition to explicit time integration, the semi-discrete ODE can be solved using 
+ \b Implicit \b methods (Eg. backward Euler (TSBEULER - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSBEULER.html),
Crank-Nicholson (TSCN - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSCN.html#TSCN), \f$\theta\f$-method (TSTHETA - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSTHETA.html#TSTHETA), 
etc.)
+ <B>Semi-implicit (IMEX) methods</B> (TSARKIMEX - http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSARKIMEX.html)

Implementation: see petscinterface.h and SolvePETSc().

\b Implicit and \b IMEX time integration: 
+ The Jacobian-free approach is used to compute the Jacobian of the implicit term (i.e., the action of the Jacobian on a vector
  is approximated using a directional derivative). Use the flag <B>-jfnk_epsilon \<value\></B> to specify the parameter \f$\epsilon\f$
  for the directional derivative computation (default: \f$10^{-6}\f$). See PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK().
+ A preconditioner can only be used for physical models that define the Jacobian (#HyPar::JFunction) 
  (for example, LinearADRJacobian(), Euler1DJacobian(), NavierStokes3DJacobian(), etc).
  The flag <B>-with_pc</B> should be specified to use a preconditioner.

<B>IMEX Time Integration</B>: For implicit-explicit (IMEX) time integration, the semi-discrete ODE can be written as follows:
\f{equation}{
  \frac {d {\bf u}} {d t} = {\bf F}\left({\bf u}\right) + {\bf G}\left({\bf u}\right)
\f}
where \f${\bf F}\left({\bf u}\right)\f$ is integrated explicitly in time and \f${\bf G}\left({\bf u}\right)\f$
is integrated implicitly in time. The following flags (specified in the command line or in the <B>.petscrc</B>
file) can be used to specify which of the hyperbolic, parabolic, or source terms are treated explicitly, and 
which are treated implicitly.

Term                                                        |  Explicit             | Implicit
------------------------------------------------------------|-----------------------|---------------------
Hyperbolic \f$\hat{\bf F}_{\rm hyp}\left({\bf u}\right)\f$  | -hyperbolic_explicit  | -hyperbolic_implicit
Parabolic \f$\hat{\bf F}_{\rm par}\left({\bf u}\right)\f$   | -parabolic_explicit   | -parabolic_implicit
Source \f$\hat{\bf F}_{\rm sou}\left({\bf u}\right)\f$      | -source_explicit      | -source_implicit

+ If contradictory flags are specified, i.e.,

        -parabolic_explicit -parabolic_implicit

  the flag specifying implicit treatment takes precedence.

+ In addition, if a partitioning of the hyperbolic flux is defined and is being used (#HyPar::SplitHyperbolicFlux),
  i.e,
  \f{equation}{
    \hat{\bf F}_{\rm hyp}\left({\bf u}\right) = \left[\hat{\bf F}_{\rm hyp}\left({\bf u}\right) - \delta\hat{\bf F}_{\rm hyp}\left({\bf u}\right)\right] + \delta\hat{\bf F}_{\rm hyp}\left({\bf u}\right)
  \f}
  the following flags can be used to specify which of these two terms are treated explicitly and which are treated implicitly.

    Term                                                                                                            |  Explicit               | Implicit
    ----------------------------------------------------------------------------------------------------------------|-------------------------|---------------------
    \f$\left[\hat{\bf F}_{\rm hyp}\left({\bf u}\right) - \delta\hat{\bf F}_{\rm hyp}\left({\bf u}\right)\right]\f$  | -hyperbolic_f_explicit  | -hyperbolic_f_implicit
    \f$\delta\hat{\bf F}_{\rm hyp}\left({\bf u}\right)\f$                                                           | -hyperbolic_df_explicit | -hyperbolic_df_implicit


