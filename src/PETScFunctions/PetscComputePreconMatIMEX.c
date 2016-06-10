/*! @file PetscComputePreconMatIMEX.c
    @brief Contains the function to assemble the preconditioning matrix
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscComputePreconMatIMEX"
/*!
  Compute and assemble the preconditioning matrix for the implicit-explicit (IMEX) time integration
  of the governing equations: The ODE, obtained after discretizing the governing PDE in space,
  is expressed as follows (for the purpose of IMEX time integration):
  \f{eqnarray}{
    \frac {d{\bf U}}{dt} &=& {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right), \\
    \Rightarrow \frac {d{\bf U}}{dt} - {\bf G}\left({\bf U}\right) &=& {\bf F}\left({\bf U}\right), 
  \f}
  where \f${\bf F}\f$ is non-stiff and integrated in time explicitly, and \f${\bf G}\f$
  is stiff and integrated in time implicitly, and \f${\bf U}\f$ represents the entire
  solution vector.

    Note that \f${\bf G}\left({\bf U}\right)\f$ represents all terms that the user has indicated to be
    integrated in time implicitly (#PETScContext::flag_hyperbolic_f, #PETScContext::flag_hyperbolic_df,
    #PETScContext::flag_hyperbolic, #PETScContext::flag_parabolic, and #PETScContext::flag_source).

  The Jacobian of the implicit part is thus given by:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \frac {\partial {\bf G}} {\partial {\bf U}} \right]
  \f}
  where \f$\alpha\f$ is the shift coefficient (#PETScContext::shift) of the time integration method.
  
  Currently, this function is implemented for the case where the implicitly-treated \f${\bf G}\f$
  is the spatially discretized part or whole of a hyperbolic term, i.e.,
  \f{equation}{
    {\bf G}\left({\bf U}\right) = \mathcal{D}\left\{{\bf g}\left({\bf u}\right)\right\} \approx \nabla \cdot {\bf g}\left({\bf u}\right),
  \f}
  with the governing PDE as
  \f{equation}{
    \frac {\partial {\bf u}} {\partial t} + \cdots + \nabla \cdot {\bf g}\left({\bf u}\right) + \cdots = \cdots,
  \f}
  and \f$\mathcal{D}\f$ representing the spatial discretization method. Thus, the Jacobian can be written as
  follows:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \mathcal{D}\left\{\frac {\partial {\bf g}} {\partial {\bf u}}\right\} \right]
  \f}
  The preconditioning matrix is usually a close approximation of the actual Jacobian matrix, where the actual
  Jacobian may be too expensive to evaluate and assemble. In this function, the preconditioning matrix is 
  the following approximation of the actual Jacobian:
  \f{equation}{
    {\bf J}_p = \left[\alpha{\bf I} - \mathcal{D}^{\left(1\right)}\left\{\frac {\partial {\bf g}} {\partial {\bf u}}\right\} \right] \approx {\bf J},
  \f}
  where \f$\mathcal{D}^{\left(1\right)}\f$ represents a 1st order upwind discretization operator. The matrix \f${\bf J}_p\f$
  is provided to the preconditioner. Note that #HyPar::JFunction is defined by the specific physics being solved, and computes 
  \f$\partial {\bf g}/ \partial {\bf u}\f$ at a grid point.

  + See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/index.html for more information on PETSc preconditioners.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (http://www.mcs.anl.gov/petsc/petsc-current/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
int PetscComputePreconMatIMEX(
                              Mat Pmat,   /*!< Preconditioning matrix to construct */
                              Vec Y,      /*!< Solution vector */
                              void *ctxt  /*!< Application context */
                             )
{
  /* Same implementation as PetscComputePreconMatImpl() */
  return(PetscComputePreconMatImpl(Pmat,Y,ctxt));
}

#endif
