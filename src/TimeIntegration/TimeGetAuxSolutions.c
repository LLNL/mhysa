/*! @file TimeGetAuxSolutions.c
    @brief Returns any "auxiliary" solutions
    @author Debojyoti Ghosh
*/
#include <stdio.h>
#include <string.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Return auxiliary solution: Some time integrators may have the concept of 
  auxiliary solutions that they evolve along with the main solution #HyPar::u
  (these may be used for error estimation, for example). This function returns
  a pointer to such an auxiliary solution. Note that the auxiliary solution has
  the same dimensions and array layout as the main solution.
  + Call with the final argument \a n less than 0 to get the total number of
    auxiliary solutions (this value will be stored in the argument \a N at 
    exit).
  + Call with the final argument \a n less than or equal to zero to get the 
    \a n-th auxiliary solution (the argument \a uaux will point to this at
    exit). 
    
  Note that auxiliary solutions are numbered in the C convention: 0,1,...,N-1.

  Time integration methods which use auxiliary solutions currently implemented:
  + General Linear Methods with Global Error Estimators (GLM-GEE - #_GLM_GEE_)
    (see TimeGLMGEE(), TimeGLMGEEInitialize() ).
*/
int TimeGetAuxSolutions(
                          int     *N,      /*!< Number of auxiliary solutions */
                          double  **uaux,  /*!< Pointer to the array holding the auxiliary solution */
                          void    *s,      /*!< Solver object of type #HyPar */
                          int     n        /*!< Index of the auxiliary solution to return */
                       )
{
  HyPar           *solver = (HyPar*) s;
  TimeIntegration *TS     = (TimeIntegration*) solver->time_integrator;

  if (n >= 0) {
    if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
      GLMGEEParameters *params = (GLMGEEParameters*) solver->msti;
      *uaux = TS->U[params->r+n];
    }
  } else {
    if (!TS) *N = 0;
    else {
      if (!strcmp(solver->time_scheme,_GLM_GEE_)) {
        GLMGEEParameters *params = (GLMGEEParameters*) solver->msti;
        *N = params->r - 1;
      } else *N = 0;
    }
  }

  return(0);
}
