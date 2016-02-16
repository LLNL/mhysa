/*! @file TimeGLMGEECleanup.c
    @brief Clean up the #_GLM_GEE_ time integrator
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <string.h>
#include <timeintegration.h>

/*! Clean up allocations for the GLM-GEE (#_GLM_GEE_) time integration
    method: This function frees the arrays used to store the Butcher
    tableaux 
*/
int TimeGLMGEECleanup(void *s /*!< Object of type #GLMGEEParameters*/)
{
  GLMGEEParameters *params = (GLMGEEParameters*) s;
  if (params->A_yyt ) free(params->A_yyt);
  if (params->B_yyt ) free(params->B_yyt);
  if (params->C_yyt ) free(params->C_yyt);
  if (params->D_yyt ) free(params->D_yyt);
  if (params->c_yyt ) free(params->c_yyt);
  if (params->A_yeps) free(params->A_yeps);
  if (params->B_yeps) free(params->B_yeps);
  if (params->C_yeps) free(params->C_yeps);
  if (params->D_yeps) free(params->D_yeps);
  if (params->c_yeps) free(params->c_yeps);
  return(0);
}
