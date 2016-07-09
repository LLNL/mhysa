/*! @file IBIdentifyMode.c
    @brief Identify 2D/3D mode for immersed body simulations
    @author Debojyoti Ghosh
*/

#include <string.h>
#include <immersedboundaries.h>

/*!
  Identify the "mode", i.e., whether the simulation is a true 3D simulation, 
  or a 2D simulation being run in 3D. If extent of the immersed body is larger
  than the grid along a particular axis (say, x), then we assume that the 
  intention is to simulate a 2D case around a 2D body in the plane normal
  to that axis ( y-z plane ).

  For example, to simulate a 2D cylinder in the x-y plane, we consider a cylinder
  whose extent along z is larger than the extent of the grid along z, i.e., it
  sticks out of the computational domain at both ends.

  If the immersed body is completely contained within the computational domain, or
  sticks out only on one end along a particular dimension, we assume it's a 3D 
  simulation.
*/
int IBIdentifyMode(
                    double *X,    /*!< Array of global spatial coordinates */
                    int    *dim,  /*!< global dimensions */
                    void   *ib    /*!< Immersed boundary object of type #ImmersedBoundary */
                  )
{
  ImmersedBoundary  *IB = (ImmersedBoundary*) ib;
  Body3D            *body   = IB->body;

  double *x = X, *y = x + dim[0], *z = y + dim[1];

  double grid_xmin = x[0],
         grid_xmax = x[dim[0]-1],
         grid_ymin = y[0],
         grid_ymax = y[dim[1]-1],
         grid_zmin = z[0],
         grid_zmax = z[dim[2]-1];

  if      ( (grid_xmin > body->xmin) && (grid_xmax < body->xmax) ) strcpy(IB->mode,_IB_YZ_);
  else if ( (grid_ymin > body->ymin) && (grid_ymax < body->ymax) ) strcpy(IB->mode,_IB_XZ_);
  else if ( (grid_zmin > body->zmin) && (grid_zmax < body->zmax) ) strcpy(IB->mode,_IB_XY_);
  else                                                             strcpy(IB->mode,_IB_3D_);

  return(0);
}

