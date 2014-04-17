#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*

  Rohde, A., "Eigenvalues and Eigenvectors of the Euler Equations
  in General Geometries", AIAA 2001-2609

*/

int NavierStokes2DLeftEigenvectors(double *u,double *L,void *p,int dir)
{
  NavierStokes2D *param  = (NavierStokes2D*)  p;
  _NavierStokes2DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

int NavierStokes2DRightEigenvectors(double *u,double *R,void *p,int dir)
{
  NavierStokes2D *param  = (NavierStokes2D*)  p;
  _NavierStokes2DRightEigenvectors_(u,R,param,dir);
  return(0);
}
