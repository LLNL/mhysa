#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

/*

  Rohde, A., "Eigenvalues and Eigenvectors of the Euler Equations
  in General Geometries", AIAA 2001-2609

*/

int NavierStokes3DLeftEigenvectors(double *u,double *L,void *p,int dir)
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  _NavierStokes3DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

int NavierStokes3DRightEigenvectors(double *u,double *R,void *p,int dir)
{
  NavierStokes3D *param  = (NavierStokes3D*)  p;
  _NavierStokes3DRightEigenvectors_(u,R,param,dir);
  return(0);
}
