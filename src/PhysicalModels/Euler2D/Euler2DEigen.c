#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler2d.h>
#include <hypar.h>

/*

  Rohde, A., "Eigenvalues and Eigenvectors of the Euler Equations
  in General Geometries", AIAA 2001-2609

*/

int Euler2DLeftEigenvectors(double *u,double *L,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  _Euler2DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

int Euler2DRightEigenvectors(double *u,double *R,void *p,int dir)
{
  Euler2D *param  = (Euler2D*)  p;
  _Euler2DRightEigenvectors_(u,R,param,dir);
  return(0);
}
