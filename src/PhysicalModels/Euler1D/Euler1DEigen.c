#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DLeftEigenvectors(double *u,double *L,void *p,int dir)
{
  Euler1D *param  = (Euler1D*)  p;
  _Euler1DLeftEigenvectors_(u,L,param,dir);
  return(0);
}

int Euler1DRightEigenvectors(double *u,double *R,void *p,int dir)
{
  Euler1D *param  = (Euler1D*)  p;
  _Euler1DRightEigenvectors_(u,R,param,dir);
  return(0);
}
