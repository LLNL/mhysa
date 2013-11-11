#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/fppowersystem.h>
#include <hypar.h>

double FPPowerSystemDriftFunction(int,void*,double,double,double);

int FPPowerSystemUpwind(double *fI,double *fL,double *fR,double *uL,double *uR,
                        double *u,int dir,void *s,double t)
{
  HyPar         *solver = (HyPar*) s;
  FPPowerSystem *params = (FPPowerSystem*)solver->physics;
  int           ierr    = 0,done,v;

  int ndims   = solver->ndims;
  int nvars   = solver->nvars;
  int ghosts  = solver->ghosts;
  int *dim    = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  ierr = ArrayCopy1D_int(dim,bounds_outer,ndims); CHECKERR(ierr); bounds_outer[dir] =  1;
  ierr = ArrayCopy1D_int(dim,bounds_inter,ndims); CHECKERR(ierr); bounds_inter[dir] += 1;

  done = 0; ierr = ArraySetValue_int(index_outer,ndims,0); CHECKERR(ierr);
  while (!done) {
    ierr = ArrayCopy1D_int(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p = ArrayIndex1D(ndims,bounds_inter,index_inter,NULL,0);
      double x = 0,y = 0; /* x,y coordinates of the interface */
      if (dir == 0) {
        /* x-interface */
        x = 0.5 * ( solver->GetCoordinate(0,index_inter[0]-1,dim,ghosts,solver->x) 
                  + solver->GetCoordinate(0,index_inter[0]  ,dim,ghosts,solver->x) );
        y = solver->GetCoordinate(1,index_inter[1],dim,ghosts,solver->x);
      } else if (dir == 1) {
        /* y-interface */
        x = solver->GetCoordinate(0,index_inter[0],dim,ghosts,solver->x);
        y = 0.5 * ( solver->GetCoordinate(1,index_inter[1]-1,dim,ghosts,solver->x) 
                  + solver->GetCoordinate(1,index_inter[1]  ,dim,ghosts,solver->x) );
      }
      double drift = FPPowerSystemDriftFunction(dir,params,x,y,t);
      for (v = 0; v < nvars; v++)  
        fI[nvars*p+v] = drift * (drift > 0 ? fL[nvars*p+v] : fR[nvars*p+v] );
    }
    done = ArrayIncrementIndex(ndims,bounds_outer,index_outer);
  }

  return(0);
}
