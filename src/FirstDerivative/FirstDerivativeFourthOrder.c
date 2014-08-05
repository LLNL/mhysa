#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <firstderivative.h>

#include <mpivars.h>
#include <hypar.h>
typedef MPIVariables  MPIContext;
typedef HyPar         SolverContext;

/* 
  Fourth order central differencing
*/

int FirstDerivativeFourthOrderCentral(double *Df,double *f,int dir,void *s,void *m)
{
  SolverContext *solver = (SolverContext*) s;
  int           i, v;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;


  if ((!Df) || (!f)) {
    fprintf(stderr, "Error in FirstDerivativeFourthOrder(): input arrays not allocated.\n");
    return(1);
  }

  static double one_twelve = 1.0/12.0;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    /* left boundary */
    for (i = -ghosts; i < -ghosts+1; i++) {
      int     qC, qp1, qp2, qp3, qp4;
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
      indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      indexC[dir] = i+3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp3);
      indexC[dir] = i+4; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp4);
      for (v=0; v<nvars; v++)  
        Df[qC*nvars+v] = (-25*f[qC*nvars+v]+48*f[qp1*nvars+v]-36*f[qp2*nvars+v]+16*f[qp3*nvars+v]-3*f[qp4*nvars+v])*one_twelve;
    }
    for (i = -ghosts+1; i < -ghosts+2; i++) {
      int qC, qm1, qp1, qp2, qp3;
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
      indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      indexC[dir] = i+3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp3);
      for (v=0; v<nvars; v++)  
        Df[qC*nvars+v] = (-3*f[qm1*nvars+v]-10*f[qC*nvars+v]+18*f[qp1*nvars+v]-6*f[qp2*nvars+v]+f[qp3*nvars+v])*one_twelve;
    }
    /* interior */
    for (i = -ghosts+2; i < dim[dir]+ghosts-2; i++) {
      int qC, qm1, qm2, qp1, qp2;
      indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
      indexC[dir] = i+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
      for (v=0; v<nvars; v++)  
        Df[qC*nvars+v] = (f[qm2*nvars+v]-8*f[qm1*nvars+v]+8*f[qp1*nvars+v]-f[qp2*nvars+v])*one_twelve;
    }
    /* right boundary */
    for (i = dim[dir]+ghosts-2; i < dim[dir]+ghosts-1; i++) {
      int qC, qm3, qm2, qm1, qp1;
      indexC[dir] = i-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
      indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
      for (v=0; v<nvars; v++)  
        Df[qC*nvars+v] = (-f[qm3*nvars+v]+6*f[qm2*nvars+v]-18*f[qm1*nvars+v]+10*f[qC*nvars+v]+3*f[qp1*nvars+v])*one_twelve;
    }
    for (i = dim[dir]+ghosts-1; i < dim[dir]+ghosts; i++) {
      int qC, qm4, qm3, qm2, qm1;
      indexC[dir] = i-4; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm4);
      indexC[dir] = i-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
      indexC[dir] = i-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
      indexC[dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
      indexC[dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      for (v=0; v<nvars; v++)  
        Df[qC*nvars+v] = (3*f[qm4*nvars+v]-16*f[qm3*nvars+v]+36*f[qm2*nvars+v]-48*f[qm1*nvars+v]+25*f[qC*nvars+v])*one_twelve;
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }
  
  return(0);
}
