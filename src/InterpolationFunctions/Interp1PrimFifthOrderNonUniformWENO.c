#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/* 
  Fifth order WENO interpolation (non-uniform grid)
*/

#undef  _MINIMUM_GHOSTS_
#define _MINIMUM_GHOSTS_ 3

static int Interp1PrimFifthOrderNonUniformWENO_N(double*,double*,double*,double*,int,int,void*,void*);

int Interp1PrimFifthOrderNonUniformWENO(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*) s;
  return(Interp1PrimFifthOrderNonUniformWENO_N(fI,fC,u,x,upw,dir,s,m));
}

int Interp1PrimFifthOrderNonUniformWENO_N(double *fI,double *fC,double *u,double *x,int upw,int dir,void *s,void *m)
{
  HyPar           *solver = (HyPar*)          s;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;
  int             i;

  if (!weno->var) weno->var = fC;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* create an array of interface x values */
  double *x_int = (double*) calloc (dim[dir]+1+2*ghosts,sizeof(double));
  for (i = 1; i < dim[dir]+2*ghosts; i++) x_int[i] = 0.5 * (x[i-1]+x[i]);
  x_int[0] = 1.5*x[0] - 0.5 * x[1];
  x_int[dim[dir]+2*ghosts] = 1.5*x[dim[dir]+2*ghosts-1] - 0.5*x[dim[dir]+2*ghosts-2];

  /* define some constants */
  static const double one_sixth          = 1.0/6.0;
  static const double thirteen_by_twelve = 13.0/12.0;
  static const double one_fourth         = 1.0/4.0;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;

  int done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,ndims);
    _ArrayCopy1D_(index_outer,indexI,ndims);
    for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
      int qm1,qm2,qm3,qp1,qp2,p;
      int im1,im2,im3,ip1,ip2,ip3;
      _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
      if (upw > 0) {
        indexC[dir] = im3 = indexI[dir]-3; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = im2 = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = im1 = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = ip1 = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = ip2 = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
                      ip3 = indexI[dir]+2;
      } else {
        indexC[dir] = im3 = indexI[dir]+2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm3);
        indexC[dir] = im2 = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = im1 = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = ip1 = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = ip2 = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
                      ip3 = indexI[dir]-3;
      }
      int v; 
      for (v=0; v<nvars; v++)  {
        /* Defining stencil points */
        double fm3, fm2, fm1, fp1, fp2;
        double  m3,  m2,  m1,  p1,  p2;
        double xm3, xm2, xm1, xp1, xp2, xp3;
        fm3 = fC[qm3*nvars+v]; xm3 = *(x_int+ghosts+im3);
        fm2 = fC[qm2*nvars+v]; xm2 = *(x_int+ghosts+im2);
        fm1 = fC[qm1*nvars+v]; xm1 = *(x_int+ghosts+im1);
        fp1 = fC[qp1*nvars+v]; xp1 = *(x_int+ghosts+ip1);
        fp2 = fC[qp2*nvars+v]; xp2 = *(x_int+ghosts+ip2);
                               xp3 = *(x_int+ghosts+ip3);
        m3 = weno->var[qm3*nvars+v];
        m2 = weno->var[qm2*nvars+v];
        m1 = weno->var[qm1*nvars+v];
        p1 = weno->var[qp1*nvars+v];
        p2 = weno->var[qp2*nvars+v];

        /* 3rd order reconstruction values */
        double f1, f2, f3, c1, c2, c3;

        c1 = (xp2-xp1)*(xp3-xp1) / ( (xp3-xm1)*(xp2-xm1) );
        c2 = (xp2-xp1)*(xp1-xm1) / ( (xp3-xm1)*(xp3-xp1) );
        f1 = fp1 + c1 * (fm1-fp1) - c2 * (fp2-fp1);

        c1 = (xp1-xm1)*(xp1-xm2) / ( (xp2-xm2)*(xp2-xm1) );
        c2 = (xp1-xm1)*(xp2-xp1) / ( (xp2-xm2)*(xp1-xm2) );
        f2 = fm1 + c1 * (fp1-fm1) - c2 * (fm2-fm1);

        c1 = (xp1-xm1)*(xp1-xm2) / ( (xm1-xm3)*(xp1-xm3) );
        c2 = 1.0 + (xp1-xm1)/(xp1-xm2) + ((xp1-xm1)/(xp1-xm3));
        f3 = fm2 + c1 * (fm3-fm2) + c2 * (fm1-fm2);

        /* Smoothness indicators */
        double fac1, fac2, b1, b2, b3;

        fac1 = m1 - 2.0*p1 + p2;
        fac2 = 3.0*m1 - 4.0*p1 + p2;
        b1 = thirteen_by_twelve*fac1*fac1 + one_fourth*fac2*fac2;

        fac1 = m2 - 2.0*m1 + p1;
        fac2 = m2 - p1;
        b2 = thirteen_by_twelve*fac1*fac1 + one_fourth*fac2*fac2;

        fac1 = m3 - 2.0*m2 + m1;
        fac2 = 3.0*m1 - 4.0*m2 + m3;
        b3 = thirteen_by_twelve*fac1*fac1 + one_fourth*fac2*fac2;

        /* optimal weights */
        c1 = (xp1-xm3)*(xp1-xm2) / ( (xp3-xm3)*(xp3-xm2) ); 
        c2 = (xp1-xm3)*(xp3-xp1)*((xp3-xm2)/(xp2-xm3)+1) / ( (xp3-xm3)*(xp3-xm2) ); 
        c3 = (xp2-xp1)*(xp3-xp1) / ( (xp3-xm3)*(xp2-xm3) );

        /* calculate WENO weights */
        double w1,w2,w3;

        if (weno->no_limiting) {
          /* fifth order polynomial interpolation */
          w1 = c1;
          w2 = c2;
          w3 = c3;
        } else {
  
          /* This parameter is needed for Borges' or Yamaleev-Carpenter
             implementation of non-linear weights  */
          double tau;
          if (weno->borges)   tau = absolute(b3 - b1);
          else if (weno->yc)  tau = (m3-4*m2+6*m1-4*p1+p2)*(m3-4*m2+6*m1-4*p1+p2);
          else                tau = 0;
  
          /* Defining the non-linear weights */
          double a1, a2, a3;
          if (weno->borges || weno->yc) {
            a1 = c1 * (1.0 + raiseto(tau/(b1+weno->eps),weno->p));
            a2 = c2 * (1.0 + raiseto(tau/(b2+weno->eps),weno->p));
            a3 = c3 * (1.0 + raiseto(tau/(b3+weno->eps),weno->p));
          } else {
            a1 = c1 / raiseto(b1+weno->eps,weno->p);
            a2 = c2 / raiseto(b2+weno->eps,weno->p);
            a3 = c3 / raiseto(b3+weno->eps,weno->p);
          }

          /* Convexity */
          double a_sum_inv;
          a_sum_inv = 1.0 / (a1 + a2 + a3);
          w1 = a1 * a_sum_inv;
          w2 = a2 * a_sum_inv;
          w3 = a3 * a_sum_inv;
  
          /* Mapping the weights, if needed */
          if (weno->mapped) {
            a1 = w1 * (c1 + c1*c1 - 3*c1*w1 + w1*w1) / (c1*c1 + w1*(1.0-2.0*c1));
            a2 = w2 * (c2 + c2*c2 - 3*c2*w2 + w2*w2) / (c2*c2 + w2*(1.0-2.0*c2));
            a3 = w3 * (c3 + c3*c3 - 3*c3*w3 + w3*w3) / (c3*c3 + w3*(1.0-2.0*c3));
            a_sum_inv = 1.0 / (a1 + a2 + a3);
            w1 = a1 * a_sum_inv;
            w2 = a2 * a_sum_inv;
            w3 = a3 * a_sum_inv;
          }
  
        }
        fI[p*nvars+v] = w1*f1 + w2*f2 + w3*f3;
      }
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  free(x_int);
  return(0);
}
