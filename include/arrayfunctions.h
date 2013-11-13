#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>

#define _ArrayIndex1D_(N,imax,i,ghost,index)                                                                 \
  {                                                                                                                 \
    index = i[N-1]+(ghost);                                                                                         \
    int arraycounter;                                                                                               \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) {                                                 \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost)));                                 \
    }                                                                                                               \
  }

#define _ArrayIndex1DWO_(N,imax,i,offset,ghost,index)                                                               \
  {                                                                                                                 \
    index = i[N-1]+(ghost)+ offset[N-1];                                                                            \
    int arraycounter;                                                                                               \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) {                                                 \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost)+offset[arraycounter]));            \
    }                                                                                                               \
  }

#define _ArrayIncrementIndex_(N,imax,i,done)                                                                        \
  {                                                                                                                 \
    int arraycounter = 0;                                                                                           \
    while (arraycounter < (N)) {                                                                                    \
      if (i[arraycounter] == imax[arraycounter]-1) {                                                                \
        i[arraycounter] = 0;                                                                                        \
        arraycounter++;                                                                                             \
      } else {                                                                                                      \
        i[arraycounter]++;                                                                                          \
        break;                                                                                                      \
      }                                                                                                             \
    }                                                                                                               \
    if (arraycounter == (N)) done = 1;                                                                              \
    else          done = 0;                                                                                         \
  }

#define _ArraySetValue_(x,size,value)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);                       \
  }

#define _ArraySubtract1D_(x,a,b,size)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] - b[arraycounter];    \
  }

#define _ArrayAdd1D_(x,a,b,size)                                                                                    \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) x[arraycounter] = a[arraycounter] + b[arraycounter];    \
  }

#define _ArrayAXPY_(x,a,y,size)                                                                                     \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter=0; arraycounter<size; arraycounter++) y[arraycounter] += a*x[arraycounter];                   \
  }

#define _ArrayCopy1D_(x,y,size)                                                                                     \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = x[arraycounter];                  \
  }

#if !defined(INLINE)
# define INLINE inline
#endif

INLINE int    ArrayCopynD            (int,double*,double*,int*,int,int,int*,int);
INLINE double ArrayMaxnD             (int,int,int*,int,int*,double*);
INLINE double ArraySumSquarenD       (int,int,int*,int,int*,double*);
INLINE double ArraySumAbsnD          (int,int,int*,int,int*,double*);

INLINE int ArrayCopynD(int ndims,double *x,double *y,int *dim,int g1,int g2,int *index,int nvars)
{
  if (!y) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"y\" not allocated.\n");
    return(1);
  }
  if (!x) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"x\" not allocated.\n");
    return(1);
  }
  int done = 0;
  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p1, p2, n;
    _ArrayIndex1D_(ndims,dim,index,g1,p1); 
    _ArrayIndex1D_(ndims,dim,index,g2,p2);
    for (n=0; n<nvars; n++) y[p2*nvars+n] = x[p1*nvars+n];
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(0);
}


INLINE double ArrayMaxnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) if (absolute(x[p*nvars+v]) > sum) sum = absolute(x[p*nvars+v]);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

INLINE double ArraySumAbsnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += absolute(x[p*nvars+v]);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

INLINE double ArraySumSquarenD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += (x[p*nvars+v]*x[p*nvars+v]);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}
