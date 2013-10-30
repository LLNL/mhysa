#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>

#if !defined(INLINE)
# define INLINE extern inline
#endif

INLINE int    ArrayAXPY              (double*,double,double*,int);

INLINE int    ArrayCopy1D_double     (double*,double*,int);
INLINE int    ArrayCopy1D_int        (int*,int*,int);
INLINE int    ArrayCopynD            (int,double*,double*,int*,int,int,int*,int);

INLINE int    ArrayIncrementIndex    (int,int*,int*);
INLINE int    ArrayIndex1D           (int,int*,int*,int*,int);

INLINE int    ArraySetValue_double   (double*,int,double);
INLINE int    ArraySetValue_int      (int*,int,int);

INLINE int    ArrayAdd1D_double      (double*,double*,double*,int);
INLINE int    ArrayAdd1D_int         (int*,int*,int*,int);
INLINE int    ArraySubtract1D_double (double*,double*,double*,int);
INLINE int    ArraySubtract1D_int    (int*,int*,int*,int);

INLINE double ArraySumSquarenD       (int,int,int*,int,int*,double*);
INLINE double ArraySumAbsnD          (int,int,int*,int,int*,double*);
INLINE double ArrayMaxnD             (int,int,int*,int,int*,double*);

INLINE int ArrayAdd1D_int(int *x,int *a,int *b,int size)
{
  if (!x || !a || !b) return(1);
  int i;
  for (i=0; i<size; i++) x[i] = a[i] + b[i];
  return(0);
}

INLINE int ArrayAdd1D_double(double *x,double *a,double *b,int size)
{
  if (!x || !a || !b) return(1);
  int i;
  for (i=0; i<size; i++) x[i] = a[i] + b[i];
  return(0);
}

INLINE int ArrayAXPY(double *x,double a,double *y,int size)
{
  if (!x || !y) return(1);
  int i;
  for (i=0; i<size; i++) y[i] += a*x[i];
  return(0);
}

INLINE int ArrayCopy1D_double(double *x,double *y,int size)
{
  if (!x || !y) return(1);
  int i;
  for (i = 0; i < size; i++) y[i] = x[i];
  return(0);
}

INLINE int ArrayCopy1D_int(int *x,int *y,int size)
{
  if (!x || !y) return(1);
  int i;
  for (i = 0; i < size; i++) y[i] = x[i];
  return(0);
}

INLINE int ArrayCopynD(int ndims,double *x,double *y,int *dim,int g1,int g2,int *index,int nvars)
{
  int ierr = 0;
  if (!y) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"y\" not allocated.\n");
    return(1);
  }
  if (!x) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"x\" not allocated.\n");
    return(1);
  }
  int done = 0;
  ierr = ArraySetValue_int(index,ndims,0); CHECKERR(ierr);
  while (!done) {
    int p1 = ArrayIndex1D(ndims,dim,index,NULL,g1); 
    int p2 = ArrayIndex1D(ndims,dim,index,NULL,g2);
    int n;
    for (n=0; n<nvars; n++) y[p2*nvars+n] = x[p1*nvars+n];
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(0);
}

INLINE int ArrayIncrementIndex(int N,int *imax,int *i)
{
  int n = 0;
  while (n < N) {
    if (i[n] == imax[n]-1) {
      i[n] = 0;
      n++;
    } else {
      i[n]++;
      break;
    }
  }
  if (n == N) return(1);
  else        return(0);
}

INLINE int ArrayIndex1D(int N,int *imax,int *i,int *offset,int ghost)
{
  int n,index = i[N-1]+ghost+(offset?offset[N-1]:0);
  for (n = N-2; n > -1; n--) index = ((index*(imax[n]+2*ghost)) + (i[n]+ghost+(offset?offset[n]:0)));
  return(index);
}

INLINE double ArrayMaxnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; ArraySetValue_int(index,ndims,0);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    int v; 
    for (v=0; v<nvars; v++) {
      if (absolute(x[p*nvars+v]) > sum) sum = absolute(x[p*nvars+v]);
    }
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(sum);
}

INLINE int ArraySetValue_double(double *x,int size,double value)
{
  int n;
  for (n = 0; n < size; n++)  x[n] = value;
  return(0);
}

INLINE int ArraySetValue_int(int *x,int size,int value)
{
  int n;
  for (n = 0; n < size; n++)  x[n] = value;
  return(0);
}

INLINE int ArraySubtract1D_int(int *x,int *a,int *b,int size)
{
  if (!x || !a || !b) return(1);
  int i;
  for (i=0; i<size; i++) x[i] = a[i] - b[i];
  return(0);
}

INLINE int ArraySubtract1D_double(double *x,double *a,double *b,int size)
{
  if (!x || !a || !b) return(1);
  int i;
  for (i=0; i<size; i++) x[i] = a[i] - b[i];
  return(0);
}

INLINE double ArraySumAbsnD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; ArraySetValue_int(index,ndims,0);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    int v; for (v=0; v<nvars; v++) sum += absolute(x[p*nvars+v]);
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(sum);
}

INLINE double ArraySumSquarenD(int nvars,int ndims,int *dim,int ghosts,int *index,double* x)
{
  double sum = 0;
  int done = 0; ArraySetValue_int(index,ndims,0);
  while (!done) {
    int p = ArrayIndex1D(ndims,dim,index,NULL,ghosts);
    int v; for (v=0; v<nvars; v++) sum += (x[p*nvars+v]*x[p*nvars+v]);
    done = ArrayIncrementIndex(ndims,dim,index);
  }
  return(sum);
}
