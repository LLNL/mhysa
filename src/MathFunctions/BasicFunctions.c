#include <math.h>
#include <mathfunctions.h>

inline double min(double a, double b)
{
  if (a < b)  return(a);
  else        return(b);
}

inline double max(double a, double b)
{
  if (a > b)  return(a);
  else        return(b);
}

inline double min3(double a, double b, double c)
{
  return (min(min(a,b),min(a,c)));
}

inline double max3(double a, double b, double c)
{
  return (max(max(a,b),max(a,c)));
}

inline double absolute(double a)
{
  if (a < 0) return(-a);
  else       return( a);
}

inline double raiseto(double x,double a)
{
  return (exp(a*log(x)));
}

inline double sign(double a)
{
  return (a < 0 ? -1.0 : 1.0);
}
