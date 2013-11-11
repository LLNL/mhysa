#include <math.h>
#include <mathfunctions.h>

double min(double a, double b)
{
  if (a < b)  return(a);
  else        return(b);
}

double max(double a, double b)
{
  if (a > b)  return(a);
  else        return(b);
}

double min3(double a, double b, double c)
{
  return (min(min(a,b),min(a,c)));
}

double max3(double a, double b, double c)
{
  return (max(max(a,b),max(a,c)));
}

double absolute(double a)
{
  if (a < 0) return(-a);
  else       return( a);
}

double raiseto(double x,double a)
{
  return (exp(a*log(x)));
}

double sign(double a)
{
  return (a < 0 ? -1.0 : 1.0);
}
