/* Power series evaluation of incomplete beta function.  */
/* Joel Heinrich   1 July 2004                           */

#include<math.h>
#include<float.h>
#include<assert.h>
#include "bayesianlimit.h"

static double ib(double a,double b,double x) {
  double n=a+b, d=a+1, prod=1, sum=1;
  while (prod>sum*DBL_EPSILON)
    sum += (prod *= n++*x/d++);
  return exp(a*log(x)+b*log(1-x)+lgamma(a+b)-lgamma(a+1)-lgamma(b))*sum;
}

double incompletebeta(double a,double b,double x) {
  assert(x>=0);assert(x<=1);assert(a>0);assert(b>0);
  return (x==0||x==1) ? x : (x*(a+b+2)<a+1) ? ib(a,b,x) : 1-ib(b,a,1-x);
}
