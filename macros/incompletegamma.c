/*  incomplete gamma function P(a,x)   Joel Heinrich 28 September 2004 */
#include<math.h>
#include<assert.h>
#include<float.h>

#if DBL_MANT_DIG <= 53
#define XCUT 20
#else
#error "appropriate XCUT value unknown for this precision"
#endif

#include "bayesianlimit.h"

double incompletegamma(double a,double x) {
  assert(x>=0); assert(a>=0);
  if(a==0) {
    assert( !(a==0&&x==0) );
    return 1;
  } else if (x==0) {
    return 0;
  } else if (x<XCUT || x<a) {
    const double f=exp(a*log(x)-x-lgamma(a+1));
    double p=1, s=p, den=a;
    if(f<DBL_MIN) return 0;
    while ( p > DBL_EPSILON*s || den<x )
      s += (p*=x/(++den));
    return s*f;
  } else {
    const double rx=-1/x;
    double p=exp((a-1)*log(x)-x-lgamma(a)), s=1-p, num=-a;
    while ( fabs(p) > DBL_EPSILON*s  )
      s -= (p*=(++num)*rx);
    return s;
  }
}
