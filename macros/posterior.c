/*            Joel Heinrich    27 Oct 2004

double posterior(double s,int n,double e0,double esig,
                                       double b0,double bsig,double alpha);

   returns posterior pdf evaluated at s, where

n           events are observed
e0 +/- esig is the estimated efficiency
b0 +/- bsig is the estimated background
s^(alpha-1) is the prior p.d.f. for s                         */




#include<math.h>
#include<float.h>
#include<assert.h>

#include "bayesianlimit.h"

static double pesig0(double s,int n,double e0,
		     double b0,double bsig,double alpha);

double posterior(double s,int n,double e0,double esig,
		  double b0,double bsig,double alpha) {
  assert(s>=0);assert(e0>0);assert(esig>=0);assert(b0>=0);assert(bsig>=0);
  if(bsig==0)return posteriorb0(s,n,e0,esig,b0,alpha);
  if(esig==0) {
    return pesig0(s,n,e0,b0,bsig,alpha);
  } else {
    const double kappa = e0/(esig*esig), mu=e0*kappa;
    const double omega = b0/(bsig*bsig), rho = b0*omega;
/* fixed 7 July version by substituting "omega" -> "omega+1" below */
    const double t=s*(omega+1)/(s+kappa);
    double aa=n, a2=rho, bb=n+alpha-1, prod=1, sum=prod;
    int k;
    if(alpha<=0) return (s==0)?DBL_MAX:0;
    for(k=1;k<=n;++k)
      sum += ( prod *= (aa--)*(a2++)/((bb--)*k*(omega+1)) );
    if(s==0) return (alpha<1)?DBL_MAX:(alpha>1)?0:
      (mu-1)*exp(lgamma(rho+n)-lgamma(rho)-n*log(omega+1)-lgamma(1+n))/
	       (kappa*sum);
    aa=n; a2=rho; bb=n+mu-1; prod=1/sum; sum=prod;
    for(k=1;k<=n;++k)
      sum += ( prod *= (aa--)*(a2++)/((bb--)*k*t) );
    return sum*exp(lgamma(mu+n) + (alpha+n-1)*log(s) + (mu-alpha)*log(kappa) -
		   lgamma(mu-alpha) - lgamma(alpha+n) - (mu+n)*log(s+kappa));
  }
}

static double pesig0(double s,int n,double e,
		     double b0,double bsig,double alpha) {
  const double omega = b0/(bsig*bsig), rho = b0*omega;
/* fixed 7 July version by substituting "omega" -> "omega+1" below */
  const double t=e*s*(omega+1);
  double aa=n, a2=alpha, bb=n+rho-1, prod=1, sum=prod;
  int k;
  for(k=1;k<=n;++k)
     sum += ( prod *= (aa--)*(a2++)*(omega+1)/((bb--)*k) );
  if(s==0) return (alpha<1)?DBL_MAX:(alpha>1)?0:e/sum;
  aa=n; bb=n+rho-1; prod=1/sum; sum=prod;
  for(k=1;k<=n;++k)
    sum += ( prod *= (aa--)*t/((bb--)*k) );
  return sum*exp(alpha*log(e)+(alpha-1)*log(s)-e*s-lgamma(alpha));
}

/*            Joel Heinrich    6 July 2004

double posteriorb0(double s,int n,double e0,double esig,double b,double alpha);

   returns posterior pdf evaluated at s, where

n           events are observed
e0 +/- esig is the estimated efficiency
b           is the background (known exactly)
s^(alpha-1) is the prior p.d.f. for s                         */

static double pb0(double s,int n,double mu,double kappa,double alpha);
static double pbsig0(double s,int n,double e0,double b,double alpha);

double posteriorb0(double s,int n,double e0,double esig,double b,double alpha){
  assert(s>=0);assert(e0>0);assert(esig>=0);assert(b>=0);
  if(esig==0) {
    return pbsig0(s,n,e0,b,alpha);
  } else {
    const double kappa = e0/(esig*esig), mu=e0*kappa, t=s/(s+kappa);
    double aa=n, bb=n+alpha-1, prod=1, sum=prod;
    int k;
    if(b==0)return pb0(s,n,mu,kappa,alpha);
    if(alpha<=0) return (s==0)?DBL_MAX:0;
    for(k=1;k<=n;++k)
      sum += ( prod *= (aa--)*b/((bb--)*k) );
    if(s==0) return (alpha<1)?DBL_MAX:(alpha>1)?0:
      (mu-1)*exp(n*log(b)-lgamma(1+n))/(kappa*sum);
    aa=n; bb=n+mu-1; prod=1/sum; sum=prod;
    for(k=1;k<=n;++k)
      sum += ( prod *= (aa--)*b/((bb--)*k*t) );
    return sum*pb0(s,n,mu,kappa,alpha);
  }
}

static double pb0(double s,int n,double mu,double kappa,double alpha){
  if(s==0) {
    return (alpha+n>1)?0:(alpha+n<1)?DBL_MAX:(mu+n-1)/kappa;
  } else if(mu-alpha<=0||alpha+n<=0) {
    return 0;
  }
  return exp(lgamma(mu+n) + (alpha+n-1)*log(s) + (mu-alpha)*log(kappa) -
	     lgamma(mu-alpha) - lgamma(alpha+n) - (mu+n)*log(s+kappa));
}

static double pbsig0(double s,int n,double e,double b,double alpha) {
  if(b==0) {
    if(s==0) {
      return (alpha+n>1)?0:(alpha+n<1)?DBL_MAX:pow(e,alpha+n);
    } else if(alpha+n<=0) {
      return 0;
    }
    return exp((alpha+n)*log(e) - e*s + (alpha+n-1)*log(s) -
	       lgamma(alpha+n));
  } else {
    double aa=n, bb=n+alpha-1, prod=1, sum=prod;
    int k;
    for(k=1;k<=n;++k)
      sum += ( prod *= (aa--)*b/((bb--)*k) );
    return exp(alpha*log(e) - e*s + n*log(e*s+b) + (alpha-1)*log(s) -
	       lgamma(alpha+n))/sum;
  }
}
