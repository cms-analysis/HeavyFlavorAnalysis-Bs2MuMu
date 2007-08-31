/*
                 Joel Heinrich   27 Oct 2004

Integrates the Bayesian posterior p.d.f. from 0 to su where:

n           events are observed
e0 +/- esig is the estimated efficiency
b0 +/- bsig is the background
s^(alpha-1) is the prior p.d.f. for s

*/


#include<math.h>
#include<assert.h>
#include "bayesianlimit.h"

double postint(double su,int n,double e0,double esig,
		double b0,double bsig,double alpha) {
  if(bsig==0) {
    return postintb0(su,n,e0,esig,b0,alpha);
  } else {
    const double kappa = ( assert(e0>0) , assert(esig>0) , e0/(esig*esig) );
    const double omega = ( assert(b0>0) , assert(bsig>0) , b0/(bsig*bsig) );
    const double mu=e0*kappa, x = ( assert(su>=0) , su/(su+kappa) );
    double ix = ( assert(alpha+n>0) , assert(mu>alpha) ,
		  incompletebeta(alpha+n,mu-alpha,x) );
    double sum1=ix, sum2=1, rho=b0*omega;
    if(x==0)return 0;
    if(x==1)return 1;
    if(n>0) {
      double prod1 = exp(lgamma(mu+n) + (alpha+n)*log(x) +
			 (mu-alpha)*log(1-x) -
			 lgamma(alpha+n+1) - lgamma(mu-alpha));
      double prod2=1;
      int k; 
      assert(ix>0); assert(alpha>0);
      for(k=1;k<=n;++k) {
	ix += (prod1 *= (alpha+1+n-k)/(x*(mu+n-k)));
/* Jun 30 2004 version had "omega" instead of "(omega+1)" in the next line */
	prod2 *= (n+1-k)*(rho++)/((omega+1)*k*(alpha+n-k));
	sum1 += ix*prod2;
	sum2 += prod2;
      }
    }
    return sum1/sum2;
  }
}


/*
                 Joel Heinrich   30 Jun 2004

Integrates the Bayesian posterior p.d.f. from 0 to su where:

n           events are observed
e0 +/- esig is the estimated efficiency
b           is the background (known exactly)
s^(alpha-1) is the prior p.d.f. for s

*/


double postintb0(double su,int n,double e0,double esig,double b,double alpha) {
  if(esig==0) {
    return postinte0b0(su,n,e0,b,alpha);
  } else {
    const double kappa = ( assert(e0>0) , assert(esig>0) , e0/(esig*esig) );
    const double mu=e0*kappa, x = ( assert(su>=0) , su/(su+kappa) );
    double ix = ( assert(alpha+n>0) , assert(mu>alpha) ,
		  incompletebeta(alpha+n,mu-alpha,x) );
    double sum1=ix, sum2=1;
    if(x==0)return 0;
    if(x==1)return 1;
    if(b>0 && n>0) {
      double prod1 = exp(lgamma(mu+n) + (alpha+n)*log(x) +
			 (mu-alpha)*log(1-x) -
			 lgamma(alpha+n+1) - lgamma(mu-alpha));
      double prod2=1;
      int k; 
      assert(ix>0); assert(alpha>0);
      for(k=1;k<=n;++k) {
	ix += (prod1 *= (alpha+1+n-k)/(x*(mu+n-k)));
	prod2 *= (n+1-k)*b/(k*(alpha+n-k));
	sum1 += ix*prod2;
	sum2 += prod2;
      }
    }
    return sum1/sum2;
  }
}


/*
                 Joel Heinrich   30 Jun 2004

Integrates the Bayesian posterior p.d.f. from 0 to su where:

n           events are observed
e           is the efficiency (known exacyly)
b           is the background (known exactly)
s^(alpha-1) is the prior p.d.f. for s

*/

double postinte0b0(double su,int n,double e,double b,double alpha) {
  const double es=e*su, res=1/es;
  double Px = ( assert(alpha+n>0) , incompletegamma(alpha+n,es) );
  double sum1=Px, sum2=1;
  if(b>0 && n>0) {
    double prod1 = exp( (alpha+n)*log(es) - es - lgamma(alpha+n+1) );
    double prod2=1;
    int k; 
    assert(Px>0); assert(alpha>0);
    for(k=1;k<=n;++k) {
      Px += (prod1 *= (alpha+1+n-k)*res);
      prod2 *= (n+1-k)*b/(k*(alpha+n-k));
      sum1 += Px*prod2;
      sum2 += prod2;
    }
  }
  return sum1/sum2;
}
