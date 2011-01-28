/* Power series evaluation of incomplete beta function.  */
/* Joel Heinrich   1 July 2004                           */

#include<math.h>
#include<float.h>
#include<assert.h>
#include "bayesianlimit.hh"

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


/*  incomplete gamma function P(a,x)   Joel Heinrich 28 September 2004 */

#if DBL_MANT_DIG <= 53
#define XCUT 20
#else
#error "appropriate XCUT value unknown for this precision"
#endif

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


/*            Joel Heinrich    27 Oct 2004

double posterior(double s,int n,double e0,double esig,
                                       double b0,double bsig,double alpha);

   returns posterior pdf evaluated at s, where

n           events are observed
e0 +/- esig is the estimated efficiency
b0 +/- bsig is the estimated background
s^(alpha-1) is the prior p.d.f. for s                         */


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


/*
                 Joel Heinrich   27 Oct 2004

Integrates the Bayesian posterior p.d.f. from 0 to su where:

n           events are observed
e0 +/- esig is the estimated efficiency
b0 +/- bsig is the background
s^(alpha-1) is the prior p.d.f. for s

*/


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



/*
                 Joel Heinrich   30 Jun 2004

Bayesian limit calculator:
Returns su such that

             beta == postint(su,n,e0,esig,b0,bsig,alpha);

*/

double blimit(double beta,int n,double e0,double esig,
	      double b0,double bsig,double alpha) {
  const double kappa = (esig!=0) ? e0/(esig*esig) : 1;
  double x0=0, x1=1, dx=1, dxp=2;
  assert(beta>=0);assert(beta<1);
  if(beta==0) return 0;
  while( (dx=x1-x0)>DBL_EPSILON && dx<dxp ) {
    const double xmid = 0.5*(x0+x1);
    const double ymid = postint(xmid*kappa/(1-xmid),n,e0,esig,b0,bsig,alpha);
    dxp = dx;
    (ymid<beta) ? (x0=xmid) : (x1=xmid);
  }
  return x0*kappa/(1-x0);
}


/*
                 Joel Heinrich   30 Jun 2004

Bayesian limit calculator:
Returns su such that

             beta == postintb0(su,n,e0,esig,b,alpha);

*/

double blimitb0(double beta,int n,double e0,double esig,double b,double alpha){
  const double kappa = (esig!=0) ? e0/(esig*esig) : 1;
  double x0=0, x1=1, dx=1, dxp=2;
  assert(beta>=0);assert(beta<1);
  if(beta==0) return 0;
  while( (dx=x1-x0)>DBL_EPSILON && dx<dxp ) {
    const double xmid = 0.5*(x0+x1);
    const double ymid = postintb0(xmid*kappa/(1-xmid),n,e0,esig,b,alpha);
    dxp = dx;
    (ymid<beta) ? (x0=xmid) : (x1=xmid);
  }
  return x0*kappa/(1-x0);
}


/*
                 Joel Heinrich   30 Jun 2004

Bayesian limit calculator:
Returns su such that

             beta == postinte0b0(su,n,e0,b,alpha);

*/


double blimite0b0(double beta,int n,double e,double b,double alpha) {
  double x0=0, x1=1, dx=1, dxp=2;
  assert(beta>=0);assert(beta<1);
  if(beta==0) return 0;
  while( (dx=x1-x0)>DBL_EPSILON && dx<dxp ) {
    const double xmid = 0.5*(x0+x1);
    const double ymid = postinte0b0(xmid/(1-xmid),n,e,b,alpha);
    dxp = dx;
    (ymid<beta) ? (x0=xmid) : (x1=xmid);
  }
  return x0/(1-x0);
}
