/*
                 Joel Heinrich   30 Jun 2004

Bayesian limit calculator:
Returns su such that

             beta == postint(su,n,e0,esig,b0,bsig,alpha);

*/



#include <assert.h>
#include <float.h>
#include "bayesianlimit.h"

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
