#include <iostream>
#include <iomanip>
//////////////////////////////////////////////////////////
//// Root version of Scp
//// Dominique Mangeol based on S.Bityukov scpf_fast.f code
//// 05-25-2006
//// Calculates poisson probability with TMath::Poisson instead of original dlpois
//// takes few seconds to return a result no matter the background size
//// Scp works for significance up to 6.25 (limited by dgausn and TMath::Poisson precision)
//// if higher will return Max(Scp,Sc12) to avoid weird scale down
////
////

double scpfor(double amus,double amub, double dmb,double dtun);
double  dalunc(double amus,double amub,double dmb,double dtun);
double dalpha(double amus,double amub,double dtun);
double dgausn(double P);
double  gfunk(double zxlam,double x);
double revgam(double z);

double scp(double sig, double bkg, double sigma_b, double delta_b, int verb=0)
{
double fac2  = sqrt(bkg+delta_b)/sqrt(bkg+sigma_b*sigma_b+delta_b);
double Sc12_sys = 2*(sqrt(sig+bkg)-sqrt(bkg+delta_b))*fac2;

double Scp = scpfor(sig,bkg, sigma_b,delta_b);
 if (verb>0) cout << " Scp " << Scp << " Sc12 " << Sc12_sys << endl; 
 if (Scp>6.25) {return TMath::Max(Sc12_sys,Scp);} else {return Scp;}
}
double scpfor(double amus,double amub, double dmb,double dtun)
{
              double dbeta=0;
              double alpha=0.;
	double dsgnf = 0.;
	if(dmb > 0.000000001) {
//            calculation with uncertainties in amub
              alpha = dalunc(amus,amub,dmb,dtun);
   	 dbeta = 1. - alpha;} else {
//             calculation without uncertainties in amub
	 alpha = dalpha(amus,amub,dtun);
  	 dbeta = 1. - alpha;}
 	 dsgnf = dgausn(alpha);
	 return dsgnf;
}

      double dalpha(double amus,double amub,double dtun)
{
//---------------------------------------
//                       S.Bityukov, S.Erofeeva, June, 2005
//---------------------------------------

	int n;

	     double alpha = 0.;
	     n = amus + amub - 1.;
	     double x = amus + amub  - 1.; 
	     double zxlam = amub + dtun;

	     double prob;
	     for (int k=0;k<=n;k++) {
	     prob = TMath::Poisson(k,amub+dtun);
	     alpha = alpha + prob;
	     }
 	     double prob1=  TMath::Poisson(n+1,amub + dtun);
                   double delta2 = gfunk(x,zxlam);
	     double delta = 0.;
	     if (prob1!=prob) {
                 delta = prob1*(prob-delta2)/(prob-prob1);
	     alpha = alpha + delta;}
	 return alpha;
}
//---------------------------------------

      double  dalunc(double amus,double amub,double dmb,double dtun) 
{
//---------------------------------------
//                       S.Bityukov, S.Erofeeva, June, 2005
//---------------------------------------
double pi = 3.141592653589793238;
	double     alpha = 0.;
	     int n = amus + amub - 1;
	     double z = amus + amub  - 1.; 
	     double sum = 0.;
	double d = dmb;
	double b = amub+dtun;
	double xmin = amub+dtun - 3.*d;
	if (xmin<0.) {xmin = 0.;}
	double xmax = amub+dtun + 3.*d;
              double step =0.;
              double dx = (xmax - xmin)*0.01;
              
	while (xmin<xmax) {
	double p = exp(-(xmin-b)*(xmin-b)/(2.*d*d));
	p = p/(d*sqrt(2.*pi));
	sum = sum + p;
	double x = xmin;
              double lastprob=0.;
	     for (int k=0;k<=n;k++) {
	 double prob = TMath::Poisson(k,x);
	 alpha = alpha + prob*p;
               lastprob = prob;
	   }
 	double   prob1= TMath::Poisson(n+1,x);
              double  xz = x;
              double delta2 = gfunk(z,xz);
	double     delta = 0.;
	  if(prob1!=lastprob) {
          delta = prob1*(lastprob-delta2)/(lastprob-prob1);
          alpha = alpha + delta*p;}
	xmin = xmin + dx;}

	alpha = alpha/sum;
	return alpha;
}
//---------------------------------------


	double  gfunk(double zxlam,double x)
{
	double zfact = revgam(zxlam);
	double zprob = zxlam*log(x) - x + zfact;
	double gfunk = exp(zprob);
      return gfunk;
}

	double revgam(double z)
{

//             calculation of  ln(1./x!)   if x<0 then 
//  
	double x = z+1.;
	double zprob = -1.;
	if(x<0.) {return zprob;}
	if(x<3.) {
	zprob = TMath::Gamma(x);
	zprob = -log(zprob);
	return zprob;}

	int i = (int) x;
	double zk = x - (float) i; 
	double zi = TMath::Gamma(zk+2.);
	zprob = -log(zi);

      for (int k = 3;k<=i;k++) {
	zi = zk + float(k-1);
        zprob = zprob - log(zi);}
      return zprob;
}
 double dgausn(double P)
{
//     Computes a "Normal Deviate"
//     Based on G.W. Hill & A.W. Davis, Algorithm 442 Normal Deviate
//     Collected Algorithms from CACM
//     Based on fortran cernlib routine

      double C = 2.50662827463100050;
      double Z1 = 1., HF = Z1/2., C1 = 3.*Z1/4., C2 = 7.*Z1/8., C3 = Z1/3.;
      double h=0.;
      if (P <= 0 || P> 1) {cout << " problem proba not between 0-1 " << endl; h=0.;}
      else if (P==HF) {
       h=0.;}
      else {
       double X=P;
       if (P > HF) {X=1.-P;}
       X=sqrt(-2*log(X));
       X=X-((7.47395*X+494.877)*X+1637.720)/(((X+117.9407)*X+908.401)*X+659.935);
       if(P < HF) {X=-X;}
       double S=X*X;
       double Z=C*(P-TMath::Freq(X))*exp(HF*S);
       h=(((((C1*S+C2)*Z+X)*X+HF)*C3*Z+HF*X)*Z+1)*Z+X;
      }
      return h;
}
