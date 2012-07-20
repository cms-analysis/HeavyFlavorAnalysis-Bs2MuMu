/*
 *  NCFuntions.cp
 *  NCRootUtils
 *
 *  Copied by Christoph on 4.3.10 from Url Langenegger.
 *
 */

#include "NCFunctions.h"
#include <TMath.h>
#include <cmath>

// ----------------------------------------------------------------------
double f_expo(double *x, double *par) {
  return par[0]*TMath::Exp(-x[0]*par[1]);
}

// ----------------------------------------------------------------------
double f_gauss(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}

// ----------------------------------------------------------------------
double f_Gauss(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma

  double sqrt2pi = 2.506628275;

  if (par[2] > 0.) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t fitval =  (par[0]/(sqrt2pi*par[2])) * TMath::Exp(-0.5*arg*arg);
    return fitval;
  }
  else {
    return -1.;
  }
}

// ----------------------------------------------------------------------
double f_2G(double *x, double *par) {
  // par[0] -> area
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean
  // par[5] -> sigma

  double sqrt2pi = 2.506628275;

  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  (par[0]/(sqrt2pi*par[2]))*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  (par[3]*par[0]/(sqrt2pi*par[2]))*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}

// ----------------------------------------------------------------------
double f_2g(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean
  // par[5] -> sigma
  Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
  }
  Double_t fitval = fitval1 + fitval2;
  return fitval;
}

// ----------------------------------------------------------------------
double f_3g(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> const
  // par[4] -> mean
  // par[5] -> sigma
  // par[6] -> const
  // par[7] -> mean
  // par[8] -> sigma
  Double_t arg1(0.), arg2(0.), arg3(0.), fitval1(0.), fitval2(0.), fitval3(0.); 
  if (par[2] > 0) {
    arg1 = (x[0] - par[1]) / par[2];
    fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
  }
  if (par[5] > 0.) {
    arg2 = (x[0] - par[4]) / par[5];
    fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
  }
  if (par[8] > 0.) {
    arg3 = (x[0] - par[7]) / par[8];
    fitval3 =  par[6]*par[0]*TMath::Exp(-0.5*arg3*arg3);
  }
  Double_t fitval = fitval1 + fitval2 + fitval3;
  return fitval;
}

// double gauss...
double f_double_gauss(double *x, double *par)
{
	return f_gauss(x,par) + f_gauss(x, par+3);
}

// ----------------------------------------------------------------------
double f_cb(double *x, double *par) {
  // par[0]:  mean
  // par[1]:  sigma
  // par[2]:  alpha, crossover point
  // par[3]:  n, length of tail
  // par[4]:  N, normalization

  Double_t cb = 0.0;
  Double_t exponent = 0.0;

  if (x[0] > par[0] - par[2]*par[1]) {
    exponent = (x[0] - par[0])/par[1];
    cb = TMath::Exp(-exponent*exponent/2.);
  } else {
    double nenner  = TMath::Power(par[3]/par[2], par[3])*TMath::Exp(-par[2]*par[2]/2.);
    double zaehler = (par[0] - x[0])/par[1] + par[3]/par[2] - par[2];
    zaehler = TMath::Power(zaehler, par[3]);
    cb = nenner/zaehler;
  }

  if (par[4] > 0.) {
    cb *= par[4];
  }

  return cb;
}



// ----------------------------------------------------------------------
double f_fnov(double *x, double *par) {
  //   par[0] = normalization
  //   par[1] = mean
  //   par[2] = FWHM/2.36
  //   par[3] = tail
  //

  // -- If tail is small then Gauss
  double qa=0,qb=0,qc=0,qx=0,qy=0;
  double result=0;

  if(fabs(par[3]) < 1.e-7) 
    qc = 0.5*pow(((x[0]-par[1])/par[2]),2);
  else {
    qa = par[3]*sqrt(log(4.));
    qb = sinh(qa)/qa;
    qx = (x[0]-par[1])/par[2]*qb;
    qy = 1.+par[3]*qx;
 
    // -- Cutting curve from right side
    if( qy > 1.E-7) 
      qc = 0.5*(pow((log(qy)/par[3]),2) + par[3]*par[3]);
    else
      qc = 15.;
  }
  result =  par[0] * exp(-qc);
  return result;
}


// ----------------------------------------------------------------------
// Argus only
double f_argus(double *x, double *par)
	// par[0] = normalization
	// par[1] = cutoff
	// par[2] = curvature
{
	double result = x[0] / par[1];
	if (result >= 1) return 0;
	
	result = 1 - result*result;
	result = par[0] * x[0] * TMath::Sqrt(result)*exp(par[2]*result);
	
	return result;
} // f_argus()

// ----------------------------------------------------------------------
// Argus and Gauss 
double f_aag(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = normalization of argus
  //   par[4] = exponential factor of argus

//    double ebeam = 10.58/2;
//    double signal = 0.;
//    double background = 0.;
//    double result=0.;
//    if (par[2] > 0.)  signal     = par[0] * exp(-(x[0]-par[1]) * (x[0]-par[1]) / (2*par[2]*par[2]));
//    background = par[3] * x[0] * sqrt(1 - (x[0]*x[0])/(ebeam*ebeam)) * exp(par[4] * (1 - (x[0]*x[0])/(ebeam*ebeam))); 
//    result = signal + background;
//    return result;

  return  (f_argus(x, &par[3]) + f_gauss(x, &par[0]));

}

// exponential + argus + gaus
double f_exparggau(double *x, double *par) {
	// exp:		par[0] - par[1]
	// argus:	par[2] - par[4]
	// gaus:	par[5] - par[7]
	return (f_expo(x, &par[0]) + f_argus(x, &par[2]) + f_gauss(x, &par[5]));
}

// ----------------------------------------------------------------------
double f_aacb(double *x, double *par) {
  //   par[0] = mean of cb
  //   par[1] = sigma of cb
  //   par[2] = alpha
  //   par[3] = n
  //   par[4] = N
  //   par[5] = normalization of argus
  //   par[6] = exponential factor of argus
  return  (f_argus(x, &par[5]) + f_cb(x, &par[0]));
}


// ----------------------------------------------------------------------
double f_p1(double *x, double *par) {
  return par[0] + par[1]*x[0]; 
}

// ----------------------------------------------------------------------
double f_p2(double *x, double *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]; 
}


// ----------------------------------------------------------------------
double f_p1acb(double *x, double *par) {
  //   par[0] = mean of cb
  //   par[1] = sigma of cb
  //   par[2] = alpha
  //   par[3] = n
  //   par[4] = N
  //   par[5] = par 0 of pol1
  //   par[6] = par 1 of pol1
  return  (f_p1(x, &par[5]) + f_cb(x, &par[0]));
}

// ----------------------------------------------------------------------
double f_p2acb(double *x, double *par) {
  //   par[0] = mean of cb
  //   par[1] = sigma of cb
  //   par[2] = alpha
  //   par[3] = n
  //   par[4] = N
  //   par[5] = par 0 of pol2
  //   par[6] = par 1 of pol2
  //   par[7] = par 2 of pol2
  return  (f_p2(x, &par[5]) + f_cb(x, &par[0]));
}


// ----------------------------------------------------------------------
// pol1 and Gauss 
double f_p1ag(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = par 0 of pol1
  //   par[4] = par 1 of pol1
  return  (f_p1(x, &par[3]) + f_gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// pol2 and Gauss 
double f_p2ag(double *x, double *par) {
  //   par[0] = normalization of gaussian
  //   par[1] = mean of gaussian
  //   par[2] = sigma of gaussian
  //   par[3] = par 0 of pol2
  //   par[4] = par 1 of pol2
  //   par[5] = par 2 of pol2
  return  (f_p2(x, &par[3]) + f_gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// pol0 and Gauss 
double f_p0ag(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] = par 0 of pol0

  return  (par[3] + f_gauss(x, &par[0]));
}

// ----------------------------------------------------------------------
// pol0 and double Gauss 
double f_p0a2g(double *x, double *par) {
  // par[0] -> const
  // par[1] -> mean
  // par[2] -> sigma
  // par[3] -> fraction in second gaussian
  // par[4] -> mean
  // par[5] -> sigma
  // par[6] = par 0 of pol0

  return  (par[6] + f_gauss(x, &par[0]) + f_gauss(x, &par[3]));
}

// ----------------------------------------------------------------------
// chi2 ellipsis
double f2_chi2ellipsis(double *x, double *par) {
   // x[0]   x
   // x[1]   y
   // par[0] x0
   // par[1] sigma(x)
   // par[2] y0 
   // par[3] sigma(y)
   // par[4] correlation
   // par[5] chi2

  double result = (x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]
    + (x[1]-par[2])*(x[1]-par[2])/par[3]/par[3]
    - 2*par[4]*(x[0]-par[0])*(x[1]-par[2])/par[1]/par[3];

  if (result < par[5]) {
    return result; 
  } else {
    return 0.;
  }
}

/* Computes the sum of a Gaussian with an exponential function
 *	par[0]..par[2] arguments to the function f_gauss
 *	par[3]..par[4] arguments to the function f_expo
 */
double f_gauss_expo(double *x, double *par)
{	
	return f_gauss(x, par) + f_expo(x, par+3);
} // f_gauss_expo()

/* Computes the sum of a Gaussian with a linear function
 *  par[0]..par[2] arguments to the function f_gauss
 *  par[3]..par[4] arguments to the function f_p1
 */
double f_gauss_linear(double *x, double *par)
{
	return f_gauss(x, par) + f_p1(x, par+3);
} // f_gauss_linear()

/* Double gauss with linear background
 *	par[0]..par[5] arguments of function f_double_gauss
 *	par[6]..par[7] arguments to function f_p1
 */
double f_double_gauss_linear(double *x, double *par)
{
	return f_double_gauss(x,par) + f_p1(x, par+6);
} // f_double_gauss_linear()

/* Boltzmann-like function
 *	f(x) = [0]*x^[1]*e^(-[2]x)
 */
double f_boltzmann(double *x, double *par)
{
	double exponent = par[1]*log(x[0])-par[2]*x[0];
	return par[0]*exp(exponent);
} // f_boltzmann()

/* Characteristic function.
 *	f(x) = 1_[[0],[1]](x) */
double f_charact(double *x, double *par)
{
	double result;
	if (par[0] < x[0] && x[0] < par[1])
		result = 1.0;
	else
		result = 0.0;
	
	return result;
} // f_charact()

/* Skew Normal distribution (c.f. Wikipedia.org)
 *	par[0] = xi		location parameter
 *	par[1] = omega	scale parameter
 *	par[2] = alpha	skewness parameter
 */
double f_skewnormal(double *x, double *par)
{
	double result = 2 / par[1];
	double parm = (x[0] - par[0])/par[1];
	
	result *= TMath::Gaus(parm, 0, 1, kTRUE);
	result *= (1 + TMath::Erf(par[2]*parm/sqrt(2))) / 2;
	
	return result;
} // f_skewnormal()
