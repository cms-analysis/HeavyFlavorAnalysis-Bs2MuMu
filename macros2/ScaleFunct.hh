/**
 * Scale function classes
 * Author M. De Mattia - 18/11/2008
 * Author S. Casasso   - 25/10/2012
 * Author E. Migliore  - 25/10/2012
 */

#ifndef SCALEFUNCTIONS_H
#define SCALEFUNCTIONS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TRandom.h"

/**
 * Used to define parameters inside the functions.
 */
struct ParameterSet
{
  ParameterSet() {}
  ParameterSet(const TString & inputName, const double & inputStep, const double & inputMini, const double & inputMaxi) :
    step(inputStep),
    mini(inputMini),
    maxi(inputMaxi)
  {
    std::cout << "setting name = " << inputName << std::endl;
    name = inputName;
  }
  TString name;
  double step, mini, maxi;
};



template <class T>
class scaleFunctBase {
public:
  virtual double scale(const double & pt, const double & eta, const double & phi, const int chg, const T & parScale) const = 0;
  virtual ~scaleFunctBase() = 0;
  virtual int parNum() const { return parNum_; }
protected:
  int parNum_;
  virtual void setPar(double* Start, double* Step, double* Mini, double* Maxi, int* ind,
		      TString* parname, const T & parResol, const std::vector<int> & parResolOrder, const std::vector<ParameterSet> & parSet );
};

//
// Curvature: (linear eta + sinusoidal in phi (both in 5 eta bins)) * global scale 
// ------------------------------------------------------------
template <class T>
class scaleFunct50 : public scaleFunctBase<T> {
public:
  scaleFunct50() { this->parNum_ = 27; }
  virtual double scale(const double & pt, const double & eta, const double & phi, const int chg, const T & parScale) const;
};

/// Service to build the scale functor corresponding to the passed identifier                                                                               
scaleFunctBase<double * > * scaleFunctService( const int identifier );

#endif
