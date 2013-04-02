/**
 * Smear function classes
 * Author M. De Mattia - 18/11/2008
 * Author S. Casasso   - 25/10/2012
 * Author E. Migliore  - 25/10/2012
 */


#ifndef RESOLUTIONFUNCTION_H
#define RESOLUTIONFUNCTION_H

#include <iostream>
#include <vector>
#include <cmath>
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TRandom.h"


template <class T>
class resolFunctBase {
 public:
  virtual double sigmaPt(const double & pt, const double & eta, const T & parval) = 0;

  resolFunctBase() {}
  virtual ~resolFunctBase() = 0;
  virtual int parNum() const { return parNum_; }

 protected:
  int parNum_;
};

template <class T>
class resolFunct45 : public resolFunctBase<T> {
 public:
  resolFunct45() { this->parNum_ = 13; }

  virtual double sigmaPt(const double & pt, const double & eta, const T & parval);
};


/// Service to build the scale functor corresponding to the passed identifier
resolFunctBase<double * > * resolutionFunctService( const int identifier );


#endif


