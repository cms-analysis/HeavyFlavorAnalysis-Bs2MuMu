#include "ResolFunct.hh"


template <class T> inline resolFunctBase<T>::~resolFunctBase() { }  // defined even though it's pure virtual; should be faster this way.

template <class T>
double resolFunct45<T>::sigmaPt(const double & pt, const double & eta, const T & parval) {
  
  if( eta < -2.0 ) return( parval[0]*pt + parval[1] );
  if( eta < -1.8 ) return( parval[0]*pt + parval[2] );
  if( eta < -1.6 ) return( parval[0]*pt + parval[3] );
  if( eta < -1.2 ) return( parval[0]*pt + parval[4] );
  if( eta < -0.8 ) return( parval[0]*pt + parval[5] );
  if( eta < 0. )   return( parval[0]*pt + parval[6] );
  if( eta < 0.8 )  return( parval[0]*pt + parval[7] );
  if( eta < 1.2 )  return( parval[0]*pt + parval[8] );
  if( eta < 1.6 )  return( parval[0]*pt + parval[9] );
  if( eta < 1.8 )  return( parval[0]*pt + parval[10] );
  if( eta < 2.0 )  return( parval[0]*pt + parval[11] );
  return( parval[0]*pt + parval[12] );
}



/// Service to build the scale functor corresponding to the passed identifier       
resolFunctBase<double * > * resolutionFunctService( const int identifier ){
  switch ( identifier ) {
  case ( 45 ): return ( new resolFunct45<double * > ); break;
  default: std::cout << "resolutionFunctService error: wrong identifier = " << identifier << std::endl; exit(1);
  }
}
