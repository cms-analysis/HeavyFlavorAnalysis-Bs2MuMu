#ifndef CSM_H
#define CSM_H


// mclimit_csm -- use a chisquared minimized over nuisance parameters
// to compute cls -- pseudoexperiment loops, Bayesian limits and cross sections
// in 1D and 2D, and Markov Chain Bayesian calculations included

// Class to run the TMinuit minimization of T. Devlin's chisquared
// defined in CDF 3126, minimized over the nuisance parameters.

// version dated Nov 18, 2010
// Author:  Tom Junk, Fermilab.  trj@fnal.gov
// Contributions from Joel Heinrich, Nils Krumnack, Tom Wright, and Kevin Lannon

#include <iostream>
#include <vector>
#include <map>
#include <TH1.h>
#include <assert.h>
#include <memory>

//#include "FastTH1.hh"
//For convenience (so we only have one header file and one source file, put all the FastTH1 stuff in here)

#ifndef NK_MCLIMIT_FAST_TH1_HH
#define NK_MCLIMIT_FAST_TH1_HH

//          Copyright Nils Krumnack 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Please feel free to contact me (nils@fnal.gov) for bug reports,
// feature suggestions, praise and complaints.

// This class provides a fast implementation of a histogram based on
// the std::vector class, providing the functions mclimit needs.


// to switch between using fast TH1's and slow TH1's in mclimit,
// set the typedef of TH1Type later.

// uncomment the following to turn on extra debugging checks for this
// class
//#define NK_MCLIMIT_FAST_TH1_DEBUG

// the class defined in this module
class FastTH1;

// effects: make a copy/clone of the given histogram
// returns: the newly created object
// guarantee: strong
// throws: std::bad_alloc
// postcondition: result.get() != NULL
// rationale: this template function allows to interface quickly and
//   easily between fast and regular TH1 objects (if the interface
//   involves a copy operation).
template <class Hist> std::auto_ptr<Hist>
copy_TH1 (const TH1& that, const std::string& name = "");
template <class Hist> std::auto_ptr<Hist>
copy_TH1 (const FastTH1& that, const std::string& name = "");


// effects: copy the contents from histogram b to histogram a
// guarantee: no-throw
template<class Hist> void
copy_TH1_content (Hist& a, const Hist& b);

class TH1Input
{
  // this class is meant to be used for passing lists of histograms
  // that are either TH1 or FastTH1

  //
  // public interface
  //

  // effects: test the invariant of the class
  // guarantee: no-throw
public:
  void test_invariant () const;


  // effects: create an input list with traditional TH1 objects
  // guarantee: no-throw
public:
  TH1Input (TH1 *list);


  // effects: create an input list with fast TH1 objects
  // guarantee: no-throw
public:
  TH1Input (FastTH1 *list);


  // effects: makes a copy of the histogram with the given type
  // guarantee: strong
  // throws: std::bad_alloc
public:
  template<class Hist> std::auto_ptr<Hist> copy_TH1 () const;



  //
  // private interface
  //

  // description: the slow histogram (if used)
private:
  TH1 *slow_;


  // description: the fast histogram (if used)
private:
  FastTH1 *fast_;
};





class TH1InputList
{
  // this class is meant to be used for passing lists of histograms
  // that are either TH1 or FastTH1

  //
  // public interface
  //

  // effects: test the invariant of the class
  // guarantee: no-throw
public:
  void test_invariant () const;


  // effects: create an input list with traditional TH1 objects
  // guarantee: no-throw
public:
  TH1InputList (TH1 **list);


  // effects: create an input list with fast TH1 objects
  // guarantee: no-throw
public:
  TH1InputList (FastTH1 **list);


  // returns: whether the histogram at the given index exists
  // guarantee: no-throw
public:
  bool has (unsigned index) const;


  // effects: makes a copy of the histogram at the given index with
  //   the given type
  // guarantee: strong
  // throws: std::bad_alloc
public:
  template<class Hist> std::auto_ptr<Hist> copy_TH1 (unsigned index) const;



  //
  // private interface
  //

  // description: the slow histogram (if used)
private:
  TH1 **slow_;


  // description: the fast histogram (if used)
private:
  FastTH1 **fast_;
};





class FastTH1
{
  //
  // public interface
  //

  // description: the type I use for storing values
public:
  typedef float value_type;


  // description: the maximum number of dimensions supported
public:
  static const unsigned max_dim = 3;


  // effects: test the invariant of the class
  // guarantee: no-throw
public:
  void test_invariant () const;


  // effects: create a copy of the given histogram
  // guarantee: strong
  // throws: std::bad_alloc
  // requires: unsigned (that.GetDimension()) <= max_dim
public:
  explicit FastTH1 (const TH1& that);


  // returns: a TH1 with the same parameters.  if !name.empty(), the
  //   histogram will be given the passed in name
  // guarantee: strong
  // throws: std::bad_alloc
public:
  std::auto_ptr<TH1> make_TH1 (std::string name = "") const;


  // effects: set the content of this histogram to the contents of
  //   that histogram
  // guarantee: no-throw
  // requires: GetNbinsX() == that.GetNbinsX()
  // requires: GetNbinsY() == that.GetNbinsY()
  // requires: GetNbinsZ() == that.GetNbinsZ()
public:
  void set_content (const FastTH1& that);



  //
  // public TH1 interface
  //

  // returns: the name of the histogram
  // guarantee: no-throw
  // postcondition: result != NULL
public:
  const char *GetName () const;


  // effects: sets the property returned by GetName()
  // guarantee: strong
  // throws: std::bad_alloc
  // requires: name != NULL
  // postcondition: std::string (name) == this->GetName()
public:
  void SetName (const char *name);



  // returns: the title of the histogram
  // guarantee: no-throw
  // postcondition: result != NULL
public:
  const char *GetTitle () const;


  // effects: sets the property returned by GetTitle()
  // guarantee: strong
  // throws: std::bad_alloc
  // requires: title != NULL
  // postcondition: std::string (title) == this->GetTitle()
public:
  void SetTitle (const char *title);


  // returns: the number of dimension
  // guarantee: no-throw
public:
  unsigned GetDimension () const;


  // returns: the number of bins in x or 1 if GetDimension() < 1
  // guarantee: no-throw
public:
  unsigned GetNbinsX () const;


  // returns: the number of bins in y or 1 if GetDimension() < 2
  // guarantee: no-throw
public:
  unsigned GetNbinsY () const;


  // returns: the number of bins in z or 1 if GetDimension() < 3
  // guarantee: no-throw
public:
  unsigned GetNbinsZ () const;


  // returns: the content of bin binx
  // guarantee: no-throw
  // requires: GetDimension() == 1
  // requires: binx >= 1 && binx <= GetNbinsX()
public:
  value_type GetBinContent (unsigned binx) const;


  // returns: the error of bin binx
  // guarantee: no-throw
  // requires: GetDimension() == 1
  // requires: binx >= 1 && binx <= GetNbinsX()
public:
  value_type GetBinError (unsigned binx) const;


  // returns: the content of bin (binx,biny)
  // guarantee: no-throw
  // requires: GetDimension() <= 2
  // requires: binx >= 1 && binx <= GetNbinsX()
  // requires: biny >= 1 && biny <= GetNbinsY()
public:
  value_type GetBinContent (unsigned binx, unsigned biny) const;


  // returns: the error of bin (binx,biny)
  // guarantee: no-throw
  // requires: GetDimension() <= 2
  // requires: binx >= 1 && binx <= GetNbinsX()
  // requires: biny >= 1 && biny <= GetNbinsY()
public:
  value_type GetBinError (unsigned binx, unsigned biny) const;


  // returns: the content of bin binx
  // guarantee: no-throw
  // requires: GetDimension() == 1
  // requires: binx >= 1 && binx <= GetNbinsX()
public:
  void SetBinContent (unsigned binx, value_type content);


  // returns: the error of bin binx
  // guarantee: no-throw
  // requires: GetDimension() == 1
  // requires: binx >= 1 && binx <= GetNbinsX()
  // requires: error >= 0
public:
  void SetBinError (unsigned binx, value_type error);


  // returns: the content of bin (binx,biny)
  // guarantee: no-throw
  // requires: GetDimension() <= 2
  // requires: binx >= 1 && binx <= GetNbinsX()
  // requires: biny >= 1 && biny <= GetNbinsY()
public:
  void SetBinContent (unsigned binx, unsigned biny, value_type content);


  // returns: the error of bin (binx,biny)
  // guarantee: no-throw
  // requires: GetDimension() <= 2
  // requires: binx >= 1 && binx <= GetNbinsX()
  // requires: biny >= 1 && biny <= GetNbinsY()
  // requires: error >= 0
public:
  void SetBinError (unsigned binx, unsigned biny, value_type error);


  // returns: the sum of all bins
  // guarantee: no-throw
public:
  float Integral () const;



  //
  // private interface
  //

  // description: the name and title of the histogram
private:
  std::string name_, title_;


  // description: the dimension of the histogram
private:
  unsigned dimension_;


  // description: the binning used
private:
  unsigned nbins_[max_dim];
  float low_[max_dim];
  float high_[max_dim];


  // description: a multiplication factor to apply for bins in the
  //   given dimension
private:
  unsigned mult_[max_dim];


  // description: the vector containing the bin contents
private:
  typedef std::vector<value_type>::const_iterator bin_contents_iter;
  typedef std::vector<value_type>::iterator bin_contents_miter;
  std::vector<value_type> bin_contents_;


  // description: the vector containing the bin errors
private:
  typedef std::vector<value_type>::const_iterator bin_errors_iter;
  typedef std::vector<value_type>::iterator bin_errors_miter;
  std::vector<value_type> bin_errors_;
};





//
// inline methods
//

template<> inline std::auto_ptr<TH1>
copy_TH1 (const TH1& that, const std::string& name)
{
  std::auto_ptr<TH1> result;

  result.reset (dynamic_cast<TH1*>(that.Clone (name.c_str())));

  assert (result.get() != NULL && "postcondition failed"); //spec
  return result;
};



template<> inline std::auto_ptr<FastTH1>
copy_TH1 (const TH1& that, const std::string& name)
{
  std::auto_ptr<FastTH1> result;

  result.reset (new FastTH1 (that));
  if (!name.empty())
    result->SetName (name.c_str());

  assert (result.get() != NULL && "postcondition failed"); //spec
  return result;
};



template<> inline std::auto_ptr<TH1>
copy_TH1 (const FastTH1& that, const std::string& name)
{
  std::auto_ptr<TH1> result;

  result = that.make_TH1 (name);

  assert (result.get() != NULL && "postcondition failed"); //spec
  return result;
};



template<> inline std::auto_ptr<FastTH1>
copy_TH1 (const FastTH1& that, const std::string& name)
{
  std::auto_ptr<FastTH1> result;

  result.reset (new FastTH1 (that));
  if (!name.empty())
    result->SetName (name.c_str());

  assert (result.get() != NULL && "postcondition failed"); //spec
  return result;
};



template<> inline void
copy_TH1_content (TH1& a, const TH1& b)
{
  a.Reset ();
  a.Add (&b, 1);
};



template<> inline void
copy_TH1_content (FastTH1& a, const FastTH1& b)
{
  a.set_content (b);
};



inline void TH1Input ::
test_invariant () const
{
#ifndef NDEBUG
  assert (!(slow_ && fast_) && "invariant violated"); //spec
#endif
};



inline TH1Input ::
TH1Input (TH1 *list)
  : slow_ (list), fast_ (NULL)
{
  test_invariant ();
};



inline TH1Input ::
TH1Input (FastTH1 *list)
  : slow_ (NULL), fast_ (list)
{
  test_invariant ();
};



template<class Hist> inline std::auto_ptr<Hist> TH1Input ::
copy_TH1 () const
{
  test_invariant ();

  std::auto_ptr<Hist> result;

  if (slow_)
    result = ::copy_TH1<Hist> (*slow_);
  if (fast_)
    result = ::copy_TH1<Hist> (*fast_);

  return result;
};



inline void TH1InputList ::
test_invariant () const
{
#ifndef NDEBUG
  assert (!(slow_ && fast_) && "invariant violated"); //spec
#endif
};



inline TH1InputList ::
TH1InputList (TH1 **list)
  : slow_ (list), fast_ (NULL)
{
  test_invariant ();
};



inline TH1InputList ::
TH1InputList (FastTH1 **list)
  : slow_ (NULL), fast_ (list)
{
  test_invariant ();
};



inline bool TH1InputList ::
has (unsigned index) const
{
  test_invariant ();

  bool result = false;
  if (slow_)
    result = slow_[index];
  if (fast_)
    result = fast_[index];
  return result;
};



template<class Hist> inline std::auto_ptr<Hist> TH1InputList ::
copy_TH1 (unsigned index) const
{
  test_invariant ();

  std::auto_ptr<Hist> result;
  if (slow_ && slow_[index])
    result = ::copy_TH1<Hist> (*slow_[index]);
  if (fast_ && fast_[index])
    result = ::copy_TH1<Hist> (*fast_[index]);
  return result;
};



inline void FastTH1 ::
test_invariant () const
{
#ifndef NDEBUG
  // this is required to make a TH1
  assert (dimension_ >= 1 && dimension_ <= 3 && "invariant violated"); //spec
  // this is required to store the binning
  assert (dimension_ < max_dim && "invariant violated"); //spec
  unsigned nbin = 1;
  for (unsigned dim = dimension_; dim != max_dim; ++ dim)
    assert (nbins_[dim] == 1 && "invariant violated"); //spec
  for (unsigned dim = 0; dim != dimension_; ++ dim)
  {
    assert (nbins_[dim] > 0 && "invariant violated"); //spec
    assert (low_[dim] < high_[dim] && "invariant violated"); //spec
    nbin *= nbins_[dim];
  };
  assert (nbin == bin_contents_.size() && "invariant violated"); //spec
  assert (nbin == bin_errors_.size() && "invariant violated"); //spec
  for (bin_errors_iter error = bin_errors_.begin();
       error != bin_errors_.end(); ++ error)
  {
    assert (*error >= 0 && "invariant violated"); //spec
  };
#endif
};



inline void FastTH1 ::
set_content (const FastTH1& that)
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetNbinsX() == that.GetNbinsX() && "requirement failed"); //spec
  assert (GetNbinsY() == that.GetNbinsY() && "requirement failed"); //spec
  assert (GetNbinsZ() == that.GetNbinsZ() && "requirement failed"); //spec
#endif

  bin_contents_ = that.bin_contents_;
  bin_errors_ = that.bin_errors_;

#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif
};



inline const char *FastTH1 ::
GetName () const
{
  test_invariant ();

  const char *result = name_.c_str();

  assert (result != NULL && "postcondition failed");
  return result;
};



inline void FastTH1 ::
SetName (const char *name)
{
  test_invariant ();
  assert (name != NULL && "requirement failed"); //spec

  name_ = name;

  test_invariant ();
  assert (std::string (name) == this->GetName() && "postcondition failed"); //spec
};



inline const char *FastTH1 ::
GetTitle () const
{
  test_invariant ();

  const char *result = title_.c_str();

  assert (result != NULL && "postcondition failed");
  return result;
};



inline void FastTH1 ::
SetTitle (const char *title)
{
  test_invariant ();
  assert (title != NULL && "requirement failed"); //spec

  title_ = title;

  test_invariant ();
  assert (std::string (title) == this->GetTitle() && "postcondition failed"); //spec
};



inline unsigned FastTH1 ::
GetDimension () const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif

  return dimension_;
};



inline unsigned FastTH1 ::
GetNbinsX () const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif

  return nbins_[0];
};



inline unsigned FastTH1 ::
GetNbinsY () const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif

  return nbins_[1];
};



inline unsigned FastTH1 ::
GetNbinsZ () const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif

  return nbins_[2];
};



inline FastTH1::value_type FastTH1 ::
GetBinContent (unsigned binx) const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 1 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
#endif

  return bin_contents_[binx-1];
};



inline FastTH1::value_type FastTH1 ::
GetBinError (unsigned binx) const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 1 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
#endif

  return bin_errors_[binx-1];
};



inline FastTH1::value_type FastTH1 ::
GetBinContent (unsigned binx, unsigned biny) const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() <= 2 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
  assert (biny >= 1 && biny <= GetNbinsY() && "requirement failed"); //spec
#endif

  return bin_contents_[(binx-1) * mult_[0] + (biny-1) * mult_[1]];
};



inline FastTH1::value_type FastTH1 ::
GetBinError (unsigned binx, unsigned biny) const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 2 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
  assert (biny >= 1 && biny <= GetNbinsY() && "requirement failed"); //spec
#endif

  return bin_errors_[(binx-1) * mult_[0] + (biny-1) * mult_[1]];
};



inline void FastTH1 ::
SetBinContent (unsigned binx, value_type content)
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 1 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
#endif

  bin_contents_[binx-1] = content;

#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif
};



inline void FastTH1 ::
SetBinError (unsigned binx, value_type error)
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 1 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
  assert (error >= 0 && "requirement failed"); //spec
#endif

  bin_errors_[binx-1] = error;

#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif
};



inline void FastTH1 ::
SetBinContent (unsigned binx, unsigned biny, value_type content)
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 2 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
  assert (biny >= 1 && biny <= GetNbinsY() && "requirement failed"); //spec
#endif

  bin_contents_[(binx-1) * mult_[0] + (biny-1) * mult_[1]] = content;

#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif
};



inline void FastTH1 ::
SetBinError (unsigned binx, unsigned biny, value_type error)
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
  assert (GetDimension() == 2 && "requirement failed"); //spec
  assert (binx >= 1 && binx <= GetNbinsX() && "requirement failed"); //spec
  assert (biny >= 1 && biny <= GetNbinsY() && "requirement failed"); //spec
  assert (error >= 0 && "requirement failed"); //spec
#endif

  bin_errors_[(binx-1) * mult_[0] + (biny-1) * mult_[1]] = error;

#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif
};



inline float FastTH1 ::
Integral () const
{
#ifdef NK_MCLIMIT_FAST_TH1_DEBUG
  test_invariant ();
#endif

  float result = 0;
  for (bin_contents_iter content = bin_contents_.begin();
       content != bin_contents_.end(); ++ content)
    result += *content;
  return result;
};

#endif

//typedef TH1 TH1Type;
typedef FastTH1 TH1Type;

namespace NK { namespace McLimit { class ExpMcLimit; }; };

// for the map below -- from the SGI example for a map

struct csm_ltstr
{
  bool operator()(char* s1, char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};


typedef enum {
  flat=10,
  corr=20
} PRIOR;

struct svstruct_s
{
  Int_t itemplate;   // which template histo this syst. variation applies to
  char *sysname;     // name of nuisance parameter this corresponds to
  Double_t sysfracl;
  Double_t sysfrach;
  // if there is no shape uncertainty associated with this dependence of a particular
  // template on a nuisance parameter, then these pointers should be set to zero.
  TH1Type *lowshape;     // for shape uncertainty -- low histogram shape histo id
  TH1Type *highshape;    // for shape uncertainty -- high histogram shape
  Double_t xsiglow;  // how many sigma low lowshape corresponds to (should be a negative number)
  Double_t xsighigh; // how many sigma high highshape corresponds to. (should be a positive number)
};

typedef struct svstruct_s svstruct;

// Constraint equations between nuisance parameters.  Sometimes there are fewer degrees of freedom
// than there are nuisance parameters.  Example:  MET vs. Iso evaluation of the non-W
// background in a signal region, using the 4-sector method.  A*C/B=D relates four
// nuisance parameters down to three: D is computed from A, B, and C.  Note:  multiple constraints between nuisance
// parameters will not be solved for, they will just be evaluated in order.  In pseudexperiments,
// nuisance parameters that are constrained are not randomly chosen, but rather are calculated
// based on the other nuisance parameters.  The arbitrary function here means that no checking
// can be done to make sure that the computed nuisance parameters are physical.

struct npcstruct_s
{
  Int_t ninput;       // number of nuisance paramters input to the calculation of a constrained one
  char **pnameinput;  // the names of the input nuisance parameters
  char *pnameoutput;  // the name of the output parameter
  Double_t (*f)(Double_t*);  // a function taking an array of all the input parameters and computing the output parameter
};

typedef struct npcstruct_s npcstruct;

// one-sided or two-sided 3-sigma or 5-sigma

#define MCLIMIT_CSM_TWOSIDED

// what to do about 1-sided or 2-sided 2-sigmas?

#ifdef MCLIMIT_CSM_TWOSIDED
#define MCLIMIT_CSM_2S 0.02275
#define MCLIMIT_CSM_3S 1.349898E-3
#define MCLIMIT_CSM_5S 2.866516E-7
#else
#define MCLIMIT_CSM_2S 0.0455
#define MCLIMIT_CSM_3S 2.6998E-3
#define MCLIMIT_CSM_5S 5.7330E-7
#endif

// cumulative probabilities for defining bands on test statistic
// and CL plots

#define MCLIMIT_CSM_MCLM2S 0.02275
#define MCLIMIT_CSM_MCLM1S 0.16
#define MCLIMIT_CSM_MCLMED 0.5
#define MCLIMIT_CSM_MCLP1S 0.84
#define MCLIMIT_CSM_MCLP2S 0.97725

// some messages to pass around inside for the s95 calculator

#define MCLIMIT_CSM_CLS 1
#define MCLIMIT_CSM_CLSM2 2
#define MCLIMIT_CSM_CLSM1 3
#define MCLIMIT_CSM_CLSMED 4
#define MCLIMIT_CSM_CLSP1 5
#define MCLIMIT_CSM_CLSP2 6

#define MCLIMIT_CSM_LUMI95 1
#define MCLIMIT_CSM_LUMI3S 2
#define MCLIMIT_CSM_LUMI5S 3

// horizontal is the interpolation style which calls csm_pvmorph (or csm_pvmorph_2d)
// vertical interpolation just does a linear interpolation bin-by-bin, but not
// letting any bin go below zero

typedef enum {
  CSM_INTERP_HORIZONTAL,
  CSM_INTERP_VERTICAL,
  CSM_INTERP_HORIZONTAL_EXTRAP,
  CSM_INTERP_VERTICAL_EXTRAP
} INTERPSTYLE;

// use to steer the handling of bin-to-bin uncertainties in templates

#define CSM_NOBINERR 0
#define CSM_POISSON_BINERR 1
#define CSM_GAUSSIAN_BINERR 2

// use to steer the choice of integration method in bayes_heinrich and bayes_heinrich_withexpect

#define CSM_BAYESINTEGRAL_JOEL 0
#define CSM_BAYESINTEGRAL_QUICK 1

/* a full model for a particular channel.  csm_model is just a collection
   of channel models and associated names */

class csm_channel_model
{
 public:
   csm_channel_model();
   ~csm_channel_model();
   void add_template( const TH1Input&,      //template histogram
	   	      Double_t,   //scale factor to multiply template by to compare w/ data
                      Int_t,      // number of nuisance parameters (Gaussian of unit width)
                      char *[],   // nuisance parameter names 
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- low side
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
                      const TH1InputList&,    // array of low histogram shapes, one for each nuisance param
                      Double_t *, // number of sigma low for each nuisance parameter shape variation
                      const TH1InputList&,    // array of high histogram shapes, one for each nuisance param
		      Double_t *, // number of sigma high for each shape variation
                      Int_t,      // MC statistical error steering (inupt).    Values and meaning:
                                  //    2:  Pay attention to error bars in each bin
                                  //    1:  Poisson errors with unit-weight entries (see below)
                                  //    0:  No bin-by-bin errors.
                                  //  Notes on Poisson (value=1):  There is a split between the
                                  // Poisson handling when the model is an ensemble model (for pseudoexperiment
                                  // generation) and when it is a test model (for fitting to data and pseudodata).
                                  //  In an ensemble model, the bin contents are treated
                                  // as Poisson means and the values are fluctuated according to Poisson statistics.
                                  // These random components are then used with their scale factors as components
                                  // of another Poisson mean for generating pseudodata. 
                                  // In the model used to fit to the data, these values are treated as integer
                                  // measurements from subsidiary experiments.
                      Int_t);     // Scale flag -- set to 1 if you want this template to be scaled in the s95
                                  // calculation, 0 if you don't.  It's intended that applications should set this to 1 for
                                  // signal histograms, and 0 for background histograms.


   csm_channel_model* Clone();   // make an exact copy of this channel model -- all the internal
                                 // cloned histograms are cloned again.  Better than just
                                 // assigning a new model to this one because the
                                 // destructors won't delete the same memory twice.
   void nuisance_response(Int_t, char *[], Double_t []); // update the internal copies of the varied histograms
                                                         // and scale factors according to the nuisance parameters supplied
   void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values

   void print();                 // print out some details of the model
   void plotwithdata(TH1*);      // compare data with a model.
   void candcheck(TH1*);         // print out high s/b candidates
   double kstest(TH1*);          // get ROOT's raw KS prob
   double kstest_px(TH1*);       // get ROOT's raw KS px prob
   // adding two models together, and multiplying a model by a scalar

   // note that all of these operations on models create a new model (which must be cleaned up later)

   csm_channel_model* add(csm_channel_model&); // adds the components of two models together
   csm_channel_model* scale(Double_t);        // scales all parts of the model
   csm_channel_model* scalesignal(Double_t);  // scales only those parts of the model called "signal"
   csm_channel_model* scale_err(Double_t);  // scales rates up with scale factor, and scales
                                            // systematic errors down with scale factor.  n.b. --
                                            // MC statistical errors on model histos cannot be scaled
   csm_channel_model* interpolate(csm_channel_model *b, double frac); // interpolates frac of the way between this channel model
                                                      // and the channel model pointed to by b.  Returns a pointer to a newly
                                                      // created channel model which is the interpolation.

   Double_t chisquared1(TH1 *);  // Inputs a pointer to data histogram -- computes the chisquared of this model
                                 // compared to the supplied data histogram, without minimizing over nuisance
                                 // parameters, a la Tom Devlin's note

   Int_t checkneg();             // check for negative bins

   void set_interpolation_style(INTERPSTYLE);  // either CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL
                                               // horizontal (csm_pvmorph) is the default if this is not called.

   friend class csm_model;       // so we can access systematic error names and limits
   friend class mclimit_csm;

 private:
  std::vector<TH1Type*> histotemplate;
  std::vector<TH1Type*> histotemplate_varied;
  std::vector<Double_t> sft;
  std::vector<Double_t> sft_varied;
  std::vector<Int_t> poissflag;
  std::vector<Int_t> scaleflag; // 1 if signal, 0 if background
  std::vector<svstruct> syserr;
  std::map<char*, std::vector<int>, csm_ltstr> semap;  // a list of entries in the syserr vector indexed by name

  INTERPSTYLE chan_istyle;


  // description: this is a temporary object only to be used by
  //   nuisance_response
  // rationale: I hope this speeds up nuisance_response
private:
  std::auto_ptr<TH1Type> hcl_nuisance_response;

  // rationale: I made my wrapper a friend to avoid figuring out
  //   accessors
  friend class NK::McLimit::ExpMcLimit;
};

class csm_model
{
 public:
   csm_model();
   ~csm_model();
   void add_template( TH1 *,      //template histogram
	   	      Double_t,   //scale factor to multiply template by to compare w/ data
                      Int_t,      // number of nuisance parameters (Gaussian of unit width)
                      char *[],   // nuisance parameter names 
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- low side
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
                      TH1 *[],    // array of low histogram shapes, one for each nuisance param
                      Double_t *, // number of sigma low for each nuisance parameter shape variation (should be negative!)
                      TH1 *[],    // array of high histogram shapes, one for each nuisance param
		      Double_t *, // number of sigma high for each shape variation (should be positive!)
                      Int_t,      // Poisson flag -- 1 if Poisson, 0 of not.  There is a split between the
                                  // interpretation of this flag when the model is an ensemble model and when
                                  // it is a test model.  In an ensemble model, the bin contents are treated
                                  // as Poisson means and the values are fluctuated according to Poisson statistics.
                                  // These random components are then used with their scale factors as components
                                  // of another Poisson mean for generating pseudodata. 
                                  // In the model used to fit to the data, these values are treated as integer
                                  // measurements from subsidiary experiments.
                      Int_t,      // Scale flag -- set to 1 if you want this template to be scaled in the s95
                                  // calculation, 0 if you don't.  It's intended that applications should set this to 1 for
                                  // signal histograms, and 0 for background histograms.
                      char *);    // Channel name

   void set_interpolation_style(char *,INTERPSTYLE);  //  sets the interpolation style 
            //  for a particlar channel -- first arg: channel name
            // of channel -- second arg:  interpolation style:  CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL

   void add_chanmodel(csm_channel_model*,
                      char*); // instead of adding a template at a time, let's add a whole channel's model

   void add_npcons(Int_t, char**, char*, Double_t (*f)(Double_t*)); // for a nuisance parameter which
	    // can be computed given the values of other nuisance parameters, this allows
	    // the user to specify such a funtion.   Arguments;  number of n.p.'s the one to calculate
            // depends on, their names, the name of the n.p. to calculate, and a pointer to the
            // function which does the job.

   void add_npbounds(char *npname, Double_t lowbound, Double_t highbound);

   void plotwithdata(char*,TH1*); // plot a named channel's model with the supplied data histogram
   void candcheck(char*,TH1*); // check candidates
   double kstest(char*,TH1*); // get raw ROOT ks prob
   double kstest_px(char*,TH1*); // get ROOT's ks test prob with px simulation (stat only)

   csm_model* Clone();   // make an exact copy of this model -- all the internal
                                 // cloned histograms are cloned again.  Better than just
                                 // assigning a new model to this one because the
                                 // destructors won't delete the same memory twice.

   void print();                 // print out some details of the model
   void print(char *channame);   // just print out the piece corresponding to channame
   // adding two models together, and multiplying a model by a scalar

   // note that all of these operations on models create a new model (which must be cleaned up later)

   csm_model* add(csm_model&); // adds the components of two models together
   csm_model* scale(Double_t);        // scales all parts of the model
   csm_model* scalesignal(Double_t);  // scales only those parts of the model called "signal"
   csm_model* scale_err(Double_t);  // scales rates up with scale factor, and scales
                                            // systematic errors down with scale factor.  n.b. --
                                            // MC statistical errors on model histos cannot be scaled
   csm_model* interpolate(csm_model *b, double frac); // interpolates frac of the way between this model
                                                      // and the model pointed to by b.  Returns a pointer to a newly
                                                      // created model which is the interpolation.
   void nuisance_response(Int_t, char *[], Double_t []); // updates the fluctuated version of the histogram templates
                                 // and scale factors inside the channel model
                                 // according to the nuisance parameters provided.  Inputs: number of nuisance
                                 // parameters and their names
   void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values
   void varysyst();   // randomly choose nuisance parameter values and call nuisance_response
   void print_nuisance_params();  // after calling varysyst, can print out the randomly chosen parameters and names
   void single_pseudoexperiment(TH1 *[]); // generate pseudodata for all the channels in this model.
   void list_nparams(std::vector<char *> *, std::vector<Double_t> *, std::vector<Double_t> *); // get a list of
      // unique nuisance parameter names for all the channels in this model and their most
      // restrictive lower and upper bounds.
   Double_t chisquared1(TH1**); // calls the chisquared routine for each channel model.

   friend class mclimit_csm;
   friend class csm;

 private:
  std::vector<char*> channame;
  std::vector<csm_channel_model*> chanmodel;
  Int_t lookup_add_channame(char *);  // look up the channel name in the channame std::vector.  If it's
                                       // not there, add it.
  std::vector<npcstruct> npcm;  // constraint equations between nuisance parameters

  std::vector<char*> npbname;   // upper and lower bounds
  std::vector<Double_t> npbhigh;
  std::vector<Double_t> npblow;
  std::vector<Double_t> npvalp; // handles to stored lists of nuisance parameter names and values for
  std::vector<char*> npnp;

  // rationale: I made my wrapper a friend to avoid figuring out
  //   accessors
  friend class NK::McLimit::ExpMcLimit;
};


class csm
{
  public:
    	csm();
	~csm();
        void set_htofit(TH1*,char*); //histogram and channel name
        void set_modeltofit(csm_model*); // a set of template histograms to fit to the data histograms
	Double_t chisquared();           // calculates chisquared
        Int_t ndof();                    // calculates an (approximate!) NDOF
	Int_t    getnparams();
        Double_t getparam(Int_t);
        Double_t getperror(Int_t);
        Double_t getcov(Int_t, Int_t);  // get an entry out of the covariance matrix.
        char* getpname(Int_t);
        csm_model* getbestmodel(); // a model with the template histograms all normalized to the best fit values.
                                   // suitable for plotting.  This varies the input model (given in set_modeltofit)
                                   //  with the best fit nuisance parameters and returns a pointer to the same model
                                   // given in set_modeltofit.
        void plotcompare(char *); // make a stacked plot of data in the named channel against the central
                                              // value model provided.
        void setminuitmaxcalls(Int_t);  // Maximum number of function calls MINUIT is allowed to do per minimization
                                        // default: 500
        Int_t getminuitmaxcalls(); 
        void setminosmaxcalls(Int_t);  // Maximum number of function calls MINOS is allowed to do per parameter
                                        // default: 500
        Int_t getminosmaxcalls(); 
        void setminuitstepsize(Double_t);  // Initial step size for MINUIT fit parameters default: 0.1
        Double_t getminuitstepsize(); 
        void setminosflag(bool);   // true: call MINOS. Best to have printing set too if 
                                   // you're going to run MINOS.  False:  Do not call MINOS (default) 
        bool getminosflag(); 
	void setprintflag(bool);  // true:  let MINUIT print stuff out;  FALSE -- turn off MINUIT printing (default)
	bool getprintflag();

 private:
        std::vector<Double_t> fitparam;   // parameters of the fit
	std::vector<Double_t> fiterror;   // errors on fit parameters
        std::vector<char*> fitparamname; //  names of fit parameters
	Double_t *fitcov;           //  pointer to covariance matrix
	Int_t nfitcov;              //  size of covariance matrix
        Int_t minuitmaxcalls;       //  how many calls to the function Minuit is allowed to make
        Int_t minosmaxcalls;       //  how many calls to the function Minos is allowed to make
	Double_t minuitstepsize;    // initial step size for MINUIT parameters -- default 0.1
        bool minuitprintflag;       // true if we let MINUIT print stuff out.
        bool minosflag;             // true if we are to call MINOS (most useful when the printflag is on too)
};

// interpolate histogram with a pvmorph-style procedure, or with vertical interpolation.
// csm_interpolate_histogram interpolates the errors too, while csm_interpolate_histogram_noerr
// is a speedup which interpolates only the bin contents and not the errors

template<class Hist> void csm_interpolate_histogram(Hist*,Double_t,Hist*,Double_t,Hist*,Double_t,INTERPSTYLE);

template<class Hist> void csm_interpolate_histogram_noerr(Hist*,Double_t,Hist*,Double_t,Hist*,Double_t,INTERPSTYLE);

// version to be used with cascading shape errors -- needs a central shape, a varied shape,
// and a shape to apply the variations to (which may not be either of the above, but the result
// of a different shape variation) startshape.  The output is outshape.

// csm_interpolate histogram2 calls csm_interpolate_histogram3 twice, once for the
// bin contents, once for the errors.

template<class Hist> void
csm_interpolate_histogram2(Hist* central, Double_t paramcentral,
                                Hist* varied, Double_t paramvaried,
                                Hist* startshape, 
                                Hist* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

// here's a version that just interpolates the bin contents but not the errors
//  (twice as fast)

template<class Hist> void
csm_interpolate_histogram2_noerr(Hist* central, Double_t paramcentral,
                                Hist* varied, Double_t paramvaried,
                                Hist* startshape, 
                                Hist* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

// this routine just interpolates the histogram contents

template<class Hist> void
csm_interpolate_histogram3(Hist* central, Double_t paramcentral,
			   Hist* varied, Double_t paramvaried,
			   Hist* startshape, 
			   Hist* outshape,
			   Double_t param,
			   INTERPSTYLE istyle);

void csm_pvmc(Int_t nb, Double_t *dist1, Double_t *dist2, Double_t *dist3, Double_t *distn,
	      Double_t par1, Double_t par2, Double_t parn);

void csm_pvmc2d(Int_t nx, Int_t ny, Double_t *xydist1, 
                Double_t *xydist2, Double_t *xydist3, Double_t *xydistn,
                Double_t par1, Double_t par2, Double_t parn);

void csm_yproj(Int_t nx, Int_t ny, Double_t *xydist, Double_t *ydist);

void csm_ycont(Int_t ny, Double_t *ydist1, Double_t *ydist2,
               Double_t *ydist3, Double_t *ydistn,
               Double_t *alpha1, Double_t *alpha2, Double_t *alpha3);

void csm_ycontaux(Int_t ny, Double_t *y, Double_t *yn,
                  Double_t *alpha);

void csm_acnvec2(Double_t *vec, Int_t n);

//maximum number of iterations in the quadratic system solver
#define CSM_MAXITER 100
//if all rates change fractionally by less than this, then we declare
//the system to be solved.
#define PREC1 1.0e-8

#define PVMORPH_MAXBINS 5000
#define CSM_DEBUGPRINT 1

class mclimit_csm
{
 public:
   mclimit_csm();
   ~mclimit_csm();

   void print_version();

   void set_datahist(TH1 *,char *); /* data histogram and channel name */
   // the hypothesis set routines do not make clones of their inputs, they just
   // store pointers to the models. Best practice -- fully define a model (that
   // is, add all templates, before calling these routines and do not update
   // the models before the pseudoexperiments are run.

   void set_null_hypothesis(csm_model *);
   void set_test_hypothesis(csm_model *);
   void set_null_hypothesis_pe(csm_model *);
   void set_test_hypothesis_pe(csm_model *);
   void set_npe(Int_t);      // sets the number of pseudoexperiments to do.  
                             // The default is set in the constructor to 10000
   Int_t get_npe();          // returns the value set in set_npe

   void set_chisquarehistos(TH1 *,TH1 *,TH1 *,TH1 *);
   // set pointers to histograms to accumulate chisquare distributions for
   // null hyp chisquare distrib in null hyp pseudoexperiments
   // test hyp chisquare distrib in null hyp pseudoexperiments
   // null hyp chisquare distrib in test hyp pseudoexperiments
   // test hyp chisquare distrib in test hyp pseudoexperiments

   // these accessor methods are used to control the behavior of MINUIT in the calculation
   // of the test statistic.  Their values are just passed to the csm class instance when
   // the test statistic is computed.

   void setminuitmaxcalls(Int_t);  // maximum number of calls to MINUIT's minimization function
                                   // default: 500
   Int_t getminuitmaxcalls();
   void setminosmaxcalls(Int_t);  // maximum number of calls MINOS is allowed to use per parameter
                                   // default: 500
   Int_t getminosmaxcalls();
   void setminuitstepsize(Double_t);  // initial step size for MINUIT parameters (default=0.1)
   Double_t getminuitstepsize();
   void setminosflag(bool);   // true: call MINOS. Best to have printing set too if 
                              // you're going to run MINOS.  False:  Do not call MINOS (default) 
   bool getminosflag(); 
   void setprintflag(bool);  // true:  let MINUIT print stuff out;  FALSE -- turn off MINUIT printing (default)
   bool getprintflag();
   void setpxprintflag(bool);    // print out pseudoexperiment results
                             // this flag affects printing test statistic outputs in the
                             // run_pseudoexperiment loop, as well as bayes_heinrich_withexpect
                             // and bayes_heinrich_coverage check
   bool getpxprintflag();
   void setprintnpextremeflag(bool); // turn on printing of nuisance parameters if the null hyp has
                                     // a very signal-like outcome.  True: print, False: don't print.  Default: F
   void setprintextremevalue(Double_t);  // cut on how negative -2lnQ(bg) has to be before triggering
                                         // a nuisance parameter printout

   void run_pseudoexperiments(); // see set_npe() above to determine how many pseudoexperiments to run

   // output retrieval:

   Double_t cls();
   Double_t clsw();  // cls+b computed with null hyp px's reweighted -- good for small cls's
   Double_t clsb();
   Double_t clsbw();  // cls+b computed with null hyp px's reweighted -- good for small cls's
   Double_t clb();
   Double_t omclb(); // "1-clb" computed as a p-value, including the probability of the exact outcome
                     // observed in the data.  Computed with null hypothesis pseudoexperiments
   Double_t omclbw(); // Same as above, but using test hypothesis pseudoexperiments, reweighted
   // with the likelihood ratio to approximate the null hypothesis distribution

   Double_t ts();  // test statistic -- a delta chisquared -- computed for the observed data

   Double_t tsbm2(); // distributions of test statistic in null hyp pseudoexperiments 2 sigma low edge
   Double_t tsbm1(); // 1 sigma low edge
   Double_t tsbmed(); // median test statistic in null hyp pseudoexperiments
   Double_t tsbp1();  // 1 sigma upper edge
   Double_t tsbp2();  // 2 sigma upper edge

   Double_t tssm2(); // distributions of test statistic in test hyp pseudoexperiments 2 sigma low edge
   Double_t tssm1(); // 1 sigma low edge
   Double_t tssmed(); // median test statistic in null hyp pseudoexperiments
   Double_t tssp1();  // 1 sigma upper edge
   Double_t tssp2();  // 2 sigma upper edge
 
   Double_t clsexpbm2(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbm1(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbmed(); // Expected cls in null hypothesis -- median
   Double_t clsexpbp1(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t clsexpbp2(); // Expected cls in null hypothesis -- 2 sigma upper edge

   // generalize the above a bit in case CLs is not a monotonic function of -2lnQ

   Double_t gclsexpbm2(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t gclsexpbm1(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t gclsexpbmed(); // Expected cls in null hypothesis -- median
   Double_t gclsexpbp1(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t gclsexpbp2(); // Expected cls in null hypothesis -- 2 sigma upper edge

   Double_t clsbexpbm2(); // Expected clsb in null hypothesis -- 2 sigma low edge
   Double_t clsbexpbm1(); // Expected clsb in null hypothesis -- 2 sigma low edge
   Double_t clsbexpbmed(); // Expected clsb in null hypothesis -- median
   Double_t clsbexpbp1(); // Expected clsb in null hypothesis -- 1 sigma upper edge
   Double_t clsbexpbp2(); // Expected clsb in null hypothesis -- 2 sigma upper edge


   // computed with null hypothesis px's reweighted -- good for small expected CLs's

   Double_t clsexpbm2w(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbm1w(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbmedw(); // Expected cls in null hypothesis -- median
   Double_t clsexpbp1w(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t clsexpbp2w(); // Expected cls in null hypothesis -- 2 sigma upper edge

   Double_t clsexpsm2(); // Expected cls in test hypothesis -- 2 sigma low edge
   Double_t clsexpsm1(); // Expected cls in test hypothesis -- 2 sigma low edge
   Double_t clsexpsmed(); // Expected cls in test hypothesis -- median
   Double_t clsexpsp1(); // Expected cls in test hypothesis -- 1 sigma upper edge
   Double_t clsexpsp2(); // Expected cls in test hypothesis -- 2 sigma upper edge
 
   // these accessors below use the CLs definition of CLb which includes the 
   // probability of observing exactly the data outcome
   // (subtracting it from 1 makes 1-CLb computed with these routines omit the
   // probability of observing exactly the data outcome)  Not to be used
   // for discovery significance!
   Double_t clbexpsm2(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t clbexpsm1(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t clbexpsmed(); // Expected clb in test hypothesis -- median
   Double_t clbexpsp1(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t clbexpsp2(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // these accessors below use the p-value definition of 1-CLb which includes the
   // probability of observing exactly the data outcome.  These are computed
   // using null hypothesis px's to compute 1-CLb's.
   Double_t omclbexpsm2(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsm1(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsmed(); // Expected clb in test hypothesis -- median
   Double_t omclbexpsp1(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t omclbexpsp2(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // Same as above, but use the test hypothesis pseudoexperiments, reweighted with
   // the likelihood ratio, to approximate the distribution of the null hypothesis
   // distribution of the test statistic
   Double_t omclbexpsm2w(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsm1w(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsmedw(); // Expected clb in test hypothesis -- median
   Double_t omclbexpsp1w(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t omclbexpsp2w(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // these accessors below use the CLs definition of CLb which includes the 
   // probability of observing exactly the data outcome
   // (subtracting it from 1 makes 1-CLb computed with these routines omit the
   // probability of observing exactly the data outcome)
   Double_t clbexpbm2(); // Expected clb in null hypothesis -- 2 sigma low edge
   Double_t clbexpbm1(); // Expected clb in null hypothesis -- 2 sigma low edge
   Double_t clbexpbmed(); // Expected clb in null hypothesis -- median
   Double_t clbexpbp1(); // Expected clb in null hypothesis -- 1 sigma upper edge
   Double_t clbexpbp2(); // Expected clb in null hypothesis -- 2 sigma upper edge

   // these accessors below use the p-value definition of 1-CLb which includes the
   // probability of observing exactly the data outcome
   Double_t omclbexpbm2(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbm1(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbmed(); // Expected 1-clb in null hypothesis -- median
   Double_t omclbexpbp1(); // Expected 1-clb in null hypothesis -- 1 sigma upper edge
   Double_t omclbexpbp2(); // Expected 1-clb in null hypothesis -- 2 sigma upper edge

   // Same as above, but use the test hypothesis pseudoexperiments, reweighted with
   // the likelihood ratio, to approximate the distribution of the null hypothesis
   // distribution of the test statistic
   Double_t omclbexpbm2w(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbm1w(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbmedw(); // Expected 1-clb in null hypothesis -- median
   Double_t omclbexpbp1w(); // Expected 1-clb in null hypothesis -- 1 sigma upper edge
   Double_t omclbexpbp2w(); // Expected 1-clb in null hypothesis -- 2 sigma upper edge

   //these three below are computed with test hyp. px's reweighted to get the small
   //tails of the null hyp distribution better modeled with low px stats
   Double_t p2sigmat(); // Probability of a 2-sigma evidence assuming test hyp. is true
   Double_t p3sigmat(); // Probability of a 3-sigma evidence assuming test hyp. is true
   Double_t p5sigmat(); // Probability of a 5-sigma discovery assuming test hyp. is true

   //these are computed with the null hyp px's -- can be unreliable for low px stats.
   //hopefully the 1-CLb p-value is uniformly distributed between 0 and 1, but
   // these routines, as well as omclbexpb* are designed to quantify just how
   // uniform that is.
   Double_t p2sigman(); // Probability of a 2-sigma evidence assuming null hyp. is true
   Double_t p3sigman(); // Probability of a 3-sigma evidence assuming null hyp. is true
   Double_t p5sigman(); // Probability of a 5-sigma discovery assuming null hyp. is true

   Double_t calc_chi2(csm_model *, TH1 *[]);  // interface to chisquare calculator
                                                      // using our models

   Double_t weightratio(csm_model *testhyp, csm_model *nullhyp, TH1 *[]); // weight factor
   // for reweighting MC px's using the varied templates.

   // rate limit calculators
   Double_t s95();      // scale factor on signal which is excluded at exactly 95% CL
   Double_t s95m2();   // variation around the median expected s95 in the bg hypothesis, -2 sigma
   Double_t s95m1();   // variation around the median expected s95 in the bg hypothesis, -1 sigma
   Double_t s95med();  // median expected s95 in the background hypothesis
   Double_t s95p1();   // variation around the median expected s95 in the bg hypothesis, +1 sigma
   Double_t s95p2();   // variation around the median expected s95 in the bg hypothesis, +2 sigma

   Double_t lumi95(); // calculates the lumi needed for a median experiment to exclude at 95% CL
                      // what's returned is a multiplicative factor on whatever luminosity was used
                      // to construct the test and null hypotheses.
   Double_t lumi3s(); // calculates the lumi needed for a median experiment to discover at 3 sigma
   Double_t lumi5s(); // calculates the lumi needed for a median experiment to discover at 5 sigma

   void tshists(TH1*,TH1*);  // fills histograms with test statisic values in the pseudoexperiments
                             // (you define the binning).  First histo: test hypothesis, second histo:
                             // null hypothesis

   void plotlnsb(TH1 *mcb_hist, TH1 *mcs_hist, TH1 *data_hist); // make a summed plot of ln(1+s/b) given the 
                             // input histogram pointers.  They are filled in with summed MC 
                             // and data (they are reset first).

   // Call Joel Heinrich's (CDF 7587) genlimit
   //Bayesian limit calculator.  First arg:  credibility level:  e.g., 0.95.  Second arg, 
   //scale factor on signal to produce the limit.  Third arg, uncertainty on the limit scale factor.
   //as with the s95 routines above, this assumes that the signal adds incoherently to the
   //background.  Requires set_test_hypothesis_pe to be called first in order to make the
   //"Prior ensemble".  The size of the prior ensemble is set with set_npe()
   //also call the relevant set_datahist methods too.

   void bayes_heinrich(Double_t beta, Double_t* sflimit, Double_t* unc);

   // Routine to call Joel Heinrich's (CDF 7587) genlimit, a Bayesian limit calculator,
   // but to repeat the calculation for pseudoexpeirments drawn from the null hypothesis
   // (be sure to call both set_test_hypothesis_pe and set_null_hypothesis_pe before using
   // this).  This computes the observed and expected limits.
   // Arguments:  1:  beta (credibility level, e.g. 0.95)
   // Argument 2:  observed limit
   // Agrument 3:  error on observed limit
   // Argument 4:  npx pseudoexperiments to run to compute expected limits
   // Arguments 5-9: Expected limits.  Median +-1, 2 sigma expectations

   void bayes_heinrich_withexpect(Double_t beta, Double_t* sflimit, Double_t* unc,
				  Int_t npx,
                                  Double_t* sm2, Double_t* sm1, Double_t* smed, Double_t* sp1,
                                  Double_t* sp2);  // Call Joel Heinrich's (CDF 7587) genlimit

   // Some extra things that bayes_heinrich and bayes_heinrich_withexpect will compute,
   // if requested.  You can get a plot of the posterior PDF by specifying the range
   // over which it is to be evaluated and the point sample density  -- These are initialized
   // to zero by the constructor.  Just specify the beginning and the end of the interval and
   // the step size, and the bayes_posterior vector will be filled in when bayes_heinrich and 
   // bayes_heinrich_withexpect are called.  bayes_interval_end > bayes_interval_begin and
   // bayes_inteval_step > 0 for bayes_posterior to be filled in.

   Double_t bayes_interval_begin;
   Double_t bayes_interval_end;
   Double_t bayes_interval_step;

   // these two vectors are filled if the three paramters above are specified
   // and bayes_heinrich or bayes_heinrich_withexpect is called.  Plot the
   // contents of the first vector versus the values of the second vector.

   std::vector<Double_t> bayes_posterior;
   std::vector<Double_t> bayes_posterior_points;

   // Default method is Joel's (CDF 7587), which is option 0 (CSM_BAYESINTEGRAL_JOEL)
   // Option 1 (CSM_BAYESINTEGRAL_QUICK) uses the interval above and performs a quick and
   // dirty numerical sum integral to speed up slow calculations.

   void set_bayes_integration_method(int imethod);
   int get_bayes_integration_method();

   // With bayes_heinrich_withexpect, you also can get a histogram of expected limits on the
   // background-only pseudoexperiments.  Just specify the histogram to receive the entries
   // (it is reset on entry into bayes_heinrich_withexpect and filled in) -- the pointer is
   // initially null, and no filling occurs unless this pointer is set.

   TH1F *bayes_pseudoexperiment_limits;

   void bayes_heinrich_coverage_check(Double_t beta,
                                      Double_t sfscale,
				      Int_t npx,
                                      Double_t* falsex);  // checks the false exclusion rate (coverage)

// accessor methods (rather cobbled onto a method of getting access to a 
// value which is global to mclimit_csm.C, since the Bayesian routines are 
// written in C and not C++, this value wasn't put in a nicer place.

// dlcsn is an exponential factor on the posterior.  If you want the sum of the likelihoods,
// you need to take the posterior and multiply it by exp(-dlcsn)/bnorm, where dlcsn is the
// output of getdlcsn, and bnorm is the output of getbnorm

    double getdlcsn();

    double getdlcsn2d();

    double getbhnorm();

// Markov chain limit calculator -- all arguments are optional
// Beta is the credibility level desired for the limits -- the default is 0.95.
// prior is flat or correlated -- See Joel's CDF 7587, sec. 6

    double bayeslimit_mcmc1(double beta=0.95, PRIOR prior=corr, TString histoutfile="");

    void bayeslimit_mcmc1_expect(double beta, PRIOR prior, int npx, 
					      double *sm2, double *sm1, double *smed, 
					      double *sp1, double *sp2);


// log likelihood function with Gauss constraints on nuisance parameters with the optional
// correlated signal prior
// note -- for efficiency reasons, this function will not call nuisance_response, assuming it has
// already been called for the model.  So only the bin-by-bin Poisson uncertainties and the signal
// scale factor are used to adjust predictions.

    double bayesmcmcllf(csm_model *model, int nnptot, double *npvalues, 
                        int nperrtot, double *pnpvalues, double ssf, PRIOR prior);

   //---------------------------------------------------------------

   // fit a cross section using Joel Heinrich's routines  -- only testhyp_pe is used here.
   // cross sections are fit as scale factors of the cross section put in the templates
   // uses the data histograms supplied.  This also fills in the vectors
   // bayes_posterior and bayees_posterior_points
   // for now it needs a suggestion of what interval the maximum is expected
   // (to set the scale), so set bayes_interval_begin and bayes_interval_end and bayes_interval_step
   // make bayes_interval_step very small (otherwise there will be some discretization in your answer)

   void bh_xsfit(Double_t *xsfit, Double_t *downerr, Double_t *uperr);

   // draw pseudoexperiments from the fluctuated testhyp and find the cross section
   // from each one using testhyp_pe  Get expected values.

   void bh_xsfit_expect(Int_t npx, Double_t *xsfitavg, Double_t *m2s, 
     Double_t *m1s, Double_t *med, Double_t *p1s, Double_t *p2s);

   // scan over two signals and print out the marginalized posterior for later analysis
   // always assume the two signals are declared in the same order in all channels.
   // If more than two signals are present, the first one in  each channel is called signal 1,
   // and the sum of all others is called signal 2
   // Additional optional default arguments -- these are so bh_2d_scan_expect can use bh_2d_scan.

   void bh_2d_scan(Double_t s1low, Double_t s1high, Double_t ds1,
                   Double_t s2low, Double_t s2high, Double_t ds2,
		   // optional arguments -- can be omitted when just calling bh_2d_scan
		   Double_t *s1fit=0, Double_t *s2fit=0,
		   Int_t pflag2d=1, Double_t *s1test=0, Double_t *s2test=0,
		   Int_t *in68 = 0, Int_t *in95=0);

// This routine below runs pseudoexperiments to get the expected distributions.  It does all the
// calculations that bh_2d_scan does for the real data, but for the pseudoexperiments instead, and prints
// out the 2D posterior distributions, cumulative integrals, and best-fit points for each pseudoexperiment.
// As an added feature, you put in the signal scale factors used in the pseudoexperiment generation
//  (testhyp is used for pseudoexperiments, and testhyp_pe -- yes, I know that seems backwards -- is the
//  model used to compute the posterior), relative to the model used to test as s1true and s2true, and it
// will compute 68% and 95% coverage fractions based on the nearest grid point.  Verbose printout is
// supplied.

   void bh_2d_scan_expect(Int_t npx, Double_t s1true, Double_t s2true,
                   Double_t s1low, Double_t s1high, Double_t ds1,
                   Double_t s2low, Double_t s2high, Double_t ds2);

   TH1* get_datahist(Int_t);

   void printsbd();              // dump all signal, background, candidates

   void systsumplot(char *channame, int nsamples, TH1 **ehm2, TH1 **ehm1, TH1 **ehmed, TH1 **ehp1, TH1 **ehp2);

 private:
   csm_model *null_hypothesis;
   csm_model *test_hypothesis;
   csm_model *null_hypothesis_pe;
   csm_model *test_hypothesis_pe;

   Int_t nmc;  // number of pseudoexperiments which have been done
               // this is set to zero by the constructor and by
               // anything that modifies the hypotheses, which is the
               // indication that the pseudoexperiments need to be redone.

   Int_t nmc_req;  // number of pseudoexperiments to do

   Int_t recalctsflag; // 1 if we need to redo the data test statistic

   Double_t *tsb;  // test statistic in null hypothesis -- one per pseudoexperiment
   Double_t *tss;  // test statistic in test hypothesis -- one per pseudoexperiment
   Double_t *wtsb;  // weight for converting the corresponding tsb into a s+b px probability
   Double_t *wtss;  // weight for converting the corresponding tss into a b px probability
   Int_t *itsb; // sort order for tsb
   Int_t *itss; // sort order for tss
   Double_t tsd;                     // test statistic for data
   std::vector<TH1*> datahist;
   std::vector<char*> dhname;            // channel names for each data histogram -- must match
                                    // the channel names used in the models.

   Double_t s95aux(Int_t); // s95 calculator using a function you pass in.
   Double_t lumipaux(Int_t); // luminosity threshold calculator
   Double_t clsaux(Double_t); // cls, clsb and clb for an arbitrary test statistic
   Double_t clsauxw(Double_t); // cls, clsb and clb for an arbitrary test statistic
   Double_t clsbaux(Double_t);
   Double_t clsbauxw(Double_t);
   Double_t clbaux(Double_t);
   Double_t omclbaux(Double_t);
   Double_t omclbauxw(Double_t);
   Double_t gclsaux(Double_t);

   TH1 *nullnullchisquare;
   TH1 *nulltestchisquare;
   TH1 *testnullchisquare;
   TH1 *testtestchisquare;

   Int_t minuitmaxcalls;
   Int_t minosmaxcalls;
   Double_t minuitstepsize;
   bool minuitprintflag;       // true if we let MINUIT print stuff out.
   bool minosflag;             // true if we are to call MINOS
   bool pxprintflag;           // print out results of pseudoexperiments
   bool extremenpprintflag;    // print out nuisance parameters on extreme null hypothesis pseudoexperiments
   Double_t extremenpprintvalue;  // cut on how negative -2lnQ(bg) has to be before triggering a printout

   int bayesintegralmethod;
   Double_t quickbint(Double_t); // quick and dirty integrator

   double bayesmcmcllf(csm_model *model, int nnptot, std::vector<char*> npnames,
			double *npvalues, int nperrtot, double *pnpvalues, double ssf, PRIOR prior);

  // rationale: I made my wrapper a friend to avoid figuring out
  //   accessors
  friend class NK::McLimit::ExpMcLimit;
};

#endif

#ifndef GENLIMIT

typedef struct {
  float e,b;
}EB;

typedef struct {
  float e1,e2,b;
}EB2D;

double cspdf(double s,double norm,
             int nchan,int nens,const int nobs[],const EB* ens,PRIOR prior);

double cspdf2d(double s1, double s2, double norm,
             int nchan,int nens,const int nobs[],const EB2D* ens, PRIOR prior);

double csint(double s0,int nchan,int nens,const int nobs[],const EB* ens,
	     int* ngl,double xgl[],double lwgl[],
	     PRIOR prior,double* uncertainty);

void csint2(double s1,double s2,
	    int nchan,int nens,const int nobs[],const EB* ens,
	    int* ngl,double xgl[],double lwgl[],PRIOR prior,
	    double* int1,double* int2,
	    double* v11,double* v12,double* v22);

void csint2cut(double s1,double s2,double shi,
	       int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],PRIOR prior,
	       double* int1,double* int2,
	       double* v11,double* v12,double* v22);

void setdlcsn(int nchan, int nens, int nobs[], const EB* ens);

void setdlcsn2d(int nchan, int nens, int nobs[], const EB2D* ens);

double cslimit(double beta,int nchan,int nens,const int nobs[],const EB* ens,
               int* ngl,double xgl[],double lwgl[],
	       PRIOR prior,double* uncertainty);

double cscutlimit(double beta,double smax,
		  int nchan,int nens,const int nobs[],const EB* ens,
		  int* ngl,double xgl[],double lwgl[],
		  PRIOR prior,double* uncertainty);

double galim(double beta,int nchan,int nens,const int nobs[],const EB* ens);

void gausslaguerre(double x[],double lw[],int n,double alpha);

#define GENLIMIT 1

/*------------------------------------------------------------------------*/

// interpolate 1D histograms and 2D histograms
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output

// new version -- rely on the more general interpolator with three inputs, but reduce the argument
// count for backward compatibility

// csm_interpolate_histogram interpolates the bin contents and errors

template<class Hist> inline void
csm_interpolate_histogram(Hist* a, Double_t xa, 
			  Hist* b, Double_t xb,
			  Hist* c, Double_t xc,
			  INTERPSTYLE istyle)
{
  csm_interpolate_histogram2(a,xa,b,xb,a,c,xc,istyle);
}

// csm_interpolate_histogram_noerr interpolates just the bin contents but not the errors

template<class Hist> void
csm_interpolate_histogram_noerr(Hist* a, Double_t xa, 
				Hist* b, Double_t xb,
				Hist* c, Double_t xc,
				INTERPSTYLE istyle)
{
  csm_interpolate_histogram2_noerr(a,xa,b,xb,a,c,xc,istyle);
}

// interpolate 1D histograms and 2D histograms 
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output
// d is the histogram to apply the shift given by a and b to, for compounded interpolations.

// approximate attempt to interpolate the uncertainties too.  Problem is, an interpolated
// histogram is a long sum of pieces interpolated from the same central value histogram,
// and thus the errors are correlated in interesting ways.
// A subterfuge -- jut linearly interpolate the errors in the same way that the
// bin contents are linearly interpolated.  It's not a full error propagation.  Halfway interpolations
// really are averages of statistically uncertain histograms, and thus the error on the average should
// be a bit better than the error on either one.  But itnterpolate again, and correlations have to be
// taken into account to do it right.
// we've also lost at this point whether the errors need to be interpolated, but let's
// do them for all histograms.
// speedup 9 Dec 2007 -- avoid cloning TH1's as this is slow

template<class Hist> inline void
csm_interpolate_histogram2(Hist* a, Double_t xa, 
			   Hist* b, Double_t xb,
			   Hist* d,
			   Hist* c, Double_t xc,
			   INTERPSTYLE istyle)
{
  Int_t i,j;
  Double_t xtmp;
  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      std::cout << "nbins mismatch1 in csm_interpolate_histogram2: " << nbinsa << " " << nbinsb << std::endl;
    }
  if (nbinsb != nbinsc)
    {
      std::cout << "nbins mismatch2 in csm_interpolate_histogram2: " << nbinsb << " " << nbinsc << std::endl;
    }
  if (nbinsc != nbinsd)
    {
      std::cout << "nbins mismatch3 in csm_interpolate_histogram2: " << nbinsc << " " << nbinsd << std::endl;
    }
  if (nbinsya != nbinsyb)
    {
      std::cout << "nbinsy mismatch1 in csm_interpolate_histogram2: " << nbinsya << " " << nbinsyb << std::endl;
    }
  if (nbinsyb != nbinsyc)
    {
      std::cout << "nbinsy mismatch2 in csm_interpolate_histogram2 " << nbinsyb << " " << nbinsyc << std::endl;
    }
  if (nbinsyc != nbinsyd)
    {
      std::cout << "nbinsy mismatch3 in csm_interpolate_histogram2: " << nbinsyc << " " << nbinsyd << std::endl;
    }

  if (xb == xa)
    {
      std::cout << "xb == xa in csm_interpolate_histogram2 " << xa << std::endl;
      std::cout << "fatal error -- exiting." << std::endl;
      throw("nla");
      //      exit(0);
    }

  // interpolate contents

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

  // swap errors and contents and interpolate again  (approximate method for evaluating
  // errors on interpolated histograms)
  // be careful to swap only once, even if some pointers are repeated.

  if (nbinsya == 1)
    {
      for (i=1;i<=nbinsa;i++)
	{
	  xtmp = a->GetBinContent(i);
	  a->SetBinContent(i,a->GetBinError(i));
	  a->SetBinError(i,xtmp);
	  if (a != b)
	    {
	      xtmp = b->GetBinContent(i);
	      b->SetBinContent(i,b->GetBinError(i));
	      b->SetBinError(i,xtmp);
	    }
	  // c is the output histogram -- hopefully it is not the same as one of the input histograms
          xtmp = c->GetBinContent(i);
          // c->SetBinContent(i,c->GetBinError(i));
	  c->SetBinError(i,xtmp);

	  if (a != d && b != d)
	    {
	      xtmp = d->GetBinContent(i);
	      d->SetBinContent(i,d->GetBinError(i));
	      d->SetBinError(i,xtmp);
	    }
	}
    }
  else
    {
      for (i=1;i<=nbinsa;i++)
	{
	  for (j=1;j<=nbinsya;j++)
	    {
	       xtmp = a->GetBinContent(i,j);
	       a->SetBinContent(i,j,a->GetBinError(i,j));
	       a->SetBinError(i,j,xtmp);
	       if (a != b)
		 {
	           xtmp = b->GetBinContent(i,j);
	           b->SetBinContent(i,j,b->GetBinError(i,j));
	           b->SetBinError(i,j,xtmp);
		 }
	       xtmp = c->GetBinContent(i,j);
	       //c->SetBinContent(i,j,c->GetBinError(i,j));
	       c->SetBinError(i,j,xtmp);
	       if (a != d && b != d)
		 {
	           xtmp = d->GetBinContent(i,j);
	           d->SetBinContent(i,j,d->GetBinError(i,j));
	           d->SetBinError(i,j,xtmp);
		 }
	    }
	}
    }

  // interpolate the errors now and swap them back -- put the
  // original histograms back together again too

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

  if (nbinsya == 1)
    {
      for (i=1;i<=nbinsa;i++)
	{
	  xtmp = a->GetBinContent(i);
	  a->SetBinContent(i,a->GetBinError(i));
	  a->SetBinError(i,xtmp);
	  if (a != b)
	    {
	      xtmp = b->GetBinContent(i);
	      b->SetBinContent(i,b->GetBinError(i));
	      b->SetBinError(i,xtmp);
	    }
	  xtmp = c->GetBinContent(i);
	  c->SetBinContent(i,c->GetBinError(i));
	  c->SetBinError(i,xtmp);
	  if (a != d && b != d)
	    {
	      xtmp = d->GetBinContent(i);
	      d->SetBinContent(i,d->GetBinError(i));
	      d->SetBinError(i,xtmp);
	    }
	}
    }
  else
    {
      for (i=1;i<=nbinsa;i++)
	{
	  for (j=1;j<=nbinsya;j++)
	    {
	       xtmp = a->GetBinContent(i,j);
	       a->SetBinContent(i,j,a->GetBinError(i,j));
	       a->SetBinError(i,j,xtmp);
	       if (a != b)
		 {
	           xtmp = b->GetBinContent(i,j);
	           b->SetBinContent(i,j,b->GetBinError(i,j));
	           b->SetBinError(i,j,xtmp);
		 }
	       xtmp = c->GetBinContent(i,j);
	       c->SetBinContent(i,j,c->GetBinError(i,j));
	       c->SetBinError(i,j,xtmp);
	       if (a != d && b != d)
		 {
	           xtmp = d->GetBinContent(i,j);
	           d->SetBinContent(i,j,d->GetBinError(i,j));
	           d->SetBinError(i,j,xtmp);
		 }
	    }
	}
    }
}

template<class Hist> void
csm_interpolate_histogram2_noerr(Hist* a, Double_t xa, 
				 Hist* b, Double_t xb,
				 Hist* d,
				 Hist* c, Double_t xc,
				 INTERPSTYLE istyle)
{
  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      std::cout << "nbins mismatch1 in csm_interpolate_histogram2_noerr: " << nbinsa << " " << nbinsb << std::endl;
    }
  if (nbinsb != nbinsc)
    {
      std::cout << "nbins mismatch2 in csm_interpolate_histogram2_noerr: " << nbinsb << " " << nbinsc << std::endl;
    }
  if (nbinsc != nbinsd)
    {
      std::cout << "nbins mismatch3 in csm_interpolate_histogram2_noerr: " << nbinsc << " " << nbinsd << std::endl;
    }
  if (nbinsya != nbinsyb)
    {
      std::cout << "nbinsy mismatch1 in csm_interpolate_histogram2_noerr: " << nbinsya << " " << nbinsyb << std::endl;
    }
  if (nbinsyb != nbinsyc)
    {
      std::cout << "nbinsy mismatch2 in csm_interpolate_histogram2_noerr " << nbinsyb << " " << nbinsyc << std::endl;
    }
  if (nbinsyc != nbinsyd)
    {
      std::cout << "nbinsy mismatch3 in csm_interpolate_histogram2_noerr: " << nbinsyc << " " << nbinsyd << std::endl;
    }

  if (xb == xa)
    {
      std::cout << "xb == xa in csm_interpolate_histogram2_noerr " << xa << std::endl;
      std::cout << "fatal error -- exiting." << std::endl;
      //      exit(0);
      throw("nla");
    }

  // interpolate just the bin contents

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

}


template<class Hist> inline void
csm_interpolate_histogram3(Hist* a, Double_t xa, 
			   Hist* b, Double_t xb,
			   Hist* d,
			   Hist* c, Double_t xc,
			   INTERPSTYLE istyle)
{
  Double_t hnorma,hnormb,hnormc,hnormd,hnormci;
  Int_t i,j;
  Double_t gbc;

  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (a->Integral()<=0 || b->Integral()<=0)
    { 
      for (i=1;i<=nbinsc;i++)
	{
	   for (j=1;j<=nbinsyc;j++)
	     { c->SetBinContent(i,j,0);
	     }
	} 
      //c->Reset();
      return;
    }
    
  if (nbinsya == 1)
    {
      Double_t *dista = new Double_t[nbinsa];
      Double_t *distb = new Double_t[nbinsb];
      Double_t *distc = new Double_t[nbinsc];
      Double_t *distd = new Double_t[nbinsd];

      hnorma = 0;
      hnormb = 0;
      hnormd = 0;
      for (i=0;i<nbinsa;i++)
        { dista[i] = a->GetBinContent(i+1); hnorma += dista[i]; }
      for (i=0;i<nbinsb;i++)
        { distb[i] = b->GetBinContent(i+1); hnormb += distb[i]; }
      for (i=0;i<nbinsd;i++)
        { distd[i] = d->GetBinContent(i+1); hnormd += distd[i]; }

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);
      // linearly interpolate the normalization between the central value and
      // the varied template.
      hnormc = hnormd*(hnormc/hnorma); // scale the normalization with the new template

      if (istyle == CSM_INTERP_HORIZONTAL || istyle == CSM_INTERP_HORIZONTAL_EXTRAP)
	{
           csm_pvmc(nbinsa,dista,distb,distd,distc,xa,xb,xc);
           hnormci = 0;
           for (i=0;i<nbinsc;i++) { hnormci += distc[i]; }

           for (i=0;i<nbinsc;i++)
           {
             c->SetBinContent(i+1,distc[i]*hnormc/hnormci);
           }
	}
      else if (istyle == CSM_INTERP_VERTICAL || istyle == CSM_INTERP_VERTICAL_EXTRAP)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = distd[i] + ((xc-xa)/(xb-xa))*(distb[i]-dista[i]);
	      if (gbc < 0) 
		{
		  gbc = 0;
		}

	      c->SetBinContent(i+1,gbc);
	    }
	}
      else
	{
	  std::cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << std::endl;
	  throw("nla");
	  //	  exit(0);
	}

      //std::cout << xa << " " << xb << " " << xc << std::endl;

      delete[] dista;
      delete[] distb;
      delete[] distc;
      delete[] distd;
    }
  else         // 2d case
    {
      Double_t *distxya = new Double_t[nbinsa*nbinsya];
      Double_t *distxyb = new Double_t[nbinsb*nbinsyb];
      Double_t *distxyc = new Double_t[nbinsc*nbinsyc];
      Double_t *distxyd = new Double_t[nbinsd*nbinsyd];

      hnorma = 0;
      for (j=0;j<nbinsya;j++)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = a->GetBinContent(i+1,j+1);
	      distxya[i+nbinsa*j] = gbc;
	      hnorma += gbc;
	    }
	}

      hnormb = 0;
      for (j=0;j<nbinsyb;j++)
	{
	  for (i=0;i<nbinsb;i++)
	    {
	      gbc = b->GetBinContent(i+1,j+1);
	      distxyb[i+nbinsb*j] = gbc;
	      hnormb += gbc;
	    }
	}

      hnormd = 0;
      for (j=0;j<nbinsyd;j++)
	{
	  for (i=0;i<nbinsd;i++)
	    {
	      gbc = d->GetBinContent(i+1,j+1);
	      distxyd[i+nbinsb*j] = gbc;
	      hnormd += gbc;
	    }
	}

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);
      // linearly interpolate the normalization between the central value and
      // the varied template.
      hnormc = hnormd*(hnormc/hnorma); // scale the normalization with the new template

      if (istyle == CSM_INTERP_HORIZONTAL || istyle == CSM_INTERP_HORIZONTAL_EXTRAP)
	{
          csm_pvmc2d(nbinsa,nbinsya,
                     distxya,
                     distxyb,
                     distxyd,
		     distxyc,
                     xa, xb, xc);

          hnormci = 0;
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          hnormci += distxyc[i+nbinsc*j];
	        }
	    }
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          c->SetBinContent(i+1,j+1,distxyc[i+nbinsc*j]*hnormc/hnormci);
	        }
	    }
	}
      else if (istyle == CSM_INTERP_VERTICAL || istyle == CSM_INTERP_VERTICAL_EXTRAP)
	{
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
		  gbc = distxyd[i+nbinsc*j] + ((xc-xa)/(xb-xa))*(distxyb[i+nbinsc*j]-distxya[i+nbinsc*j]);
		  if (gbc < 0)
		    {
		      gbc = 0;
		    }
	          c->SetBinContent(i+1,j+1,gbc);
	        }
	    }
	}
      else
	{
	  std::cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << std::endl;
	  throw("nla");
	  //	  exit(0);
	}

      delete[] distxya;
      delete[] distxyb;
      delete[] distxyc;
      delete[] distxyd;
    }

}

#endif
