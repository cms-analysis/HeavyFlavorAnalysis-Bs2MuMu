/* This is a set of routines that works like mclimit.f but uses csm.c to compute
   the chisquared of the data histogram compared to a sum of models, each of which
   may (or may not) have sensitivities to nuisance parameters in their shapes,
   normalizations, and even statistical uncertainty in each bin */

//    25 May 2007 fix MC statistical error problem in Bayesian limit calc.
//    30 May 2007 add a version printout method to the mclimit_csm class
//                add a histogram of ln(1+s/b) in the spirit of plotwithdata, but which
//                combines all channels together.
//    30 July 2007 Work around a problem in Joel's genlimit involving double-precision
//                underflow computing small likelihood functions
//    26 Sep 2007  Add an access method for the MINUIT covariance matrix
//                 put in protection against zero background in a systematic variation in
//                 the Bayesian limit calculators
//     2 Oct 2007 Update the coverage checker to allow for an arbitrary "true" signal rate for
//                 which we'd like to check the coverage.
//     3 Oct 2007 Found another place where the max minuit calls was set to 500 -- raise to 5000
//     5 Oct 2007 Put in access methods for setting max minuit calls, printout switches for pseudoexperiments
//                and whether to call MINOS or not.
//    23 Oct 2007 Put in access methods to control MINUIT's initial step size
//     8 Nov 2007 Fix memory leak -- fitcov not deleted after csm is deleted.  Other fixes from valgrind output --
//                (output should not change). thanks to Kevin Lannon!
//    28 Nov 2007 Scale horizontally interpolated histograms like the vertically interpolated ones.
//     8 Dec 2007 Add a new MC statistical error mode which pays attention to the error bars supplied with
//                each template.  Treatment is approximate -- meant to ease the cases where MC contributions
//                come from a sample of differently-weighted events.
//     9 Dec 2007 Minor change to speed up the new histogram interpolation with errors.
//                Avoid cloning TH1's as that's a very slow operation.
//    11 Dec 2007 Change to formatting of pseudoexperiment printout in the CLs calculation.  MINOS
//                has some unavoidable printf statements which can corrupt the printout, so label
//                the pseudoexperiment data and make sure they both go on the same line.
//    12 Dec 2007 Small changes to the error in each bin handling -- bugfix to option=2
//    18 Dec 2007 Another bugfix to option=2 -- swapping contents and errors was buggy if two of the
//                histogram pointers were identical
//    24 Dec 2007 Speedup -- make a histogram interpolator that does not interpolate
//                errors for use with options 0 and 1 for poissflag
//     7 Feb 2008 Have a separate parameter controlling MINOS's maxcalls
//    11 Feb 2008 Always call IMPROVE after MINIMIZE in the fits
//    12 Feb 2008 Protect against divide by zero in MC stat fitter if a channel's scale factor is zero
//    14 Feb 2008 Always clear bayes_posterior_points when clearing the bayes_posterior vector
//    14 Feb 2008 Add on Marginalized posterior methods for measuring cross sections bh_xsfit and
//                bh_xsfit_withexpect
//    17 Feb 2008 Minor bugfix in bh_xsfit and bh_xsfit_withexpect -- make sure not to drop the
//                largest point in the marginalized posterior when integrating it!
//    26 Feb 2008 Update the interpolator style option to allow or disallow shape extrapolations
//    27 Feb 2008 Add a feature to put bounds on nuisance parameters
//    27 Feb 2008 Fix two bugs in the 2D interpolation with error bars -- one cleared the histogram contents
//                if all errors were zero, and the other set the wrong contents.
//    28 Feb 2008 Add in gclsexpb* routines -- Tom Wright's analysis with unconstrianed fits makes for
//                a case where CLs is not a monotonic function of -2lnQ, so need to scan over all
//                possible outcomes.
//     7 Mar 2008 Add in a quicker Bayesian limit integrator
//                with a cutoff in finite steps.  Sometimes bayes_heinrich and bayes_heinrich_withexpect
//                take too long, particularly for channels with lots of candidates.
//    26 Mar 2008 Fix a minor bug in the gclsexp* routines -- it didn't deal with discrete outcomes
//                as well as it could.
//    21 Apr 2008 Make gclsexp* the default in s95aux.
//    21 Apr 2008 Add a 2D cspdf2d and a function to print it out.  No 2D limits yet -- just print them out
//                for now for further plotting and analysis.  Also wrap a pseudoexperiment loop
//                around it.  New functions:  bh_2d_scan and bh_2d_scan_expect
//    22 Apr 2008 Add to bh_2d_scan_expect the same cross section fits and integrals as in bh_2d_scan.  Also
//                add to the argument list the input cross section scale factors so that coverages can be
//                computed.
//    23 Apr 2008 Add a printout of nuisance parameters for extreme null-hypothesis pseudoexperiments
//     8 May 2008 Fix bugs in the nuisance parameter constraint equations.
//    15 May 2008 2D cross section fit bugfix -- wrong low bound chosen.
//    18 Jun 2008 Add a printout of all bins' s, b, d
//    18 Sep 2009 Add some accessors to get a hold of normalization factors applied to the Bayesian posteriors.
//                Useful when parallelizing a big computation in order to add up many posteriors
//    26 Feb 2009 setdclscn in bh_xsfit
//    20 May 2009 break up cross section fit into pieces so that we can do many sytematic samples with many
//                channels without running out of memory
//    18 Sep 2009 re-order systematic samples and integration to reduce memory usage in the Bayesian integrators --
//                add in fast TH1's from Nils
//     2 Dec 2009 Nils had put the C-style routines into the .h file but that makes people building shared objects
//                get into trouble.  Put them back in this source file
//     6 Apr 2010 Include native adaptive Markov Chain code in this code
//    27 Apr 2010 Adjust signal proposal function and also Wade's asymmetric nuisance parameter handling
//    28 Apr 2010 Adjust signal proposal function again
//     7 May 2010 Add a lookup map to systematic tables so we don't have to go digging through the name list
//                in nuisance_response
//    11 May 2010 do a find in the map instead of [] to keep from adding null vectors for parameters not
//                on the list
//    18 Nov 2010 Small off-by-one bugfix in clsbaux and clbaux -- thanks to Pasha Murat

#define MCLIMIT_CSM_VERSION_NUMBER 4.16
#define MCLIMIT_CSM_VERSION_DATE "Nov 18, 2010"

// Author:  Tom Junk, Fermilab.  trj@fnal.gov
// Contributions from Joel Heinrich, Nils Krumnack, Tom Wright, and Kevin Lannon

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include "TRandom.h"
#include "TMinuit.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "mclimit_csm.hh"
#include "TString.h"
#include "TFile.h"

using namespace std;

#define max(max_a,max_b) (((max_a)>(max_b)) ? (max_a) : (max_b))
#define min(min_a,min_b) (((min_a)>(min_b)) ? (min_b) : (min_a))

// Minuit ugliness -- data communication with the function to fit is either
// via member functions of a new class inherited from TObject (not done here),
// or in global data (at least global to this source file) storage, which is
// the method chosen here because it's easier.

static vector<TH1*> datatofit;
static vector<char*> datatofitname;
static csm_model *modeltofit;
static vector<Int_t> constrainedfitparam;
static vector<char*> npfitname;

void csm_minuit_fcn(Int_t &npar, double *gin, double &f,
		    double *par, Int_t iflag);

double csint0(double xlo,double logscale,
	      int nchan,const int nobs[],const EB chan[],
	      int ngl,const double xgl[],const double lwgl[],
	      PRIOR prior);

void csint02cut(double xlo1,double xlo2,double xhi,double logscale,
		int nchan,const int nobs[],const EB chan[],
		int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		double* int1,double* int2);

void gameansigma(double *mean,double *sigma,
		 int nchan,int nens,const int nobs[],const EB* ens);
double arcfreq(double y);
#define freq(x) (0.5*erfc(-0.707106781186547524*(x)))

void csint02(double xlo1,double xlo2,double logscale,
	     int nchan,const int nobs[],const EB chan[],
	     int ngl,const double xgl[],const double lwgl[],PRIOR prior,
	     double* int1,double* int2);

// some globals which really should be put into classes, but the Bayesian routines are
// written in C and not C++

double dlcsn=0;
double dlcsn2d=0;
double bhnorm=0;


/*----------------------------------------------------------------------------*/

// constructor

mclimit_csm::mclimit_csm()
{
  nmc = 0;
  nmc_req = 10000;
  recalctsflag = 1;
  // set null pointers to our cumulative histograms -- if we
  // don't get any from the user, don't bother filling them.
  nullnullchisquare = 0;
  nulltestchisquare = 0;
  testnullchisquare = 0;
  testtestchisquare = 0;
  // null pointers to the test statistic arrays -- need to allocate memory when
  // we know how many to do
  tss = 0;
  tsb = 0;
  wtss = 0;
  wtsb = 0;
  itss = 0;
  itsb = 0;

  bayes_interval_begin = 0;
  bayes_interval_end = 0;
  bayes_interval_step = 0;
  bayes_posterior.clear();
  bayes_posterior_points.clear();
  bayes_pseudoexperiment_limits = 0;

  minuitmaxcalls = 500;
  minosmaxcalls = 500;
  minuitstepsize = 0.1;
  minuitprintflag = 0;
  minosflag = 0;
  pxprintflag = 0;
  bayesintegralmethod = CSM_BAYESINTEGRAL_JOEL;

  extremenpprintflag = 0;
  extremenpprintvalue = 0;
}

/*----------------------------------------------------------------------------*/

// destructor

mclimit_csm::~mclimit_csm()
{
  Int_t i;

  // deallocate cloned input data histograms and their names.

  for (i=0; i<(Int_t) datahist.size(); i++)
    {
      delete datahist[i];
      delete[] dhname[i]; 
    }
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::print_version()
{
  cout << "Version information for mclimit_csm.C: " << endl;
  cout << "Version number: " << MCLIMIT_CSM_VERSION_NUMBER << endl;
  cout << "Version date:   " << MCLIMIT_CSM_VERSION_DATE << endl;
}

void mclimit_csm::setminuitmaxcalls(Int_t maxcalls)
{
  minuitmaxcalls = maxcalls;
}
Int_t mclimit_csm::getminuitmaxcalls()
{
  return(minuitmaxcalls);
}

void mclimit_csm::setminosmaxcalls(Int_t maxcalls)
{
  minosmaxcalls = maxcalls;
}
Int_t mclimit_csm::getminosmaxcalls()
{
  return(minosmaxcalls);
}

void mclimit_csm::setminuitstepsize(Double_t stepsize)
{
  minuitstepsize = stepsize;
}
Double_t mclimit_csm::getminuitstepsize()
{
  return(minuitstepsize);
}

void mclimit_csm::setprintflag(bool pf)
{
  minuitprintflag = pf;
}
bool mclimit_csm::getprintflag()
{
  return(minuitprintflag);
}

void mclimit_csm::setminosflag(bool mf)
{
  minosflag = mf;
}
bool mclimit_csm::getminosflag()
{
  return(minosflag);
}

void mclimit_csm::setpxprintflag(bool pf)
{
  pxprintflag = pf;
}
bool mclimit_csm::getpxprintflag()
{
  return(pxprintflag);
}

void mclimit_csm::set_bayes_integration_method(int imethod)
{
  bayesintegralmethod = imethod;
}

int mclimit_csm::get_bayes_integration_method()
{
  return(bayesintegralmethod);
}

/*----------------------------------------------------------------------------*/
/* Build the list of channel data histograms and channel names that is sorted by channel */
/* name during the building process.  Use the same sorting procedure as used in */
/* csm_model::lookup_add_channame so that our data histograms are stored in */
/* the same order as our models */

void mclimit_csm::set_datahist(TH1 *h, char *cname)
{
  Int_t i,ifound,j,jfound;
  char *s;
  vector<char*>::iterator nhi;
  vector<TH1*>::iterator dhi;

  recalctsflag = 1;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) dhname.size(); i++)
    {
      j = (Int_t) strcmp(cname,dhname[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings.  If the name is on the 
     list, replace the existing data histogram with a clone of the one supplied. */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      strcpy(s,cname);
      if (jfound == -1)
	{
          dhname.push_back(s);
          datahist.push_back(copy_TH1<TH1> (*h).release());
	}
      else
	{
	  nhi = dhname.begin() + jfound;
	  dhname.insert(nhi,s);
	  dhi = datahist.begin() + jfound;
	  datahist.insert(dhi, copy_TH1<TH1> (*h).release());
	}
    }
  else
    {
      delete datahist[ifound];
      datahist[ifound] = copy_TH1<TH1> (*h).release();
    }
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_npe(Int_t nperequest)
{
  if (nperequest < 0)
    {
      cout << "mclimit_csm::set_npe: Invalid pseudoexperiment request: " << nperequest << endl;
      exit(0);
    }
  nmc_req = nperequest;
}

/*----------------------------------------------------------------------------*/

Int_t mclimit_csm::get_npe()
{
  return(nmc_req);
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_null_hypothesis(csm_model *model)
{
  null_hypothesis = model;
  recalctsflag = 1;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_test_hypothesis(csm_model *model)
{
  test_hypothesis = model;
  recalctsflag = 1;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_null_hypothesis_pe(csm_model *model)
{
  null_hypothesis_pe = model;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_test_hypothesis_pe(csm_model *model)
{
  test_hypothesis_pe = model;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

// test statistic of observed data

Double_t mclimit_csm::ts()
{
  Int_t i;

  if (recalctsflag) 
    {
      // copy the data histogram pointers into a flat array
      // for the chisquare calculator.
      TH1** darray = new TH1*[datahist.size()];
      for (i=0;i<(Int_t)datahist.size();i++)
	{
	  darray[i] = datahist[i];
	}
      tsd = calc_chi2(test_hypothesis,darray) -
	calc_chi2(null_hypothesis,darray);
      recalctsflag = 0;
      delete[] darray;
    }
  return(tsd);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbm2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbm1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbmed()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbp1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbp2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssm2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssm1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssmed()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssp1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssp2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/
// confidence levels for an arbitrary test statistic.  Used for computing
// actual and expected confidence levels.
// may need to work on this a bit because of finite MC statistics -- it's hard to
// compute the confidence levels on the tails unless the histograms are properly
// filled out there.  Particularly values of clb very close to 1 need to be computed
// with extreme care.
// This is addressed with the reweighting procedure suggested by Alex Read --
// the likelihood ratio can be used to reweight test hypothesis pseudoexperiments
// to model the null hypotheis background distribution, and vice versa.

Double_t mclimit_csm::clsbaux(Double_t tsaux)
{
  Int_t i;
  Double_t clsbloc;
  if (nmc == 0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clsbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tss[itss[i]] < tsaux)
	{ 
	  clsbloc = ((Double_t) (i+1))/((Double_t) nmc);
	}
      else break;
    }
  return(1-clsbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbauxw(Double_t tsaux)
{
  Int_t i;
  Double_t clsbloc;
  if (nmc == 0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clsbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] >= tsaux)
	{ 
          clsbloc += wtsb[itsb[i]];
	}
    }
  clsbloc /= ((Double_t) nmc);
  return(clsbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbaux(Double_t tsaux)
{
  Int_t i;
  Double_t clbloc;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] < tsaux)
	{ 
	  clbloc = ((Double_t) (i+1))/((Double_t) nmc);
	}
      else break;
    }
  return(1-clbloc);
}

/*----------------------------------------------------------------------------*/

// compute 1-CLb (as a p-value, including the outcomes with exactly the -2lnQ
// observed) using null hypothesis px's.

Double_t mclimit_csm::omclbaux(Double_t tsaux)
{
  Int_t i;
  Double_t omclbloc;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  omclbloc = 1;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] <= tsaux)
	{ 
	  omclbloc = ((Double_t) (i+1))/((Double_t) nmc);
	}
      else break;
    }
  return(omclbloc);
}

/*----------------------------------------------------------------------------*/

/*  Compute 1-CLb using reweighted test hypothesis pseudoexperiments.  Reweight
    using the inverese of the likelihood ratio, 1/(p(data|test)/p(data|null)) */

Double_t mclimit_csm::omclbauxw(Double_t tsaux)
{
  Int_t i;
  Double_t omclbloc = 0;

  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  omclbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tss[itss[i]] <= tsaux)
	{ 
          omclbloc += wtss[itss[i]];
	}
    }
  omclbloc /= ((Double_t) nmc);
  return(omclbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsaux(Double_t tsaux)
{
  Double_t clbloc,clsloc;
  clbloc = clbaux(tsaux);
  if (clbloc > 0)
    {
      clsloc = clsbaux(tsaux)/clbloc;
    }
  else
    {
      clsloc = 1;
    }
  return(clsloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsauxw(Double_t tsaux)
{
  Double_t clbloc,clsloc;
  clbloc = clbaux(tsaux);
  if (clbloc > 0)
    {
      clsloc = clsbauxw(tsaux)/clbloc;
    }
  else
    {
      clsloc = 1;
    }
  return(clsloc);
}

/*----------------------------------------------------------------------------*/
// confidence levels using the data test statistic.  Recompute the data test
// statistic if need be.

Double_t mclimit_csm::cls()
{
  return(clsaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsb()
{
  return(clsbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsw()
{
  return(clsauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbw()
{
  return(clsbauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clb()
{
  return(clbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclb()
{
  return(omclbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbw()
{
  return(omclbauxw(ts()));
}

/*----------------------------------------------------------------------------*/
// keep the convention of the other clsexp routines

Double_t mclimit_csm::gclsexpbm2()  { return(gclsaux(MCLIMIT_CSM_MCLP2S)); }
Double_t mclimit_csm::gclsexpbm1()  { return(gclsaux(MCLIMIT_CSM_MCLP1S)); }
Double_t mclimit_csm::gclsexpbmed() { return(gclsaux(MCLIMIT_CSM_MCLMED)); }
Double_t mclimit_csm::gclsexpbp1()  { return(gclsaux(MCLIMIT_CSM_MCLM1S)); }
Double_t mclimit_csm::gclsexpbp2()  { return(gclsaux(MCLIMIT_CSM_MCLM2S)); }

Double_t mclimit_csm::gclsaux(Double_t thresh)
{

  vector<Double_t> clslist;
  if (nmc == 0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  // get CLs for each background outcome, sort them (they could be out of order as a function of -2lnQ)
  // tss and tsb are already sorted in ascending order, using the sort indices
  // fix up by counting duplicates properly

  int k=0;
  int j=0;
  for (int i=0;i<nmc;i++)
    {
      double tstmp = tsb[itsb[i]];
      while (tsb[itsb[k+1]] < tstmp && k<nmc-2) k++;
      while (tss[itss[j+1]] < tstmp && j<nmc-2) j++;
      clslist.push_back(  ((Double_t) (nmc-j)) / ((Double_t) (nmc-k)) );
    }
  std::sort(clslist.begin(),clslist.end());
  int i =  (int) nearbyint( ((Double_t) nmc)*thresh);
  return(clslist[i]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm2()
{
  return(clsaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm1()
{
  return(clsaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbmed()
{
  return(clsaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp1()
{
  return(clsaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp2()
{
  return(clsaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm2w()
{
  return(clsauxw(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm1w()
{
  return(clsauxw(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbmedw()
{
  return(clsauxw(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp1w()
{
  return(clsauxw(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp2w()
{
  return(clsauxw(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbm2()
{
  return(clsbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbm1()
{
  return(clsbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbmed()
{
  return(clsbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbp1()
{
  return(clsbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbp2()
{
  return(clsbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsm2()
{
  return(clsaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsm1()
{
  return(clsaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsmed()
{
  return(clsaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsp1()
{
  return(clsaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsp2()
{
  return(clsaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbm2()
{
  return(clbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbm1()
{
  return(clbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbmed()
{
  return(clbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbp1()
{
  return(clbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbp2()
{
  return(clbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsm2()
{
  return(clbaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsm1()
{
  return(clbaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsmed()
{
  return(clbaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsp1()
{
  return(clbaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsp2()
{
  return(clbaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm2()
{
  return(omclbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm1()
{
  return(omclbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbmed()
{
  return(omclbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp1()
{
  return(omclbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp2()
{
  return(omclbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm2()
{
  return(omclbaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm1()
{
  return(omclbaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsmed()
{
  return(omclbaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp1()
{
  return(omclbaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp2()
{
  return(omclbaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm2w()
{
  return(omclbauxw(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm1w()
{
  return(omclbauxw(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbmedw()
{
  return(omclbauxw(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp1w()
{
  return(omclbauxw(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp2w()
{
  return(omclbauxw(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm2w()
{
  return(omclbauxw(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm1w()
{
  return(omclbauxw(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsmedw()
{
  return(omclbauxw(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp1w()
{
  return(omclbauxw(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp2w()
{
  return(omclbauxw(tssp2()));
}



/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p2sigmat()
{
  Int_t i;
  Double_t p2s;
  p2s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p2sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_2S)
	{
	  p2s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p2s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p3sigmat()
{
  Int_t i;
  Double_t p3s;
  p3s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p3sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_3S)
	{
	  p3s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p3s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p5sigmat()
{
  Int_t i;
  Double_t p5s;
  p5s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p5sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_5S)
	{
	  p5s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p5s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p2sigman()
{
  Int_t i;
  Double_t p2s;
  p2s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p2sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_2S)
	{
	  p2s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p2s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p3sigman()
{
  Int_t i;
  Double_t p3s;
  p3s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p3sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_3S)
	{
	  p3s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p3s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p5sigman()
{
  Int_t i;
  Double_t p5s;
  p5s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p5sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_5S)
	{
	  p5s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p5s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95()
{
  return(s95aux(MCLIMIT_CSM_CLS));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95m2()
{
  return(s95aux(MCLIMIT_CSM_CLSM2));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95m1()
{
  return(s95aux(MCLIMIT_CSM_CLSM1));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95med()
{
  return(s95aux(MCLIMIT_CSM_CLSMED));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95p1()
{
  return(s95aux(MCLIMIT_CSM_CLSP1));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95p2()
{
  return(s95aux(MCLIMIT_CSM_CLSP2));
}

Double_t mclimit_csm::s95aux(Int_t itype)
{
  Double_t sf,sfl,sfh,cltest,cll,cla,clh;
  Int_t j,foundit;
  csm_model *testhypsave,*testhyppesave;

  sfl = 0;
  sfh = 0;
  cltest = 0;

  // this hypothesis never gets destroyed, but we lose the pointer to it

  testhypsave = test_hypothesis;
  testhyppesave = test_hypothesis_pe;

  sf = 1.0;
  cll = 0.0;
  cla = 0.0;
  clh = 0.0; 
  foundit = 0;

  for (j=0;j<32;j++)
    {
      //cout << "in s95aux seek " << j << " scale: " << sf << endl;
      if (foundit == 0)
	{
	  csm_model* scaledsignal = testhypsave->scalesignal(sf);
	  csm_model* scaledsignalpe = testhyppesave->scalesignal(sf);
	  test_hypothesis = scaledsignal;
	  test_hypothesis_pe = scaledsignalpe;
	  run_pseudoexperiments();
	  recalctsflag = 1;
	  if (itype == MCLIMIT_CSM_CLS)
	    {
	      cltest = cls(); 
	    }
	  else if (itype == MCLIMIT_CSM_CLSM2)
	    {
	      cltest = gclsexpbm2();
	    }
	  else if (itype == MCLIMIT_CSM_CLSM1)
	    {
	      cltest = gclsexpbm1();
	    }
	  else if (itype == MCLIMIT_CSM_CLSMED)
	    {
	      cltest = gclsexpbmed();
	    }
	  else if (itype == MCLIMIT_CSM_CLSP1)
	    {
	      cltest = gclsexpbp1();
	    }
	  else if (itype == MCLIMIT_CSM_CLSP2)
	    {
	      cltest = gclsexpbp2();
	    }
	  delete scaledsignal;
	  delete scaledsignalpe;
	  if (j==0)
	    {
	      cla = cltest;
	    }
	  if (cltest<0.05)
	    {
	      if (cla>0.05)
		{
		  sfh = sf;
		  clh = cltest;
		  sfl = sf/2.0;
		  cll = cla;
		  foundit = 1;
		}
	      sf /= 2.0;
	    }
	  else if (cltest>0.05)
	    {
	      if (cla<0.05)
		{
		  sfl = sf;
		  cll = cltest;
		  sfh = sf*2.0;
		  clh = cla;
		  foundit = 1;
		}
	      sf *= 2.0;
	    }
	  else
	    {
	      test_hypothesis = testhypsave;
	      test_hypothesis_pe = testhyppesave;
	      return(sf);
	    }
	  cla = cltest;
	}
    } // end of loop over 32 powers of 2 in search of a signal scale factor which brackets
  // 95% CL exclusion

  //cout << "done with seek loop " << sf << endl;
  sf = sfh;
  if (foundit == 0)
    { 
      cout << "mclimit_csm::s95** could not find s95 within 2**32 of original guess" << endl;
      sf = 0;
    }
  else
    {
      // From Tom Wright -- speed up by doing a deterministic five more
      // calcs of CL and a linear fit of log(CL) vs. sf.

      // find error on 0.05 CL for number of PEs
      double dcl=sqrt(0.05*0.95/nmc);

      // put in some protection against logarithms of negative numbers
      // makes sure -5*dcl + 0.05 is not negative.

      dcl = min(dcl,0.0099);

      // try +6sigma, +3sigma, 0sigma, -3sigma, -6sigma regions
      // increment stuff used for linear fit of ln(CL) vs sf
      double lf_a=0, lf_b=0, lf_c=0, lf_d=0, lf_e=0, lf_f=0;
      for( int j=-5; j<6; j+=2 )
	{
	  sf = sfl + (log(0.05+j*dcl) - log(cll))*
	    (sfl-sfh)/(log(cll)-log(clh));

	  double clsbtest=0;
	  double clbtest=0;

	  // calculate CL for this sf
	  csm_model* scaledsignal = testhypsave->scalesignal(sf);
          csm_model* scaledsignalpe = testhyppesave->scalesignal(sf);
          test_hypothesis = scaledsignal;
          test_hypothesis_pe = scaledsignalpe;
          run_pseudoexperiments();
          recalctsflag = 1;
          if (itype == MCLIMIT_CSM_CLS)
	    {
	      cltest = cls();
	      clbtest = clb();
	      clsbtest = clsb();
  	    }
          else if (itype == MCLIMIT_CSM_CLSM2)
	    {
	      cltest = clsexpbm2();
	      clbtest = clbexpbm2();
	      clsbtest = clsbexpbm2();
	    }
          else if (itype == MCLIMIT_CSM_CLSM1)
	    {
	      cltest = clsexpbm1();
	      clbtest = clbexpbm1();
	      clsbtest = clsbexpbm1();
	    }
          else if (itype == MCLIMIT_CSM_CLSMED)
	    {
	      cltest = clsexpbmed();
	      clbtest = clbexpbmed();
	      clsbtest = clsbexpbmed();
	    }
          else if (itype == MCLIMIT_CSM_CLSP1)
	    {
	      cltest = clsexpbp1();
	      clbtest = clbexpbp1();
	      clsbtest = clsbexpbp1();
	    }
          else if (itype == MCLIMIT_CSM_CLSP2)
	    {
	      cltest = clsexpbp2();
	      clbtest = clbexpbp2();
	      clsbtest = clsbexpbp2();
	    }
	  delete scaledsignal;
	  delete scaledsignalpe;

	  // double dcltest=sqrt(cltest*(1-cltest)/nmc);
	  // 	  double dcltest = cltest*sqrt((1-clbtest)/clbtest/nmc +
	  // 				       (1-clsbtest)/clsbtest/nmc);
	  double dcltest=sqrt(clsbtest*(1-clsbtest)/nmc)/clbtest;

          //  printf("%f %f %f %f %f\n",sf,clbtest,clsbtest,cltest,dcltest);

	  double lcl = log(cltest);
	  double dlcl = dcltest/cltest;

	  lf_a += sf/dlcl/dlcl;
	  lf_b += 1/dlcl/dlcl;
	  lf_c += lcl/dlcl/dlcl;
	  lf_d += sf*sf/dlcl/dlcl;
	  lf_e += sf*lcl/dlcl/dlcl;
	  lf_f += lcl*lcl/dlcl/dlcl;
	}

      // Find fit parameters for log(CL)=p1+p2*sf
      double lf_p1 = (lf_d*lf_c-lf_e*lf_a)/(lf_d*lf_b-lf_a*lf_a);
      double lf_p2 = (lf_e*lf_b-lf_c*lf_a)/(lf_d*lf_b-lf_a*lf_a);

      //double lf_dp1 = sqrt(lf_d/(lf_b*lf_d-lf_a*lf_a));
      //double lf_dp2 = sqrt(lf_b/(lf_b*lf_d-lf_a*lf_a));
      //double lf_rho = -lf_a/(lf_b*lf_d-lf_a*lf_a)/lf_dp1/lf_dp2;

      //printf("fit results %f %f %f %f %f\n",lf_p1,lf_dp1,lf_p2,lf_dp2,lf_rho);

      //double lf_x2 = lf_f-2*lf_p2*lf_e-2*lf_p1*lf_c+lf_p2*lf_p2*lf_d+
      //	2*lf_p1*lf_p2*lf_a+lf_p1*lf_p1*lf_b;
      //printf("chisuare/dof: %f\n",lf_x2/4);

      // invert to get sf at 0.05 and its error
      // assuming 100% anticorrelation
      double lf_sf = (log(0.05)-lf_p1)/lf_p2;
      //printf("CL variation at %f: %f %f %f\n",lf_sf,
      //     exp(lf_p1-lf_dp1+(lf_p2+lf_dp2)*lf_sf),
      //     exp(lf_p1+lf_p2*lf_sf),
      //     exp(lf_p1+lf_dp1+(lf_p2-lf_dp2)*lf_sf));

      //double lf_dsf1 = (log(0.05)-lf_p1-lf_dp1)/(lf_p2-lf_dp2);
      //double lf_dsf2 = (log(0.05)-lf_p1+lf_dp1)/(lf_p2+lf_dp2);

      //printf("SF variation: %f %f %f\n",lf_dsf1,lf_sf,lf_dsf2);
      sf = lf_sf;

    }

  test_hypothesis = testhypsave;
  test_hypothesis_pe = testhyppesave;
  recalctsflag = 1;
  return(sf);
}


/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi95()
{
  return(lumipaux(MCLIMIT_CSM_LUMI95));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi3s()
{
  return(lumipaux(MCLIMIT_CSM_LUMI3S));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi5s()
{
  return(lumipaux(MCLIMIT_CSM_LUMI5S));
}

/*----------------------------------------------------------------------------*/
// compute median amounts of luminosity needed for 95% CL exclusion, 3 sigma
// evidence, or 5 sigma discovery -- scale the systematic errors with 1/sqrt(lumi/lumi_0)

Double_t mclimit_csm::lumipaux(Int_t itype)
{
  Double_t sf,sfl,sfh,cltest,cll,cla,clh;
  Int_t foundit,j;
  csm_model *testhypsave,*nullhypsave;
  csm_model *testhyppesave,*nullhyppesave;
  Double_t resdes;

  resdes = 0.5; // do median luminosity thresholds

  sfl = 0;
  sfh = 0;
  cltest = 0;

  sf = 1.0;
  cll = 0.0;
  cla = 0.0;
  clh = 0.0; 
  foundit = 0;

  testhypsave = test_hypothesis;
  nullhypsave = null_hypothesis;
  testhyppesave = test_hypothesis_pe;
  nullhyppesave = null_hypothesis_pe;
  
  for (j=0;j<32;j++)
    {
      if (foundit == 0)
	{
	  csm_model* scaledtest = testhypsave->scale_err(sf);
	  csm_model* scalednull = nullhypsave->scale_err(sf);
	  csm_model* scaledtestpe = testhyppesave->scale_err(sf);
	  csm_model* scalednullpe = nullhyppesave->scale_err(sf);
	  test_hypothesis = scaledtest;
	  null_hypothesis = scalednull;
	  test_hypothesis_pe = scaledtestpe;
	  null_hypothesis_pe = scalednullpe;
          run_pseudoexperiments();
          recalctsflag = 1;

          if (itype == MCLIMIT_CSM_LUMI95)
	    {
	      cltest = p2sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI3S)
	    {
	      cltest = p3sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI5S)
	    {
	      cltest = p5sigmat();
	    }
	  delete scaledtest;
	  delete scalednull;
	  delete scaledtestpe;
	  delete scalednullpe;

	  if (j==0)
	    {
	      cla = cltest;
	    }
	  if (cltest < resdes)
	    {
	      if (cla > resdes)
	        {
		  sfl = sf;
	  	  cll = cltest;
		  sfh = sf*2.0;
		  clh = cla;
		  foundit = 1;
	        }
	      sf *= 2.0;
	    }
  	  else if (cltest > resdes)
	    {
	      if (cla < resdes)
	        {
 		  sfh = sf;
		  clh = cltest;
		  sfl = sf/2.0;
		  cll = cla;
		  foundit = 1;
	        }
	      sf /= 2.0;
	    }
	  else
	    {
	      test_hypothesis = testhypsave;
	      null_hypothesis = nullhypsave;
	      test_hypothesis_pe = testhyppesave;
	      null_hypothesis_pe = nullhyppesave;
	      return(sf);
	    }
	  cla = cltest;
        }
    } // end of loop over 32 powers of 2 in search of a luminosity scale factor which
      // brackets the desired sensitvity

  sf = sfl;
  if (foundit == 0)
    { 
      cout << "mclimit_csm::lumipaux** could not find s95 within 2**32 of original guess" << endl;
      sf = 0;
    }
  else
    {

      // From Tom Wright -- speed up by doing a deterministic five more
      // calcs of CL and a linear fit of log(CL) vs. sf.

      // find error on resdes CL for number of PEs
      double dcl=sqrt(resdes*(1.0-resdes)/nmc);

      // put in some protection against logarithms of negative numbers
      // makes sure -5*dcl + resdes is not negative.

      dcl = min(dcl,resdes/5 - 0.0001);

      // try +6sigma, +3sigma, 0sigma, -3sigma, -6sigma regions
      // increment stuff used for linear fit of ln(CL) vs sf
      double lf_a=0, lf_b=0, lf_c=0, lf_d=0, lf_e=0, lf_f=0;
      for( int j=-5; j<6; j+=2 )
	{
	  sf = sfl + (log(resdes+j*dcl) - log(cll))*
	    (sfl-sfh)/(log(cll)-log(clh));

	  //double clsbtest, clbtest;

	  // calculate CL for this sf
	  csm_model* scaledtest = testhypsave->scale_err(sf);
	  csm_model* scalednull = nullhypsave->scale_err(sf);
	  csm_model* scaledtestpe = testhyppesave->scale_err(sf);
	  csm_model* scalednullpe = nullhyppesave->scale_err(sf);
	  test_hypothesis = scaledtest;
	  null_hypothesis = scalednull;
	  test_hypothesis_pe = scaledtestpe;
	  null_hypothesis_pe = scalednullpe;
          run_pseudoexperiments();
          recalctsflag = 1;

          if (itype == MCLIMIT_CSM_LUMI95)
	    {
	      cltest = p2sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI3S)
	    {
	      cltest = p3sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI5S)
	    {
	      cltest = p5sigmat();
	    }

	  delete scaledtest;
	  delete scalednull;
	  delete scaledtestpe;
	  delete scalednullpe;

	  // double dcltest=sqrt(cltest*(1-cltest)/nmc);
	  // 	  double dcltest = cltest*sqrt((1-clbtest)/clbtest/nmc +
	  // 				       (1-clsbtest)/clsbtest/nmc);
	  double dcltest=sqrt(cltest*(1-cltest)/nmc)/cltest;

          //  printf("%f %f %f %f %f\n",sf,clbtest,clsbtest,cltest,dcltest);

	  double lcl = log(cltest);
	  double dlcl = dcltest/cltest;

	  lf_a += sf/dlcl/dlcl;
	  lf_b += 1/dlcl/dlcl;
	  lf_c += lcl/dlcl/dlcl;
	  lf_d += sf*sf/dlcl/dlcl;
	  lf_e += sf*lcl/dlcl/dlcl;
	  lf_f += lcl*lcl/dlcl/dlcl;
	}

      // Find fit parameters for log(CL)=p1+p2*sf
      double lf_p1 = (lf_d*lf_c-lf_e*lf_a)/(lf_d*lf_b-lf_a*lf_a);
      double lf_p2 = (lf_e*lf_b-lf_c*lf_a)/(lf_d*lf_b-lf_a*lf_a);

      //double lf_dp1 = sqrt(lf_d/(lf_b*lf_d-lf_a*lf_a));
      //double lf_dp2 = sqrt(lf_b/(lf_b*lf_d-lf_a*lf_a));
      //double lf_rho = -lf_a/(lf_b*lf_d-lf_a*lf_a)/lf_dp1/lf_dp2;

      //printf("fit results %f %f %f %f %f\n",lf_p1,lf_dp1,lf_p2,lf_dp2,lf_rho);

      //double lf_x2 = lf_f-2*lf_p2*lf_e-2*lf_p1*lf_c+lf_p2*lf_p2*lf_d+
      //	2*lf_p1*lf_p2*lf_a+lf_p1*lf_p1*lf_b;
      //printf("chisuare/dof: %f\n",lf_x2/4);

      // invert to get sf at resdes and its error
      // assuming 100% correlation
      double lf_sf = (log(resdes)-lf_p1)/lf_p2;
      //printf("CL variation at %f: %f %f %f\n",lf_sf,
      //     exp(lf_p1-lf_dp1+(lf_p2+lf_dp2)*lf_sf),
      //     exp(lf_p1+lf_p2*lf_sf),
      //     exp(lf_p1+lf_dp1+(lf_p2-lf_dp2)*lf_sf));

      //double lf_dsf1 = (log(resdes)-lf_p1-lf_dp1)/(lf_p2-lf_dp2);
      //double lf_dsf2 = (log(resdes)-lf_p1+lf_dp1)/(lf_p2+lf_dp2);

      //printf("SF variation: %f %f %f\n",lf_dsf1,lf_sf,lf_dsf2);
      sf = lf_sf;

    }
   
  test_hypothesis = testhypsave;
  null_hypothesis = nullhypsave;
  test_hypothesis_pe = testhyppesave;
  null_hypothesis_pe = nullhyppesave;
  return(sf);
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::tshists(TH1* testhypts, TH1* nullhypts)
{
  int i;
  testhypts->Reset();
  nullhypts->Reset();
  for (i=0;i<nmc;i++)
    { 
      testhypts->Fill(tss[i]);
      nullhypts->Fill(tsb[i]);
    }
}

/*----------------------------------------------------------------------------*/
// makes and overlaid plot of ln(1+s/b) in the user's binning
// assumes a canvas and pad are already set up and plot options are set up

void mclimit_csm::plotlnsb(TH1 *mcb_hist, TH1 *mcs_hist, TH1 *data_hist)
{
  Int_t i,ibinx,ibiny,nbinsx,nbinsy,ic,nc;
  Double_t s,b,dtb,gbc,sbln;
  csm_channel_model *cm;
  TH1 *dhp;

  mcb_hist->Reset();
  mcs_hist->Reset();
  data_hist->Reset();

  for (i=0;i<(Int_t)(test_hypothesis_pe->chanmodel.size());i++)
    {
      // dhp is the data histogram we're using to fill in the ln(1+s/b) histo
      // cm is the channel model for this data histogram
      dhp = datahist[i];
      cm = test_hypothesis_pe->chanmodel[i];
      nc = (Int_t) cm->histotemplate.size();
      nbinsx = dhp->GetNbinsX();
      nbinsy = dhp->GetNbinsY();
      for (ibinx=0;ibinx<nbinsx;ibinx++)
	{
	  for (ibiny=0;ibiny<nbinsy;ibiny++)
	    {
	      s = 0;
	      b = 0;
	      if (nbinsy == 1)
		{ dtb = dhp->GetBinContent(ibinx+1); }
	      else
		{ dtb = dhp->GetBinContent(ibinx+1,ibiny+1); }
	      for (ic=0;ic<nc;ic++)
		{
		  if (nbinsy == 1)
		    { gbc = cm->histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { gbc = cm->histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
		  gbc *= cm->sft_varied[ic];
		  if (cm->scaleflag[ic])
		    { s += gbc; }
		  else
		    { b += gbc; }
		}
	      if (b>0 && s>= 0) 
		{ 
		  sbln = log(1.0+s/b);
		  mcs_hist->Fill(sbln,s);
		  mcb_hist->Fill(sbln,b);
		  data_hist->Fill(sbln,dtb);
		}
	    }
	}
    }

  THStack *hs = new THStack("lnsbstack",data_hist->GetTitle());
  mcs_hist->SetFillColor(2);
  mcb_hist->SetFillColor(3);
  hs->Add(mcb_hist);
  hs->Add(mcs_hist);
  data_hist->GetXaxis()->SetTitle("ln(1+s/b)");
  TLegend *slegend = new TLegend(0.7,0.6,0.89,0.89);
  slegend->AddEntry(mcs_hist,"Signal","F");
  slegend->AddEntry(mcb_hist,"Background","F");
  slegend->AddEntry(data_hist,data_hist->GetName(),"P");
  //slegend->SetHeader(data_hist->GetTitle());

  Double_t mcmax,datamax,plotmax;
  mcmax = hs->GetMaximum();
  datamax = data_hist->GetMaximum();
  datamax += sqrt(datamax);
  plotmax = max(datamax,mcmax);
  hs->SetMaximum(plotmax);
  data_hist->SetMaximum(plotmax);
  hs->Draw("HIST");
  data_hist->SetMarkerStyle(20);
  data_hist->SetMarkerColor(kBlack);
  data_hist->DrawCopy("E0SAME");
  slegend->Draw();
}


/*----------------------------------------------------------------------------*/

void mclimit_csm::set_chisquarehistos(TH1 *nn, TH1 *nt, TH1 *tn, TH1 *tt)
{
  nullnullchisquare = nn;
  nulltestchisquare = nt;
  testnullchisquare = tn;
  testtestchisquare = tt;

}

// steer printing of nuisance parameters for extreme pseudoexperiments

void mclimit_csm::setprintnpextremeflag(bool prf)
{
  extremenpprintflag = prf;
}

void mclimit_csm::setprintextremevalue(Double_t pexval)
{
  extremenpprintvalue = pexval;
}

/*----------------------------------------------------------------------------*/
// run pseudoexperiments, allocating and filling the tss and tsb arrays
// also fill in wtss and wtsb, for use in reweighting pseudoexperiments, using
// the varied nuisance parameters.  Sort tss and tsb, and keep the sort order
// arrays itss and itsb around so the corresponding wtss and wtsb arrays can
// be used with them

void mclimit_csm::run_pseudoexperiments()
{
  Int_t i;
  char *pdname;
  // Double_t tmp;
  Double_t csnull,cstest;

  // make some histograms to store the pseudodata.

  TH1** pdarray = new TH1*[null_hypothesis_pe->channame.size()];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      pdname = new char[strlen(null_hypothesis_pe->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,null_hypothesis_pe->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = copy_TH1<TH1> (*null_hypothesis_pe->chanmodel[i]->histotemplate[0], pdname).release();
      delete[] pdname;
    }

  // allocate memory for test statistic and weight and sort order storage
  if (tss != 0) delete[] tss;
  if (tsb != 0) delete[] tsb;
  if (wtss != 0) delete[] wtss;
  if (wtsb != 0) delete[] wtsb;
  if (itss != 0) delete[] itss;
  if (itsb != 0) delete[] itsb;
  tss = new Double_t[nmc_req];
  tsb = new Double_t[nmc_req];
  wtss = new Double_t[nmc_req];
  wtsb = new Double_t[nmc_req];
  itss = new Int_t[nmc_req];
  itsb = new Int_t[nmc_req];

  for(i=0;i<nmc_req;i++)
    {
      null_hypothesis_pe->single_pseudoexperiment(pdarray);
      wtsb[i] = weightratio(test_hypothesis_pe,null_hypothesis_pe,pdarray);
      csnull = calc_chi2(null_hypothesis,pdarray);
      cstest = calc_chi2(test_hypothesis,pdarray);
      if (nullnullchisquare != 0)
	{ nullnullchisquare->Fill(csnull); }
      if (nulltestchisquare != 0)
	{ nulltestchisquare->Fill(cstest); }
      tsb[i] = cstest - csnull;
      // cout << "null hyp chisquared: " << csnull << " test hyp chisquared: " << cstest << endl;
      // diagnostic code below to print out all the nuisance parameter values in the case that the null hypothesis
      // has fluctuated into a very signal-like region of -2lnQ.
      if (extremenpprintflag && tsb[i]<extremenpprintvalue)
	{
	  cout << "Extreme value of -2lnQ for a null hypothesis pseudoexperiment: " << tsb[i] << " < " << extremenpprintvalue << endl;
	  null_hypothesis_pe->print_nuisance_params();
	}

      test_hypothesis_pe->single_pseudoexperiment(pdarray);
      wtss[i] = weightratio(null_hypothesis_pe,test_hypothesis_pe,pdarray);
      csnull = calc_chi2(null_hypothesis,pdarray);
      cstest = calc_chi2(test_hypothesis,pdarray);
      if (testnullchisquare != 0)
	{ testnullchisquare->Fill(csnull); }
      if (testtestchisquare != 0)
	{ testtestchisquare->Fill(cstest); }
      tss[i] = cstest - csnull;
      // cout << "null hyp chisquared: " << csnull << " test hyp chisquared: " << cstest << endl;

      if (pxprintflag)
	{
          cout << " Null, Test hyp px: " << tsb[i] << " " << tss[i] << endl;
	}
    }

  TMath::Sort(nmc_req,tss,itss,0);
  TMath::Sort(nmc_req,tsb,itsb,0);

  for (i=0;i<(Int_t)null_hypothesis_pe->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  nmc = nmc_req;
}

/*----------------------------------------------------------------------------*/
// update the internal representation of the nuisance response of the model.
// Use the method for doing this with each
// channel separately.

void csm_model::nuisance_response(Int_t nparams,
                                  char *paramname[],
                                  Double_t paramvalue[])
{
  Int_t i,nchans;
  Int_t ipar,icons,j,k,ifound,jfound;
  Double_t *cinput;

  /*
    cout << "in model nuisance response: " << endl;
    for (i=0;i < (Int_t) nparams; i++)
    {
    cout << "param: " << i << " name: " << paramname[i] << endl;
    }
  */

  // compute the nuisance parameters which are functions of the others
  // loop over constraints and replace the parameter values with the constrained ones
  // do not overwrite the input parameters but make our own copy

  Double_t *parloc = new Double_t[nparams];
  for (i=0;i<nparams;i++)
    {
      parloc[i] = paramvalue[i];
    }
  for (icons=0;icons<(Int_t)npcm.size();icons++)
    {
      jfound = 0;
      for (ipar=0;ipar<nparams;ipar++)
	{
	  //cout << "Comparing: " << npcm[icons].pnameoutput << " with " << paramname[ipar] << endl;
          if (strcmp(npcm[icons].pnameoutput,paramname[ipar])==0)
	    {
	      jfound = 1;
	      cinput = new Double_t[npcm[icons].ninput];
	      for (j=0;j<(Int_t)npcm[icons].ninput;j++)
	        {
	          ifound = 0;
	          for (k=0;k<nparams;k++)
		    {
		      if (strcmp(npcm[icons].pnameinput[j],paramname[k])==0)
		        {
		          cinput[j] = parloc[k];
		          ifound = 1;
		        }
		    }
	          if (ifound == 0)
		    {
		      cout << "Didn't find parameter name: " << 
			npcm[icons].pnameinput[j] << 
			" in the list of nuisance parameters" << endl;
		      exit(0);
		    }
	        }
	      parloc[ipar] = npcm[icons].f(cinput);
	      delete[] cinput;
	    }
	}
      // Not a disaster if we didn't find a parameter on the list.  Sometimes we
      // only fit for a subset of parameters and they aren't on the list.
      //  if (jfound == 0) 
      //	{
      //  cout << "Constraint equation found for nuisance parameter not on our list: " 
      //       << npcm[icons].pnameoutput << endl;
      //  exit(0);
      //  }
    }

  nchans = (Int_t) channame.size();
  for (i=0;i< nchans; i++)
    {
      chanmodel[i]->nuisance_response(nparams,paramname,parloc);
    }
  delete[] parloc;
}

/*----------------------------------------------------------------------------*/

void csm_model::undo_nuisance_response()
{
  Int_t i,nchans;
  nchans = (Int_t) channame.size();
  for (i=0;i< nchans; i++)
    {
      chanmodel[i]->undo_nuisance_response();
    }
}

/*----------------------------------------------------------------------------*/

// Create a fluctuated channel model -- input a list of nuisance parameter names
// and values, and return a pointer to a new channel model which has 
// responded to those nuisance parameters.  Be sure to delete it when done.
// todo -- make sure that the shape errors accumulate.  Suggestion of John
// Zhou: average all shape error interpolations.  --  Better compounded interpolations
// introduced in Spring 2007

void csm_channel_model::nuisance_response(Int_t nparams,
                                          char *paramname[],
                                          Double_t paramvalue[])
{
  Int_t i,j,itpl,ntemplates;

  /*
    cout << "in channel nuisance response: " << endl;
    for (i=0;i < (Int_t) syserr.size(); i++)
    {
    cout << "error source: " << i << " name: " << syserr[i].sysname << endl;
    }
    for (i=0;i < (Int_t) nparams; i++)
    {
    cout << "param: " << i << " name: " << paramname[i] << endl;
    }
  */

  undo_nuisance_response();
  ntemplates = (Int_t) histotemplate.size();

  // add the rate contributions linearly as Wade does.  I like to multiply them
  // but we go for consistency

  double rsum[ntemplates];
  for (itpl=0;itpl<ntemplates;itpl++) rsum[itpl] = 0.0;

  // NK: check performance
  std::auto_ptr<TH1Type>& hcl = hcl_nuisance_response;

  std::vector<int> *vpt;
  std::map<char*, std::vector<int>, csm_ltstr>::iterator vpi;

  for (j=0; j<nparams; j++)
    {
      vpi = semap.find(paramname[j]);
      if (vpi != semap.end())
	{
          vpt = &(vpi->second);
	  int nsloc = vpt->size();
	  for (int k1=0;k1<nsloc;k1++)
	    {
	      i = (*vpt)[k1];
	      itpl = syserr[i].itemplate;

	      // use Wade's notation here -- we keep around the minus sign on the negative
	      // variations.

	      double sigmaN = -syserr[i].sysfracl;
	      double sigmaP = syserr[i].sysfrach;
	      double r = paramvalue[j];
	      double sig = sigmaN;
	      if (r>0) sig = sigmaP;
	      double quadMatch = r*(sigmaP+sigmaN)/2.0 + r*r*(sigmaP-sigmaN)/2.0;

	      // Tom's old notation
	      //double quadMatch = 
	      //	(syserr[i].sysfrach+syserr[i].sysfracl)*paramvalue[j]*paramvalue[j]/2.0 +
	      //	(syserr[i].sysfrach-syserr[i].sysfracl)*paramvalue[j]/2.0 + 1.0;

	      double rf = 1.0/(1.0+3.0*fabs(r));
	      double bridge = r*sig*(1.0-rf) + rf*quadMatch;
	      //    double lnB = exp(bridge); // cannot go below 0 by constuction
	      double lnB = 0;
	      if (bridge<0)
		{ lnB = exp(bridge); } // cannot go below 0 by constuction
	      else
		{ lnB = bridge + 1.0; } // regulate large excursions -- prevent large variation

	      // the additive version as Wade does it
	      rsum[itpl] += lnB - 1.0;
	      // multiplicative version
	      // sft_varied[itpl] *= lnB;

	      if (paramvalue[j]>0)
		{
		  if (syserr[i].highshape != 0)
		    {
		      if (hcl.get() == 0)
			{
			  hcl = copy_TH1<TH1Type> (*histotemplate[itpl]);
			}
		      else
			{
			  copy_TH1_content (*hcl, *histotemplate_varied[itpl]);
			}
		      if (poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{
		          csm_interpolate_histogram2(histotemplate[itpl],0.0,
						     syserr[i].highshape,syserr[i].xsighigh,
						     hcl.get(),histotemplate_varied[itpl],paramvalue[j],chan_istyle);
			}
		      else
			{
		          csm_interpolate_histogram2_noerr(histotemplate[itpl],0.0,
							   syserr[i].highshape,syserr[i].xsighigh,
							   hcl.get(),histotemplate_varied[itpl],paramvalue[j],chan_istyle);
			}

		      //cout << "did a +interpolation " << i << " " << j << " param: " << paramvalue[j] << endl;
		    }
		}
	      else
		{
		  if (syserr[i].lowshape != 0)
		    {
		      if (hcl.get() == 0)
			{
			  hcl = copy_TH1<TH1Type> (*histotemplate[itpl]);
			}
		      else
			{
			  copy_TH1_content (*hcl, *histotemplate_varied[itpl]);
			}
		      if (poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{
		          csm_interpolate_histogram2(histotemplate[itpl],0.0,
						     syserr[i].lowshape, -syserr[i].xsiglow,
						     hcl.get(),histotemplate_varied[itpl],-paramvalue[j],chan_istyle);
			}
		      else
			{
		          csm_interpolate_histogram2_noerr(histotemplate[itpl],0.0,
							   syserr[i].lowshape, -syserr[i].xsiglow,
							   hcl.get(),histotemplate_varied[itpl],-paramvalue[j],chan_istyle);
			}
		      //cout << "did a -interpolation " << i << " " << j << " param: " << paramvalue[j] <<  endl;
		    }
		}
	    }
	}
    }
  for (itpl=0;itpl<ntemplates;itpl++) sft_varied[itpl] *= max(1E-8,rsum[itpl]+1.0);

}

/*----------------------------------------------------------------------------*/
// resets all the varied histograms and scales to their unvaried states.

void csm_channel_model::undo_nuisance_response()
{
  Int_t i,ntemplates,nbinsx,nbinsy,ix,iy;

  ntemplates = (Int_t) histotemplate.size();
  for (i=0;i<ntemplates;i++)
    {
      sft_varied[i] = sft[i];
      nbinsx = histotemplate[i]->GetNbinsX();
      nbinsy = histotemplate[i]->GetNbinsY();
      if (nbinsy==1)
	{
          for (ix=1;ix<=nbinsx;ix++)
	    {
	      histotemplate_varied[i]->SetBinContent(ix,histotemplate[i]->GetBinContent(ix));
	      histotemplate_varied[i]->SetBinError(ix,histotemplate[i]->GetBinError(ix));
	    }
	}
      else
	{
          for (ix=1;ix<=nbinsx;ix++)
	    {
	      for (iy=1;iy<=nbinsy;iy++)
		{
	          histotemplate_varied[i]->SetBinContent(ix,iy,histotemplate[i]->GetBinContent(ix,iy));
	          histotemplate_varied[i]->SetBinError(ix,iy,histotemplate[i]->GetBinError(ix,iy));
		}
	    }
	}
    }
}

/*----------------------------------------------------------------------------*/
// check to see if any bin has a total negative prediction in this channel

Int_t csm_channel_model::checkneg()
{
  cout << "csm_channel_model::checkneg() to be written" << endl;
  return(0);
}

/*----------------------------------------------------------------------------*/
// collect all nuisance parameter names and upper and lower bounds for this model
// Where the bounds come from -- do not allow extrapolation on histogram shapes.  
// (all histogram extrapolation should be done and verified by the user)
// Also do not allow any contribution to signal or background to go negative.

void csm_model::list_nparams(vector<char *> *npn, vector<Double_t> *nplb, vector<Double_t> *nphb)
{
  Int_t i,j,k,ifound;
  csm_channel_model* cm;
  Double_t nplb_tmp,nphb_tmp;
  //Double_t a,b,c,disc,xp,xm,xht,xlt;
  npn->clear();
  nplb->clear();
  nphb->clear();

  for (i=0;i<(Int_t) channame.size();i++)
    {
      cm = chanmodel[i];
      for (j=0;j<(Int_t) cm->syserr.size();j++)
	{
	  //cout << "sys error item channel: " << i << 
          //" error index: " << j << " " << cm->syserr[j].sysname << " " <<
	  //(cm->syserr[j].highshape != 0) << " " <<  
          //(cm->syserr[j].lowshape != 0) << " " <<
          //cm->syserr[j].xsiglow << " " << cm->syserr[j].xsighigh << endl;  

          // question -- do we need to consider nuisance parameter variations beyond 20 sigma?
          // probably not if we only need 5-sigma discovery significance.

          nplb_tmp = -20;
	  nphb_tmp = 20;

          // Require the user to supply shape variations out to the number of sigma
          // we will investigate here.  This program won't do shape extrapolations internally,
          // but the csm_pvmorph subroutine supplied will in fact extrapolate.  Users should
          // look at what they get when extrapolating histograms, though -- check and validate.

	  // 26 Feb 2008 -- allow shape extrapolations of templates beyond the provided ranges.

	  if (cm->chan_istyle != CSM_INTERP_HORIZONTAL_EXTRAP && cm->chan_istyle != CSM_INTERP_VERTICAL_EXTRAP)
	    {
              if (cm->syserr[j].lowshape != 0)
	        { nplb_tmp = max(nplb_tmp,cm->syserr[j].xsiglow); }
              if (cm->syserr[j].highshape != 0)
	        { nphb_tmp = min(nphb_tmp,cm->syserr[j].xsighigh); }
	    }

          // limit the nuisance paramters also so that individual scale factors do not go negative.
          // There's protection in the fit function, but we need the pseudoexperiments also to
          // be sensible -- this is the equivalent (using the asymmetric errors supplied) of the
          // truncated Gaussian

	  // December 7, 2009 -- with Wade's new lognormal priors we don't need this any more

	  //a = (cm->syserr[j].sysfrach + cm->syserr[j].sysfracl)/2.0;
	  //b = (cm->syserr[j].sysfrach - cm->syserr[j].sysfracl)/2.0;
          //c = 1;
	  //if (a == 0)
	  //  {
	  //    if (b > 0)
	  //	{ nplb_tmp = max(nplb_tmp,-1.0/b); }
	  //    if (b < 0)
	  //	{ nphb_tmp = min(nphb_tmp,-1.0/b); }
	  //  }
	  //else
	  //  {
	  //    disc = b*b - 4.0*a*c;
	  //    if (disc > 0)
	  //	{ 
	  //	  xp = (-b + sqrt(disc))/(2.0*a);
	  //	  xm = (-b - sqrt(disc))/(2.0*a);
	  //	  xht = max(xp,xm);
	  //	  xlt = min(xp,xm); 
	  //	  // we know that a nuisance parameter value of 0 has a non-negative prediction,
	  //       // but the choice of which of these two solutions to a quadratic to take depends
	  //	  // on which side of zero they are on.
	  //	  if (xht < 0)
	  //	    {
	  //	      nplb_tmp = max(nplb_tmp,xht);
	  //	    }
	  //	  else if (xlt > 0) 
	  //	    {
	  //	      nphb_tmp = min(nphb_tmp,xlt);
	  //	    }
	  //	  else
	  //	    {
	  //	      nphb_tmp = min(nphb_tmp,xht);
	  //	      nplb_tmp = max(nplb_tmp,xlt);
	  //	    }
	  //	}
	  // }

	  ifound = -1;
	  for (k=0;k<(Int_t) npn->size();k++)
	    {
	      if (strcmp(cm->syserr[j].sysname,(*npn)[k]) == 0) { ifound = k; }
	    }
	  if (ifound == -1)
	    {
	      npn->push_back(cm->syserr[j].sysname);
	      nplb->push_back(nplb_tmp);
	      nphb->push_back(nphb_tmp);
	      //cout << "sysname: " << cm->syserr[j].sysname << " assigned ranges: " << nplb_tmp << " " << nphb_tmp << endl;
	    }
	  else
	    {
	      (*nplb)[ifound] = max((*nplb)[ifound],nplb_tmp);
	      (*nphb)[ifound] = min((*nphb)[ifound],nphb_tmp); 
	      //cout << "sysname: " << cm->syserr[j].sysname << " reassigned ranges: " << nplb_tmp << " " << nphb_tmp <<  " " <<
	      // (*nplb)[ifound] << " " << (*nphb)[ifound] << endl;
	    }


	}
    }
  // add in user-specified bounds (27 Feb 2008)

  for (k=0;k<(Int_t) npbname.size();k++)
    {
      for (j=0;j<(Int_t) npn->size();j++)
	{
	  if (strcmp(npbname[k],(*npn)[j])==0)
	    {
	      (*nplb)[j] = max((*nplb)[j],npblow[k]);
	      (*nphb)[j] = min((*nphb)[j],npbhigh[k]); 
	    }
	}
    }
}

/*----------------------------------------------------------------------------*/
/* a splitoff from single_pseudoexperiment -- just vary the templates but do  */
/* not generate pseudodata.   Useful for interfacing with Joel's program      */
/*----------------------------------------------------------------------------*/

void csm_model::varysyst()
{
  vector<Double_t> nplb;
  vector<Double_t> nphb;
  Int_t i;
  Double_t xval;

  // systematically fluctuate our model

  list_nparams(&npnp, &nplb, &nphb);   // this clears and fills the vectors of name pointers and bounds
  //cout << " in pe: npn.size " << npn.size() << endl;
  npvalp.clear();

  for (i=0;i<(Int_t)npnp.size();i++)
    {
      do 
	{
          xval = gRandom->Gaus(0,1);
          //cout << i << " " << nplb[i] << " " << xval << " " << nphb[i] << endl;
	}
      while (xval < nplb[i] || xval > nphb[i]);
      npvalp.push_back(xval);
    }
  nuisance_response(npnp.size(),&(npnp[0]),&(npvalp[0]));
}

// prints out names and values of nuisance parameters generated in the latest call to varysyst.

void csm_model::print_nuisance_params()
{
  cout << "Nuisance parameter listing: " << endl;
  for (int i=0;i<(int) npnp.size();i++)
    {
      cout << i << " " << npnp[i] << " " << npvalp[i] << endl; 
    }
}

/*----------------------------------------------------------------------------*/

/* Generate a single pseudoexperiment from a model -- fluctuate all nuisance parameters
   with their uncertainties -- the pseudodata histograms are in the same order as 
   the channels in the model description, with the same binning assumed.  The psuedodata
   histograms should be allocated in the calling routine.  That way the histograms don't
   have to be continually created and destroyed for each pseudoexperiment but can be
   re-used.*/

void csm_model::single_pseudoexperiment(TH1 *pseudodata[])
{
  Int_t ichan,itpl,ibinx,ibiny,nbinsx,nbinsy,nchans,ntemplates;
  csm_channel_model* cm;
  Double_t bintot;
  TH1Type* ht;
  Double_t r;

  // call nuisance_response with random nuisance parameters

  varysyst();
 
  // generate random pseudodata.  Randomly fluctuate the Poisson subsidiary
  // experiments (a "systematic effect") to figure out what the proper mean
  // is for the main experiment.

  nchans = (Int_t) channame.size();
  for (ichan=0;ichan<nchans;ichan++)
    {
      cm = chanmodel[ichan];
      nbinsx = cm->histotemplate[0]->GetNbinsX();
      nbinsy = cm->histotemplate[0]->GetNbinsY(); 
      for (ibinx=0;ibinx<nbinsx;ibinx++)
	{
	  for (ibiny=0;ibiny<nbinsy;ibiny++)
	    {
	      bintot = 0;
	      ntemplates = (Int_t) cm->histotemplate.size();
              for (itpl=0;itpl<ntemplates;itpl++)
	        {
	          ht = cm->histotemplate_varied[itpl];
		  if (nbinsy == 1) 
		    { r = ht->GetBinContent(ibinx+1); }
		  else
		    { r = ht->GetBinContent(ibinx+1,ibiny+1); }
	          if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
		    { r = gRandom->Poisson(r); }
		  else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
		    { 
		      double histerr,edraw;
		      if (nbinsy==1)
		        { histerr = ht->GetBinError(ibinx+1);}
		      else
		        { histerr = ht->GetBinError(ibinx+1,ibiny+1);}
		      do
		        { edraw = gRandom->Gaus(0,histerr); }
		      while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
		      r += edraw;
		    }
	          r *= cm->sft_varied[itpl];
	          bintot += r;
		}
              r = gRandom->Poisson(bintot);
	      if (nbinsy == 1)
		{ pseudodata[ichan]->SetBinContent(ibinx+1,r); }
	      else
		{ pseudodata[ichan]->SetBinContent(ibinx+1,ibiny+1,r); }
	    } // end loop over binsy
	} // end loop over binsx
    } // end loop over channels
}

/*----------------------------------------------------------------------------*/

// calculate the ratio p(data|nmodel)/p(data|dmodel)

Double_t mclimit_csm::weightratio(csm_model *nmodel, csm_model *dmodel, TH1 *hist[])
{
  int nchans = nmodel->channame.size();
  int ichan;
  csm_channel_model *ncm;
  csm_channel_model *dcm;
  int nbinsx,nbinsy,ibin,jbin,ic;
  double wr,pn,pd;
  int dtb,ncc,dcc;

  wr = 0;
  for (ichan=0;ichan<nchans;ichan++)
    {
      ncm = nmodel->chanmodel[ichan];
      ncc = ncm->histotemplate.size();
      dcm = dmodel->chanmodel[ichan];
      dcc = dcm->histotemplate.size();

      nbinsx = hist[ichan]->GetNbinsX();
      nbinsy = hist[ichan]->GetNbinsY();
      for (ibin=0;ibin<nbinsx;ibin++)
	{
	  for (jbin=0;jbin<nbinsy;jbin++)
	    {
	      dtb = 0;
	      if (nbinsy == 1)
		{
		  dtb = max(0,(int) nearbyint(hist[ichan]->GetBinContent(ibin+1)));
		}
	      else
		{
		  dtb = max(0,(int) nearbyint(hist[ichan]->GetBinContent(ibin+1,jbin+1)));
		}

	      pn = 0; // prediction for numerator model
	      for (ic=0;ic<ncc;ic++)
		{
		  if (nbinsy == 1)
		    {
		      pn += ncm->histotemplate_varied[ic]->GetBinContent(ibin+1)*
			ncm->sft_varied[ic];
		    }
		  else
		    {
		      pn += ncm->histotemplate_varied[ic]->GetBinContent(ibin+1,jbin+1)*
			ncm->sft_varied[ic];
		    }
		}
	      pd = 0; // prediction for denominator model
	      for (ic=0;ic<dcc;ic++)
		{
		  if (nbinsy == 1)
		    {
		      pd += dcm->histotemplate_varied[ic]->GetBinContent(ibin+1)*
			dcm->sft_varied[ic];
		    }
		  else
		    {
		      pd += dcm->histotemplate_varied[ic]->GetBinContent(ibin+1,jbin+1)*
			dcm->sft_varied[ic];
		    }
		}
	      if (pd > 0)
		{
	          wr += dtb*log(pn/pd) - pn + pd;
		}
	    }
	}
    }
  // about the limit of exponentials
  if (wr>680.) 
    { wr = 680.; }
  wr = exp(wr);
  return(wr);
}

/*----------------------------------------------------------------------------*/

// minimize the chisquared function over the nuisance parameters

Double_t mclimit_csm::calc_chi2(csm_model *model, TH1 *hist[])
{
  Int_t i;
  Double_t chisquared;
  csm* mycsm = new csm;
  mycsm->setminuitmaxcalls(minuitmaxcalls);
  mycsm->setminosmaxcalls(minosmaxcalls);
  mycsm->setminuitstepsize(minuitstepsize);
  mycsm->setprintflag(minuitprintflag);
  mycsm->setminosflag(minosflag);
  mycsm->set_modeltofit(model);

  // assume (check?) that the model channel names match up with the data channel names

  for (i=0;i<(Int_t)model->channame.size();i++)
    {
      mycsm->set_htofit(hist[i],model->channame[i]);
    }
  chisquared = mycsm->chisquared();
  //  cout << "total number of nuisance parameters: " << mycsm->getnparams() << endl;
  //for (int i=0;i<mycsm->getnparams();i++)
  //  {
  //    cout << mycsm->getpname(i) << endl;
  //  }
  delete mycsm;
  return(chisquared);
}

Double_t csm::chisquared()
{
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;
  Double_t arglist[20];
  Int_t ierflag = 0;
  Int_t i,j,icons;
  Double_t cresult;
  Double_t param,paramerror;

  modeltofit->list_nparams(&npn, &nplb, &nphb);
  Int_t npns = npn.size();

  if (npns > 0)
    {
      
      TMinuit *mnp = new TMinuit(npns+1);

      mnp->SetFCN(csm_minuit_fcn);
      if (!minuitprintflag)
	{
          arglist[0] = -1;
          mnp->mnexcm("SET PRINT", arglist, 1, ierflag);
          mnp->mnexcm("SET NOW",arglist,1,ierflag);
	}

      arglist[0] = 2; 
      mnp->mnexcm("SET STRATEGY", arglist, 1, ierflag);	
      mnp->mnexcm("SET NOG",arglist,1,ierflag);  // no gradiants required of FCN


      mnp->SetMaxIterations(minuitmaxcalls); // doesn't seem to do anything in TMinuit
      char npname[10];

      for (i=0;i < (Int_t) npns;i++)
        {
          sprintf(npname,"np%d",i);
	  //cout << "setting minuit parameter: " << npname << " " << nplb[i] << " " << nphb[i] << endl;
	  TString npname2 = npname;
          mnp->mnparm(i,npname2,0.0,minuitstepsize,nplb[i],nphb[i],ierflag);
          icons = 1;
	  for (j=0;j<(Int_t) modeltofit->npcm.size();j++)
	    {
	      if (strcmp(modeltofit->npcm[j].pnameoutput,npn[i])==0)
		{
		  ierflag = mnp->FixParameter(i);
		  icons = 0;
		}
	    }
          char *s = new char[strlen(npn[i])+1];
          strcpy(s,npn[i]);
          npfitname.push_back(s);  // this copy is in static global storage so the minuit function knows about it
          if (strstr(npn[i],"UNCONSTRAINED") != 0)
  	    {
	      icons = 0;
	    }
          constrainedfitparam.push_back(icons);
        }

      arglist[0] = 1;
      mnp->mnexcm("SET ERR", arglist ,1,ierflag); 
      ierflag = 0;
      arglist[0] = minuitmaxcalls;  // here's where maxcalls makes a difference
      arglist[1] = 1.;

      //      mnp->mnexcm("SIMPLEX", arglist ,2,ierflag);
      //      mnp->mnexcm("MIGRAD", arglist ,2,ierflag);

      mnp->mnexcm("MINI", arglist ,2,ierflag);
      mnp->mnexcm("IMPROVE", arglist ,2,ierflag);

      if (minosflag) 
	{
          arglist[0] = minosmaxcalls;
          mnp->mnexcm("MINOS",arglist,1,ierflag);
	}

      //cout << "Number of function calls in Minuit: " << mnp->fNfcn << endl;

      // copy best fit parameters for outside use

      cresult = mnp->fAmin;
      fitparam.clear();
      fiterror.clear();
      for (i=0;i<(Int_t) fitparamname.size();i++)
        {
          delete[] fitparamname[i];
        }
      fitparamname.clear();

      // allocate memory for the covariance matrix only if we have to.
      // (re-use the old memory if it has the right size)

      if (nfitcov != npns)
	{
	  if (fitcov != 0)
	    { 
	      delete[] fitcov;
	    }
	  nfitcov = 0;
	}
      if (nfitcov == 0)
        { 
	  fitcov = new Double_t[npns*npns];
	  nfitcov = npns;
	}

      for (i=0;i < (Int_t) npns;i++)
        {
          mnp->GetParameter(i,param,paramerror);
          fitparam.push_back(param);
          fiterror.push_back(paramerror);
          char *s = new char[strlen(npn[i])+1];
          strcpy(s,npn[i]);
          fitparamname.push_back(s); // this copy's part of the class private members
	  mnp->mnemat(fitcov,npns);
        }

      delete mnp;
    }
  else
    {
      i = 0;
      csm_minuit_fcn(i,0,cresult,0,0);
    }
  if (cresult < 0)
    { 
      //cout << "chisquared less than zero: " << cresult << " setting it to zero" << endl;
      cresult = 0;
    }

  for (i=0;i<(Int_t) npfitname.size();i++)
    {
      delete[] npfitname[i];
    }
  npfitname.clear();
  constrainedfitparam.clear();

  return(cresult);
}

/*----------------------------------------------------------------------------*/

// A model is a collection of channel models and names

csm_model::csm_model()
{
}

/*----------------------------------------------------------------------------*/

csm_model::~csm_model()
{
  Int_t i,j;
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      delete[] channame[i];
      delete chanmodel[i];
    }
  for (i=0;i<(Int_t) npcm.size();i++)
    {
      for (j=0;j<npcm[i].ninput;j++)
	{
	  delete[] npcm[i].pnameinput[j];
	}
      delete[] npcm[i].pnameinput;
      delete[] npcm[i].pnameoutput;
    }
  for (i=0;i<(Int_t) npbname.size();i++)
    {
      delete[] npbname[i];
    }

  /* The vectors themselves are deleted when the class instance is deleted */
}

/*----------------------------------------------------------------------------*/

void csm_model::add_template(TH1 *template_hist, //Poisson or non-Poisson histogram
			     Double_t sf,        //scale factor to multiply template by to compare w/ data 
			     //(e.g., (data_lum/MC_lum) for a MC Poisson histogram
			     Int_t nnp,          // number of nuisance parameters (Gaussian of unit width)
			     char* npname[],     // nuisance parameter names 
			     Double_t *nps_low,  // fractional uncertainty on sf due to each nuisance parameter -- low side
			     Double_t *nps_high, // fractional uncertainty on sf due to each nuisance parameter -- high side
			     // typically nps_low and nps_high are input with opposite signs -- if opposite
			     // variations of the nuisance parameter create opposite changes in sf.  The sign
			     // is retained in the calculation in case both variations of a nuisance parameter
			     // shift the normalization in the same way (either both + or both -)
			     TH1 *lowshape[],    // array of low hisogram shapes, one for each nuisance param (null if no shape error)
			     Double_t *lowsigma, // number of sigma low for each nuisance parameter shape variation
			     TH1 *highshape[],   // array of high histogram shapes, one for each nuisance param (null if no shape error)
			     Double_t *highsigma, // number of sigma high for each shape variation
			     Int_t pflag,         // Poisson flag -- 1 if Poisson, 0 of not.  2 if Gaussian error from the histo contents
			     Int_t sflag,         // scale flag -- 1 if signal, 0 if background (for use with s95 calculator)
			     char *cname)
{
  Int_t i;
  i = lookup_add_channame(cname);
  chanmodel[i]->add_template(template_hist,sf,nnp,npname,nps_low,
                             nps_high,lowshape,lowsigma,highshape,
                             highsigma,pflag,sflag);
}

/*----------------------------------------------------------------------------*/
// add a whole channel's model to the total set of models.

void csm_model::add_chanmodel(csm_channel_model *cm, char *cname)
{
  Int_t ichan;

  ichan = lookup_add_channame(cname);
  chanmodel[ichan] = cm->Clone();
}

void csm_model::add_npbounds(char *pname, Double_t lowbound, Double_t highbound)
{
  char *s = new char[strlen(pname)+1];
  strcpy(s,pname);
  npbname.push_back(s);
  npblow.push_back(lowbound);
  npbhigh.push_back(highbound);
}

/*----------------------------------------------------------------------------*/
// add a constraint function between nuisance parameters.  Make our own copies of all
// the names and the function pointer.

void csm_model::add_npcons(Int_t nparin,char **parin, char *parout, Double_t (*f)(Double_t*))
{
  Int_t i,j;
  npcstruct npc;
  char *s;

  for (i=0;i<(Int_t) npcm.size();i++)
    {
      if (strcmp(parout,npcm[i].pnameoutput)==0)
	{
	  cout << "Warning: Two constraint functions for the same nuisance parameter: " << parout << " defined" << endl;
          exit(0); // bad enough to crash
	}
      for (j=0;j<npcm[i].ninput;j++)
	{
	  if (strcmp(parout,npcm[i].pnameinput[j]) == 0)
	    {
	      cout << "Warning: nuisance parameter: " << npcm[i].pnameinput[j] << " depends on " << parout << endl;
              cout << "but " << parout << "is computed itsef by a constraint after " << npcm[i].pnameinput[j] << "is computed." << endl;
              exit(0); // bad enough to crash
	    }
	}
    }

  npc.ninput = nparin;
  npc.pnameinput = new char*[nparin];
  for (i=0;i<nparin;i++)
    {
      s = new char[strlen(parin[i])+1];
      strcpy(s,parin[i]);
      if (strcmp(s,parout)==0)
	{
	  cout << "Constraint function for nuisance parameter: " << s 
               << " depends on nuisance parameter " << s << endl;
	  exit(0);
	}
      npc.pnameinput[i] = s;
    }
  s = new char[strlen(parout)+1];
  strcpy(s,parout);
  npc.pnameoutput = s;
  npc.f = f;
  npcm.push_back(npc);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::Clone()
{
  Int_t i;
  csm_model* mclone = new csm_model;

  for (i=0;i < (Int_t) channame.size(); i++)
    {
      mclone->add_chanmodel(chanmodel[i],channame[i]);
    } 
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      mclone->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      mclone->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(mclone);
}

/*----------------------------------------------------------------------------*/
// includes all the new contributing histograms, and also collects together all constraint
// relationships between nuisance parameters.

csm_model* csm_model::add(csm_model &a)
{
  Int_t i;
  csm_model* mclone = a.Clone();
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      mclone->add_chanmodel(chanmodel[i],channame[i]);
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      mclone->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      mclone->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(mclone);
}

/*----------------------------------------------------------------------------*/

// returns a new model, interpolated frac of the way from a to b.  frac=0: makes a clone of a (= this model),
// frac=1, makes a clone of b.  No real checking here that a and b are commensurate -- they should
// have an identical list of channels. nuisance parameter
// bounds and nuisance parameter constraints are gotten from this model.

csm_model* csm_model::interpolate(csm_model *b, double frac)
{
  int i,j,nchans;
  csm_channel_model *bc;
  csm_model *imodel = new csm_model();

  nchans = channame.size();
  for (i=0;i<nchans;i++)
    {
      for (j=0;j<nchans;j++)
	{
	  if (strcmp(channame[i],b->channame[j])==0)
	    {
	      bc = chanmodel[i]->interpolate(b->chanmodel[j],frac);
	      imodel->add_chanmodel(bc,channame[i]);
	      delete bc;
	    }
	}
    }

  for (i=0;i < (int) npcm.size(); i++)
    {
      imodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }

  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      imodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }

  return(imodel);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::scale(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scale(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      smodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::scalesignal(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scalesignal(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      smodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::scale_err(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scale_err(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      smodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/
/* Build the list of channel models and channel names that is sorted by channel */
/* name during the building process.  Return the vector index to use to refer */
/* to this particular channel. */

Int_t csm_model::lookup_add_channame(char *cname)
{
  Int_t i,ifound,j,jfound;
  char *s;
  csm_channel_model *cm; 
  vector<char*>::iterator cni;
  vector<csm_channel_model*>::iterator cmi;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      j = (Int_t) strcmp(cname,channame[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      cm = new csm_channel_model;
      strcpy(s,cname);
      if (jfound == -1)
	{
          ifound = channame.size();
          channame.push_back(s);
          chanmodel.push_back(cm);
	}
      else
	{
	  ifound = jfound;
	  cni = channame.begin() + jfound;
	  channame.insert(cni,s);
	  cmi = chanmodel.begin() + jfound;
	  chanmodel.insert(cmi,cm);
	}
    }

  return(ifound);
}

/*----------------------------------------------------------------------------*/

/* A channel model is a sum of template histograms along with systematic errors */

// constructor

csm_channel_model::csm_channel_model()
{
  chan_istyle = CSM_INTERP_HORIZONTAL;  //  defaults to csm_pvmorph interpolation
}

/*----------------------------------------------------------------------------*/

// destructor

csm_channel_model::~csm_channel_model()
{
  Int_t i;

  //cout << "Called csm_channel_model destructor" << endl;

  // deallocate memory used to save sytematic error names

  for (i=0;i < (Int_t) syserr.size();i++)
    {
      delete[] syserr[i].sysname;
    }
  // deallocate cloned input histograms
  for (i=0;i < (Int_t) histotemplate.size();i++)
    {
      delete histotemplate[i];
      delete histotemplate_varied[i];
    }
  for (i=0;i < (Int_t) syserr.size();i++)
    {
      if (syserr[i].lowshape !=0)
	{
	  delete syserr[i].lowshape;
	}
      if (syserr[i].highshape !=0)
	{
	  delete syserr[i].highshape;
	}
    }
}

TH1* mclimit_csm::get_datahist(Int_t i)
{
  return(datahist[i]);
}

// Print out signal, background and data for all bins in all channels.
// use test_hypothesis_pe which has all components.

void mclimit_csm::printsbd()
{
  Int_t i,j,k;
  Int_t nbinsx,nbinsy;

  Int_t nchans = (Int_t) test_hypothesis_pe->channame.size();

  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      Int_t nobs;
	      Double_t nsig = 0;
	      Double_t nbkg = 0;
	      if (nbinsy==1)
		{ nobs = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      csm_channel_model *cm = test_hypothesis_pe->chanmodel[i];
	      Int_t ntemplates = (Int_t) cm->histotemplate.size();
	      for(Int_t itpl=0;itpl<ntemplates;itpl++)
		{
		  Double_t r;
		  if (nbinsy==1)
		    { r = cm->histotemplate_varied[itpl]->GetBinContent(j+1); }
		  else
		    { r = cm->histotemplate_varied[itpl]->GetBinContent(j+1,k+1); }
                  r *= cm->sft_varied[itpl];
                  if (cm->scaleflag[itpl] != 0)
		    { nsig += r; }
		  else
		    { nbkg += r; }
		}
	      cout << "DumpSBD: " << nsig << " " << nbkg << " " << nobs << endl;
	    }
	}
    }  
}


void csm_model::print()
{
  Int_t i,j;

  cout << "csm_model::print  -- printing out model information" << endl;
  for (i=0;i<(Int_t) channame.size();i++)
    {
      cout << "Channel: " << i << " Name: " << channame[i] << endl;
      chanmodel[i]->print();
    }
  for (i=0;i<(Int_t) npcm.size();i++)
    {
      cout << "-------------------" << endl;
      cout << "Constraint equation:  " << npcm[i].pnameoutput << " is computed from " << endl;
      for (j=0;j<npcm[i].ninput;j++)
	{
	  cout << npcm[i].pnameinput[j] << endl;
	}
    }
  cout << "-------------------" << endl;
  for (i=0;i<(Int_t) npbname.size(); i++)
    {
      cout << "NP Bounds.  Name, lowbound highbound: " << npbname[i] << " " << npblow[i] << " " << npbhigh[i] << endl; 
    }
}

// print out just a piece of a model, indexed by channel name.

void csm_model::print(char *channame)
{
  Int_t i = lookup_add_channame(channame);
  cout << "Printing One Channel: " << channame << endl;
  chanmodel[i]->print();
}

void csm_channel_model::print()
{
  Int_t i,j;
  Double_t ssum = 0;
  Double_t bsum = 0;
  Double_t ssumv = 0;
  Double_t bsumv = 0;
  Double_t central_integral=0;

  cout << "Begin-----------------csm_channel_model::print()------------" << endl;

  for(i=0;i < (Int_t) histotemplate.size();i++)
    {
      cout << endl;
      cout << "Template " << i << endl;
      cout << "  Histogram name: " << histotemplate[i]->GetName();
      cout << "  Histogram title: " << histotemplate[i]->GetTitle();
      cout << "  sft: " << sft[i] << endl;
      cout << "  sft_varied: " << sft_varied[i] << endl;
      cout << "  poissflag: " << poissflag[i] << endl;
      cout << "  signalflag: " << scaleflag[i] << endl;
      central_integral = histotemplate[i]->Integral();
      cout << "  Integral: " << central_integral << endl;
      cout << "  Scaled Integral: " << histotemplate[i]->Integral()*sft[i] << endl;
      cout << "  Scaled Integral with all syst: " << histotemplate_varied[i]->Integral()*sft_varied[i] << endl;
      cout << "  Template bbeta: " << (histotemplate_varied[i]->Integral()*sft_varied[i])/(histotemplate[i]->Integral()*sft[i]) << endl;
      if (scaleflag[i]) 
	{
	  ssum += histotemplate[i]->Integral()*sft[i];
	  ssumv += histotemplate_varied[i]->Integral()*sft_varied[i];
	}
      else
	{
	  bsum += histotemplate[i]->Integral()*sft[i];
	  bsumv += histotemplate_varied[i]->Integral()*sft_varied[i];
	}

      //histotemplate[i]->Print("all");
      Double_t errtotup = 0;
      Double_t errtotdown = 0;
      Double_t uperrloc,downerrloc,hsi,hsl,fracerr;

      for (j=0;j< (Int_t) syserr.size();j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      cout << "Syst: " << syserr[j].sysname << endl;
	      uperrloc = syserr[j].sysfrach;
	      downerrloc = syserr[j].sysfracl;
	      cout << "  Up rate error: " << uperrloc << endl;
	      cout << "  Down rate error: " << downerrloc << endl;

	      if (syserr[j].highshape != 0)
		{
		  hsi = syserr[j].highshape->Integral();
		  cout << "  Up shape error provided sigma: " << syserr[j].xsighigh << 
		    " integral: " << hsi << endl;
		  fracerr = (hsi-central_integral)/central_integral;
		  cout << "    Fractional error due to the shape: " << fracerr << endl;
		  uperrloc += fracerr;
		  cout << "    Total up rate error incl. shape: " << uperrloc << endl;
		  //syserr[j].highshape->Print("all");
		}
	      if (syserr[j].lowshape != 0)
		{
		  hsl = syserr[j].lowshape->Integral();
		  cout << "  Down shape error provided sigma: " << syserr[j].xsiglow << 
		    " integral: " << hsl << endl;

		  fracerr = (hsl-central_integral)/central_integral;
		  cout << "    Fractional error due to the shape: " << fracerr << endl;
		  downerrloc += fracerr;
		  cout << "    Total down rate error incl. shape: " << downerrloc << endl;

		  //syserr[j].lowshape->Print("all");
		}

	      errtotup += uperrloc*uperrloc;
	      errtotdown += downerrloc*downerrloc;

	    }
	}
      errtotup = sqrt(errtotup);
      errtotdown = sqrt(errtotdown);
      cout << "Total relative error on this template (up): " << errtotup << endl;
      cout << "Total relative error on this template (down): " << errtotdown << endl;
    }

  cout << "Total signal this channel in this model: " << ssum << endl;
  cout << "Total background this channel in this model: " << bsum << endl;
  cout << "Syst. Varied Total signal this channel in this model: " << ssumv << endl;
  cout << "Syst. Varied Total background this channel in this model: " << bsumv << endl;
  cout << "End-------------------csm_channel_model::print()------------" << endl;
}
/*----------------------------------------------------------------------------*/

void csm_channel_model::add_template
(const TH1Input& template_hist, //Poisson or non-Poisson histogram
 Double_t sf,        //scale factor to multiply template by to compare w/ data 
 //(e.g., (data_lum/MC_lum) for a MC Poisson histogram
 Int_t nnp,          // number of nuisance parameters (Gaussian of unit width)
 char* npname[],     // nuisance parameter names 
 Double_t *nps_low,  // fractional uncertainty on sf due to each nuisance parameter -- low side
 Double_t *nps_high, // fractional uncertainty on sf due to each nuisance parameter -- high side
 // typically nps_low and nps_high are input with opposite signs -- if opposite
 // variations of the nuisance parameter create opposite changes in sf.  The sign
 // is retained in the calculation in case both variations of a nuisance parameter
 // shift the normalization in the same way (either both + or both -)
 const TH1InputList& lowshape,   // array of low hisogram shapes, one for each nuisance param (null if no shape error)
 Double_t *lowsigma,  // number of sigma low for each nuisance parameter shape variation
 const TH1InputList& highshape,   // array of high histogram shapes, one for each nuisance param (null if no shape error)
 Double_t *highsigma, // number of sigma high for each shape variation
 Int_t pflag,         // Poisson flag -- 1 if Poisson, 0 of not.  2 if Gaussian error from the histo contents
 Int_t sflag)        // scale flag -- 1 if signal, 0 if background (for use with s95 calculator)
{
  int i;
  svstruct ses;
  char *s;

  sft.push_back(sf);
  sft_varied.push_back(sf);
  poissflag.push_back(pflag);
  scaleflag.push_back(sflag);
  TH1Type *template_hist_arg;
  histotemplate.push_back(template_hist_arg = template_hist.copy_TH1<TH1Type>().release());
  histotemplate_varied.push_back(template_hist.copy_TH1<TH1Type>().release());
  for (i=0;i<nnp;i++)
    {
      s = new char[strlen(npname[i])+1];
      strcpy(s,npname[i]);
      ses.sysname = s;
      ses.itemplate = histotemplate.size()-1;
      ses.sysfracl = nps_low[i];
      ses.sysfrach = nps_high[i];
      if (lowshape.has(i))
	{
	  ses.lowshape = lowshape.copy_TH1<TH1Type> (i).release();
	}
      else
	{
	  ses.lowshape = 0;
	}
      if (highshape.has(i))
	{
	  ses.highshape = highshape.copy_TH1<TH1Type> (i).release();
	}
      else
	{
	  ses.highshape = 0;
	}
      ses.xsiglow = lowsigma[i];
      ses.xsighigh = highsigma[i];
      semap[ses.sysname].push_back(syserr.size());
      syserr.push_back(ses);

      if (ses.highshape != 0)
	{
          if (ses.highshape->GetNbinsX() != template_hist_arg->GetNbinsX())
            {
              cout << "Chisquared minmization:  histo template and high shape have different bin counts." << endl;
              cout << template_hist_arg->GetNbinsX() << " != " << ses.highshape->GetNbinsX() << endl;
              exit(0);
	    }
        }
      if (ses.lowshape != 0)
	{
          if (ses.lowshape->GetNbinsX() != template_hist_arg->GetNbinsX())
            {
              cout << "Chisquared minmization:  histo template and low shape have different bin counts." << endl;
              cout <<  template_hist_arg->GetNbinsX() << " != " << ses.lowshape->GetNbinsX() << endl;
              exit(0);
	    }
	}
    }
  //cout << "model::add_template: " << histotemplate.size() << endl;
  //gDirectory->ls();
}

/*----------------------------------------------------------------------------*/

// make a copy of this model by adding the templates over again.
// that way the clone can be deleted by itself, and the destructor
// won't try to delete allocated memory twice

csm_channel_model* csm_channel_model::Clone()
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1Type **lowshape = new TH1Type *[syserr.size()];
  TH1Type **highshape = new TH1Type *[syserr.size()];
  char **ename = new char *[syserr.size()];

  csm_channel_model* mclone = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      //cout << "Model clone adding template " << i << endl;
      mclone->add_template(histotemplate[i],sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }

  mclone->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(mclone);
}

/*----------------------------------------------------------------------------*/

// interpolates channel models template by template.  Also interpolate systematic error
// rates and shapes.   Assumes exact congruence between systematic error lists and orderings
// of the two channel models to be interpolated.  Returns an interpolated result frac of the way from this
// channel model to the one pointed to by b.  Put in frac=0, get this channel model back, put in
// frac=1, get *b back.  This routine assumes that the high and low shape errors correspond to the same
// number of sigma.

csm_channel_model* csm_channel_model::interpolate(csm_channel_model *b, double frac)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1Type **lowshape = new TH1Type *[syserr.size()];
  TH1Type **highshape = new TH1Type *[syserr.size()];
  char **ename = new char *[syserr.size()];

  csm_channel_model* mint = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl + frac*(b->syserr[j].sysfracl-syserr[j].sysfracl);
	      nps_high[nnp] = syserr[j].sysfrach + frac*(b->syserr[j].sysfrach-syserr[j].sysfrach);
	      if (syserr[j].lowshape != 0 && b->syserr[j].lowshape != 0)
		{
		  std::auto_ptr<TH1Type> htj;
		  htj = copy_TH1<TH1Type> (*syserr[j].lowshape);
		  if (poissflag[i] == CSM_GAUSSIAN_BINERR)
		    {
		      csm_interpolate_histogram(syserr[j].lowshape,0.0,
						b->syserr[j].lowshape,1.0,
						htj.get(),frac,chan_istyle);
		    }
                  else
                    {
		      csm_interpolate_histogram_noerr(syserr[j].lowshape,0.0,
						      b->syserr[j].lowshape,1.0,
						      htj.get(),frac,chan_istyle);
	            }
		}
	      else
		{
		  lowshape[nnp] = 0;
		}

	      if (syserr[j].highshape != 0 && b->syserr[j].highshape != 0)
		{
		  std::auto_ptr<TH1Type> htj;
		  htj = copy_TH1<TH1Type> (*syserr[j].highshape);
		  if (poissflag[i] == CSM_GAUSSIAN_BINERR)
		    {
		      csm_interpolate_histogram(syserr[j].highshape,0.0,
						b->syserr[j].highshape,1.0,
						htj.get(),frac,chan_istyle);
		    }
                  else
                    {
		      csm_interpolate_histogram_noerr(syserr[j].highshape,0.0,
						      b->syserr[j].highshape,1.0,
						      htj.get(),frac,chan_istyle);
	            }
		}
	      else
		{
		  highshape[nnp] = 0;
		}

	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      //cout << "channel model interpolation adding template " << i << endl;
      std::auto_ptr<TH1Type> hti;
      hti = copy_TH1<TH1Type> (*histotemplate[i]);
      if (poissflag[i] == CSM_GAUSSIAN_BINERR)
	{
	  csm_interpolate_histogram(histotemplate[i],0.0,
				    b->histotemplate[i],1.0,
				    hti.get(),frac,chan_istyle);
	}
      else
	{
	  csm_interpolate_histogram_noerr(histotemplate[i],0.0,
					  b->histotemplate[i],1.0,
					  hti.get(),frac,chan_istyle);
	}
      mint->add_template(hti.get(),sft[i],nnp,ename,
			 nps_low,nps_high,lowshape,lowsigma,
			 highshape,highsigma,poissflag[i],scaleflag[i]);
    }

  mint->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(mint);
}

/*----------------------------------------------------------------------------*/

// addition of two models makes a new model with the sum of the templates

csm_channel_model* csm_channel_model::add(csm_channel_model &a)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1Type **lowshape = new TH1Type*[syserr.size()];
  TH1Type **highshape = new TH1Type*[syserr.size()];
  char **ename = new char*[syserr.size()];

  csm_channel_model* mclone = a.Clone();

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      mclone->add_template(histotemplate[i],sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  mclone->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(mclone);
}

/*----------------------------------------------------------------------------*/

// multiplication of a model and a scalar

csm_channel_model* csm_channel_model::scale(Double_t coefficient)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1Type **lowshape = new TH1Type*[syserr.size()];
  TH1Type **highshape = new TH1Type*[syserr.size()];
  char **ename = new char*[syserr.size()];

  csm_channel_model* smodel = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      smodel->add_template(histotemplate[i],coefficient*sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

/*----------------------------------------------------------------------------*/

// multiplication of a model and a scalar -- scale the systematic
// errors down with 1/sqrt(coefficient)

csm_channel_model* csm_channel_model::scale_err(Double_t coefficient)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t escale;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1Type **lowshape = new TH1Type*[syserr.size()];
  TH1Type **highshape = new TH1Type*[syserr.size()];
  char **ename = new char*[syserr.size()];

  csm_channel_model* smodel = new csm_channel_model;

  escale = 1.0/sqrt(coefficient);

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i< ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j< nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl*escale;
	      nps_high[nnp] = syserr[j].sysfrach*escale;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow/escale;
	      highsigma[nnp] = syserr[j].xsighigh/escale;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      smodel->add_template(histotemplate[i],coefficient*sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

/*----------------------------------------------------------------------------*/

// multiplication of only parts of a model and a scalar

csm_channel_model* csm_channel_model::scalesignal(Double_t coefficient)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1Type **lowshape = new TH1Type*[syserr.size()];
  TH1Type **highshape = new TH1Type*[syserr.size()];
  char **ename = new char*[syserr.size()];
  Double_t sc1;

  csm_channel_model* smodel = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      sc1 = sft[i];
      if (scaleflag[i] != 0)
	{
	  sc1 = coefficient*sft[i];
	}
      smodel->add_template(histotemplate[i],sc1,nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

/*----------------------------------------------------------------------------*/


// Use TMinuit to minimize the chisquared in T. Devlin's note CDF 3126 wrt the
// nuisance parameters.

// updated here -- do a joint minimization over shared nuisance parameters
// for several histograms (channels).

// global (in this file) declarations are at the top of the source file

// constructor
csm::csm()
{
  //  cout << "Chisquared Minimizer Constructor called\n";
  datatofit.clear();
  datatofitname.clear();
  constrainedfitparam.clear();
  npfitname.clear();
  fitcov=0;  // we have not yet allocated memory for the fit error matrix.
  nfitcov=0;
  minuitmaxcalls = 500;
  minosmaxcalls = 500;
  minuitstepsize = 0.1;
  minuitprintflag = 0;
  minosflag = 0;
}

//destuctor

csm::~csm()
{
  Int_t i;
  //cout << "Chisquared Minimizer Destructor called\n";

  // clear out static global variables

  for (i=0;i<(Int_t) datatofit.size(); i++)
    {
      delete datatofit[i];
      delete[] datatofitname[i];
    }
  datatofit.clear();
  datatofitname.clear();

  // clear out allocated memory pointed to by our private members

  for (i=0;i<(Int_t) fitparamname.size();i++)
    {
      delete[] fitparamname[i];
    }
  if (fitcov) delete[] fitcov;
}

void csm::setminuitmaxcalls(Int_t maxcalls)
{
  minuitmaxcalls = maxcalls;
}
Int_t csm::getminuitmaxcalls()
{
  return(minuitmaxcalls);
}

void csm::setminosmaxcalls(Int_t maxcalls)
{
  minosmaxcalls = maxcalls;
}
Int_t csm::getminosmaxcalls()
{
  return(minosmaxcalls);
}

void csm::setminuitstepsize(Double_t stepsize)
{
  minuitstepsize = stepsize;
}
Double_t csm::getminuitstepsize()
{
  return(minuitstepsize);
}

void csm::setprintflag(bool pf)
{
  minuitprintflag = pf;
}
bool csm::getprintflag()
{
  return(minuitprintflag);
}

void csm::setminosflag(bool mf)
{
  minosflag = mf;
}
bool csm::getminosflag()
{
  return(minosflag);
}

// put in the data histograms in the same order we built up the model histograms

void csm::set_htofit(TH1 *h, char *cname)
{
  Int_t i,ifound,j,jfound;
  vector<char*>::iterator dni;
  vector<TH1*>::iterator dfi;
  char *s;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) datatofitname.size(); i++)
    {
      j = (Int_t) strcmp(cname,datatofitname[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings.  If the name is on the 
     list, replace the existing data histogram with a clone of the one supplied. */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      strcpy(s,cname);
      if (jfound == -1)
	{
          datatofitname.push_back(s);
          datatofit.push_back((TH1*) h->Clone());
	}
      else
	{
          dni = datatofitname.begin() + jfound;
	  datatofitname.insert(dni,s);
	  dfi = datatofit.begin() + jfound;
	  datatofit.insert(dfi,(TH1*) h->Clone());
	}
    }
  else
    {
      delete datatofit[ifound];
      datatofit[ifound] = (TH1*) h->Clone();
    }
}

void csm::set_modeltofit(csm_model* mtf)
{
  modeltofit = mtf;
}

csm_model* csm::getbestmodel()
{
  Int_t i;

  // make a local array of pointers to nuisance parameter names
  
  char **fpnameloc = new char *[fitparamname.size()];
  for (i=0;i<(Int_t) fitparamname.size();i++)
    {
      fpnameloc[i] = fitparamname[i];
      //cout << "in getbestmodel, paramname: " << fpnameloc[i] << endl;
    }
  Double_t *parloc = new Double_t[fitparam.size()];
  for (i=0;i<(Int_t) fitparam.size();i++)
    {
      parloc[i] = fitparam[i];
      //cout << "in getbestmodel, param: " << parloc[i] << endl;
    }

  modeltofit->nuisance_response(fitparam.size(),fpnameloc,parloc);
  delete[] fpnameloc;
  delete[] parloc;
  return(modeltofit);
}

void csm::plotcompare(char *cname)
{
  Int_t i;
  for (i=0;i<(Int_t)datatofitname.size();i++)
    {
      if (strcmp(datatofitname[i],cname)==0)
	{
	  modeltofit->plotwithdata(cname,datatofit[i]);
	}
    }
}

// Number of degrees of freedom -- this is approximately true for
// large statistics (in fact, the whole chisquared idea is only approximately
// true in cases of large statistics where distributions are approximately Gaussian)
// Degrees of freedom "freeze out" as the expected number of events gets small
// (<5 or so).  A bin with no data and no expectation shouldn't contribute either
// to the chisquared or the number of degrees of freedom, and neither really should
// a bin with 1E-6 expected and no observed events.  This routine won't draw the
// line (and even interpolated histograms can have variable numbers of bins with
// zero expectation).  This routine's very naive and just counts bins, filled or not.

Int_t csm::ndof()
{
  Int_t ndofl,i;
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;

  ndofl = 0;
  for (i=0;i<(Int_t) datatofit.size();i++)
    {
      ndofl += datatofit[i]->GetNbinsX()*datatofit[i]->GetNbinsY();
    }
  modeltofit->list_nparams(&npn, &nplb, &nphb);
  ndofl -= npn.size();
  cout << "nDOF isn't very clearly defined here... todo" << endl;
  return(ndofl);
}


Int_t csm::getnparams()
{
  return(fitparam.size());
}

Double_t csm::getparam(Int_t iparam)
{
  return(fitparam[iparam]);
}

Double_t csm::getcov(Int_t iparam, Int_t jparam)
{
  Int_t nparams = fitparam.size();
  return(fitcov[iparam+nparams*jparam]);
}

Double_t csm::getperror(Int_t iparam)
{
  return(fiterror[iparam]);
}

char* csm::getpname(Int_t iparam)
{
  return(fitparamname[iparam]);
}

// Call the individual channel chisquared calculators inside here.
// the parameters par are labeled by their names npfitname in the static global vector.

void csm_minuit_fcn(Int_t &npar, Double_t */*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
  Int_t i;

  //adjust the model according to the nuisance paramters
  
  char **fpnameloc = new char *[npfitname.size()];
  for (i=0;i<(Int_t) npfitname.size();i++)
    {
      fpnameloc[i] = npfitname[i];
      //      cout << "in minuit fit fcn: " << i << " " << npfitname[i] << endl;
    }
  modeltofit->nuisance_response(npfitname.size(),fpnameloc,par);

  TH1** dfloc = new TH1*[datatofit.size()];
  for (i=0;i<(Int_t) datatofit.size();i++)
    {
      dfloc[i] = datatofit[i];
    }

  //cout << "In minuit function: printing out the model" << endl;
  //modeltofit->print();

  f = modeltofit->chisquared1(dfloc);

  // Gaussian constraints for variables which are constrained.

  for (i=0;i<npar;i++) 
    {
      if (constrainedfitparam[i] != 0)
        {
	  //cout << "In fcn: " << i << " " << par[i] << endl;
          f += par[i]*par[i];
        }
    }

  //cout << "end of computation of f in minuit_fit_fcn: " << f << endl;
  //printf("%20.14f\n",f);

  delete[] fpnameloc;
  delete[] dfloc;
}

/*------------------------------------------------------------------------*/
// make a plot of the results, along with some data

void csm_channel_model::plotwithdata(TH1* dh)
{
  Int_t i,ntemplates,nbinsy;
  Double_t stackmax,datamax,plotmax;
  THStack *hs = new THStack("hs",dh->GetTitle());
  ntemplates = (Int_t) histotemplate.size();
  TLegend *slegend = (TLegend*) new TLegend(0.7,0.6,0.89,0.89);

  for (i=0;i<ntemplates;i++)
    {
      std::auto_ptr<TH1> htl = copy_TH1<TH1> (*histotemplate_varied[i]);
      htl->Scale(sft_varied[i]);
      htl->SetFillColor(i+40);
      htl->SetFillStyle(1001);
      htl->SetLineColor(kBlack);
      htl->SetLineWidth(1);
      hs->Add(htl.release());
    }

  TList *hlist = hs->GetHists();
  TObjLink *lnk = hlist->LastLink();          
  while (lnk)
    {  slegend->AddEntry(lnk->GetObject(),lnk->GetObject()->GetName(),"F");
    lnk = lnk->Prev();                       
    }     
  // make sure the plot is big enough to fit the data, the model stack,
  // and the data error bars with a little room to spare
  stackmax = hs->GetMaximum();
  datamax = dh->GetMaximum();
  nbinsy = dh->GetNbinsY();
  datamax += sqrt(datamax);
  plotmax = max(datamax,stackmax);
  hs->SetMaximum(plotmax);
  dh->SetMaximum(plotmax);
  hs->SetMinimum(0);
  dh->SetMinimum(0);
  if (nbinsy==1)
    {
      hs->Draw("HIST");
      dh->SetMarkerStyle(20);
      dh->SetLineColor(kBlack);
      dh->SetMarkerColor(kBlack);
      dh->DrawCopy("E0SAME");
    }
  else
    {
      hs->Draw();
      dh->SetMarkerStyle(20);
      dh->SetLineColor(kBlack);
      dh->SetMarkerColor(kBlack);
      dh->DrawCopy("LEGO,SAME");
    }
  slegend->AddEntry(dh,dh->GetName(),"P");
  slegend->SetHeader(dh->GetTitle());
  slegend->Draw();
}

void csm_channel_model::candcheck(TH1 *dh)
{
  cout << dh->GetTitle() << " Candidate Check " << endl;

  Double_t sumsb = 0;
  Double_t ssum = 0;
  Double_t bsum = 0;
  Double_t dsum = 0;
  Int_t nbinsx = histotemplate[0]->GetNbinsX();
  Int_t nbinsy = histotemplate[0]->GetNbinsY();
  if (nbinsx != dh->GetNbinsX())
    {
      cout << "data histogram and model template have different numbers of x bins: " <<
	nbinsx << " != " << dh->GetNbinsX() << endl;
      return;
    }
  if (nbinsy != dh->GetNbinsY())
    {
      cout << "data histogram and model template have different numbers of y bins: " <<
	nbinsy << " != " << dh->GetNbinsY() << endl;
      return;
    }

  Int_t ibinx,ibiny;

  Double_t hcb[nbinsx][nbinsy];
  Double_t hcs[nbinsx][nbinsy];
  for (ibinx=0;ibinx<nbinsx;ibinx++)
    {
      for (ibiny=0;ibiny<nbinsy;ibiny++)
	{
	  hcb[ibinx][ibiny] = 0;
	  hcs[ibinx][ibiny] = 0;
	}
    }

  for(Int_t ic=0;ic < (Int_t) histotemplate.size();ic++)
    {
      for (ibinx=0;ibinx<nbinsx;ibinx++)
	{
	  for (ibiny=0;ibiny<nbinsy;ibiny++)
	    {
	      if (scaleflag[ic])
		{
		  hcs[ibinx][ibiny] += sft_varied[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1);
		  if (hcs[ibinx][ibiny] < 0) 
		    {
		      cout << "Negative signal expectation (" << ibinx+1 << "," << ibiny+1 << "):" << hcs[ibinx][ibiny] << endl;
		    }
		}
	      else
		{
		  hcb[ibinx][ibiny] += sft_varied[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1);
		  if (hcb[ibinx][ibiny] < 0) 
		    {
		      cout << "Negative background expectation (" << ibinx+1 << "," << ibiny+1 << "):" << hcb[ibinx][ibiny] << endl;
		    }
		}
	    }
	}
    }

  for (ibinx=0;ibinx<nbinsx;ibinx++)
    {
      for (ibiny=0;ibiny<nbinsy;ibiny++)
	{
	  ssum += hcs[ibinx][ibiny];
	  bsum += hcb[ibinx][ibiny];
          Double_t dc = dh->GetBinContent(ibinx+1,ibiny+1);
	  //cout << "*@*" << hcs[ibinx][ibiny] << " " << hcb[ibinx][ibiny] << " " << dc << endl;
	  dsum += dc;
	  if (hcs[ibinx][ibiny]>0 && hcb[ibinx][ibiny]<=0)
	    {
	      cout << "Null background with expected signal: (" << ibinx+1 << "," << ibiny+1 
                   << ") Cands: " << dc << " Signal: " << hcs[ibinx][ibiny] << endl;
	    }
	  else
	    {
	      if (dc > 0)
		{
		  if (hcs[ibinx][ibiny]<=0 && hcb[ibinx][ibiny]<=0)
		    {
		      cout << "Null background with observed candidate(s): (" << ibinx+1 << "," << ibiny+1 
			   << ") Cands: " << dc << " Signal: " << hcs[ibinx][ibiny] << endl;
		    }
		  Double_t sbratio = hcs[ibinx][ibiny]/hcb[ibinx][ibiny];
		  sumsb += dc*sbratio;
		  if (sbratio>0.3)
		    {
		      cout << "High s/b candidate(s): (" << ibinx+1 << "," << ibiny+1 << ") cands: " << dc << " s/b: " << sbratio << 
			" s: " << hcs[ibinx][ibiny] << " b: " << hcb[ibinx][ibiny] << endl;
		    }
		}
	    }
	}
    }
  cout << "S/B sum over all candidates: " << sumsb << endl;
  cout << "S sum over all bins: " << ssum << endl;
  cout << "B sum over all bins: " << bsum << endl;
  cout << "D sum over all bins: " << dsum << endl;
}

double csm_channel_model::kstest(TH1* dh)
{
  Int_t i,ntemplates;
  ntemplates = (Int_t) histotemplate.size();
  double tout;

  std::auto_ptr<TH1> hsum = copy_TH1<TH1> (*histotemplate_varied[0]);
  hsum->Sumw2();
  hsum->Reset();


  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i];
      hsum->Add(htl,sft_varied[i]);
    }

  tout = hsum->KolmogorovTest(dh);
  return(tout);
}

double csm_channel_model::kstest_px(TH1* dh)
{
  Int_t i,ntemplates;
  ntemplates = (Int_t) histotemplate.size();
  double tout;

  std::auto_ptr<TH1> hsum = copy_TH1<TH1> (*histotemplate_varied[0]);
  hsum->Sumw2();
  hsum->Reset();


  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i];
      hsum->Add(htl,sft_varied[i]);
    }

  tout = hsum->KolmogorovTest(dh,"X");
  return(tout);
}

/*------------------------------------------------------------------------*/

// and a method to allow an object of type csm_model to plot up one of its
// channels with some data compared.

void csm_model::plotwithdata(char* cname, TH1* dh)
{
  Int_t i;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  chanmodel[i]->plotwithdata(dh);
	}
    }
}
 
/*------------------------------------------------------------------------*/

// check candidates

void csm_model::candcheck(char* cname, TH1* dh)
{
  Int_t i;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  chanmodel[i]->candcheck(dh);
	}
    }
}
 
/*------------------------------------------------------------------------*/

double csm_model::kstest(char* cname, TH1* dh)
{
  Int_t i;
  double ksresult=0; 
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  ksresult = chanmodel[i]->kstest(dh);
	}
    }
  return(ksresult);
}
 
/*------------------------------------------------------------------------*/

double csm_model::kstest_px(char* cname, TH1* dh)
{
  Int_t i;
  double ksresult=0;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  ksresult = chanmodel[i]->kstest_px(dh);
	}
    }
  return(ksresult);
}
 
/*------------------------------------------------------------------------*/

/* chisquared1 evaluates a chisquared function in the style of T. Devlin's CDF 3126, eq's 9 and 10
   The signal is a sum of signal contributions and the background is a sum of
   background contributions.  This chisquared function is meant to be minimized 
   with respect to the free nuisance parameters 

   This version does one 1D or 2D histogram at a time.

   This version allows for multiple sources of signal and multiple sources of background,
   some of each of which are estimated using finite MC or data statistics in each bin.
   This routine does not distinguish between a signal source and a background source --
   finding the chisquared of a data distribution to a sum of models does not need a distinction
   at this level.  Instead, one may compare the chisquared of the same data against collections
   of models that include signals and those that do not include signals, calling this routine
   twice (or more times).

   CDF 3126 describes how to minimize the chisquared function over each bin's uncertain
   Poisson-constrained rates.  When multiple sources are allowed to be estimated from Poisson
   distributed subsidiary measurements, the quadratic polynomial to be solved for turns
   into a system of quadratic equations which is solved here iteratively.

   This function is meant to be part of a MINUIT minimization over the nuisance
   parameters.


   input: TH1 *dh -- data histogram to compare the channel's model against
            
   output:  chi squared, the function value.

   Update 5 July 2006 -- Reading Barlow and Beeston about bins with zero MC prediction in one
   or more source.  Take the one with the strongest contribution (here taken from the normalization
   scale factors), and set the others to zero when solving the n coupled quadratic equations.

   Update 8 Dec, 2007 -- put in the error bars on the template histograms when the flag is
   CSM_GAUSSIAN_BINERR (new feature) as if they were Poisson (i.e., no different terms in
   the likelihood function -- they're probably Poisson underneath anyhow, but more often, they
   are a complicated mixture of MC or data events with different weights, and all treatments
   of them are approximations).  

*/

Double_t csm_channel_model::chisquared1(TH1 *dh)
{
  Double_t chi2;
  Int_t ip1,ic,ibinx,ibiny,iter,iprec;
  Double_t A,B,C,D;
  Double_t csum,cpsum,gbc,gbe;
  Int_t nbinsx,nbinsy;
  Double_t nsubs;
  Int_t dtb;  // data observed in a single bin
  Int_t nc;

  // number of template histograms
  nc = (Int_t) histotemplate.size();

  Int_t lpoissflag[nc];  // local Poisson flag -- if error is zero for a template in a bin
                         // reclassify it as no bin error

  // allocate rho1 and rho2 for all templates, even though we're only going to need
  // them for the Poisson-distributed ones
  // push them on the stack -- 8 dec 2007

  Double_t rho1[nc];
  Double_t rho2[nc];
  Int_t zlist[nc];
  Double_t sfgp[nc];  // scale factor needed to approximate Gaussian errors as Poisson

  nbinsx = dh->GetNbinsX();
  nbinsy = dh->GetNbinsY();

  chi2 = 0;

  for (ibinx=0;ibinx<nbinsx;ibinx++)
    {
      for (ibiny=0;ibiny<nbinsy;ibiny++)
        {
	  if (nbinsy==1)
	    { dtb = (Int_t) dh->GetBinContent(ibinx+1); }
	  else
	    { dtb = (Int_t) dh->GetBinContent(ibinx+1,ibiny+1); }
          //cout << "In chi2calc: " << ibinx << " " << ibiny << " " << dtb << endl;

	  // if we are told to pay attention to the error but it's zero, reclassify it locally
	  // as a zero-error bin in this template

	  for (ic=0;ic<nc;ic++)
	    {
	      sfgp[ic] = 1.0;
	      lpoissflag[ic] = poissflag[ic];
	      if (nbinsy == 1)
	        { 
		  gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); 
		  gbe = histotemplate_varied[ic]->GetBinError(ibinx+1); 
		}
	      else
	        {
		  gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); 
		  gbe = histotemplate_varied[ic]->GetBinError(ibinx+1,ibiny+1); 
		}
	      if (poissflag[ic] == CSM_GAUSSIAN_BINERR && gbe == 0)
		{ lpoissflag[ic] = CSM_NOBINERR; }
	      if (lpoissflag[ic] == CSM_GAUSSIAN_BINERR)
		{
		  sfgp[ic] = gbc/(gbe*gbe);
		  lpoissflag[ic] = CSM_POISSON_BINERR;
		}
	    }

	  /* the sum of zero-bin-error contributions, varied by the nuisance parameters */
          csum = 0;
          for (ic=0;ic<nc;ic++)
 	    {
	      if (lpoissflag[ic] == CSM_NOBINERR)
	        {
		  if (nbinsy == 1)
		    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
	          csum += gbc*sft_varied[ic];
	        }
	    }
	  if (csum < 0) 
	    { chi2 += 1E10; }

	  /* solve for the rho's in each bin for each source of Poisson-estimated model rate
	     rho1 is the current estimate used to compute the rho2's.  On each iteration,
	     copy the previous iteration's rho2 into the rho1 array and re-solve for rho2.
	     start with nominal central values from the subsidiary measurements */

	  int haszero = 0;
          int im1=-1;
          double xm1=0;
          for (ic=0;ic<nc;ic++) 
	    { 
	      rho1[ic] = 0;
  	      rho2[ic] = 0;
	      if (lpoissflag[ic] == CSM_POISSON_BINERR)
		{
		  if (nbinsy == 1)
		    { gbc = max(0,histotemplate_varied[ic]->GetBinContent(ibinx+1)); }
		  else
		    { gbc = max(0,histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1)); }
		  rho2[ic] = gbc*sft_varied[ic];
		  rho1[ic] = rho2[ic];
		  if (gbc == 0 || sft_varied[ic] == 0)
		    { 
		      haszero = 1;
		      zlist[ic] = 1;
		      if (sft_varied[ic]>xm1)
			{ 
			  xm1 = sft_varied[ic]; 
			  im1 = ic;
			}
		    }
		  else
		    {
		      zlist[ic] = 0;
		    }
		}
	    }
	  if (haszero != 0 && im1 > -1)
	    {
	      zlist[im1] = 0;
	    }

          for (iter=0;iter<CSM_MAXITER;iter++)
	    {
              for (ic=0;ic<nc;ic++) 
                { 
		  if (lpoissflag[ic] == CSM_POISSON_BINERR)
		    {
		      rho1[ic] = rho2[ic];
		      if (zlist[ic] == 1)
			{ 
			  rho1[ic] = 0;
			}
		    }
                }

	      for (ic=0;ic<nc;ic++)
	        {
	          if (lpoissflag[ic] == CSM_POISSON_BINERR)
	            {
		      if (zlist[ic] == 0)
			{
	                  D = csum;
                          for (ip1=0;ip1<nc;ip1++)
		            { if ( (lpoissflag[ip1]==CSM_POISSON_BINERR) && ip1 != ic) D += rho1[ip1]; }
	                  A = 1.0 + sfgp[ic]/sft_varied[ic];
		          if (nbinsy == 1)
			    { gbc = sfgp[ic]*max(0,histotemplate_varied[ic]->GetBinContent(ibinx+1)); }
		          else
			    { gbc = sfgp[ic]*max(0,histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1)); }
		          B = A*D - dtb - gbc;
		          C = -gbc*D;
	                  rho2[ic] = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
			  //cout << "ABC: " << A << " " << B << " " << C << endl;
			}
		      else // a la Barlow and Beeston, set only one prediction to nonzero
			{                    // if we have zero MC -- the "strongest" one among all the contributions
			  rho2[ic] = 0;     // with zero MC prediction
			}
		    }
	        }
	      iprec = 0;

	      for (ic=0;ic<nc;ic++)
	        {
	          if (lpoissflag[ic] == CSM_POISSON_BINERR)
	            {
	              if (fabs(rho1[ic]) < PREC1)
		        { 
		          if (fabs(rho2[ic]-rho1[ic]) > PREC1)
		            { 
                              iprec = 1;
		              break;
		            }
		        }
	              else
		        {
		          if ( fabs((rho2[ic]-rho1[ic])/rho1[ic])>PREC1 )
		            {
		              iprec = 1;
		              break;
		            }
		        }
	            }
	        }
              if (iprec == 0) break;
  	    }  /* end loop over iterations to compute the rho's.  rho2 is the computed array */

	  /*
	    if (CSM_DEBUGPRINT >0 && iprec ==1)
	    {
	    // cout << "csm_chisquared1: iterations failed to converge " << endl;
	    cout << "In chi2calc: " << ibinx << " " << ibiny << " " << dtb << endl;
	    for (ic=0;ic<nc;ic++)
	    {
	    if (nbinsy == 1)
	    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); }
	    else
	    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
	    if (lpoissflag[ic] == CSM_POISSON_BINERR)
	    {
	    cout << "Poisson contrib " << ic << " " << gbc << " " << sft_varied[ic] << endl;
	    }
	    else
	    {
	    cout << "Non-poisson contrib " << ic << " " << gbc << " " << sft_varied[ic] << endl;
	    }
	    }
	    }
	  */

	  // When the iterations fail to converge, it is usually an oscillatory
	  // solution.  Pick the rho1 or the rho2 array which minimizes chisquare
	  // first compute the chisquare using the rho2 array, and if we need to,
          // redo it with the rho1 array, and pick the smaller of the two.

	  Double_t chi2a = chi2;

	  cpsum = csum;
	  for (ic=0;ic<nc;ic++)
	    {
	      if (lpoissflag[ic] == CSM_POISSON_BINERR)
		{
		  cpsum += rho2[ic];
		}
	    }

	  if (cpsum>0)
	    {
	      chi2a += cpsum;
	      chi2a -= dtb;
	      if (dtb>0) {chi2a -= dtb*log(cpsum/((Double_t) dtb));}
	    }
	  else if (dtb>0)
	    { chi2a += 1E10; }

	  for (ic=0;ic<nc;ic++)
	    {
	      if ( (lpoissflag[ic] == CSM_POISSON_BINERR) && sft_varied[ic] > 0)
		{
		  if (nbinsy == 1)
		    { nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }

#ifdef DEBUGPRINTC2
		  Double_t c2cont = ( rho2[ic]*sfgp[ic]/sft_varied[ic] - nsubs);
		  if (nsubs > 0 && rho2[ic] > 0)
		    c2cont -= ((Double_t) nsubs)*log(rho2[ic]*sfgp[ic]/(sft_varied[ic]*((Double_t) nsubs)));
		  if (c2cont > 0.0)
		    //  if (ibinx == 9 && ibiny == 4)
		    {
		      cout << "in chi2calc: " << ibinx << " " << ibiny << " " << ic << " " << 
			rho2[ic] << " " << sft_varied[ic] << " " << nsubs << " " << c2cont;
		      if (c2cont>0.01) 
			{
			  cout << "*" << endl;
			}
		      else
			{
			  cout << endl;
			}
		    }
#endif
		  chi2a += (rho2[ic]*sfgp[ic]/sft_varied[ic] - nsubs);
		  if (nsubs > 0 && rho2[ic] > 0)
		    chi2a -= ((Double_t) nsubs)*log(rho2[ic]*sfgp[ic]/(sft_varied[ic]*((Double_t) nsubs)));
		}
	    }

	  Double_t chi2b = chi2;

	  if (iprec == 1)
	    {
	      cpsum = csum;
	      for (ic=0;ic<nc;ic++)
		{
		  if (lpoissflag[ic] == CSM_POISSON_BINERR)
		    {
		      cpsum += rho1[ic];
		    }
		}

	      if (cpsum > 0)
		{
		  chi2b += cpsum;
		  chi2b -= dtb;
		  if (dtb>0) {chi2b -= dtb*log(cpsum/((Double_t) dtb));}
		}
	      else if (dtb>0)
		{ chi2b += 1E10; }

	      for (ic=0;ic<nc;ic++)
		{
		  if ( (lpoissflag[ic] == CSM_POISSON_BINERR) && sft_varied[ic] > 0)
		    {
		      if (nbinsy == 1)
			{ nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		      else
			{ nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }

		      /*
			cout << "in chi2calc: " << ibinx << " " << ibiny << " " << ic << " " << 
			rho1[ic] << " " << sft_varied[ic] << " " << nsubs << " " sfgp[ic] << endl;
		      */
		      chi2b += (rho1[ic]*sfgp[ic]/sft_varied[ic] - nsubs);
		      if (nsubs > 0 && rho1[ic] > 0)
			chi2b -= ((Double_t) nsubs)*log(rho1[ic]*sfgp[ic]/(sft_varied[ic]*((Double_t) nsubs)));
		    }
		}
	    }
          else
	    {
	      chi2b = chi2a;
	    }
	  chi2 = min(chi2a,chi2b);
        } /* end loop over binsy */
    } /* end loop over binsx */

  chi2 *= 2.0;

  //cout << "chisquared calc: " << chi2 << endl;
  //if ( !(chi2<0 || chi2>=0))
  //  {
  //    cout << "Bad chi2: " << endl;
  //    exit(0);
  //  }

  return(chi2);

}

Double_t csm_model::chisquared1(TH1 **dh)
{
  Int_t i;
  Double_t cs;
  cs = 0;
  for (i=0;i<(Int_t)chanmodel.size();i++)
    {
      cs += chanmodel[i]->chisquared1(dh[i]);
    }
  return(cs);
}

/*------------------------------------------------------------------------*/
//  Set the interpolation style for a particular channel.  Two methods    
//  one for channel models, and one if you just have a pointer to a csm_model
/*------------------------------------------------------------------------*/


void csm_model::set_interpolation_style(char *cname, INTERPSTYLE istyle)
{
  Int_t i;
  for (i=0;i<(Int_t) channame.size(); i++)
    {
      if (strcmp(channame[i],cname)==0)
	{
	  chanmodel[i]->set_interpolation_style(istyle);
	}
    }
}

void csm_channel_model::set_interpolation_style(INTERPSTYLE istyle)
{
  chan_istyle = istyle;
}

/*------------------------------------------------------------------------*/

/* compounded interpolation -- dist1 = central value shape, dist2 = syst. varied shape,
   dist3 = shape to distort (may be the result of previous distortions for compounded shape
   variations), distn = resultant shape.  par1 = value of parameter for dist1.  par2 = value of
   parameter (like # of sigma) for dist2.  parn = value of parameter for the output histogram
   Built on the idea of d_pvmorph, but generalized a bit.  Returns a null histogram if any of
   the three input histograms has zero or negative total sum.
*/

//#define DEBUGPVMC

void csm_pvmc(Int_t nb, Double_t *dist1, Double_t *dist2, Double_t *dist3, Double_t *distn,
	      Double_t par1, Double_t par2, Double_t parn)
{
  Int_t nb3 = nb*3 + 3;
  Int_t i,j,k,k1,k2;
  Double_t total;
  Double_t wta,wtb;
  Double_t yd[nb3];
  Double_t id[nb3];
  Double_t xd[nb3];
  Double_t xdis[nb3];
  Double_t ydis[nb3];
  Double_t ydi[nb + 1];
  Int_t idx[nb3];
  Int_t ifirst;
  Double_t x1l,x2l,x3l,y1l,y2l,y3l;
  Double_t x1,x2,x3,y1,y2,y3;
  Double_t xloc,yloc,x1i,x2i,x3i;

  // default output -- empty distribution

  for (i=0;i<nb;i++)
    {
      distn[i] = 0.0;
    }

  // default index list

  for (i=0;i<nb3;i++)
    { idx[i] = i;}

  // parameter weights

  if (par2 != par1) 
    {
      wta = 1. - (parn-par1)/(par2-par1);
      wtb = 1. + (parn-par2)/(par2-par1);
    }
  else
    {
      wta = 0.5;
      wtb = 0.5;
    }

  // suppress warning messages in case of extrapolations

  //  if ( (parn>par1 && parn>par2) || (parn<par1 && parn<par2) )
  //  {
  //    cout << "CSM_PVMC: Histogram Extrapolation: " << parn << 
  //            " is not between " << par1 << " and " << par2 << endl;
  //  }

  // Fill cumulative distribution arrays -- squeeze out repeated entries
  // due to empty bins

  // The first point in the cumulative distributions has zero integral and
  // starts at the left-hand edge of the first bin with any value in it.
  // The id array says which distribution it came from, and the
  // xd array gives the x value at which the cumulative distribution is evaluated
  // (at the right-hand edge of the bin)

  j = 0;
  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist1[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist1 in csm_pvmc" << endl;
	  cout << i << " " << dist1[i] << endl;
	  exit(0);
	}
      total += dist1[i];
    }
  if (total <= 0) return;

  yd[j] = 0;
  id[j] = 1;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist1[i] > 0)
	{ 
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
	  yd[j] = yd[j-1] + dist1[i]/total;
	  id[j] = 1;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist2[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist2 in csm_pvmc" << endl;
	  cout << i << " " << dist2[i] << endl;
	  exit(0);
	}
      total += dist2[i];
    }
  if (total <= 0) return;
  yd[j] = 0;
  id[j] = 2;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist2[i]>0)
	{
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
          yd[j] = yd[j-1] + dist2[i]/total;
	  id[j] = 2;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist3[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist3 in csm_pvmc" << endl;
	  cout << i << " " << dist3[i] << endl;
	  exit(0);
	}
      total += dist3[i];
    }
  if (total <= 0) return;
  yd[j] = 0;
  id[j] = 3;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist3[i]>0)
	{
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
          yd[j] = yd[j-1] + dist3[i]/total;
	  id[j] = 3;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  // Sort all of the edges of the cumulative distribution functions
  // j is the number of entries in the yd, xd and id arrays

  TMath::Sort(j,yd,idx,0);

#ifdef DEBUGPVMC
  for (i=0;i<j;i++)
    {
      cout << i << " " << xd[i] << " " << " " << yd[i] << " " << id[i] << endl;
    }
  cout << "Sort index" << endl;
  for (i=0;i<j;i++)
    {
      cout << i << " " << idx[i] << endl;
    }
#endif

  x1l = 0;
  x2l = 0;
  x3l = 0;

  y1l = 0;
  y2l = 0;
  y3l = 0;

  x1 = 0;
  x2 = 0;
  x3 = 0;

  y1 = 0;
  y2 = 0;
  y3 = 0;

  // the three lowest points in the sort should all have zero integral --
  // interpolate the x's of these

  for (i=0;i<3;i++)
    {
      if ( id[idx[i]] == 1 )
	{
	  x1 = xd[idx[i]];
	  y1 = yd[idx[i]]; // should be zero
	}
      else if ( id[idx[i]] == 2 )
	{
	  x2 = xd[idx[i]];
	  y2 = yd[idx[i]];  // should be zero
	}
      else if ( id[idx[i]] == 3 )
	{
	  x3 = xd[idx[i]];
	  y3 = yd[idx[i]];  // should be zero
	}
    }
  // don't have the other ends of the line segments yet -- find them as we go along.

#ifdef DEBUGPVMC
  cout << "first bins: " << x1 << " " << x2 << " " << x3 << endl;
#endif
  y1l = y1;
  y2l = y2;
  y3l = y3;
  x1l = x1;
  x2l = x2;
  x3l = x3;

  // first point on interpolated curve -- zero integral.

  k = 0;
  xdis[k] = wta*x1l + wtb*x2l - x1l + x3l;
  xdis[k] = min((Double_t) (nb+1),max(0.0,xdis[k]));
  ydis[k] = 0;

#ifdef DEBUGPVMC
  cout << "first point: " << xdis[0] << " " << ydis[0] << endl;
#endif

  for (i=3;i<j;i++)
    {
      xloc = xd[idx[i]];
      yloc = yd[idx[i]];

      if (id[idx[i]] == 1 )
	{
	  x1l = x1;
	  y1l = y1;
	  x1 = xloc;
	  y1 = yloc;

	  if (yloc>y2)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==2)
		    {
		      y2l = y2;
		      x2l = x2;
		      y2 = yd[idx[k1]];
		      x2 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y3)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==3)
		    {
		      y3l = y3;
		      x3l = x3;
		      y3 = yd[idx[k1]];
		      x3 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}
      else if (id[idx[i]] == 2 )
	{
	  x2l = x2;
	  y2l = y2;
	  x2 = xloc;
	  y2 = yloc;

	  if (yloc>y1)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==1)
		    {
		      y1l = y1;
		      x1l = x1;
		      y1 = yd[idx[k1]];
		      x1 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y3)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==3)
		    {
		      y3l = y3;
		      x3l = x3;
		      y3 = yd[idx[k1]];
		      x3 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}
      else if (id[idx[i]] == 3 )
	{
	  x3l = x3;
	  y3l = y3;
	  x3 = xloc;
	  y3 = yloc;

	  if (yloc>y2)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==2)
		    {
		      y2l = y2;
		      x2l = x2;
		      y2 = yd[idx[k1]];
		      x2 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y1)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==1)
		    {
		      y1l = y1;
		      x1l = x1;
		      y1 = yd[idx[k1]];
		      x1 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}

      if (yloc>ydis[k] && ydis[k] < 0.999999999)
        {

#ifdef DEBUGPVMC
	  cout << "Interpolating: " << x1 << " " << x2 << " " << x3 << endl;
	  cout << "Interpolating: " << x1l << " " << x2l << " " << x3l << endl;
	  cout << "Interpolating: " << y1 << " " << y2 << " " << y3 << endl;
	  cout << "Interpolating: " << y1l << " " << y2l << " " << y3l << endl;
#endif

	  k++;
	  ydis[k] = yloc;
	  if (yloc == y1l)
	    {
	      x1i = x1l;
	    }
	  else
	    {
	      x1i = x1l + (yloc-y1l)*(x1-x1l)/(y1-y1l);
	    }
	  if (yloc == y2l)
	    {
	      x2i = x2l;
	    }
	  else
	    {
	      x2i = x2l + (yloc-y2l)*(x2-x2l)/(y2-y2l);
	    }
	  if (yloc == y3l)
	    {
	      x3i = x3l;
	    }
	  else
	    {
	      x3i = x3l + (yloc-y3l)*(x3-x3l)/(y3-y3l);
	    }
	  xdis[k] = x3i + wta*x1i + wtb*x2i - x1i;
          xdis[k] = min((Double_t) (nb+1),max(0.0,xdis[k]));
	  if (xdis[k]<xdis[k-1])
	    {
	      k--;
	      ydis[k] = yloc;
	    }

#ifdef DEBUGPVMC
	  cout << "point: " << k << endl;
	  cout << "x1i, x2i, x3i: " << x1i << " " << x2i << " " << x3i << endl;
	  cout << "Interpolated: " << xdis[k] << " " << yloc << endl;
#endif
	}
    }

#ifdef DEBUGPVMC
  for (i=0;i<=k;i++)
    {
      cout << "IC before bin: " << i << " " << xdis[i] << " " << ydis[i] << endl;
    }
#endif


  // k is the index of the last entry in the xdis, ydis interpolated array.
  // find the places where the piecewise linear interpolated cumulative distribution
  // crosses the bin edges  the index on ydi is the low bin edge.

  ydi[0] = 0.0;
  for (i=1;i<(1 + nb);i++)
    {
      ydi[i] = 1.0;
    }

  Int_t k2last;
  k1 = 0;
  for (i=0;i<=k;i++)
    {
      k2last = k1;
      for (k2=k1+1;k2<(nb+1);k2++)
	{
	  if ( (Double_t) k2 < xdis[i] )
	    {
	      if (i==0)
		{
		  ydi[k2] = 0;
		}
	      else
		{
		  ydi[k2] = ydis[i-1] + ( (Double_t) k2  - xdis[i-1] )*
		    (ydis[i]-ydis[i-1])/(xdis[i]-xdis[i-1]);
#ifdef DEBUGPVMC
		  cout << "filling bins: " << i << " " << k1 << " " << k2 << " " << ydi[k2] << endl; 
#endif

		}
	      k2last = k2;
	    } 
	  if ( (Double_t) k2 > xdis[i] ) break;
	}
      k1 = k2last;
    }

#ifdef DEBUGPVMC
  for (i=0;i<(nb+1);i++)
    {
      cout << "interp. cumulative: " << i << " " << ydi[i] << endl;
    }
#endif

  // differentiate to get the output distn

  for (i=0;i<(nb);i++)
    {
      distn[i] = ydi[i+1] - ydi[i]; 
    }
}

/*------------------------------------------------------------------------*/


/*  
    Re-coded version of d_pvmorph_2d from Alex Read.  C version from Tom Junk
    Added feature of compounding shape variations as systematic errors.
    February 2007

    ......Do a linear interpolation between three two-dimensional
    probability distributions (scatterplots) as a function
    of the characteristic parameter of the distribution, for use
    in both interpolation and in application of systematic uncertainties.
    xydist1 is the "central value" histogram
    xydist2 is the "systematically varied" histogram
    xydist3 is the histogram to apply the variation to
    xydistn is the output histogram.  See csm_pvmc
    for compounded systematic variation applicaiton in 1D
    (used repeatedly in here).

    This is a generalization of csm_pvmc. The 2d distribution
    can be move around and be stretched or squeezed in two
    dimenions but finite rotations (changes in the correlation)
    are poorly approximated.

    nx        : Number of x-bins in the input and output distributions.
    ny        : Number of y-bins in the input and output distributions.
    xydist1,xydist2,xydist3,xydistn
    : Bin contents of scatterplots. The arrays should be
    packed with the index running fastest over the x
    dimension.
    Contents are in xydist[ix+nx*iy]
    par1,par2,parn     : Values of the linear parameters that characterise the
    histograms (e.g. the Higgs mass).

    Output: xydistn.  Same binning as xydist1,xydist2,xydist3.
    Its memory must be allocated outside of the routine

*/

//#define DEBUGPVMC2D

void csm_pvmc2d(Int_t nx, Int_t ny, Double_t *xydist1, 
                Double_t *xydist2, Double_t *xydist3, Double_t *xydistn,
                Double_t par1, Double_t par2, Double_t parn)
{
  Double_t ydist1[ny],ydist2[ny],ydist3[ny],ydistn[ny];
  Double_t xtemp1[nx],xtemp2[nx],xtemp3[nx],xtempn[nx];
  Double_t alpha1[ny*ny],alpha2[ny*ny],alpha3[ny*ny];
  Int_t i,j,k;

  // Project xydist1,2,3 onto the y axis and normalize

  csm_yproj(nx,ny,xydist1,ydist1);
  csm_yproj(nx,ny,xydist2,ydist2);
  csm_yproj(nx,ny,xydist3,ydist3);

  // Interpolate the y-projections 

  csm_pvmc(ny,ydist1,ydist2,ydist3,ydistn,par1,par2,parn);

#ifdef DEBUGPVMC2D
  for (i=0;i<ny;i++)
    {
      cout << "iy: " << i << " " << ydist1[i] << " " << ydist2[i] << " " << ydistn[i] << endl;
    }
#endif

  // Find out which y bins of histograms 1,2,3 contribute
  // to each bin of the interpolated distribution ydistn

  csm_ycont(ny,ydist1,ydist2,ydist3,ydistn,alpha1,alpha2,alpha3);

  // Extract the x-distributions in the y-slice determined above
  // and interpolate them

  for (i=0;i<ny;i++)  // loop over resulting bins
    {
      for (k=0;k<nx;k++) xtemp1[k] = 0;
      for (k=0;k<nx;k++) xtemp2[k] = 0;
      for (k=0;k<nx;k++) xtemp3[k] = 0;
      for (j=0;j<ny;j++) // loop over contributing bins
	{
	  for (k=0;k<nx;k++)
	    {
	      xtemp1[k] += alpha1[j+ny*i]*xydist1[k+nx*j];
	      xtemp2[k] += alpha2[j+ny*i]*xydist2[k+nx*j];
	      xtemp3[k] += alpha3[j+ny*i]*xydist3[k+nx*j];
	    }
	}
      // Interpolate the x distributions
      csm_pvmc(nx,xtemp1,xtemp2,xtemp3,xtempn,par1,par2,parn);

      // Insert the interpolated x distribution into the final output dist
      for (k=0;k<nx;k++) xydistn[k+nx*i] = xtempn[k]*ydistn[i];
    }
}

/*
  Re-coded version of d_ypvscat from Alex Read.  C version from Tom Junk
  Project a scatterplot onto the y-axis.The
  projection is normalized so that the sum of the bin contents is 1.0.

  nx,ny    : Number of bins in the scatterplot for the x and y coordindates.
  The projection is done onto <ny> bins.
  xydist   : The 2-dimensional array of the probabilities
  ydist    : The 1-dimensional array of the 2d probabilities projected onto
  the y-axis.

  Inputs : nx,ny,xydist
  Outputs: ydist (ny is unchanged from input to output)
*/

void csm_yproj(Int_t nx, Int_t ny, Double_t *xydist, Double_t *ydist)
{
  Int_t i,j;
  Double_t total;

  for (i=0;i<ny;i++) ydist[i] = 0;
  total = 0;

  for (i=0;i<ny;i++)
    {
      for (j=0;j<nx;j++)
	{
          ydist[i] += xydist[j+nx*i];
	}
      total += ydist[i]; 
    }

  if (total>0)
    {
      for (i=0;i<ny;i++) ydist[i] /= total;
    }
}

/*
  Recoded d_getycont -- original by Alex Read, recoded by Tom Junk
  February 2007

  <ydist1> and <ydist2> and <ydist3>
  are the projections on the y-axis of three
  scatterplots which are going to be interpolated. <ydistn> is
  the interpolated 1d distribution which represent the projection
  of the interpolated scatterplot on the y-axis. This routine determines
  which bins of <ydist1> and <ydist2> and <ydist3>
  contribute and by what amount to
  each bin of <ydistn>. This information is used in csm_pvmc2d to 
  determine the input distributions in the x-direction of each
  y-bin: these are then interpolated and accumulated in the interpolated
  2d distribution.

  Inputs : ny,ydist1,ydist2,ydist3,ydistn
  Outputs: alpha1,alpha2,alpha3

  alpha1[iyc+ny*iy] encodes the contribution of bin iyc in ydist1
  to to bin iy in ydistn

*/

void csm_ycont(Int_t ny, Double_t *ydist1, Double_t *ydist2,
               Double_t *ydist3, Double_t *ydistn,
               Double_t *alpha1, Double_t *alpha2, Double_t *alpha3)
{
  Double_t y[ny+1];
  Double_t yn[ny+1];
  Int_t i;

  // Make arrays to describe the straight-line approximations
  // to the four cumulative distributions, y1,y2,y3,yn 
  // Make sure to start out with a point at 0

  yn[0] = 0;
  for (i=0;i<ny;i++) yn[i+1] = ydistn[i];
  csm_acnvec2(yn,ny+1);

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist1[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha1);
#ifdef DEBUGPVMC2D
  for (i=0;i<ny+1;i++)
    {
      cout << "getting alpha1: " << i << " " << y[i] << " " << yn[i] << endl;
    }
  Int_t j; 
  for (i=0;i<ny;i++)
    {
      for (j=0;j<ny;j++)
	{
	  cout << i << " " << j << " " << alpha1[i+ny*j] << endl;
	}
    }
  
#endif

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist2[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha2);

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist3[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha3);
}

void csm_ycontaux(Int_t ny, Double_t *y, Double_t *yn,
                  Double_t *alpha)
{
  Int_t ny2;
  Int_t i,j;

  ny2 = ny*ny;

  // clear out the alpha array

  for (i=0;i<ny2;i++) alpha[i] = 0;

  // loop over bins and see what fraction each contributes

  for (i=0;i<ny;i++) // interpolated histogram
    {
      for (j=0;j<ny;j++) // contributing histogram bin
	{
          if (y[j+1]-y[j]>0)
	    {
	      // first case -- contributing bin entirely contained
	      // within the interpolated output bin.
	      if (y[j]>=yn[i] && y[j]<yn[i+1] && 
		  y[j+1]>=yn[i] && y[j+1]<yn[i+1])
		{
		  alpha[j+ny*i] = 1;
		}
	      // second case -- interpolated output bin is entirely
	      // contained within the contributing bin
	      else if (y[j]<yn[i] && y[j+1] >= yn[i+1])
		{
		  alpha[j+ny*i] = (yn[i+1]-yn[i])/(y[j+1]-y[j]);
		}
	      // third case -- contributing bin straddles the
	      // left edge of the interpolated output bin but ends inside
	      // the bin
	      else if (y[j]<yn[i] && y[j+1]>=yn[i] && y[j+1]<yn[i+1])
		{
		  alpha[j+ny*i] = (y[j+1]-yn[i])/(y[j+1]-y[j]);
		}
	      // fourth case -- contributing bin straddles the
	      // right edge of the interpolated output bin but starts inside
	      // the output bin
	      else if (y[j]>=yn[i] && y[j]<yn[i+1] && y[j+1]>=yn[i+1])
		{
		  alpha[j+ny*i] = (yn[i+1]-y[j])/(y[j+1]-y[j]);
		}
	      // non-overlapping case -- do nothing.
	      // save some time if we're beyond the edge
	      if (y[j]>yn[i+1]) break; 
	    }
	}
    } 
}

// Integrate and normalize a vector -- pretty much the same as
// csm_acnvec

void csm_acnvec2(Double_t *vec, Int_t n)
{
  Int_t i;
  Double_t tot;
  for (i=1;i<n;i++)
    {
      vec[i] += vec[i-1];
    }
  tot = vec[n-1];
  if (tot > 0)
    {
      for (i=0;i<n;i++) vec[i] /= tot;
    }
}


/*-------------------------------------------------------------------------*/
/*    Interface to Joel's genlimit.c program for a Bayesian calcualtion    */
/*    of an upper limit.  Also run pseudoexperiments if need be to compute */
/*    expected limits.                                                     */
/*    Bayesian limit calculation uses test_hypothesis_pe to compute the    */
/*    "Bayesian ensemble" because it should have signal and background     */
/*    components marked, and because it should have all systematic         */
/*    errors included.  It uses null_hypothesis_pe in order to generate    */
/*    pseudoexperiments to compute expected limits however.                */
/*-------------------------------------------------------------------------*/
/* arguments:  beta: credibility level:L  0.95 for 95% CL limits
   sflimit:  the observed limit
   unc:      MC statistical unc. on observed limit
   npx:      Number of pseudoexperiments to run to compute expected limits
   sm2:      -2 sigma expected limit      *put in null pointers for all
   sm1:      -1 sigma expected limit      *five of these to skip the
   smed:     median expected limit        *calculation and speed it up.
   sp1:      +1 sigma expected limit
   sp2:      +2 sigma expected limit
*/

void mclimit_csm::bayes_heinrich_withexpect(Double_t beta,
                                            Double_t* sflimit,
                                            Double_t* unc,
					    Int_t npx,
                                            Double_t* sm2,
                                            Double_t* sm1,
                                            Double_t* smed,
                                            Double_t* sp1,
                                            Double_t* sp2)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1Type* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;
  double *xgl;
  double *lwgl;

  vector<Double_t> bpploc;
  vector<Double_t> bploc;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      //cout << ibin << "    " << nobs[ibin] << endl;
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = gRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = gRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
                                                  // if r is already zero this won't get stuck in a loop
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  //cout << iens << " " << ens[iens].b << " " << ens[iens].e << endl;
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.  We will enlarge these arrays
  //later if the need arises.

  nglmax = nobstot;
  if (nglmax<10000) {nglmax = 10000;}
  xgl = new double[nglmax];
  lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;

  *sflimit = 0;
  if (bayesintegralmethod == CSM_BAYESINTEGRAL_JOEL)
    {
      *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);
    }
  else
    {
      setdlcsn(nbinstot,nens,nobs,ens);
    }

  // make a vector of the posterior likelihood function
  if (bayes_interval_step > 0)
    {
      if ( (bayes_interval_end-bayes_interval_begin)>0 )
	{
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  Double_t b,p;
	  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
	    {
	      p = cspdf(b,1.0,nbinstot,nens,nobs,ens,prior);
	      bayes_posterior.push_back(p);
	      bploc.push_back(p);
	      bayes_posterior_points.push_back(b);
	      bpploc.push_back(b);
	    }
	  int jsiz = bayes_posterior.size();
	  int ic;
	  Double_t bptot = 0;
	  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	  if (bptot>0)
	    {
    	      Double_t scaleb = 1.0/(bptot*bayes_interval_step);
	      bhnorm = scaleb;
	      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
	    }
	  if (bayesintegralmethod == CSM_BAYESINTEGRAL_QUICK)
	    {
	      *sflimit = quickbint(beta);
	    }
	}
    }


  // compute expected limits

  cslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      pdname = new char[strlen(test_hypothesis_pe->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,null_hypothesis_pe->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = copy_TH1<TH1> (*null_hypothesis_pe->chanmodel[i]->histotemplate[0], pdname).release();
      delete[] pdname;
    }
  for (ipx=0;ipx<npx;ipx++)
    {
      null_hypothesis_pe->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
		  if (nbinsy==1)
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
		  else
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		  nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
		  ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  if (nglmax < nobstotlist[nobsindex[0]]/2 + 1)
    {
      nglmax = nobstotlist[nobsindex[0]]/2 + 1; 
      delete[] xgl;
      delete[] lwgl;
      xgl = new double[nglmax];
      lwgl = new double[nglmax];
    }

  ngl = 0;
  for (ipx=0;ipx<npx;ipx++)
    {
      Double_t p=0;
      if (bayesintegralmethod == CSM_BAYESINTEGRAL_JOEL)
	{
          if (ipx>0)
	    { if (nobstotlist[nobsindex[ipx]] != nobstotlist[nobsindex[ipx-1]])
	      { 
	        ngl = 0;
	      }
	    }
	  p = cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,unc);
	}
      else
	{

	  setdlcsn(nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens);

	  // make a vector of the posterior likelihood function
	  if (bayes_interval_step > 0)
	    {
	      if ( (bayes_interval_end-bayes_interval_begin)>0 )
		{
		  bayes_posterior.clear();
		  bayes_posterior_points.clear();
		  Double_t b,p1;
		  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
		    {
		      p1 = cspdf(b,1.0,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,prior);
		      bayes_posterior.push_back(p1);
		      bayes_posterior_points.push_back(b);
		    }
		  int jsiz = bayes_posterior.size();
		  int ic;
		  Double_t bptot = 0;
		  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
		  if (bptot>0)
		    {
		      Double_t scaleb = 1.0/(bptot*bayes_interval_step);
	              bhnorm = scaleb;
		      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
		    }
		}
	      p = quickbint(beta);
	    }

	}


      if (pxprintflag)
	{
	  cout << "bayespx: " << p << endl;
	}
      if (bayes_pseudoexperiment_limits != 0) bayes_pseudoexperiment_limits->Fill(p);
      cslist.push_back(p);
    }
  std::sort(cslist.begin(),cslist.end());
  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm2 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm1 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *smed = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp1 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp2 = cslist[i];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;

  // copy the observed posterior curve into the externally visible vectors

  bayes_posterior.clear();
  bayes_posterior_points.clear();
  for (int k=0;k< (int) bpploc.size();k++)
    {
      bayes_posterior.push_back(bploc[k]);
      bayes_posterior_points.push_back(bpploc[k]);
    }

}

/*-------------------------------------------------------------------------*/
//  Same thing as above, but only compute observed limit (much quicker)
/*-------------------------------------------------------------------------*/

void mclimit_csm::bayes_heinrich(Double_t beta,
                                 Double_t* sflimit,
                                 Double_t* unc)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1Type* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  Int_t nobstot;
  int nglmax;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = gRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = gRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.

  nglmax = nobstot;
  if (nglmax<10000) {nglmax = 10000;}
  double* xgl = new double[nglmax];
  double* lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;

  *sflimit = 0;
  if (bayesintegralmethod == CSM_BAYESINTEGRAL_JOEL)
    {
      *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);
    }
  else
    {
      setdlcsn(nbinstot,nens,nobs,ens);
    }

  //cout << "ngl: " << ngl << " " << nobstot << endl;

  // make a vector of the posterior likelihood function
  if (bayes_interval_step > 0)
    {
      if ( (bayes_interval_end-bayes_interval_begin)>0 )
	{
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  Double_t b,p;
	  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
	    {
	      p = cspdf(b,1.0,nbinstot,nens,nobs,ens,prior);
	      bayes_posterior.push_back(p);
	      bayes_posterior_points.push_back(b);
	    }
	  int jsiz = bayes_posterior.size();
	  int ic;
	  Double_t bptot = 0;
	  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	  if (bptot>0)
	    {
    	      Double_t scaleb = 1.0/(bptot*bayes_interval_step);
	      bhnorm = scaleb;
	      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
	    }
	  if (bayesintegralmethod == CSM_BAYESINTEGRAL_QUICK)
	    {
	      *sflimit = quickbint(beta);
	    }
	}
    }

  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
}

// expected limits with the Markov Chain MC limit calculator
// run pseudoexperiments based on null_hypothesis_pe and call
// bayeslimit_mcmc1 for each one

void mclimit_csm::bayeslimit_mcmc1_expect(double beta, PRIOR prior, int npx, 
					  double *sm2, double *sm1, double *smed, 
					  double *sp1, double *sp2)
{
  vector<double> clist;
  int nchans = null_hypothesis_pe->channame.size();
  TH1** pdarray = new TH1*[nchans];     // data and pseudodata stay as TH1's
  TH1** dhsavearray = new TH1*[nchans];
  char *pdname;

  for (int i=0;i<nchans; i++)
    {
      pdname = new char[strlen(null_hypothesis_pe->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,null_hypothesis_pe->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) datahist[i]->Clone(pdname);
      delete[] pdname;
      dhsavearray[i] = datahist[i];
      datahist[i] = pdarray[i];
    }

  // caclculate a limit for each pseudoexperiment

  for (int ipx=0; ipx<npx; ipx++)
    {
      null_hypothesis_pe->single_pseudoexperiment(pdarray);
      double slimit = bayeslimit_mcmc1(beta, prior);
      if (pxprintflag) cout << "bayespx: " << slimit << endl; 
      clist.push_back(slimit);
    }

  std::sort(clist.begin(),clist.end());
  int i =  (int) nearbyint(npx*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm2 = clist[i];

  i =  (int) nearbyint(npx*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm1 = clist[i];

  i =  (int) nearbyint(npx*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *smed = clist[i];

  i =  (int) nearbyint(npx*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp1 = clist[i];

  i =  (int) nearbyint(npx*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp2 = clist[i];

  for (int i=0;i<nchans; i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  for (int i=0;i<nchans; i++)
    {
      datahist[i] = dhsavearray[i];
    }
  delete[] dhsavearray;

}

//-------------------------------------------------------------------------
//  Adaptive Markov Chain limits
//  Uses testhyp_pe as the model to test
//-------------------------------------------------------------------------
// all arguments are optional
// Beta is the credibility level desired for the limits -- the default is 0.95.

double mclimit_csm::bayeslimit_mcmc1(double beta, PRIOR prior, TString histoutfile)
{
  csm_channel_model* cm;

  // count up all the nuisance parameters -- add in the bin-by-bin errors (Poisson flag = 1 or 2)
  // take sqrt(contents) for poisson flag = 1 and the error supplied for poisson flag = 2

  vector<char*> npnames;
  vector<double> nplb;
  vector<double> nphb;

  Bool_t addStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); // cloned histograms go in memory, and aren't deleted when files are closed.
                             // be sure to restore this state when we're out of the routine.

  test_hypothesis_pe->list_nparams(&npnames, &nplb, &nphb);   // this clears and fills the vectors of name pointers and bounds

  int nbinstot = 0;
  int nperrtot = 0;

  // count up channels, bins, and independent Poisson uncertainties

  int nchans = test_hypothesis_pe->channame.size();
  for (int i=0;i<nchans;i++)
    {
      cm = test_hypothesis_pe->chanmodel[i];
      int nbtmp = cm->histotemplate[0]->GetNbinsX() * cm->histotemplate[0]->GetNbinsY();
      nbinstot += nbtmp;
      int ntemplates = cm->histotemplate.size();
      for (int itpl = 0; itpl < ntemplates; itpl++)
	{
          if (cm->poissflag[itpl] != 0) nperrtot += nbtmp;
	}
    }

  // allocate storage for nuisance parameter values and set initial values to zero

  // rate and shape errors
  int nnptot = npnames.size();
  double npvalues[nnptot];
  double proposed_npvalues[nnptot];
  double npresolutions[nnptot];
  double npf1[nnptot];

  // bin by bin errors (Poisson)
  double pnpvalues[nperrtot];
  double proposed_pnpvalues[nperrtot];
  double pnpresolutions[nperrtot];
  double pnpf1[nperrtot];

  double proposed_ssf=1.0;
  for (int i=0;i<nnptot;i++) npvalues[i] = 0;
  for (int i=0;i<nperrtot;i++) pnpvalues[i] = 0;
  double ssf = 1.0;
  double slog = 0.0;
  int phistflag = 0;
  TFile *histfile = 0;
  vector<TH1*> nphist;
  vector<double> ssflist;
  TH1D *totcount = 0;
  const int nburnin = 500; // tuning parameter -- do not like...
  double stotdim = sqrt(nnptot + nperrtot + 1.0)/2.0;
  double stotdims = stotdim/2.0;       // see code below -- this gets overwritten.

  const int accept_interval = 5000;    // every 5000 samples re-evalulate the acceptance rate.
  const double accept_low = 0.3;       // desire an acceptance ratio between low and high values.
  const double accept_high = 0.5;
  const double dimscale = 1.2;         // scale the step by this up or down if we don't have the right acceptance rate.

  // set up posterior histograms if requested
  if (histoutfile != "")
    {
      phistflag = 1;
      for (int i=0;i<nnptot;i++)
	{
	  TH1 *h = new TH1F(npnames[i],npnames[i],100,-5,5);
	  nphist.push_back(h);
	}
      for (int i=0;i<nperrtot;i++)
	{
	  char pnpntmp[15];
	  sprintf(pnpntmp,"Poiss_%d",i);
	  TH1 *h = new TH1F(pnpntmp,pnpntmp,100,-5,5);
	  nphist.push_back(h);
	}
      totcount = new TH1D("totcount","totcount",3,-0.5,2.5);
    }

  double llf = 0;

  // first-time through initialization
  // first time through use default nuisance parameter values, and a signal scale factor of 1.0

  test_hypothesis_pe->undo_nuisance_response();  // since bayesmcmcllf no longer calls nuisance_response
  // for the non-Poisson nuisance parameters.
  llf = bayesmcmcllf(test_hypothesis_pe,nnptot,npvalues,
		     nperrtot,pnpvalues,ssf,prior);
  slog = llf;

  // at the default values of the nuisance parameters, explore along the signal scale factor axis
  // so we can set the signal scale appropriately.
  // do a scan over six orders of magnitude, in 80 steps.  

  const double sfbase = 1.1885;   // the 40th root of 1000

  double lmx = 0;
  int imx = -40;
  double llflist[81];
  for (int isf=-40;isf<41;isf++)
    {
      double esf = pow(sfbase,(double) isf);
      double llftest = bayesmcmcllf(test_hypothesis_pe,nnptot,npvalues,
				    nperrtot,pnpvalues,esf,prior);
      if (isf == -40 || llftest > lmx) 
	{
	  imx = isf;
	  lmx = llftest;
	}
      llflist[isf+40] = llftest;
    }

  // start out the signal scale factor at the peak value
  ssf = pow(sfbase, (double) imx);
  // figure out how far we have to go in order to drop by a factor of two

  const double g2 = 0.693147181;  // log(2)
  int isfl = imx;
  int isfh = imx;
  for (int i=1;i<81;i++)
    {
      if (imx-i >= -40)
	{
	  if ( (lmx - llflist[imx-i+40]) < g2 ) isfl = imx-i;
	}
      if (imx+i <= 40)
	{
	  if ( (lmx - llflist[imx+i+40]) < g2 ) isfh = imx+i;
	}
    }
  double siz = pow(sfbase, (double) isfh) - pow(sfbase, (double) isfl);
  if (siz>0.0)
    {
      stotdims = 1.0/siz;
    }
  else
    {
      stotdims = stotdim/2.0;
    }

  // explore the space of nuisance parameters along the axes and see which ones are constrained 

  for (int i=0;i<nnptot;i++)
    {
      double llft[3],slf[3];
      int k=0;
      for (double v=-4.0; v<4.1; v += 0.5)  // vary over +- 4 sigma in steps of one sigma and find the
	// three highest points
	{
	  if (v>nplb[i] && v<nphb[i])
	    {
	      npvalues[i] = v;
	      test_hypothesis_pe->nuisance_response(nnptot,&(npnames[0]),npvalues);
	      double llft1 = bayesmcmcllf(test_hypothesis_pe,nnptot,npvalues,
					  nperrtot,pnpvalues,ssf,prior) - slog;
	      npvalues[i] = 0;  // udno this so we explore along the axes
	      if (k<3)
		{
		  llft[k] = llft1;
		  slf[k] = v;
		  k++;
		}
	      else
		{
		  int imin=0;  // find the smallest and drop it off the list if we found a bigger one
		  for (int k1=0;k1<3;k1++)
		    {
		      if (llft[k1]<llft[imin]) imin=k1;
		    }
		  if (llft1>llft[imin]) 
		    {
		      llft[imin] = llft1;
		      slf[imin] = v;
		    }
		}
	    }
	}
      // draw a parabola through the top three points -- find the maximum, require it to be in the
      // parameter rage, and find the width, called the resolution.

      //cout << "llft: " << llft[0] << " " << llft[1] << " " << llft[2] << endl;
      //cout << "slf : " << slf[0] << " " << slf[1] << " " << slf[2] << endl;

      double x1 = slf[0];
      double x2 = slf[1];
      double x3 = slf[2];
      double xi[3][3];
      double det = x1*x1*x2 - x1*x1*x3 -x2*x2*x1 + x3*x3*x1 + x2*x2*x3 - x2*x3*x3;

      if (det != 0)
	{
	  xi[0][0] = (x2-x3)/det;
	  xi[0][1] = (x3-x1)/det;
	  xi[0][2] = (x1-x2)/det;
	  xi[1][0] = (x3*x3-x2*x2)/det;
	  xi[1][1] = (x1*x1-x3*x3)/det;
	  xi[1][2] = (x2*x2-x1*x1)/det;
	  xi[2][0] = (x2*x2*x3-x2*x3*x3)/det;
	  xi[2][1] = (x1*x3*x3-x1*x1*x3)/det;
	  xi[2][2] = (x1*x1*x2-x1*x2*x2)/det;
	  double a=0;
	  double b=0;
	  double c=0;
	  for (int j=0;j<3;j++)
	    {
	      a += xi[0][j]*llft[j];
	      b += xi[1][j]*llft[j];
	      c += xi[2][j]*llft[j];
	    }
          //cout << "abc: " << a << " " << b << " " << c << endl;
	  //cout << a*x1*x1 + b*x1 + c - llft[0] << endl;
	  //cout << a*x2*x2 + b*x2 + c - llft[1] << endl;
	  //cout << a*x3*x3 + b*x3 + c - llft[2] << endl;

	  if (a == 0)  // flat parabola
	    {
	      npf1[i] = 0;
	      npresolutions[i] = 1.0;
	    }
	  else
	    {
	      npf1[i] = -b/(2.0*a);
	      npresolutions[i] = -0.5/a;
	    }
	}
      else // case of zero determinant
	{
	  npf1[i] = 0;
	  npresolutions[i] = 1.0;
	}
      //cout << "resolution: " << npresolutions[i] << " " << i << endl;
      npf1[i] = min(nphb[i],max(nphb[i],npf1[i]));  // bound it in the min, max range
      npresolutions[i] = max(0.1,min(1.0,npresolutions[i]));  // bound the resolutions too
    }

  // for now let's set the resolutions of all bin by bin parameters to 1.0
  for (int i=0;i<nperrtot;i++)
    {
      pnpf1[i] = 0.0;
      pnpresolutions[i] = 1.0;
    }

  int iadcount=0;
  int iaccept_count=0;

  for (int isample=0; isample<nmc_req; isample++)
    {
      // the Metropolis-Hastings proposal function -- Gaussians for the nuisance parameters,
      // uniform for the signal scale factor.

      int itest=0;
      if (totcount) totcount->Fill(0);
      for (int i=0;i<nnptot;i++)
	{
	  do 
	    {
	      proposed_npvalues[i] = npvalues[i] + gRandom->Gaus(0,1)*npresolutions[i]/stotdim;
	      //cout << "proposed npvalue: " << proposed_npvalues[i] << endl;
	    }
	  while (proposed_npvalues[i] > nphb[i] || proposed_npvalues[i] < nplb[i]); // never go out of bounds
	}

      test_hypothesis_pe->nuisance_response(nnptot,&(npnames[0]),proposed_npvalues);

      // all the work below is to determine if a bin's content would go negative if we throw
      // the Poisson nuisance parameter -- rethrow it if it does go negative.  This depends on
      // the systematically varied bin contents from varying the other nuisance parameters so we
      // have to apply nuisance_response here in the proposal function
      // don't apply sft_varied as they scale both the contents and the errors.  Do it in the likelihood function though
      // also compute the llf asymmetry factor for the bin by bin errors as they are trucated to keep the
      // predictions from going negative

      int inp = 0;
      double proposed_bbcorr = 0;
      for (int i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  int nbinsx = cm->histotemplate[0]->GetNbinsX();
	  int nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (int j=0;j<nbinsx;j++)
	    {
	      for (int k=0;k<nbinsy;k++)
		{
		  int ntemplates = cm->histotemplate.size();
		  for (int itpl = 0; itpl < ntemplates; itpl++)
		    { 
                      if (cm->poissflag[itpl] != 0) 
			{
			  double rtmp = 0;
			  double dr=0;
			  if (nbinsy == 1) 
			    {
			      rtmp = cm->histotemplate_varied[itpl]->GetBinContent(j+1);
			    }
			  else
			    {
			      rtmp = cm->histotemplate_varied[itpl]->GetBinContent(j+1,k+1);
			    }
			  if (cm->poissflag[itpl] == CSM_POISSON_BINERR) 
			    {
			      if (nbinsy == 1) 
				{
				  dr = cm->histotemplate_varied[itpl]->GetBinContent(j+1);
				}
			      else
				{
				  dr = cm->histotemplate_varied[itpl]->GetBinContent(j+1,k+1);
				}
			      if (dr>0) dr = sqrt(dr);
			    }
			  if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR) 
			    {
			      if (nbinsy == 1) 
				{
				  dr = cm->histotemplate_varied[itpl]->GetBinError(j+1);
				}
			      else
				{
				  dr = cm->histotemplate_varied[itpl]->GetBinError(j+1,k+1);
				}
			    }
			  double rte=rtmp;
			  do 
			    {
			      proposed_pnpvalues[inp] = pnpvalues[inp] + gRandom->Gaus(0,1)/stotdim;
			      rte = rtmp + dr*proposed_pnpvalues[inp];
			    }
			  while (rte<0);
			  if (dr>0)
			    {
			      proposed_bbcorr -= log( 0.5*(1.0 + erf(rtmp/(dr*sqrt(2.0)))));  
			    }

			  inp++;
			}
		    }
		}
	    }
	}
      do
	{
	  double ur=0;
	  do 
	    {
	      ur = gRandom->Uniform(1.0);
	    }
	  while (ur==1.0); // rethrow if we get exactly 1.0 to make it completely symmetric
	  proposed_ssf = ssf + (ur - 0.5)/stotdims;
	}
      while (proposed_ssf < 0);

      double proposed_llf = bayesmcmcllf(test_hypothesis_pe,nnptot,proposed_npvalues,
					 nperrtot,proposed_pnpvalues,proposed_ssf,prior);

      // since the parameter space is truncated, put in a factor for the asymmetry in the
      // step from here to another place and back again.  

      for (int i=0;i<nnptot;i++)
	{
	  proposed_llf -= log( 0.5*( erf((nphb[i]-proposed_npvalues[i])/(sqrt(2.0)/(stotdim/npresolutions[i]))) - 
	  		             erf((nplb[i]-proposed_npvalues[i])/(sqrt(2.0)/(stotdim/npresolutions[i]))) ) );
	}
      if (proposed_ssf < 0.5/stotdims)
	{
	  proposed_llf -= log(0.5 + proposed_ssf*stotdims);
	}
      proposed_llf += proposed_bbcorr;

      double dllf = proposed_llf - llf;
      double pratio = 0;
      if (dllf>-150 && dllf < 150)
	{
	  pratio = exp(dllf);
	}
      double ur = gRandom->Uniform(1.0);
      if (ur < pratio) 
	{
	  itest = 1;
	  ssf = proposed_ssf;
	  for (int i=0;i<nnptot;i++) npvalues[i] = proposed_npvalues[i];
	  for (int i=0;i<nperrtot;i++) pnpvalues[i] = proposed_pnpvalues[i];
	  llf = proposed_llf;
	}
      if (isample>=nburnin)
	{
          if (totcount) 
	    { 
	      totcount->Fill(1);
	      if (itest==1) totcount->Fill(2);
	    }

	  ssflist.push_back(ssf);
	  if (phistflag)
	    {
	      int k=0;
              for (int i=0;i<nnptot;i++)
		{
		  nphist[k]->Fill(npvalues[i]);
		  k++;
		}
              for (int i=0;i<nperrtot;i++)
		{
		  nphist[k]->Fill(pnpvalues[i]);
	          k++;
		}
	    }
	  // adapt the acceptance ratio

	  iadcount++;
	  if (itest==1) iaccept_count++;
	  if (iadcount >= accept_interval)
	    {
	      double acceptrate = ((double) iaccept_count)/((double) iadcount);
	      if (acceptrate>accept_high) 
		{
		  stotdim /= dimscale;
		  stotdims /= dimscale;
		}
	      if (acceptrate<accept_low)
		{
		  stotdim *= dimscale;
		  stotdims *= dimscale;
		}
	      iadcount = 0;
	      iaccept_count = 0;
	    }

	}
    }

  std::sort(ssflist.begin(),ssflist.end());

  if (phistflag) 
    {
      histfile = new TFile(histoutfile,"RECREATE");
      if (histfile==0)
	{
	  cout << "Markov Chain output histogram file open failed " << histoutfile << endl;
	  exit(0);
	}
      if (histfile->IsZombie())
	{ 
	  cout << "Markov Chain output histogram file open failed " << histoutfile << endl;
	  exit(0);
	}

      histfile->cd();
      TH1F ssfhist("mcmcssf","Signal Scale Factor",100,0,ssflist.back());
      int nssflist = ssflist.size();
      for (int i=0; i<nssflist; i++) ssfhist.Fill(ssflist[i]);

      int nnphist=nphist.size();
      for (int i=0;i<nnphist; i++)
	{
	  nphist[i]->Write();
	}
      ssfhist.Write();
      totcount->Write();
      histfile->Close();
    }

  // free up memory used by allocated histograms
  int nnphist=nphist.size();
  for (int i=0; i<nnphist; i++)
    {
      delete nphist[i];
    }

  TH1::AddDirectory(addStatus);

  // find the beta quantile of ssflist (already sorted)

  int idx = (int) nearbyint(beta*ssflist.size());
  idx = min(ssflist.size()-1,max(0,idx));
  return(ssflist[idx]);

}

// logarithm of the product of probabilities to observe the data and the Gaussian constraints
// on the nuisance parameters -- still need the nuisance parameter values to compute the Gaussian constraints
// assumes nuisance_response has already been called for *model

double mclimit_csm::bayesmcmcllf(csm_model *model, int nnptot, double *npvalues, 
                                 int nperrtot, double *pnpvalues, double ssf, PRIOR prior)
{
  csm_channel_model *cm;

  double llf = 0;
  double stot = 0;

  // put in constraint for nuisance parameters.  Treat bin by bin errors as Gaussian

  for (int i=0;i<nnptot;i++)
    {
      llf -= npvalues[i]*npvalues[i]/2.0;
    }
  for (int i=0;i<nperrtot;i++)
    {
      llf -= pnpvalues[i]*pnpvalues[i]/2.0;
    }

  int nchans = model->channame.size();
  int inp = 0;
  for (int i=0;i<nchans;i++)
    {
      cm = model->chanmodel[i];
      int nbinsx = cm->histotemplate[0]->GetNbinsX();
      int nbinsy = cm->histotemplate[0]->GetNbinsY();
      for (int j=0;j<nbinsx;j++)
	{
	  for (int k=0;k<nbinsy;k++)
	    {
              double r = 0;
	      int ntemplates = cm->histotemplate.size();
              for (int itpl = 0; itpl < ntemplates; itpl++)
	        { 
		  double rtmp = 0;
	          double dr=0;
		  if (nbinsy==1) 
		    {
		      rtmp = cm->histotemplate_varied[itpl]->GetBinContent(j+1)*cm->sft_varied[itpl];
		    }
		  else
		    {
		      rtmp = cm->histotemplate_varied[itpl]->GetBinContent(j+1,k+1)*cm->sft_varied[itpl];
		    }
                  if (cm->poissflag[itpl] == CSM_POISSON_BINERR) 
	            {
		      if (nbinsy==1) 
		        {
			  dr = cm->histotemplate_varied[itpl]->GetBinContent(j+1);
		        }
		      else
		        {
			  dr = cm->histotemplate_varied[itpl]->GetBinContent(j+1,k+1);
		        }
		      if (dr>0) dr = sqrt(dr)*cm->sft_varied[itpl];
		    }
                  if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR) 
	            {
		      if (nbinsy==1) 
		        {
			  dr = cm->histotemplate_varied[itpl]->GetBinError(j+1)*cm->sft_varied[itpl];
		        }
		      else
		        {
			  dr = cm->histotemplate_varied[itpl]->GetBinError(j+1,k+1)*cm->sft_varied[itpl];
		        }
		    }
		  double rte = rtmp;
		  if (cm->poissflag[itpl] != 0) 
                    { 
		      rte += dr*pnpvalues[inp];
		      inp ++;
		    }
		  if (rte < 0) 
		    { 
		      rte = 0; // truncate bin by bin errors by template at zero
		      cout << "calling routine of bayesmcmcllf gives a poisson parameter resulting in a negative prediction" << endl;
		      exit(0); // should never happen.
		    }
		  if (cm->scaleflag[itpl] != 0) 
                    { 
		      stot += rte;  // accumulate the unscaled signal and use that for the correlated prior
		      rte*=ssf;
		    }
		  r += rte;
	        } // end loop on template
	      int nobs = 0;
	      if (nbinsy==1)
		{ nobs = (int) nearbyint(datahist[i]->GetBinContent(j+1));}
	      else
		{ nobs = (int) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }

	      if (r<1E-10) r=1e-10; // protect against zero expectation

	      llf += nobs*log(r) - r - lgamma(nobs+1);

	    } // end loop on ybins
	} // end loop on xbins
    } // end loop on channels
  assert(prior==flat || prior==corr);
  if (prior == corr)
    {
      if (stot > 0.0) 
        { llf += log(stot); }
      else
	{ llf = -1.0E40; } // to exponentiate to zero.
    }
  return(llf);
}


// Make channel error band plots.   This routine will create the histograms and fill in pointers to them.

void mclimit_csm::systsumplot(char *channame, int nsamples, TH1 **ehm2, TH1 **ehm1, TH1 **ehmed, TH1 **ehp1, TH1 **ehp2)
{
  int isample,i,j;
  int ichan=0;
  int nchans = test_hypothesis_pe->channame.size();
  int nbinsx = 0;
  int nbinsy = 0;
  int nbinstot = 0;
  csm_channel_model *cm;
  TH1Type* ht;
  double r;
  // Nils likes std_auto_ptrs but we'll have to clone them at the end
  // in order for the histograms to persist after the function ends.
  std::auto_ptr<TH1> ahm2;

  for (i=0;i<nchans;i++)
    {
      if (strcmp(channame,test_hypothesis_pe->channame[i])==0) 
        { 
	  ichan=i;
          cm = test_hypothesis_pe->chanmodel[ichan];
          nbinsx = cm->histotemplate[0]->GetNbinsX();
          nbinsy = cm->histotemplate[0]->GetNbinsY();
          nbinstot = nbinsx*nbinsy;
          ahm2  = copy_TH1<TH1> (*cm->histotemplate[0]);
	}
    }

  // persistent versions of these histograms.

  TH1 *hm2 = (TH1*) ahm2->Clone("ehm2");
  TH1 *hm1 = (TH1*) ahm2->Clone("ehm1");
  TH1 *hmed = (TH1*) ahm2->Clone("ehmed");
  TH1 *hp1 = (TH1*) ahm2->Clone("ehp1");
  TH1 *hp2 = (TH1*) ahm2->Clone("ehp2");

  *ehm2 = hm2;
  *ehm1 = hm1;
  *ehmed = hmed;
  *ehp1 = hp1;
  *ehp2 = hp2;

  hm2->Reset();
  hm1->Reset();
  hmed->Reset();
  hp1->Reset();
  hp2->Reset();

  double bv[nbinstot*nsamples];
  int idx = 0;

  // fill up an array of fluctuated total yields for each bin for each systematic variation

  for (isample=0;isample<nsamples;isample++)
    {
      test_hypothesis_pe->varysyst();
      cm = test_hypothesis_pe->chanmodel[ichan];
      int ntemplates = cm->histotemplate.size();
      for (i=0;i<nbinsx;i++)
	{
	  for (j=0;j<nbinsy;j++)
	    {
	      bv[idx] = 0;
	      for (int itpl=0;itpl<ntemplates;itpl++)
		{
		  ht = cm->histotemplate_varied[itpl];
		  if (nbinsy==1) 
		    { r = ht->GetBinContent(i+1); }
		  else
		    { r = ht->GetBinContent(i+1,j+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = gRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(i+1);}
			  else
			    { histerr = ht->GetBinError(i+1,j+1);}
			  do
			    { edraw = gRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      bv[idx] += r;
		}
	      idx++;
	    }
	}
    }

  // find the quantiles of each bin's distributions

  int im2s = (int) nearbyint(nsamples*MCLIMIT_CSM_MCLM2S);
  int im1s = (int) nearbyint(nsamples*MCLIMIT_CSM_MCLM1S);
  int imed = (int) nearbyint(nsamples*MCLIMIT_CSM_MCLMED);
  int ip1s = (int) nearbyint(nsamples*MCLIMIT_CSM_MCLP1S);
  int ip2s = (int) nearbyint(nsamples*MCLIMIT_CSM_MCLP2S);

  for (i=0;i<nbinsx;i++)
    {
      for (j=0;j<nbinsy;j++)
	{
	  vector<double> bd;
	  for (isample=0;isample<nsamples;isample++)
	    {
	      bd.push_back(bv[isample*nbinstot + i*nbinsy + j]);
	    }
	  std::sort(bd.begin(),bd.end());
	  if (nbinsy == 1)
	    {
  	      hm2->SetBinContent(i+1,bd[im2s]);
  	      hm1->SetBinContent(i+1,bd[im1s]);
  	      hmed->SetBinContent(i+1,bd[imed]);
  	      hp1->SetBinContent(i+1,bd[ip1s]);
  	      hp2->SetBinContent(i+1,bd[ip2s]);
	    }
	  else
	    {
  	      hm2->SetBinContent(i+1,j+1,bd[im2s]);
  	      hm1->SetBinContent(i+1,j+1,bd[im1s]);
  	      hmed->SetBinContent(i+1,j+1,bd[imed]);
  	      hp1->SetBinContent(i+1,j+1,bd[ip1s]);
  	      hp2->SetBinContent(i+1,j+1,bd[ip2s]);
	    }
	}
    }

}

// Fit a cross section using Joel Heinrich's marginalized posterior
// use testhyp_pe for this fit (has all nuisance parameters defined)

// reverse order of summation -- do only one systematic sample at a time; saves memory

void mclimit_csm::bh_xsfit(Double_t *xsfit, Double_t *downerr, Double_t *uperr)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1Type* ht;
  Double_t r;
  const PRIOR prior=flat;
  Int_t nobstot;


  *xsfit = 0.0;
  *downerr = 0.0;
  *uperr = 0.0;

  if (bayes_interval_step < 0 || (bayes_interval_end-bayes_interval_begin)<0 )
    {
      cout << "bh_xsfit: invalid integration region or step." << endl;
      cout << "bh_xsfit:  begin: " << bayes_interval_begin << " end: " << bayes_interval_end 
           << " step: " << bayes_interval_step << endl;
      cout << "Not doing cross section fit." << endl;
    }


  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  // do one systematic sample at a time and accumulate the posterior

  for (ipx=0;ipx<nmc_req;ipx++)
    {
      iens = 0;
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = gRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = gRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    
      nens = 1;
      if (ipx == 0)
	{
	  setdlcsn(nbinstot,nens,nobs,ens);
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  for (double b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
            {
              bayes_posterior.push_back(0);
	      bayes_posterior_points.push_back(b);
            }
	}
      int itmp = 0;
      for (double b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
	{
	  double p = cspdf(b,1.0,nbinstot,nens,nobs,ens,prior);
	  bayes_posterior[itmp] += p;
	  itmp++;
	}
    }


  // find the fit cross section and errors

  double lmax = 0;
  int imax=0;
  int itmp=0;

  for (double b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
    {
      double p = bayes_posterior[itmp];
      if (p>lmax)
	{ 
	  lmax=p; 
	  imax = itmp;
	  *xsfit = b;
	}
      itmp++;
    }

  // normalize to unit area, using the trapezoidal rule

  int jsiz = bayes_posterior.size();
  int ic;
  Double_t bptot = 0;
  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
  bptot -= 0.5*bayes_posterior[0];
  bptot -= 0.5*bayes_posterior.back();
  if (bptot>0)
    {
      Double_t scaleb = 1.0/(bptot);
      bhnorm = scaleb;
      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
    }

  // start at the maximum, and integrate left and right until we get 68%

  int ilow=imax;
  int ihigh=imax;
  double psum = bayes_posterior[imax];
  double psumtrap = 0.0;
  do 
    {
      int ibest = ilow;
      double bpbest = 0.0;
      if (ilow>0)
	{
	  ibest = ilow-1;
	  bpbest = bayes_posterior[ibest];
	}
      if (ihigh<jsiz-1)
	{
	  if (bayes_posterior[ihigh+1] > bpbest)
	    {
	      ibest = ihigh+1;
	      bpbest = bayes_posterior[ibest];
	    }
	}
      psum += bpbest;
      if (ibest == ilow-1) 
	{ ilow --; }
      else
	{ ihigh ++; }
      psumtrap = psum - 0.5*bayes_posterior[ilow] - 0.5*bayes_posterior[ihigh];
    }
  while (psumtrap < 0.68);
  *downerr = *xsfit-bayes_posterior_points[ilow];
  *uperr = bayes_posterior_points[ihigh] - *xsfit;

  delete[] nobs;
  delete[] ens;
}


// quick and dirty posterior integrator -- assumes bayes_posterior and bayes_posterior_points
// have been filled.

Double_t mclimit_csm::quickbint(Double_t beta)
{
  // normalize to unit area, using the trapezoidal rule

  int jsiz = bayes_posterior.size();
  int ic;
  Double_t bptot = 0;
  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
  bptot -= 0.5*bayes_posterior[0];
  bptot -= 0.5*bayes_posterior.back();
  if (bptot>0)
    {
      Double_t scaleb = 1.0/(bptot);
      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
    }

  // count bacwards from the end until we get 1-beta of the integral
  // do a little trapezoidal trick to add back part of the interval if we overstep

  bptot = 0;
  for (ic=jsiz-1;ic>=0;ic--)
    {
      bptot += bayes_posterior[ic];
      if (bptot >= 1-beta) return(bayes_posterior_points[ic]+
				  bayes_interval_step*((bptot-(1-beta))/bayes_posterior[ic]));
    }
  return(0);
}

// expected fitted cross section distribuitons -- run pseudoexperiments from the test hypothesis
// and give distributions of fitted values.   Use bh_xsfit above to get the cross section in each one
// in order to save memory and make the code more maintainable.

void mclimit_csm::bh_xsfit_expect(Int_t npx, Double_t *xsfitavg, Double_t *m2s, 
				  Double_t *m1s, Double_t *med, Double_t *p1s, Double_t *p2s)
{
  int i,ipx,nchans;
  vector<double> xslist;
  double xsfit=0;
  double downerr=0;
  double uperr=0;

  // run pseudoexperiments and fit the cross section for each one. -- Use test_hypothesis
  // to generate the pseudoexperiments.

  nchans = test_hypothesis->channame.size();
  xslist.clear();
  TH1** pdarray = new TH1*[nchans];     // data and pseudodata stay as TH1's
  TH1** dhsavearray = new TH1*[nchans];
  char *pdname;

  for (i=0;i<nchans; i++)
    {
      pdname = new char[strlen(test_hypothesis->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,test_hypothesis->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) datahist[i]->Clone(pdname);
      delete[] pdname;
      dhsavearray[i] = datahist[i];
      datahist[i] = pdarray[i];
    }

  for (ipx=0;ipx<npx;ipx++)
    {
      test_hypothesis->single_pseudoexperiment(pdarray);
      bh_xsfit(&xsfit,&downerr,&uperr);

      if (pxprintflag)
	{
	  cout << "Marginalized fit px: " << xsfit << " + " << uperr << " - " << downerr << endl;
	}
      xslist.push_back(xsfit);
      *xsfitavg += xsfit;
    }
  if (xslist.size()>0) 
    {
      *xsfitavg /= xslist.size();
    }

  std::sort(xslist.begin(),xslist.end());
  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *m2s = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *m1s = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *med = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *p1s = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *p2s = xslist[i];

  for (i=0;i<nchans; i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  for (i=0;i<nchans; i++)
    {
      datahist[i] = dhsavearray[i];
    }
  delete[] dhsavearray;

}

// scan over two signals -- always assume they are in the same order in all channels.  Assume
// a flat prior in the two signals and print out the marginalized posterior.
// If more than two signals are present, the first one in  each channel is called signal 1,
// and the sum of all others is called signal 2

// update March 20, 2009 -- do one systematic sample at a time.

void mclimit_csm::bh_2d_scan(Double_t s1low, Double_t s1high, Double_t ds1,
			     Double_t s2low, Double_t s2high, Double_t ds2,
			     // arguments from here are given defaults in mclimit_csm.h 
                             // for backwards compatibility
			     Double_t *s1fit, Double_t *s2fit,
                             Int_t pflag2d, Double_t *s1test, Double_t *s2test,
                             Int_t *in68, Int_t *in95)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1Type* ht;
  Double_t r;
  Int_t nobstot;
  const PRIOR prior=flat;

  double xsig1,xsig2;
  vector <double> xsig1v;
  vector <double> xsig2v;
  vector <double> postv;
  vector <double> postint;

  int ibest = -1;
  double rtest;
  double rbest = 0;

  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB2D* ens = new EB2D[nbinstot];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments
  // but split up the first and subsequent signals.

  for (ipx=0;ipx<nmc_req;ipx++)
    {
      iens = 0;
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e1 = 0;
                  ens[iens].e2 = 0;
                  ens[iens].b = 0;
		  int isc = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = gRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = gRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ 
			  if (isc == 0)
			    { 
			      ens[iens].e1 += r;
			      isc++;
			    }
			  else
			    {
			      ens[iens].e2 += r;
			    }
			}
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}

      nens = 1;
      if (ipx == 0)
	{
	  setdlcsn2d(nbinstot,nens,nobs,ens);
	  for (xsig1=s1low;xsig1<=s1high;xsig1+=ds1)
	    {
	      for (xsig2=s2low;xsig2<=s2high;xsig2+=ds2)
		{
		  xsig1v.push_back(xsig1);
		  xsig2v.push_back(xsig2);
		  postv.push_back(0);
		  postint.push_back(0);

		  // identify which entry in the posterior grid is closest to the
		  // (s1test,s2test) point.  This will be used to determine if
		  // (s1test,s2test) is inside 68 or 95% CL.

		  if (s1test != 0 && s2test != 0)
		    {
                      rtest = sqrt((xsig1 - *s1test)*(xsig1 - *s1test) + (xsig2 - *s2test)*(xsig2 - *s2test));
		      if (ibest == -1 || (rtest < rbest))
			{
			  ibest = xsig1v.size() - 1;
			  rbest = rtest;
			}
		    }
		}
	    }
	}
      int itmp = 0;
      for (xsig1=s1low;xsig1<=s1high;xsig1+=ds1)
	{
	  for (xsig2=s2low;xsig2<=s2high;xsig2+=ds2)
	    {
	      postv[itmp] += cspdf2d(xsig1,xsig2,1.0,nbinstot,nens,nobs,ens,prior);
	      itmp ++;
	    }
	}
    }

  // to prevent a segfault in case s1test and s2test weren't given but in68 and in95 was given.
  if (ibest<0) ibest = 0;

  // for now just print out the posterior (normalized to unit sum over the interval chosen)

  double psum = 0;
  int itot = postv.size();
 
  for (i=0;i<itot;i++) psum += postv[i];
  if (psum > 0)
    {
      for (i=0;i<itot;i++) postv[i] = postv[i]/psum;
    }
  else
    {
      if (pflag2d) cout << "bh_2d_scan -- normalization failed.  Not printing posterior" << endl;
    }

  // find the best fit, and the integrals needed to compute the
  // 68% and 95% regions -- but only if it's normalized.

  if (in68) *in68 = 0;
  if (in95) *in95 = 0;

  if (psum>0)
    {
      int *idx = new int[itot];
      TMath::Sort(itot,&(postv[0]),idx,kTRUE);
      if (pflag2d) cout << "Two-Dimensional Maximum Posterior: " << xsig1v[idx[0]] << " " << xsig2v[idx[0]] << endl;
      if (s1fit) *s1fit = xsig1v[idx[0]];
      if (s2fit) *s2fit = xsig2v[idx[0]];
      double ps2 = 0;
      for (i=0;i<itot;i++)
	{
	  ps2 += postv[idx[i]];
	  postint[idx[i]] = ps2;
	}
      if (in68)
	{
	  if (postint[ibest]<=0.68) *in68=1;
	}
      if (in95)
	{
	  if (postint[ibest]<=0.95) *in95=1;
	}
      delete[] idx;
    }

  if (pflag2d)
    {
       for (i=0;i<itot;i++)
         {	   
            cout << "bh_2d_scan: " << xsig1v[i] << " " << xsig2v[i] << " " 
	      << postv[i] << " " << postint[i] << endl;
         }
    }

  delete[] nobs;
  delete[] ens;
}

// scan over two signals -- always assume they are in the same order in all channels.  Assume
// a flat prior in the two signals and print out the marginalized posterior.
// If more than two signals are present, the first one in  each channel is called signal 1,
// and the sum of all others is called signal 2

// This routine below runs pseudoexperiments to get the expected distributions.  It does all the
// calculations that bh_2d_scan does for the real data, but for the pseudoexperiments instead, and prints
// out the 2D posterior distributions, cumulative integrals, and best-fit points for each pseudoexperiment.
// As an added feature, you put in the signal scale factors used in the pseudoexperiment generation
//  (testhyp is used for pseudoexperiments, and testhyp_pe -- yes, I know that seems backwards -- is the
//  model used to compute the posterior), relative to the model used to test as s1true and s2true, and it
// will compute 68% and 95% coverage fractions based on the nearest grid point.  Verbose printout is
// supplied.

// Modified April 2, 2009 -- call the bh_2d_scan routine instead, after a memory optimization.  Same
// functionality as before.

void mclimit_csm::bh_2d_scan_expect(Int_t npx, Double_t s1true, Double_t s2true, 
                                    Double_t s1low, Double_t s1high, Double_t ds1,
				    Double_t s2low, Double_t s2high, Double_t ds2)
{
  int i,ipx,nchans;

  // run pseudoexperiments and fit the cross section for each one. -- Use test_hypothesis
  // to generate the pseudoexperiments.

  nchans = test_hypothesis->channame.size();
  TH1** pdarray = new TH1*[nchans];
  TH1** dhsavearray = new TH1*[nchans];
  char *pdname;

  for (i=0;i<nchans; i++)
    {
      pdname = new char[strlen(test_hypothesis->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,test_hypothesis->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) datahist[i]->Clone(pdname);
      delete[] pdname;
      dhsavearray[i] = datahist[i];      // save data pointers (if any) so we don't clobber them with our pseudodata
      datahist[i] = pdarray[i];          // wire up new pointers to pseudodata
    }

  int n68 = 0;
  int n95 = 0;

  for (ipx=0;ipx<npx;ipx++)
    {
      test_hypothesis->single_pseudoexperiment(pdarray);
      Double_t s1fit,s2fit;
      Int_t in68,in95;
      bh_2d_scan(s1low,s1high,ds1,s2low,s2high,ds2,&s1fit,&s2fit,0,&s1true,&s2true,&in68,&in95);
      cout << "Two-Dimensional Maximum Posterior px: " << s1fit << " " << s2fit << endl;
      if (in68) {
        cout << "bh_2d_scan_expect true signal falls within 68 percent region" << endl;
	n68++;
      }
      if (in95) {
        cout << "bh_2d_scan_expect true signal falls within 95 percent region" << endl;
	n95++;
      }
      if (in68==0 && in95==0) {
	cout << "bh_2d_scan_expect true signal falls outside 95 percent region" << endl;
      }
    }

  cout << "bh_2d_scan_expect: n68, npx, Fraction of time true signal lies within 68 percent region: " <<
    n68 << " " << npx << " " << ((double) n68)/((double) npx) << endl;
  cout << "bh_2d_scan_expect: n95, npx, Fraction of time true signal lies within 95 percent region: " <<
    n95 << " " << npx << " " << ((double) n95)/((double) npx) << endl;

  delete[] pdarray;
  for (i=0;i<nchans; i++)
    {
      datahist[i] = dhsavearray[i];
    }
  delete[] dhsavearray;
}



/*-------------------------------------------------------------------------*/
/* Coverage checker for Joel's Bayesian limit calc.  Based on              */
/* bayes_heinrich_withexpect, but now uses test_hypothesis_pe scaled       */
/* so the signal is at a user-setabble desired rate (a useful test is to   */
/* test the coverage for the rate excluded by bayes_heinrich, but one can  */
/* also test the coverage for any other signal rate, scaling the           */
/* signal in testhyp_pe by sflimit. The coverage should be 100% for zero   */
/* signal, for example.    The px's are done                               */
/* assuming the signal+background is present, and the false exclusion rate */
/* is computed.                                                            */
/*-------------------------------------------------------------------------*/
/* arguments:  beta: credibility level:L  0.95 for 95% CL limits
   sfscale:  INPUT -- desired multiplier on the signal in the testhyp_pe hypothesis
   for which we'd like to check the coverage.
   npx:      Number of pseudoexperiments to run to compute fales exclusion rate
   falsex:   false exclusion rate:  Should be no more than 1-beta.

*/

void mclimit_csm::bayes_heinrich_coverage_check(Double_t beta,
                                                Double_t sfscale,
					        Int_t npx,
                                                Double_t* falsex)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1Type* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;
  double *xgl;
  double *lwgl;

  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = gRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = gRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.  We will enlarge these arrays
  //later if the need arises.

  nglmax = nobstot;
  if (nglmax<10000) {nglmax = 10000;}
  xgl = new double[nglmax];
  lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;
  // do not compute sflimit here -- input it instead as an adjustable parameter

  //  *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);
  
  // Run signal+background pseudoexperiments at the observed limit, and see what
  // the distribution of limits we get out is.  Limits are computed using the
  // same Bayesian ensemble with the unscaled test_hypothesis_pe and so the
  // limits that are more restrictive than *sfscale are false exclusions.

  csm_model* testhyppescale = test_hypothesis_pe->scalesignal(sfscale);

  cslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) testhyppescale->channame.size(); i++)
    {
      pdname = new char[strlen(testhyppescale->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,testhyppescale->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = copy_TH1<TH1> (*testhyppescale->chanmodel[i]->histotemplate[0], pdname).release();
      delete[] pdname;
    }
  for (ipx=0;ipx<npx;ipx++)
    {
      testhyppescale->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
		  if (nbinsy==1)
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
		  else
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		  nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
		  ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  if (nglmax < nobstotlist[nobsindex[0]]/2 + 1)
    {
      nglmax = nobstotlist[nobsindex[0]]/2 + 1; 
      delete[] xgl;
      delete[] lwgl;
      xgl = new double[nglmax];
      lwgl = new double[nglmax];
    }

  ngl = 0;
  for (ipx=0;ipx<npx;ipx++)
    {
      if (ipx>0)
	{ if (nobstotlist[nobsindex[ipx]] != nobstotlist[nobsindex[ipx-1]])
	  { 
	    ngl = 0;
	  }
	}
      Double_t unc;
      Double_t p = (Double_t) cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,&unc);
      if (pxprintflag)
	{
          cout << "Bayes Coverage Check px: " << p << endl;
	}
      cslist.push_back(p);
    }

  Int_t nfalse = 0;
  for (i = 0; i<npx; i++)
    {
      if (cslist[i] <= sfscale)
	{
	  nfalse++;
	}
    }

  *falsex = ((Double_t) nfalse)/((Double_t) npx);

  // clean up

  for (i=0;i<(Int_t) testhyppescale->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
  delete testhyppescale;
}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*    genlimit Bayesian code from Joel                                     */
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

/*  Joel Heinrich  8 April 2005

Returns cross section posterior p.d.f. evaluated at s.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/

void setdlcsn(int nchan, int nens, int nobs[], const EB* ens)
{
  dlcsn = 0;
  int set=0;
  double tmax=0;
  double tmin=0;
  double s=0;

  int i,k;
  double lgp = 0;
  const EB* p = ens;

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*p->e + p->b;
      const int n = nobs[k];
      esum += p->e;
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    if (!set)
      { tmax = t; tmin = t; set = 1;}
    else
      { tmax=max(tmax,t); tmin = min(tmin,t); }
    dlcsn = -tmax;
  }
}


void setdlcsn2d(int nchan, int nens, int nobs[], const EB2D* ens)
{
  dlcsn2d = 0;
  int set=0;
  double tmax=0;
  double tmin=0;
  double s=0;

  int i,k;
  double lgp = 0;
  const EB2D* p = ens;

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*(p->e1 + p->e2) + p->b;
      const int n = nobs[k];
      esum += (p->e1 + p->e2);
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    if (!set)
      { tmax = t; tmin = t; set = 1;}
    else
      { tmax=max(tmax,t); tmin = min(tmin,t); }
    dlcsn2d = -tmax;
  }
}

double cspdf(double s,double norm,
	     int nchan,int nens,const int nobs[],const EB* ens,PRIOR prior) {
  int i,k;
  double sum = 0, lgp = 0;
  const EB* p = ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*p->e + p->b;
      const int n = nobs[k];
      esum += p->e;
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    sum += (prior==flat) ? exp(t+dlcsn) : esum*exp(t+dlcsn);
  }
  //cout << "pdf: " << s << " " << sum/nens << endl;
  return sum/(norm*nens);
}

// Tom Junk 21 April 2008 -- generalize to two signals

double cspdf2d(double s1, double s2, double norm,
	       int nchan,int nens,const int nobs[],const EB2D* ens,PRIOR prior) {
  int i,k;
  double sum = 0, lgp = 0;
  const EB2D* p = ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s1*p->e1+ s2*p->e2 + p->b;
      const int n = nobs[k];
      esum += (p->e1+p->e2);
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    sum += (prior==flat) ? exp(t+dlcsn2d) : esum*exp(t+dlcsn2d);
  }
  //cout << "pdf: " << s << " " << sum/nens << endl;
  return sum/(norm*nens);
}

/*  Joel Heinrich  8 April 2005

Function returns integral from xlo to infinity.

*uncertainty (if not null pointer) returned with uncertainty of
integral due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf


*/


double csint(double xlo,int nchan,int nens,const int nobs[],const EB* ens,
	     int* ngl,double xgl[],double lwgl[],
	     PRIOR prior,double* uncertainty) {
  int i,ntot=0;
  double sum=0, sum2=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    const double t = csint0(xlo,logscale,nchan,nobs,p,*ngl,xgl,lwgl,prior);
    sum += t;
    sum2 += t*t;
    p+=nchan;
  }
  sum /= nens;
  sum2 /= nens;
  if(uncertainty)
    *uncertainty = (nens>1) ? sqrt((sum2-sum*sum)/(nens-1)) : 1.0 ;
  return sum;
}

double csint0(double xlo,double logscale,
	      int nchan,const int nobs[],const EB chan[],
	      int ngl,const double xgl[],const double lwgl[],
	      PRIOR prior) {
  int i,k;
  double sum=0, esum=0, bsum=0, resum;

  for(i=0;i<nchan;++i) {
    esum += chan[i].e;
    bsum += chan[i].b + xlo*chan[i].e;
  }
  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo;
    double t = logscale-bsum, v;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t += nobs[i] * log( xr*chan[i].e + chan[i].b );
    if (dlcsn==0)
      {
	dlcsn = -t;
	if (dlcsn==0)
	  {
	    dlcsn = 1.0;
	  }
      }
    sum += v = exp(lwgl[k]+t+dlcsn);
    if(v<DBL_EPSILON*sum) break;
  }
  if(prior==flat)
    sum *= resum;
  return sum;
}


/*      Joel Heinrich  8 April 2005

returns:
*int1 = integral from xlo1 to xhi
*int2 = integral from xlo2 to xhi
*v11  = variance of *int1
*v12  = covariance between *int1 and *int2
*v22  = variance of *int2

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/



void csint2cut(double xlo1,double xlo2,double xhi,
	       int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],PRIOR prior,
	       double* int1,double* int2,
	       double* v11,double* v12,double* v22) {
  int i,ntot=0;
  double sum1=0, sum21=0, sum2=0, sum22=0, sump=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    double t1=0, t2=0;
    csint02cut(xlo1,xlo2,xhi,logscale,nchan,nobs,p,*ngl,xgl,lwgl,
	       prior,&t1,&t2);
    sum1 += t1;
    sum21 += t1*t1;
    sum2 += t2;
    sum22 += t2*t2;
    sump += t1*t2;
    p+=nchan;
  }

  {
    const double rnens = 1.0/nens;
    sum1 *= rnens;
    sum21 *= rnens;
    sum2 *= rnens;
    sum22 *= rnens;
    sump *= rnens; 
    if(nens>1) {
      const double rn1 = 1.0/(nens-1);
      *v11 = (sum21-sum1*sum1)*rn1;
      *v22 = (sum22-sum2*sum2)*rn1;
      *v12 = (sump-sum1*sum2)*rn1;
    } else {
      *v11 = 1;
      *v22 = 1;
      *v12 = 0;
    }
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

void csint02cut(double xlo1,double xlo2,double xhi,double logscale,
		int nchan,const int nobs[],const EB chan[],
		int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		double* int1,double* int2) {
  int i,k;
  double sum1=0, sum2=0, sum3=0, esum=0, bsum1=0, bsum2=0, bsum3=0, resum;

  for(i=0;i<nchan;++i) {
    const double ee=chan[i].e, bb=chan[i].b;
    esum += ee;
    bsum1 += bb + xlo1*ee;
    bsum2 += bb + xlo2*ee;
    bsum3 += bb + xhi*ee;
  }

  if(esum==0 && prior==corr) {
    *int1 = *int2 = 0;
    return;
  }

  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo1;
    double t1 = logscale-bsum1, v1;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t1 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum1 += v1 = exp(lwgl[k]+t1+dlcsn);
    if(v1<DBL_EPSILON*sum1) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo2;
    double t2 = logscale-bsum2, v2;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t2 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum2 += v2 = exp(lwgl[k]+t2+dlcsn);
    if(v2<DBL_EPSILON*sum2) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xhi;
    double t3 = logscale-bsum3, v3;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t3 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum3 += v3 = exp(lwgl[k]+t3+dlcsn);
    if(v3<DBL_EPSILON*sum3) break;
  }

  if(prior==flat) {
    *int1 = (sum1-sum3)*resum;
    *int2 = (sum2-sum3)*resum;
  } else {
    *int1 = sum1-sum3;
    *int2 = sum2-sum3;
  }
  return;
}


/*  Joel Heinrich  8 April 2005

Returns cross section upper limit.

*uncertainty (if not null pointer) returned with uncertainty of
limit due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/


double cslimit(double beta,int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],
	       PRIOR prior,double* uncertainty) {
  const double eps=1.0e-6;

  
  /*
    cout << "nchan: " << nchan << endl;
    cout << " nens: " << nens << endl;
    cout << " prior: " << prior << endl;
    int i;
    for (i=0;i<nchan;i++)
    {
    if (nobs[i]>0)
    {  cout << "i, n, b, s: " << i << " " << nobs[i] << " " << ens[i].b << " " << ens[i].e << endl;
    //cout << "nobs(" << i << ") = " << nobs[i] << endl;
    }
    }
  */

  dlcsn = 0;

  double norm = csint(0,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  double limit = galim(beta,nchan,nens,nobs,ens);
  double dl=limit, rpdf=0;
  double lo=0, hi=1e200;

  while(fabs(dl)>1.0e-10*limit) {
    const double pbeta =
      1-csint(limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL)/norm;
    rpdf = 1/cspdf(limit,norm,nchan,nens,nobs,ens,prior);
    if(pbeta>beta) {
      hi=limit*(1+eps);
    } else {
      lo=limit*(1-eps);
    }
    dl = (pbeta-beta)*rpdf;
    if(limit-dl>=lo && limit-dl<=hi) {
      limit -= dl;
    } else {
      dl = limit - 0.5*(hi+lo);
      limit = 0.5*(hi+lo);
    }
  }

  if (uncertainty) {
    double i1=0, i2=0, v11=0, v22=0, v12=0;
    const double c = 1-beta;
    csint2(0,limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,
	   &i1,&i2,&v11,&v12,&v22);
    *uncertainty = rpdf*sqrt(v22 + v11*c*c - 2*v12*c)/norm   ;
  }


  return limit;
}


/*  Joel Heinrich  8 April 2005

Returns cross section upper limit.

*uncertainty (if not null pointer) returned with uncertainty of
limit due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/




double cscutlimit(double beta,double smax,
		  int nchan,int nens,const int nobs[],const EB* ens,
		  int* ngl,double xgl[],double lwgl[],
		  PRIOR prior,double* uncertainty) {
  const double norm = csint(0,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  const double tail = csint(smax,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  const double eps=1.0e-6;
  double limit = galim(beta,nchan,nens,nobs,ens), rpdf=0;
  double dl=limit, lo=0, hi=smax;

  if(beta<=0) return 0;
  if(beta>=1) return smax;

  if(limit>smax || limit<0) dl = limit = 0.5*smax;

  while(fabs(dl)>1.0e-10*limit) {
    double pbeta =
      1-(csint(limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL)-tail)/
      (norm-tail);
    rpdf = 1/cspdf(limit,norm-tail,nchan,nens,nobs,ens,prior);
    if(pbeta>beta) {
      hi=limit*(1+eps);
    } else {
      lo=limit*(1-eps);
    }
    dl = (pbeta-beta)*rpdf;
    if(limit-dl>=lo && limit-dl<=hi) {
      limit -= dl;
    } else {
      dl = limit - 0.5*(hi+lo);
      limit = 0.5*(hi+lo);
    }

  }


  if (uncertainty) {
    double i1=0, i2=0, v11=0, v22=0, v12=0;
    const double c = 1-beta;

    csint2cut(0,limit,smax,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,
	      &i1,&i2,&v11,&v12,&v22);
    *uncertainty = rpdf*sqrt(v22 + v11*c*c - 2*v12*c)/(norm-tail)   ;
  }


  return limit;
}


/*         Joel Heinrich  24 March 2005

returns crude Gaussian approximation to upper limit.
For use as a starting point.
*/

double galim(double beta,int nchan,int nens,const int nobs[],const EB* ens) {
  double mean=0,sigma=0;
  gameansigma(&mean,&sigma,nchan,nens,nobs,ens);
  return mean-sigma*arcfreq( (1-freq(-mean/sigma))*(1-beta) );
}



void gameansigma(double *mean,double *sigma,
		 int nchan,int nens,const int nobs[],const EB* ens) {

  double sum=0,sum2=0,vsum=0;
  const EB* p = ens;
  int i,j;
  for(i=0;i<nens;++i) {
    double s=0,s2=0;
    for(j=0;j<nchan;++j) {
      const int n = nobs[j];
      const double eps = p->e;
      s += (n-p->b)*eps/(n+1);
      s2 += eps*eps/(n+1);
      ++p;
    }
    s /= s2;
    vsum += 1/s2;
    sum += s;
    sum2 += s*s;
  }
  
  *mean = sum/nens;
  *sigma = sqrt(vsum/nens + sum2/nens - (*mean)*(*mean));
  return;
}

#define rdfreq(x) (exp(0.5*(x)*(x))*2.50662827463100050242)

#define C0 2.515517
#define C1 0.802853
#define C2 0.010328
#define D0 1.0
#define D1 1.432788
#define D2 0.189269
#define D3 0.001308

double arcfreq(double y) {
  const double yy = (y>0.5) ? 1-y : y, t = sqrt(-2*log(yy));
  double x = (C0+t*(C1+t*C2))/(D0+t*(D1+t*(D2+t*D3))) - t;
  x -= (freq(x) - yy)*rdfreq(x);
  x -= (freq(x) - yy)*rdfreq(x);
  return (y>0.5) ? -x : x;
}


/*

Joel Heinrich
February 10 2005

Returns Gauss-Laguerre quadrature abscissas and log(weights) which can
be used to approximate

integral u=0 to infinity pow(u,alpha)*exp(-u)*f(u) du
as
sum k=0 to n-1  exp(lw[k])*f(x[k])

or equivalently

sum k=0 to n-1  exp(lw[k]+log(f(x[k])))

The quadrature is exact for polynomial f of degree 2n-1 or less.

*/


void gausslaguerre(double x[],double lw[],int n,double alpha){
  const int nshift = 20;
  const double shift = 1<<nshift, rshift=1/shift;
  int i;
  double z=0;
  
  for(i=0;i<n;++i) {
    int j=0, k=2, nscale=0;
    double dz=0.0, p1=0, p2=0;
    if(i==0) {
      z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
    } else if(i==1) {
      z += (15.0+6.25*alpha)/(1.0+2.5*n+0.9*alpha);
    } else if(i==2) {
      const double ai=i-1;
      z += ( (1.0+2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1.0+3.5*ai) )*
	(z-x[i-2])/(1.0+0.3*alpha);
    } else if(i==3) {
      z = 3.0*(x[2]-x[1])+x[0];
    } else if(i==4) {
      z = 4.0*x[3] - 6.0*x[2] + 4.0*x[1] - x[0];
    } else if(i==5) {
      z = 5.0*x[4] - 10.0*x[3] + 10.0*x[2] - 5.0*x[1] + x[0];
    } else {
      z = 6.0*x[i-1] - 15.0*x[i-2] + 20.0*x[i-3] -
	15.0*x[i-4] + 6.0*x[i-5] - x[i-6];
    }
    while(k>0) {
      p1=1;
      p2=0;
      nscale=0;
      z -= dz;
      for(j=1;j<=n;++j){
	const double p3=p2;
	p2=p1;
	p1=((2*j-1+alpha-z)*p2 - (j-1+alpha)*p3)/j;
	if(fabs(p2)>shift) {
	  ++nscale;
	  p1 *= rshift;
	  p2 *= rshift;
	}
      }
      dz = p1*z/(n*p1-(n+alpha)*p2);
      if(fabs(dz)<1.0e-10*z)--k;
    }
    x[i]=z;
    lw[i] = log(z/(p2*p2)) - 2*nshift*nscale*M_LN2 ;
  }
  
  {
    double t = 0.0;
    for(i=n-1;i>=0;--i)
      t += exp(lw[i]);
    t = lgamma(alpha+1)-log(t);
    for(i=0;i<n;++i)
      lw[i] += t;
  }

  return;
}

/*    Joel Heinrich  8 April 2005

returns:
*int1 = integral from xlo1 to infinity
*int2 = integral from xlo2 to infinity
*v11  = variance of *int1
*v12  = covariance between *int1 and *int2
*v22  = variance of *int2


See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf


*/



void csint2(double xlo1,double xlo2,
	    int nchan,int nens,const int nobs[],const EB* ens,
	    int* ngl,double xgl[],double lwgl[],PRIOR prior,
	    double* int1,double* int2,
	    double* v11,double* v12,double* v22) {
  int i,ntot=0;
  double sum1=0, sum21=0, sum2=0, sum22=0, sump=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    double t1=0, t2=0;
    csint02(xlo1,xlo2,logscale,nchan,nobs,p,*ngl,xgl,lwgl,prior,&t1,&t2);
    sum1 += t1;
    sum21 += t1*t1;
    sum2 += t2;
    sum22 += t2*t2;
    sump += t1*t2;
    p+=nchan;
  }

  {
    const double rnens = 1.0/nens;
    sum1 *= rnens;
    sum21 *= rnens;
    sum2 *= rnens;
    sum22 *= rnens;
    sump *= rnens; 
    if(nens>1) {
      const double rn1 = 1.0/(nens-1);
      *v11 = (sum21-sum1*sum1)*rn1;
      *v22 = (sum22-sum2*sum2)*rn1;
      *v12 = (sump-sum1*sum2)*rn1;
    } else {
      *v11 = 1;
      *v22 = 1;
      *v12 = 0;
    }
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

void csint02(double xlo1,double xlo2,double logscale,
	     int nchan,const int nobs[],const EB chan[],
	     int ngl,const double xgl[],const double lwgl[],PRIOR prior,
	     double* int1,double* int2) {
  int i,k;
  double sum1=0, sum2=0, esum=0, bsum1=0, bsum2=0, resum;

  for(i=0;i<nchan;++i) {
    esum += chan[i].e;
    bsum1 += chan[i].b + xlo1*chan[i].e;
    bsum2 += chan[i].b + xlo2*chan[i].e;
  }
  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo1;
    double t1 = logscale-bsum1, v1;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t1 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum1 += v1 = exp(lwgl[k]+t1+dlcsn);
    if(v1<DBL_EPSILON*sum1) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo2;
    double t2 = logscale-bsum2, v2;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t2 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum2 += v2 = exp(lwgl[k]+t2+dlcsn);
    if(v2<DBL_EPSILON*sum2) break;
  }

  if(prior==flat) {
    sum1 *=resum;
    sum2 *=resum;
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

// access the globals in here with accessor methods.

double mclimit_csm::getdlcsn() { return dlcsn; }
double mclimit_csm::getdlcsn2d() { return dlcsn2d; }
double mclimit_csm::getbhnorm() { return bhnorm; }


//          Copyright Nils Krumnack 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Please feel free to contact me (nils@fnal.gov) for bug reports,
// feature suggestions, praise and complaints.

//
// includes
//

//#include "FastTH1.hh"

#include <TAxis.h>
#include <TH2.h>
#include <TH3.h>

//
// method implementations
//

FastTH1 ::
FastTH1 (const TH1& that)
  : name_ (that.GetName()), title_ (that.GetTitle()),
    dimension_ (that.GetDimension())
{
  assert (unsigned (that.GetDimension()) <= max_dim && "requirement failed"); //spec

  for (unsigned dim = 0; dim != max_dim; ++ dim)
  {
    nbins_[dim] = 1;
    mult_[dim] = 0;
  };

  for (unsigned dim = 0; dim != dimension_; ++ dim)
  {
    TAxis *axis = NULL;
    switch (dim)
    {
    case 0:
      axis = that.GetXaxis();
      break;
    case 1:
      axis = that.GetYaxis();
      break;
    case 2:
      axis = that.GetZaxis();
      break;
    };
    nbins_[dim] = axis->GetNbins();
    low_[dim] = axis->GetBinLowEdge (1);
    high_[dim] = axis->GetBinUpEdge (nbins_[dim]);

    mult_[dim] = 1;
    for (unsigned dim2 = 0; dim2 != dim; ++ dim2)
      mult_[dim2] *= nbins_[dim];
  };

  bin_contents_.reserve (nbins_[0] * nbins_[1] * nbins_[2]);
  bin_errors_.reserve (nbins_[0] * nbins_[1] * nbins_[2]);
  for (unsigned binx = 0; binx != nbins_[0]; ++ binx)
  {
    for (unsigned biny = 0; biny != nbins_[1]; ++ biny)
    {
      for (unsigned binz = 0; binz != nbins_[2]; ++ binz)
      {
	bin_contents_
	  .push_back (that.GetBinContent (1 + binx, 1 + biny, 1 + binz));
	const float error = that.GetBinError (1 + binx, 1 + biny, 1 + binz);
	assert (error >= 0 && "internal consistency error"); //spec
	bin_errors_.push_back (error);
      };
    };
  };

  test_invariant ();
};



std::auto_ptr<TH1> FastTH1 ::
make_TH1 (std::string name) const
{
  test_invariant ();

  std::auto_ptr<TH1> result;

  if (name.empty())
    name = name_;
  switch (dimension_)
  {
  case 1:
    result.reset (new TH1F (name.c_str(), title_.c_str(),
			    nbins_[0], low_[0], high_[0]));
    break;
  case 2:
    result.reset (new TH2F (name.c_str(), title_.c_str(),
			    nbins_[0], low_[0], high_[0],
			    nbins_[1], low_[1], high_[1]));
    break;
  case 3:
    result.reset (new TH3F (name.c_str(), title_.c_str(),
			    nbins_[0], low_[0], high_[0],
			    nbins_[1], low_[1], high_[1],
			    nbins_[2], low_[2], high_[2]));
    break;
  };

  assert (result.get() != NULL && "internal consistency error"); //spec

  FastTH1::bin_contents_iter content; content = bin_contents_.begin();
  bin_errors_iter error = bin_errors_.begin();
  for (unsigned binx = 0; binx != nbins_[0]; ++ binx)
  {
    for (unsigned biny = 0; biny != nbins_[1]; ++ biny)
    {
      for (unsigned binz = 0; binz != nbins_[2]; ++ binz)
      {
	result->SetBinContent (1 + binx, 1 + biny, 1 + binz, *content);
	result->SetBinError (1 + binx, 1 + biny, 1 + binz, *error);
	++ content;
	++ error;
      };
    };
  };

  assert (content == bin_contents_.end() && "internal consistency error"); //spec
  assert (error == bin_errors_.end() && "internal consistency error"); //spec

  assert (result.get() != NULL && "postcondition failed"); //spec
  return result;
};


