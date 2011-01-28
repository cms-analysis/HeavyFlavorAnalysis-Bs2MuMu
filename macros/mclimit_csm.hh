

// mclimit_csm -- use a chisquared minimized over nuisance parameters
// to compute cls.

// Class to run the TMinuit minimization of T. Devlin's chisquared
// defined in CDF 3126, minimized over the nuisance parameters.

// version dated January 14, 2009
// Author:  Tom Junk, Fermilab.  trj@fnal.gov

#ifndef CSM_H
#define CSM_H

#include "TH1.h"
#include <vector>

using std::vector;

struct svstruct_s
{
  Int_t itemplate;   // which template histo this syst. variation applies to
  char *sysname;     // name of nuisance parameter this corresponds to
  Double_t sysfracl;
  Double_t sysfrach;
  // if there is no shape uncertainty associated with this dependence of a particular
  // template on a nuisance parameter, then these pointers should be set to zero.
  TH1 *lowshape;     // for shape uncertainty -- low histogram shape histo id
  TH1 *highshape;    // for shape uncertainty -- high histogram shape
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
   void add_template( TH1 *,      //template histogram
	   	      Double_t,   //scale factor to multiply template by to compare w/ data
                      Int_t,      // number of nuisance parameters (Gaussian of unit width)
                      char *[],   // nuisance parameter names 
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- low side
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
                      TH1 *[],    // array of low histogram shapes, one for each nuisance param
                      Double_t *, // number of sigma low for each nuisance parameter shape variation
                      TH1 *[],    // array of high histogram shapes, one for each nuisance param
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
   Double_t chisquared1(TH1 *);  // Inputs a pointer to data histogram -- computes the chisquared of this model
                                 // compared to the supplied data histogram, without minimizing over nuisance
                                 // parameters, a la Tom Devlin's note

   Int_t checkneg();             // check for negative bins

   void set_interpolation_style(INTERPSTYLE);  // either CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL
                                               // horizontal (csm_pvmorph) is the default if this is not called.

   friend class csm_model;       // so we can access systematic error names and limits
   friend class mclimit_csm;

 private:
  vector<TH1*> histotemplate;
  vector<TH1*> histotemplate_varied;
  vector<Double_t> sft;
  vector<Double_t> sft_varied;
  vector<Int_t> poissflag;
  vector<Int_t> scaleflag;
  vector<svstruct> syserr;
  INTERPSTYLE chan_istyle;
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
   void nuisance_response(Int_t, char *[], Double_t []); // updates the fluctuated version of the histogram templates
                                 // and scale factors inside the channel model
                                 // according to the nuisance parameters provided.  Inputs: number of nuisance
                                 // parameters and their names
   void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values
   void varysyst();   // randomly choose nuisance parameter values and call nuisance_response
   void print_nuisance_params();  // after calling varysyst, can print out the randomly chosen parameters and names
   void single_pseudoexperiment(TH1 *[]); // generate pseudodata for all the channels in this model.
   void list_nparams(vector<char *> *, vector<Double_t> *, vector<Double_t> *); // get a list of
      // unique nuisance parameter names for all the channels in this model and their most
      // restrictive lower and upper bounds.
   Double_t chisquared1(TH1**); // calls the chisquared routine for each channel model.

   friend class mclimit_csm;
   friend class csm;

 private:
  vector<char*> channame;
  vector<csm_channel_model*> chanmodel;
  Int_t lookup_add_channame(char *);  // look up the channel name in the channame vector.  If it's
                                       // not there, add it.
  vector<npcstruct> npcm;  // constraint equations between nuisance parameters

  vector<char*> npbname;   // upper and lower bounds
  vector<Double_t> npbhigh;
  vector<Double_t> npblow;
  vector<Double_t> npvalp; // handles to stored lists of nuisance parameter names and values for
  vector<char*> npnp;
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
        vector<Double_t> fitparam;   // parameters of the fit
	vector<Double_t> fiterror;   // errors on fit parameters
        vector<char*> fitparamname; //  names of fit parameters
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

void csm_interpolate_histogram(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

void csm_interpolate_histogram_noerr(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

// version to be used with cascading shape errors -- needs a central shape, a varied shape,
// and a shape to apply the variations to (which may not be either of the above, but the result
// of a different shape variation) startshape.  The output is outshape.

// csm_interpolate histogram2 calls csm_interpolate_histogram3 twice, once for the
// bin contents, once for the errors.

void csm_interpolate_histogram2(TH1* central, Double_t paramcentral,
                                TH1* varied, Double_t paramvaried,
                                TH1* startshape, 
                                TH1* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

// here's a version that just interpolates the bin contents but not the errors
//  (twice as fast)

void csm_interpolate_histogram2_noerr(TH1* central, Double_t paramcentral,
                                TH1* varied, Double_t paramvaried,
                                TH1* startshape, 
                                TH1* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

// this routine just interpolates the histogram contents

void csm_interpolate_histogram3(TH1* central, Double_t paramcentral,
                                TH1* varied, Double_t paramvaried,
                                TH1* startshape, 
                                TH1* outshape,
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

   vector<Double_t> bayes_posterior;
   vector<Double_t> bayes_posterior_points;

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
                                      Double_t sflimit,
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

   void bh_2d_scan(Double_t s1low, Double_t s1high, Double_t ds1,
                   Double_t s2low, Double_t s2high, Double_t ds2);

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
   vector<TH1*> datahist;
   vector<char*> dhname;            // channel names for each data histogram -- must match
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
};

#endif

#ifndef GENLIMIT

typedef struct {
  float e,b;
}EB;

typedef struct {
  float e1,e2,b;
}EB2D;

typedef enum {
  flat=10,
  corr=20
} PRIOR;

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

#endif
