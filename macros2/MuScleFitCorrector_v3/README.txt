###########################################################
INSTRUCTIONS HOW TO USE MUSCLEFIT CORRECTION V1
###########################################################


- The package consists of:
	o MuScleFitCorrector.h: this is the main class implementing the corrections
	o ScaleFunct.h: this class defines the scale function and provides the flexibility needed to eventually redefine it
	o ResolFunct.h: this class defines the resolution function and provides the flexibility needed to eventually redefine it
	o MuScleFit_YEAR_DATASET_CMSSWREL.txt: this text file contains the parameters of the scale function. There exist one .txt file for every set of correction (e.g. MuScleFit_2011_DATA_42X.txt )

- Copy the classes of the package in your code (as a default it is assumed that all the files of the package are in the same directory)

- Implement the correction in your code:

	o Include the class in your code
	  #include "MuScleFitCorrector.h"

	o Define the path to the fit parameters text file
	  ...
	  TString fitParametersFile = "some/path/inmydirectories/FitParameters.txt"
	  ...

	o Instantiate the corrector object 
	  ...
	  MuScleFitCorrector* corrector_ = new MuScleFitCorrector(fitParametersFile)
	  ...

	o Apply the corrections
	  ...
	  TLorentzVector* myNegMuon = getMyNegMuonSomewhere();
	  TLorentzVector* myPosMuon = getMyPosMuonSomewhere();
	  ...
	  // The second argument of the correction is the charge of the muon (note that the correction itself does not depend on the charge, but we are scaling the curvature...)
	  corrector_->applyPtCorrection(*myNegMuon,-1);
	  corrector_->applyPtCorrection(*myPosMuon,1);
	  ...


	o Apply the smearing (!!2012 MC samples only!!)
	  ...	
	  corrector_->applyPtSmearing(*myNegMuon,-1);
	  corrector_->applyPtSmearing(*myPosMuon,1);
	  ...


ENJOY!!

If you have any question/problem with MuScleFit corrections please send us an e-mail:
- stefano.casasso@cern.ch
- migliore@to.infn.it
