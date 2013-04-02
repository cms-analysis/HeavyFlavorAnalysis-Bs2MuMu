#include "MuScleFitCorrector.hh"

/// Method to get correct pt
double MuScleFitCorrector::getCorrectPt( const TLorentzVector & lorentzVector , const int & chg) {
  
  // Loop on all the functions and apply them iteratively on the pt corrected by the previous function.
  double pt = lorentzVector.Pt();
  pt = ( scaleFunct_->scale( pt, lorentzVector.Eta(), lorentzVector.Phi(), chg, scaleParArray_) );
  return pt;
}

// Return the squared difference of relative pT resolutions data-MC
double MuScleFitCorrector::getSigmaPtDiffSquared(const double & pt, const double & eta){
  
  double sigmaPtData = resolFunct_->sigmaPt(pt, eta, resolDataParArray_);
  double sigmaPtMC = resolFunct_->sigmaPt(pt, eta, resolMCParArray_);
  return (sigmaPtData*sigmaPtData - sigmaPtMC*sigmaPtMC);
  
}


// Method to get correct pt
double MuScleFitCorrector::getSmearedPt( const TLorentzVector & lorentzVector , const int & chg) {
  
  // Loop on all the functions and apply them iteratively on the pt corrected by the previous function.
  double pt = lorentzVector.Pt();
  double eta = lorentzVector.Eta();
  double squaredDiff = getSigmaPtDiffSquared(pt,eta);
  double Cfact = 0.8;
  double sPar = Cfact*sqrt(squaredDiff);
  double curv = ((double)chg/pt);
  double smearedCurv = curv + fabs(curv)*(gRandom_->Gaus(0,sPar));
  
  return 1./((double)chg*smearedCurv);
    
}
  
// Method to apply correction to a TLorentzVector
void MuScleFitCorrector::applyPtCorrection( TLorentzVector& lorentzVector, const int & chg ){
  
  double eta = lorentzVector.Eta();
  double phi = lorentzVector.Phi();
  double m  = lorentzVector.M(); 
  double corrPt = getCorrectPt(lorentzVector, chg);
  lorentzVector.SetPtEtaPhiM(corrPt,eta,phi,m);
  
}

// Method to apply smearing to a TLorentzVector
void MuScleFitCorrector::applyPtSmearing(TLorentzVector& lorentzVector, const int & chg){
  double eta = lorentzVector.Eta();
  double phi = lorentzVector.Phi();
  double m  = lorentzVector.M(); 
  double smearedPt = getSmearedPt(lorentzVector, chg);
  lorentzVector.SetPtEtaPhiM(smearedPt,eta,phi,m);
  
}


void MuScleFitCorrector::convertToArrays()
{

  int resolParNum = resolFunct_->parNum();
  int scaleParNum = scaleFunct_->parNum();

  int resolDataParVecSize = resolDataParVec_.size();
  int resolMCParVecSize = resolMCParVec_.size();
  int scaleParVecSize = scaleParVec_.size();

  useResol_ = false;
  if (resolDataParVecSize!=0 && resolMCParVecSize!=0)  useResol_ = true;

  if( resolParNum != resolDataParVecSize && useResol_) {
    std::cout << "Error: inconsistent number of parameters: resolFunct has "<<resolParNum<<" parameters but "<<resolDataParVecSize<<" have been read from file" << std::endl;
    exit(1);
  }

  if( resolParNum != resolMCParVecSize && useResol_) {
    std::cout << "Error: inconsistent number of parameters: resolFunct has "<<resolParNum<<" parameters but "<<resolMCParVecSize<<" have been read from file" << std::endl;
    exit(1);
  }

  if( scaleParNum != scaleParVecSize ) {
    std::cout << "Error: inconsistent number of parameters: scaleFunct has "<<scaleParNum<<" parameters but "<<scaleParVecSize<<" have been read from file" << std::endl;
    exit(1);
  }

    
  std::vector<double>::const_iterator scaleParIt = scaleParVec_.begin();
  std::vector<double>::const_iterator resolDataParIt = resolDataParVec_.begin();
  std::vector<double>::const_iterator resolMCParIt = resolMCParVec_.begin();
    
  scaleParArray_ = new double[scaleParNum];
  resolDataParArray_ = new double[resolParNum];
  resolMCParArray_ = new double[resolParNum];
  
  for ( int iPar=0; iPar<scaleParNum; ++iPar) {
    double parameter = *scaleParIt;
    scaleParArray_[iPar] = parameter;
    ++scaleParIt;
  }

  if (useResol_){
    for ( int iPar=0; iPar<resolParNum; ++iPar) {
      double parameter = *resolDataParIt;
      resolDataParArray_[iPar] = parameter;
      ++resolDataParIt;
    }
    
    for ( int iPar=0; iPar<resolParNum; ++iPar) {
      double parameter = *resolMCParIt;
      resolMCParArray_[iPar] = parameter;
      ++resolMCParIt;
    }
  }
  
}


void MuScleFitCorrector::readParameters(const TString& fileName )
{
  scaleParArray_ = 0;
  resolDataParArray_ = 0;
  resolMCParArray_ = 0;

  // Read the parameters file
  ifstream parametersFile(fileName.Data());

  if( !parametersFile.is_open() ) {
    std::cout << "Error: file " << fileName << " not found. Aborting." << std::endl;
    abort();
  }

  std::string line;
  int nReadPar = 0;
  std::string iteration("Iteration ");
  // Loop on the file lines
  while (parametersFile) {

    getline( parametersFile, line );
    size_t lineInt = line.find("value");
    size_t iterationSubStr = line.find(iteration);

    if( iterationSubStr != std::string::npos ) {

      std::stringstream sLine(line);
      std::string num;
      int scaleFunctNum = 0;
      int resolFunctNum = 0;
      int wordCounter = 0;

      // Warning: this strongly depends on the parameters file structure.                                                                                   
      while( sLine >> num ) {
	++wordCounter;
	if( wordCounter == 8 ) {
	  std::stringstream in(num);
	  in >> resolFunctNum;
	}
	
	if( wordCounter == 9 ) {
	  std::stringstream in(num);
	  in >> scaleFunctNum;
	}
      }
            
      scaleFunctId_ = scaleFunctNum;
      scaleFunct_ = scaleFunctService(scaleFunctNum);
      resolFunctId_ = resolFunctNum;
      resolFunct_ = resolutionFunctService(resolFunctNum);

      std::cout<<"Function IDs: "<<std::endl;
      std::cout<<"     scale function number "<<scaleFunctId_<<std::endl;
      std::cout<<"     resolution function number "<<resolFunctId_<<std::endl;

  }
        
    int nScalePar = scaleFunct_->parNum();
    int nResolPar = resolFunct_->parNum();

    if ( (lineInt != std::string::npos) ) {
      size_t subStr1 = line.find("value");
      std::stringstream paramStr;
      double param = 0;
      paramStr << line.substr(subStr1+5);
      paramStr >> param;
      // Fill the last vector of parameters, which corresponds to this iteration.
      if (nReadPar<nScalePar) {
	scaleParVec_.push_back(param);
      }
      else if (nReadPar < (nScalePar+nResolPar) ) {
	resolDataParVec_.push_back(param);
      }
      else if (nReadPar < (nScalePar+2*nResolPar) ) {
	resolMCParVec_.push_back(param);
      }
      ++nReadPar;
    }
  }
  convertToArrays();

  std::cout<<"Scale function n. "<<scaleFunctId_<<" has "<<scaleFunct_->parNum()<<"parameters:"<<std::endl;
  for (int ii=0; ii<scaleFunct_->parNum(); ++ii){
    std::cout<<"par["<<ii<<"] = "<<scaleParArray_[ii]<<std::endl;
  }

  if (useResol_){

    std::cout<<"Resolution data function n. "<<resolFunctId_<<" has "<<resolFunct_->parNum()<<"parameters:"<<std::endl;
    for (int ii=0; ii<resolFunct_->parNum(); ++ii){
      std::cout<<"par["<<ii<<"] = "<<resolDataParArray_[ii]<<std::endl;
    }
    
    std::cout<<"Resolution MC function n. "<<resolFunctId_<<" has "<<resolFunct_->parNum()<<"parameters:"<<std::endl;
    for (int ii=0; ii<resolFunct_->parNum(); ++ii){
      std::cout<<"par["<<ii<<"] = "<<resolMCParArray_[ii]<<std::endl;
    }
    
  }
  
  
}
