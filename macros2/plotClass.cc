#include "plotClass.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"


#include "TF1.h"
#include "THStack.h"
#include "TKey.h"
#include "TTree.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVirtualPad.h"  // access to gPad
#include "TMinuit.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TEventList.h"

#include <iomanip>
#include <string>
#include <list>

using namespace std; 

ClassImp(plotClass)

// ----------------------------------------------------------------------
double fa_expo(double *x, double *par) {
  return par[0]*TMath::Exp(x[0]*par[1]);
} 

// ----------------------------------------------------------------------
double fa_err(double *x, double *par) {
  // from DK: TMath::Erf((a1-x)/a2))+a3
  return par[3]*(TMath::Erf((par[0]-x[0])/par[1])+par[2]); 
} 

// ----------------------------------------------------------------------
double fa_pol1(double *x, double *par) {
  return par[0] + par[1]*x[0]; 
}

// ----------------------------------------------------------------------
// expo and err and gauss 
double fa_expo_err(double *x, double *par) {
  //   par[0] = norm
  //   par[1] = exp
  //   par[2] = par[0] of err
  //   par[3] = par[1] of err
  //   par[4] = par[2] of err
  //   par[5] = par[3] of err
  return  (fa_err(x, &par[2]) + fa_expo(x, &par[0]));
}

// ----------------------------------------------------------------------
// expo and err and gauss 
double fa_pol1_err(double *x, double *par) {
  //   par[0] = const 
  //   par[1] = slope
  //   par[2] = par[0] of err
  //   par[3] = par[1] of err
  //   par[4] = par[2] of err
  //   par[5] = par[3] of err
  return  (fa_err(x, &par[2]) + fa_pol1(x, &par[0]));
}

// ----------------------------------------------------------------------
// This is a simplfied "landau" using an analytical formula.
// I use it because it was easiter for me to understand the integral (area) of this formula,
double fa_landausimp(double *x, double *par) {  // simplified landau
  double xx = ( x[0]- par[0] ) / par[1];
  return par[2] * 0.3989 * TMath::Exp(-0.5 * ( xx + TMath::Exp(-xx) ) );
} 

// ----------------------------------------------------------------------
plotClass::plotClass(const char *files, const char *dir, const char *cuts, int mode) { 

  fFiles = files; 

  TVirtualFitter::SetMaxIterations(20000);

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  fEpsilon = 0.00001; 

  fAccPt = 3.5;
  fAccEtaGen = 2.5;
  fAccEtaRec = 2.4;

  fRunMin = -1; 
  fRunMax = -2; 

  // PDG 2012:
  fu  = 0.401;
  fs  = 0.105;
  // -- CMS with PDG input
  fsfu = 0.282;
  fsfuE = 0.037/0.282;
  // -- CMS with LHCb input
  fsfu = 0.267;
  fsfuE = 0.021;

  fYear = 0; 

  fDoUseBDT = false; 
  fDoApplyMuonPtCuts = false; 
  fDoPrintSingleEvent = false;

  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  fInvertedIso = false; 
  fNormProcessed = false; 
  fSaveSmallTree = true;
  fSaveLargerTree = false ;

  fCutsFileName = cuts; 
  init(files, cuts, dir, mode);

  int NBINS = static_cast<int>((fMassHi - fMassLo)/0.025);

  TH1D *h; 
  TH2D *h2; 

  numbers *a(0); 
  cout << "----> fNchan = " << fNchan << " fCuts.size() = " << fCuts.size() << endl;
  for (unsigned int i = 0; i < fNchan; ++i) {
    cout << "--> set up structures for channel #" << i << endl;
    cout << "    etaMin: " << fCuts[i]->etaMin << " etaMax: " << fCuts[i]->etaMax << endl;
    // -- signal Bs2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bs2MuMu %i", i); 
    fNumbersBs.push_back(a); 

    // -- signal Bd2MuMu
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("signal Bd2MuMu %i", i); 
    fNumbersBd.push_back(a); 

    // --  normalization
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name = Form("normalization %i", i); 
    fNumbersNo.push_back(a); 

    // --  control sample
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("control sample %i", i); 
    fNumbersCs.push_back(a); 

    // -- temporary numbers (e.g. trigger efficiency)
    a = new numbers;
    initNumbers(a); 
    a->index = i; 
    a->name  = Form("Bla %i", i); 
    fNumbersBla.push_back(a); 
    
    h2 = new TH2D(Form("hAccAll%d", i), Form("hAccAll%d", i), 25, 0., 2.5, 25, 0., 50.);
    fhAccAll.push_back(h2); 
    h2 = new TH2D(Form("hAccPass%d", i), Form("hAccPass%d", i), 25, 0., 2.5, 25, 0., 50.);
    fhAccPass.push_back(h2); 

    h = new TH1D(Form("hAccPtAll%d", i), Form("hAccPtAll%d", i), 25, 0., 50.);
    fhAccPtAll.push_back(h); 
    h = new TH1D(Form("hAccPtPass%d", i), Form("hAccPtPass%d", i), 25, 0., 50.);
    fhAccPtPass.push_back(h); 

    h = new TH1D(Form("hAccEtaAll%d", i), Form("hAccEtaAll%d", i), 25, 0., 2.5);
    fhAccEtaAll.push_back(h); 
    h = new TH1D(Form("hAccEtaPass%d", i), Form("hAccEtaPass%d", i), 25, 0., 2.5);
    fhAccEtaPass.push_back(h); 


    
    h2 = new TH2D(Form("hBdtMass%d", i), Form("hBdtMass%d", i), NBINS, fMassLo, fMassHi, 100, 0., 1.0);
    fhBdtMass.push_back(h2); 

    h = new TH1D(Form("hGenAndAccNumbers%d", i), Form("hGenAndAccNumbers%d", i), 100, 0., 100.);
    fhGenAndAccNumbers.push_back(h);

    h = new TH1D(Form("hMassWithMassCuts%d", i), Form("hMassWithMassCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMassCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithMassCutsManyBins%d", i), Form("fhMassWithMassCutsManyBins%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhMassWithMassCutsManyBins.push_back(h); 
  
    h = new TH1D(Form("hMassWithAllCuts%d", i), Form("hMassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCuts.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsBlind%d", i), Form("hMassWithAllCutsBlind%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCutsBlind.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsManyBins%d", i), Form("hMassWithAllCutsManyBins%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhMassWithAllCutsManyBins.push_back(h); 


    h = new TH1D(Form("hW8MassWithAllCuts%d", i), Form("hW8MassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhW8MassWithAllCuts.push_back(h); 
    h = new TH1D(Form("hW8MassWithAllCutsBlind%d", i), Form("hW8MassWithAllCutsBlind%d", i), NBINS, fMassLo, fMassHi);
    fhW8MassWithAllCutsBlind.push_back(h); 
    h = new TH1D(Form("hW8MassWithAllCutsManyBins%d", i), Form("hW8MassWithAllCutsManyBins%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhW8MassWithAllCutsManyBins.push_back(h); 


    h = new TH1D(Form("hNorm%d", i), Form("hNorm%d", i), 100, 4.9, 5.9);
    fhNorm.push_back(h); 

    h = new TH1D(Form("hNormC%d", i), Form("hNormC%d", i), 200, 4.9, 5.9);
    fhNormC.push_back(h); 

    h = new TH1D(Form("hDstarPi%d", i), Form("hDstarPi%d", i), 20, 4.9, 5.9);
    fhDstarPi.push_back(h); 

    h = new TH1D(Form("hMassWithTriggerCuts%d", i), Form("hMassWithTriggerCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithTriggerCuts.push_back(h); 
    h = new TH1D(Form("hMassWithTriggerCutsManyBins%d", i), Form("hMassWithTriggerCutsManyBins%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhMassWithTriggerCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassWithMuonCuts%d", i), Form("hMassWithMuonCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMuonCuts.push_back(h); 
    h = new TH1D(Form("hMassWithMuonCutsManyBins%d", i), Form("hMassWithMuonCutsManyBins%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhMassWithMuonCutsManyBins.push_back(h); 

    h = new TH1D(Form("fhMassWithAnaCuts%d", i), Form("hMassChan%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAnaCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithAnaCutsManyBins%d", i), Form("hMassChan%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhMassWithAnaCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h); 
    h = new TH1D(Form("fhMassNoCutsManyBins%d", i), Form("hMassChan%d", i), 
		 static_cast<int>((fMassHi-fMassLo)*1000), fMassLo, fMassHi);
    fhMassNoCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassAbsNoCuts%d", i), Form("hMassAbsNoCuts%d", i), 100, 0, 10);
    fhMassAbsNoCuts.push_back(h); 

    h = new TH1D(Form("hMuId%d", i), Form("hMuId%d", i), 100, 0., 1.);
    fhMuId.push_back(h); 
    h = new TH1D(Form("hMuTr%d", i), Form("hMuTr%d", i), 100, 0., 1.);
    fhMuTr.push_back(h); 

    h = new TH1D(Form("hMuIdMC%d", i), Form("hMuIdMC%d", i), 100, 0., 1.);
    fhMuIdMC.push_back(h); 
    h = new TH1D(Form("hMuTrMC%d", i), Form("hMuTrMC%d", i), 100, 0., 1.);
    fhMuTrMC.push_back(h); 

  }

  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT); 
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID); 
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt); 
  fAnaCuts.addCut("fGoodMuonsEta", "#eta_{#mu}", fGoodMuonsEta); 
  fAnaCuts.addCut("fGoodTracks", "good tracks", fGoodTracks); 
  fAnaCuts.addCut("fGoodTracksPt", "p_{T,trk} [GeV]", fGoodTracksPt); 
  fAnaCuts.addCut("fGoodTracksEta", "#eta_{trk} ", fGoodTracksEta); 

  fAnaCuts.addCut("fGoodQ", "q_{1} 1_{2}", fGoodQ);   
  fAnaCuts.addCut("fGoodPvAveW8", "<w8>", fGoodPvAveW8); 
  fAnaCuts.addCut("fGoodIp", "IP", fGoodIp); 
  fAnaCuts.addCut("fGoodIpS", "IPS", fGoodIpS); 
  fAnaCuts.addCut("fGoodLip", "LIP", fGoodLip); 
  fAnaCuts.addCut("fGoodLipS", "LIPS", fGoodLipS); 
  fAnaCuts.addCut("fGoodMaxDoca", "MAXDOCA", fGoodMaxDoca); 
  fAnaCuts.addCut("fGoodPt", "p_{T,B}", fGoodPt); 
  fAnaCuts.addCut("fGoodEta", "#eta_{B}", fGoodEta); 
  fAnaCuts.addCut("fGoodAlpha", "#alpha", fGoodAlpha); 
  fAnaCuts.addCut("fGoodFLS", "l/#sigma(l)", fGoodFLS); 
  fAnaCuts.addCut("fGoodChi2", "#chi^{2}", fGoodChi2); 
  fAnaCuts.addCut("fGoodIso", "I_{trk}", fGoodIso); 
  fAnaCuts.addCut("fGoodCloseTrack", "close track veto", fGoodCloseTrack); 
  fAnaCuts.addCut("fGoodDocaTrk", "d_{ca}(trk)", fGoodDocaTrk); 
  fAnaCuts.addCut("fGoodLastCut", "lastCut", fGoodLastCut); 


  if (0) {
    fptFakePosKaons     = new PidTable("../macros/pidtables/130409/kaonFakeRate-mvaTmMuon.dat"); 
    fptFakeNegKaons     = new PidTable("../macros/pidtables/130409/kaonFakeRate-mvaTmMuon.dat"); 
    
    fptFakePosPions     = new PidTable("../macros/pidtables/130409/pionFakeRate-mvaTmMuon.dat"); 
    fptFakeNegPions     = new PidTable("../macros/pidtables/130409/pionFakeRate-mvaTmMuon.dat"); 
    
    fptFakePosProtons   = new PidTable("../macros/pidtables/130409/protonFakeRate-mvaTmMuon.dat"); 
    fptFakeNegProtons   = new PidTable("../macros/pidtables/130409/protonFakeRate-mvaTmMuon.dat"); 
  }

  if (0) {
    fptFakePosKaons     = new PidTable("../macros/pidtables/130410/kaonPlusFakeRate-mvaTmMuon.dat"); 
    fptFakeNegKaons     = new PidTable("../macros/pidtables/130410/kaonMinusFakeRate-mvaTmMuon.dat"); 
    
    fptFakePosPions     = new PidTable("../macros/pidtables/130410/pionPlusFakeRate-mvaTmMuon.dat"); 
    fptFakeNegPions     = new PidTable("../macros/pidtables/130410/pionMinusFakeRate-mvaTmMuon.dat"); 
    
    fptFakePosProtons   = new PidTable("../macros/pidtables/130410/protonPlusFakeRate-mvaTmMuon.dat"); 
    fptFakeNegProtons   = new PidTable("../macros/pidtables/130410/protonMinusFakeRate-mvaTmMuon.dat"); 
  }


  if (0) {
    fptFakePosKaons     = new PidTable("../macros/pidtables/130409/flat-kaon.dat"); 
    fptFakeNegKaons     = new PidTable("../macros/pidtables/130409/flat-kaon.dat"); 
    
    fptFakePosPions     = new PidTable("../macros/pidtables/130409/flat-pion.dat"); 
    fptFakeNegPions     = new PidTable("../macros/pidtables/130409/flat-pion.dat"); 
    
    fptFakePosProtons   = new PidTable("../macros/pidtables/130409/flat-proton.dat"); 
    fptFakeNegProtons   = new PidTable("../macros/pidtables/130409/flat-proton.dat"); 
  }

  if (0) {
    fptFakePosKaons     = new PidTable("../macros/pidtables/130419/kaonPlusFakeRate-mvaMuon.dat"); 
    fptFakeNegKaons     = new PidTable("../macros/pidtables/130419/kaonMinusFakeRate-mvaMuon.dat"); 
    
    fptFakePosPions     = new PidTable("../macros/pidtables/130419/pionPlusFakeRate-mvaMuon.dat"); 
    fptFakeNegPions     = new PidTable("../macros/pidtables/130419/pionMinusFakeRate-mvaMuon.dat"); 
    
    fptFakePosProtons   = new PidTable("../macros/pidtables/130419/protonPlusFakeRate-mvaMuon.dat"); 
    fptFakeNegProtons   = new PidTable("../macros/pidtables/130419/protonMinusFakeRate-mvaMuon.dat"); 
  }

  if (1) {
    fptFakePosKaons     = new PidTable("../macros/pidtables/130510/2012-kaonPosFakeRate-mvaMuon.dat"); 
    fptFakeNegKaons     = new PidTable("../macros/pidtables/130510/2012-kaonNegFakeRate-mvaMuon.dat"); 
    
    fptFakePosPions     = new PidTable("../macros/pidtables/130510/2012-pionPosFakeRate-mvaMuon.dat"); 
    fptFakeNegPions     = new PidTable("../macros/pidtables/130510/2012-pionNegFakeRate-mvaMuon.dat"); 
    
    fptFakePosProtons   = new PidTable("../macros/pidtables/130510/2012-protonPosFakeRate-mvaMuon.dat"); 
    fptFakeNegProtons   = new PidTable("../macros/pidtables/130510/2012-protonNegFakeRate-mvaMuon.dat"); 
  }

  if (1) { // if (2012 == fYear) {
    string name(""); 
    string directory("../macros/pidtables/130510/"); 
    name = directory + Form("%d-L1L2_data_all.dat", fYear);        fptT1     = new PidTable(name.c_str());
    name = directory + Form("%d-L3_data_all.dat", fYear);          fptT2     = new PidTable(name.c_str());
    name = directory + Form("%d-MuonID_data_all.dat", fYear);      fptM      = new PidTable(name.c_str());
    
    name = directory + Form("%d-L1L2_mc_all.dat", fYear);          fptT1MC   = new PidTable(name.c_str()); 	
    name = directory + Form("%d-L3_mc_all.dat", fYear);            fptT2MC   = new PidTable(name.c_str()); 	
    name = directory + Form("%d-MuonID_mc_all.dat", fYear);        fptMMC    = new PidTable(name.c_str());  

    name = directory + Form("%d-L1L2_data_seagulls.dat", fYear);   fptSgT1   = new PidTable(name.c_str());  
    name = directory + Form("%d-L3_data_seagulls.dat", fYear);     fptSgT2   = new PidTable(name.c_str());    
    name = directory + Form("%d-MuonID_data_seagulls.dat", fYear); fptSgM    = new PidTable(name.c_str());
    
    name = directory + Form("%d-L1L2_mc_seagulls.dat", fYear);     fptSgT1MC = new PidTable(name.c_str());    
    name = directory + Form("%d-L3_mc_seagulls.dat", fYear);       fptSgT2MC = new PidTable(name.c_str());      
    name = directory + Form("%d-MuonID_mc_seagulls.dat", fYear);   fptSgMMC  = new PidTable(name.c_str());  

    name = directory + Form("%d-L1L2_data_cowboys.dat", fYear);    fptCbT1   = new PidTable(name.c_str());  
    name = directory + Form("%d-L3_data_cowboys.dat", fYear);      fptCbT2   = new PidTable(name.c_str());    
    name = directory + Form("%d-MuonID_data_cowboys.dat", fYear);  fptCbM    = new PidTable(name.c_str());
    
    name = directory + Form("%d-L1L2_mc_cowboys.dat", fYear);      fptCbT1MC = new PidTable(name.c_str());    
    name = directory + Form("%d-L3_mc_cowboys.dat", fYear);        fptCbT2MC = new PidTable(name.c_str());      
    name = directory + Form("%d-MuonID_mc_cowboys.dat", fYear);    fptCbMMC  = new PidTable(name.c_str());  
  }


  fNumbersFileName = fDirectory + "/anaBmm.plotClass." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);
  
  dumpSamples();
  fF["SgMc"]->cd();
  dumpCutNames("candAnaMuMu/hcuts");
  fTEX.close();

}

// ----------------------------------------------------------------------
plotClass::~plotClass() {
  cout << "plotClass dtor" << endl;
  for (map<string, TFile*>::iterator imap = fF.begin(); imap != fF.end(); ++imap) {  
    cout << "    => closing " << imap->first;
    cout << ": " << fF[imap->first]->GetName() << endl;
    imap->second->Close();
  }
}


// ----------------------------------------------------------------------
void plotClass::initNumbers(numbers *a, bool initAll) {

  if (initAll) a->name = "";
  a->index          = -1; 
  a->effPtReco = 0.; a->effPtRecoE = 0.;
  a->fitYield = 0.; a->fitYieldE = 0.;
  a->fitYieldC = 0.; a->fitYieldCE = 0.;
  a->genAccFileYield = 0.; a->genYield = 0.; a->genFileYield = 0.; a->genAccYield = 0.;
  a->recoYield = 0.; a->muidYield = 0.; a->trigYield = 0.; a->chanYield = 0.; a->candYield = 0.;
  a->genFileYieldE = 0.; a->genYieldE = 0.; a->recoYieldE = 0.; a->chanYieldE = 0.; 
  a->muidYieldE = 0.; a->trigYieldE = 0.; a->candYieldE = 0.;
  a->absNoCutsYield = 0.; a->ana0Yield = 0.;  a->anaYield = 0.;  a->anaMuonYield = 0.;  a->anaTriggerYield = 0.;  a->anaWmcYield = 0.;
  a->ana0YieldE = 0.; a->anaYieldE = 0.; a->anaMuonYieldE = 0.; a->anaTriggerYieldE = 0.; a->anaWmcYieldE = 0.;
  a->acc = 0.; a->accE = 0.;
  a->effChan = 0.; a->effChanE = 0.; 
  a->effMuidMC = 0.; a->effMuidMCE = 0.; a->effTrigMC = 0.; a->effTrigMCE = 0.;
  a->effMuidTNP = 0.; a->effMuidTNPE = 0.; a->effTrigTNP = 0.; a->effTrigTNPE = 0.;
  a->effMuidTNPMC = 0.; a->effMuidTNPMCE = 0.; a->effTrigTNPMC = 0.; a->effTrigTNPMCE = 0.;
  a->effRhoMuidTNPMC = 0.; a->effRhoMuidTNPMCE = 0.; a->effRhoTrigTNPMC = 0.; a->effRhoTrigTNPMCE = 0.;
  a->effRhoMuidTNP = 0.; a->effRhoMuidTNPE = 0.; a->effRhoTrigTNP = 0.; a->effRhoTrigTNPE = 0.;
  a->rhoMuidTNPMC = 0.; a->rhoMuidTNPMCE = 0.; a->rhoTrigTNPMC = 0.; a->rhoTrigTNPMCE = 0.;
  a->effCand = 0.; a->effCandE = 0.; a->effCandTE = 0.; 
  a->effAna = 0.; a->effAnaE = 0.; 
  a->effTot = 0.; a->effTotE = 0.; a->aEffProdMC = 0.; a->aEffProdMCE = 0.; a->effProdMC = 0.; a->effProdMCE = 0.; 
  a->effProdTNP = 0.; a->effProdTNPE = 0.; a->effProdTNPMC = 0.; a->effProdTNPMCE = 0.; 
  a->effTotChan = 0.; a->effTotChanE = 0.; 
  a->prodGenYield = 0.; a->combGenYield = 0.; a->chanGenYield = 0.; 
  a->expSignal = 0.;
  a->bgObs = 0.; a->bgBsExp = 0.; a->bgBsExpE = 0.; a->bgBsExpTE = 0.; a->bgBdExp = 0.; a->bgBdExpE = 0.; a->bgBdExpTE = 0.;
  a->bsObs = 0.; a->bdObs = 0.; 
  a->bsExpObs = 0.; a->bdExpObs = 0.; a->bsExpObsE = 0.; a->bdExpObsE = 0.; a->bsExpObsTE = 0.; a->bdExpObsTE = 0.; 
  a->tauBs = 0.; a->tauBsE = 0.; a->tauBd = 0.; a->tauBdE = 0.; 
  a->offHiRare = 0.; a->offHiRareE = 0.; a->offLoRare = 0.; a->offLoRareE = 0.; 
  a->bsRare = 0.; a->bsRareE = 0.; a->bdRare = 0.; a->bdRareE = 0.; 
  a->offHi = 0.; a->offLo = 0.;
  a->bsNoScaled = 0.; a->bsNoScaledE = 0.; a->bdNoScaled = 0.; a->bdNoScaledE = 0.;
  a->pss = 0.; a->pssE = 0.; a->pssTE = 0.; a->pdd = 0.; a->pddE = 0.; a->pddTE = 0.;
  a->pls = 0.; a->plsE = 0.; a->plsTE = 0.; a->pld = 0.; a->pldE = 0.; a->pldTE = 0.;
  a->phs = 0.; a->phsE = 0.; a->phsTE = 0.; a->phd = 0.; a->phdE = 0.; a->phdTE = 0.;
  a->psd = 0.; a->psdE = 0.; a->psdTE = 0.; a->pds = 0.; a->pdsE = 0.; a->pdsTE = 0.;
  a->mBdLo = 0.; a->mBdHi = 0.; a->mBsLo = 0.; a->mBsHi = 0.;
  a->accTE = 0.; a->effAnaTE = 0.; a->effTotTE = 0.; 
  a->tauTE = 0.; 
  a->fitYieldTE = 0.;
  a->effMuidTE = 0.; a->effMuidTNPTE = 0.; a->effMuidMCTE = 0.; a->effMuidTNPMCTE = 0.; a->effTrigTE = 0.; a->effTrigMCTE = 0.;
  a->effTrigTNPTE = 0.; a->effTrigTNPMCTE = 0.; a->effRhoMuidTNPTE = 0.; 
  a->effRhoMuidTNPMCTE = 0.; a->effRhoTrigTNPTE = 0.; a->effRhoTrigTNPMCTE = 0.; 
  a->fBgPeakLo = 0.;   a->fBgPeakHi = 0.;   a->fBgPeakBs = 0.;   a->fBgPeakBd = 0.;    
  a->fBgPeakLoE1 = 0.; a->fBgPeakHiE1 = 0.; a->fBgPeakBsE1 = 0.; a->fBgPeakBdE1 = 0.;  
  a->fBgPeakLoE2 = 0.; a->fBgPeakHiE2 = 0.; a->fBgPeakBsE2 = 0.; a->fBgPeakBdE2 = 0.;  
  a->fBgRslLo = 0.;    a->fBgRslHi = 0.;    a->fBgRslBs = 0.;    a->fBgRslBd = 0.;     
  a->fBgRslLoE1 = 0.;  a->fBgRslHiE1 = 0.;  a->fBgRslBsE1 = 0.;  a->fBgRslBdE1 = 0.;   
  a->fBgRslLoE2 = 0.;  a->fBgRslHiE2 = 0.;  a->fBgRslBsE2 = 0.;  a->fBgRslBdE2 = 0.;   
  a->fBgRareLo = 0.;   a->fBgRareHi = 0.;   a->fBgRareBs = 0.;   a->fBgRareBd = 0.;    
  a->fBgRareLoE1 = 0.; a->fBgRareHiE1 = 0.; a->fBgRareBsE1 = 0.; a->fBgRareBdE1 = 0.;  
  a->fBgRareLoE2 = 0.; a->fBgRareHiE2 = 0.; a->fBgRareBsE2 = 0.; a->fBgRareBdE2 = 0.;  
  a->fBgRslsLo = 0.;   a->fBgRslsHi = 0.;   a->fBgRslsBs = 0.;   a->fBgRslsBd = 0.;    
  a->fBgRslsLoE1 = 0.; a->fBgRslsHiE1 = 0.; a->fBgRslsBsE1 = 0.; a->fBgRslsBdE1 = 0.;  
  a->fBgRslsLoE2 = 0.; a->fBgRslsHiE2 = 0.; a->fBgRslsBsE2 = 0.; a->fBgRslsBdE2 = 0.;  
  a->fBgCombLo = 0.;   a->fBgCombHi = 0.;   a->fBgCombBs = 0.;   a->fBgCombBd = 0.;    
  a->fBgCombLoE1 = 0.; a->fBgCombHiE1 = 0.; a->fBgCombBsE1 = 0.; a->fBgCombBdE1 = 0.;                            
  a->fBgCombLoE2 = 0.; a->fBgCombHiE2 = 0.; a->fBgCombBsE2 = 0.; a->fBgCombBdE2 = 0.;                            
  a->fBgNonpLo = 0.;   a->fBgNonpHi = 0.;   a->fBgNonpBs = 0.;   a->fBgNonpBd = 0.;    
  a->fBgNonpLoE1 = 0.; a->fBgNonpHiE1 = 0.; a->fBgNonpBsE1 = 0.; a->fBgNonpBdE1 = 0.;                            
  a->fBgNonpLoE2 = 0.; a->fBgNonpHiE2 = 0.; a->fBgNonpBsE2 = 0.; a->fBgNonpBdE2 = 0.;                            
  a->fBgTotLo = 0.;    a->fBgTotHi = 0.;    a->fBgTotBs = 0.;    a->fBgTotBd = 0.;    
  a->fBgTotLoE1 = 0.;  a->fBgTotHiE1 = 0.;  a->fBgTotBsE1 = 0.;  a->fBgTotBdE1 = 0.;                             
  a->fBgTotLoE2 = 0.;  a->fBgTotHiE2 = 0.;  a->fBgTotBsE2 = 0.;  a->fBgTotBdE2 = 0.;                             
  a->fSgLo = 0.;       a->fSgHi = 0.;       a->fSgBs = 0.;       a->fSgBd = 0.;       
  a->fSgLoE1 = 0.;     a->fSgHiE1 = 0.;     a->fSgBsE1 = 0.;     a->fSgBdE1 = 0.;                                
  a->fSgLoE2 = 0.;     a->fSgHiE2 = 0.;     a->fSgBsE2 = 0.;     a->fSgBdE2 = 0.;                                
  a->fBdLo = 0.;       a->fBdHi = 0.;       a->fBdBs = 0.;       a->fBdBd = 0.;       
  a->fBdLoE1 = 0.;     a->fBdHiE1 = 0.;     a->fBdBsE1 = 0.;     a->fBdBdE1 = 0.;                                
  a->fBdLoE2 = 0.;     a->fBdHiE2 = 0.;     a->fBdBsE2 = 0.;     a->fBdBdE2 = 0.;                                
  a->fSgTot = 0.;      a->fBdTot = 0.;      a->fNoTot = 0.;      a->fCsTot = 0.;      
  a->fSgTotE1 = 0.;    a->fBdTotE1 = 0.;    a->fNoTotE1 = 0.;    a->fCsTotE1 = 0.; 
  a->fSgTotE2 = 0.;    a->fBdTotE2 = 0.;    a->fNoTotE2 = 0.;    a->fCsTotE2 = 0.; 
  a->fFitSg = 0.;      a->fFitBd = 0.;      a->fFitNo = 0.;      a->fFitNoC = 0.;     a->fFitCs = 0.;    a->fFitCsC = 0.; 
  a->fFitSgE1 = 0.;    a->fFitBdE1 = 0.;    a->fFitNoE1 = 0.;    a->fFitNoCE1 = 0.;   a->fFitCsE1 = 0.;  a->fFitCsCE1 = 0.;
  a->fFitSgE2 = 0.;    a->fFitBdE2 = 0.;    a->fFitNoE2 = 0.;    a->fFitNoCE2 = 0.;   a->fFitCsE2 = 0.;  a->fFitCsCE2 = 0.;
  a->fSgAndBgLo = 0.;  a->fSgAndBgHi = 0.;  a->fSgAndBgBs = 0.;  a->fSgAndBgBd = 0.;  
  a->fSgAndBgLoE1 = 0.;a->fSgAndBgHiE1 = 0.;a->fSgAndBgBsE1 = 0.;a->fSgAndBgBdE1 = 0.;                                         
  a->fSgAndBgLoE2 = 0.;a->fSgAndBgHiE2 = 0.;a->fSgAndBgBsE2 = 0.;a->fSgAndBgBdE2 = 0.;                                         
  a->fObsLo = 0.;      a->fObsHi = 0.;      a->fObsBs = 0.;      a->fObsBd = 0.;      

}


// ----------------------------------------------------------------------
// call this on the acceptance files
void plotClass::filterEfficiency(string fname, string name) {
 
  TFile *f = fF[fname];
  if (0 == f) {
    cout << "anaBmm::filterEffciciency(" << name << "): no file " << fname << " found " << endl;
    return;
  }
  TTree *t  = (TTree*)(f->Get(Form("%s/effTree", name.c_str())));

  bool sg(false), no(false), cs(false); 

  float bg1pt, bg2pt, bg1eta, bg2eta;
  float bg3pt, bg4pt, bg3eta, bg4eta; 

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);

  if (string::npos != name.find("MuMu")) {
    cout << "anaBmm::filterEfficiency(" << name << "): SIGNAL " << endl;
    sg = true; 
  }

  if (string::npos != name.find("Bu2JpsiK")) {
    cout << "anaBmm::filterEfficiency(" << name << "): NORMALIZATION " << endl;
    no = true; 
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
  }

  if (string::npos != name.find("Bs2JpsiPhi")) {
    cout << "anaBmm::filterEfficiency(" << name << "): CONTROL SAMPLE " << endl;
    cs = true; 
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("g4pt", &bg4pt);
    t->SetBranchAddress("g4eta",&bg4eta);
  }

  int ngen(0), ngenlevel(0); 
  int nentries = Int_t(t->GetEntries());
  for (int jentry = 0; jentry < nentries; jentry++) {
    t->GetEntry(jentry);

    ++ngen;

    if (sg) {
      if ((TMath::Abs(bg1eta) < 2.5) && (TMath::Abs(bg2eta) < 2.5) 
	  && (bg1pt > 3.5) && (bg2pt > 3.5)
	  ) {
	++ngenlevel;
      }
    }

    if (no) {
      if ((TMath::Abs(bg1eta) < 2.5) && (TMath::Abs(bg2eta) < 2.5) && (TMath::Abs(bg3eta) < 2.5) 
	  && (bg1pt > 3.5) && (bg2pt > 3.5) && (bg3pt > 0.4)
	  ) {
	++ngenlevel;
      }
    }


    if (cs) {
      if ((TMath::Abs(bg1eta) < 2.5) && (TMath::Abs(bg2eta) < 2.5) && (TMath::Abs(bg3eta) < 2.5) && (TMath::Abs(bg4eta) < 2.5) 
	  && (bg1pt > 3.5) && (bg2pt > 3.5) && (bg3pt > 0.4) && (bg4pt > 0.4)
	  ) {
	++ngenlevel;
      }
    }
  }

  cout << "======================================================================" << endl;
  cout << "Filter efficiency for " << fname << " and " << name << endl;
  cout << " => " 
       << ngenlevel 
       << "/" 
       << ngen
       << " = " 
       << static_cast<double>(ngenlevel)/static_cast<double>(ngen) 
       << endl;
  cout << "======================================================================" << endl;

}


// ----------------------------------------------------------------------
void plotClass::accEffFromEffTreeBac(string fname, string dname, numbers &a, cuts &b, int proc) {

  //   TFile *f = fF[fname];
  //   if (0 == f) {
  //     cout << "anaBmm::accEffFromEffTreeBac(" << a.name << "): no file " << fname << " found " << endl;
  //     return;
  //   }
  //   TTree *t  = (TTree*)(f->Get(Form("%s/effTree", dname.c_str())));
  //   double effFilter(1.); 
  //   if (!t) {
  //     cout << "anaBmm::accEffFromEffTreeBac(" << a.name << "): no tree `effTree' found in " 
  // 	 << f->GetName() << " and dir = " << Form("%s/effTree", dname.c_str()) 
  // 	 << endl;
    
  //     return;
  //   } else {
  //     effFilter = fFilterEff[fname]; 
  //     cout << "anaBmm::accEffFromEffTreeBac(" << a.name << ")" << endl
  // 	 << " get acceptance from file " << f->GetName() << " and dir = " << Form("%s/effTree", dname.c_str()) 
  // 	 << " with filterEff = " << effFilter
  // 	 << endl;
  //   }

  //   bool sg(false), no(false), cs(false); 

  //   int   bprocid, bidx; 
  //   bool  bhlt;
  //   float bg1pt, bg2pt, bg1eta, bg2eta;
  //   float bm1pt, bm1eta, bm2pt, bm2eta;
  //   bool  bm1gt, bm2gt; 
  //   bool  bm1id, bm2id; 
  //   float bm;

  //   float bg3pt, bg4pt, bg3eta, bg4eta; 
  //   float bk1pt, bk2pt, bk1eta, bk2eta; 
  //   bool  bk1gt, bk2gt; 

  //   t->SetBranchAddress("hlt",&bhlt);
  //   t->SetBranchAddress("procid",&bprocid);
  //   t->SetBranchAddress("bidx",&bidx);

  //   t->SetBranchAddress("g1pt",&bg1pt);
  //   t->SetBranchAddress("g2pt",&bg2pt);
  //   t->SetBranchAddress("g1eta",&bg1eta);
  //   t->SetBranchAddress("g2eta",&bg2eta);

  //   t->SetBranchAddress("m1pt",&bm1pt);
  //   t->SetBranchAddress("m2pt",&bm2pt);
  //   t->SetBranchAddress("m1eta",&bm1eta);
  //   t->SetBranchAddress("m2eta",&bm2eta);

  //   t->SetBranchAddress("m1gt",&bm1gt);
  //   t->SetBranchAddress("m2gt",&bm2gt);
  //   t->SetBranchAddress("m1id",&bm1id);
  //   t->SetBranchAddress("m2id",&bm2id);


  //   if (string::npos != a.name.find("signal")) {
  //     cout << "anaBmm::accEffFromEffTreeBac(" << a.name << "): SIGNAL " << endl;
  //     sg = true; 
  //   }

  //   if (string::npos != a.name.find("normalization")) {
  //     cout << "anaBmm::accEffFromEffTreeBac(" << a.name << "): NORMALIZATION " << endl;
  //     no = true; 
  //     t->SetBranchAddress("g3pt", &bg3pt);
  //     t->SetBranchAddress("g3eta",&bg3eta);
  //     t->SetBranchAddress("k1pt", &bk1pt);
  //     t->SetBranchAddress("k1eta",&bk1eta);
  //     t->SetBranchAddress("k1gt", &bk1gt);
  //   }

  //   if (string::npos != a.name.find("control sample")) {
  //     cout << "anaBmm::accEffFromEffTreeBac(" << a.name << "): CONTROL SAMPLE " << endl;
  //     cs = true; 
  //     t->SetBranchAddress("g3pt", &bg3pt);
  //     t->SetBranchAddress("g3eta",&bg3eta);
  //     t->SetBranchAddress("k1pt", &bk1pt);
  //     t->SetBranchAddress("k1eta",&bk1eta);
  //     t->SetBranchAddress("k1gt", &bk1gt);

  //     t->SetBranchAddress("g4pt", &bg4pt);
  //     t->SetBranchAddress("g4eta",&bg4eta);
  //     t->SetBranchAddress("k2pt", &bk2pt);
  //     t->SetBranchAddress("k2eta",&bk2eta);
  //     t->SetBranchAddress("k2gt", &bk2gt);
  //   }

  //   t->SetBranchAddress("m",&bm);


  //   int nentries = Int_t(t->GetEntries());
  //   int ngen(0), nchangen(0), nreco(0), nchan(0), nmuid(0), nhlt(0), ncand(0), ncand2(0); 
  //   int chan(-1); 
  //   int recoPtA(0), recoPtB(0); 
  //   cout << "channel = " << a.index << endl;
  //   for (int jentry = 0; jentry < nentries; jentry++) {
  //     t->GetEntry(jentry);
  //     if (bidx < 0) continue;
  //     ++ngen;
  //     if (proc > 0 && bprocid != proc) continue;
  //     if (sg) {
  //       // -- Signal
  //       chan = detChan(bg1eta, bg2eta); 
  //       if (chan == a.index) {
  // 	++nchangen;
  // 	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5) {
  // 	  if (bg1pt > 1. && bg2pt > 1.) {
  // 	    if (bm1pt > 1. && bm2pt > 1.
  // 		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4
  // 		&& bm1gt && bm2gt
  // 		) {
  // 	      chan = detChan(bm1eta, bm2eta); 
  // 	      if (chan == a.index) {
  // 		++nreco;
  // 		if (bm > 0) {
  // 		  ++ncand2;
  // 		}
  // 		//if (bm1pt > b.m1pt && bm2pt > b.m2pt) {
  // 		if (bm1pt > 3.5 && bm2pt > 3.5) {
  // 		  ++nchan; 
  // 		  if (bm > 0) {
  // 		    ++ncand;
  // 		  }
  // 		  if (bm1id && bm2id) {
  // 		    ++nmuid;
  // 		    if (bhlt) {
  // 		      ++nhlt;
  // 		    }
  // 		  }
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //       }
  //     } else if (no) {
  //       // -- Normalization
  //       chan = detChan(bg1eta, bg2eta); 
  //       if (chan == a.index) {
  // 	++nchangen;
  // 	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5) {
  // 	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4) {
  // 	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.5
  // 		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4
  // 		&& bm1gt && bm2gt && bk1gt
  // 		) {
  // 	      ++recoPtA; 
  // 	      if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
  // 	      chan = detChan(bm1eta, bm2eta); 
  // 	      if (chan == a.index) {
  // 		++nreco;
  // 		if (bm1pt > b.m1pt && bm2pt > b.m2pt
  // 		    ) {
  // 		  ++nchan; 
  // 		  if (bm > 0) {
  // 		    ++ncand;
  // 		  }
  // 		  if (bm1id && bm2id) {
  // 		    ++nmuid;
  // 		    if (bhlt) {
  // 		      ++nhlt;
  // 		    }
  // 		  }
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //       }
  //     } else if (cs) {
  //       // -- Control sample
  //       chan = detChan(bg1eta, bg2eta); 
  //       if (chan == a.index) {
  // 	++nchangen;
  // 	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5 && TMath::Abs(bg4eta) < 2.5) {
  // 	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4 && bg4pt > 0.4) {
  // 	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.5 && bk2pt > 0.5
  // 		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4 && TMath::Abs(bk2eta) < 2.4
  // 		&& bm1gt && bm2gt && bk1gt && bk2gt
  // 		) {
  // 	      ++recoPtA; 
  // 	      if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
  // 	      chan = detChan(bm1eta, bm2eta); 
  // 	      if (chan == a.index) {
  // 		++nreco;
  // 		if (bm > 0) {
  // 		  ++ncand2;
  // 		}
  // 		//if (bm1pt > b.m1pt && bm2pt > b.m2pt) {
  // 		if (bm1pt > 3.5 && bm2pt > 3.5) {
  // 		  ++nchan; 
  // 		  if (bm > 0) {
  // 		    ++ncand;
  // 		  }
  // 		  if (bm1id && bm2id) {
  // 		    ++nmuid;
  // 		    if (bhlt) {
  // 		      ++nhlt;
  // 		    }
  // 		  }
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //       }
  //     }
  //   }
  //   if (recoPtA > 0) {
  //     a.effPtReco     = static_cast<double>(recoPtB)/static_cast<double>(recoPtA); 
  //     a.effPtRecoE    = dEff(recoPtB, recoPtA);
  //   } else {
  //     a.effPtReco     = 1.; 
  //     a.effPtRecoE    = 0.;
  //   }
  //   a.genAccFileYield = ngen;
  //   a.genAccYield     = a.genAccFileYield/effFilter; 
  //   a.recoYield       = nreco; // reco'ed in chan, basic global reconstruction cuts 
  //   a.muidYield       = nmuid;
  //   a.trigYield       = nhlt;
  //   a.chanYield       = nchan;
  //   a.candYield       = ncand;
  //   a.candYield       = ncand2;

  //   if (a.genAccYield > 0) {
  //     a.acc = a.recoYield/a.genAccYield;
  //     a.accE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genAccYield));
  //   }  

  //   if (a.trigYield > 0) {
  //     a.effCand  = a.candYield/a.trigYield;
  //     a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.trigYield));
  //     a.effCand  = 0.98; // estimate

  //     a.effCand  = a.candYield/a.recoYield;
  //     a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.recoYield));

  //     a.effCand  = a.candYield/a.recoYield;
  //     a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.recoYield));

  //     //     a.effCand  = 0.98; 
  //     //     a.effCandE = 0.04; 
  //   } 

  //   cout << "genAccFileYield:  " << a.genAccFileYield << endl;
  //   cout << "genAccYield:      " << a.genAccYield << endl;
  //   cout << "recoYield:        " << a.recoYield << endl;
  //   cout << "candYield:        " << a.candYield << endl;
  //   cout << "recoPtB:          " << recoPtB << endl;
  //   cout << "recoPtA:          " << recoPtA << endl;

}


// ----------------------------------------------------------------------
void plotClass::accEffFromEffTree(string fname, string dname, numbers &a, cuts &b, int proc) {

  TFile *f = fF[fname];
  if (0 == f) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no file " << fname << " found " << endl;
    return;
  }
  TTree *t  = (TTree*)(f->Get(Form("%s/effTree", dname.c_str())));
  double effFilter(1.); 
  if (!t) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no tree `effTree' found in " 
	 << f->GetName() << " and dir = " << Form("%s/effTree", dname.c_str()) 
	 << endl;
    
    return;
  } else {
    effFilter = fFilterEff[fname]; 
    cout << "anaBmm::accEffFromEffTree(" << a.name << ")" << endl
	 << " get acceptance from file " << f->GetName() << " and dir = " << Form("%s/effTree", dname.c_str()) 
	 << " with filterEff = " << effFilter
	 << endl;
  }

  bool sg(false), no(false), cs(false); 

  int   bprocid, bidx; 
  bool  bhlt;
  float bg1pt, bg2pt, bg1eta, bg2eta;
  float bm1pt, bm1eta, bm2pt, bm2eta;
  bool  bm1gt, bm2gt; 
  bool  bm1id, bm2id; 
  float bm;

  float bg3pt, bg4pt, bg3eta, bg4eta; 
  float bk1pt, bk2pt, bk1eta, bk2eta; 
  bool  bk1gt, bk2gt; 

  t->SetBranchAddress("hlt",&bhlt);
  t->SetBranchAddress("procid",&bprocid);
  t->SetBranchAddress("bidx",&bidx);

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);

  t->SetBranchAddress("m1pt",&bm1pt);
  t->SetBranchAddress("m2pt",&bm2pt);
  t->SetBranchAddress("m1eta",&bm1eta);
  t->SetBranchAddress("m2eta",&bm2eta);

  t->SetBranchAddress("m1gt",&bm1gt);
  t->SetBranchAddress("m2gt",&bm2gt);
  t->SetBranchAddress("m1id",&bm1id);
  t->SetBranchAddress("m2id",&bm2id);


  if (string::npos != a.name.find("signal")) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): SIGNAL " << endl;
    sg = true; 
  }

  if (string::npos != a.name.find("normalization")) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): NORMALIZATION " << endl;
    no = true; 
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("k1pt", &bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k1gt", &bk1gt);
  }

  if (string::npos != a.name.find("control sample")) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): CONTROL SAMPLE " << endl;
    cs = true; 
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("k1pt", &bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k1gt", &bk1gt);

    t->SetBranchAddress("g4pt", &bg4pt);
    t->SetBranchAddress("g4eta",&bg4eta);
    t->SetBranchAddress("k2pt", &bk2pt);
    t->SetBranchAddress("k2eta",&bk2eta);
    t->SetBranchAddress("k2gt", &bk2gt);
  }

  t->SetBranchAddress("m",&bm);


  int nentries = Int_t(t->GetEntries());
  int ngen(0), nreco(0), ncand(0); 
  int chan(-1); 
  int recoPtA(0), recoPtB(0); 
  cout << "channel = " << a.index << endl;
  for (int jentry = 0; jentry < nentries; jentry++) {
    t->GetEntry(jentry);
    if (bidx < 0) continue;
    ++ngen;
    if (proc > 0 && bprocid != proc) continue;
    if (sg) {
      // -- Signal
      if (TMath::Abs(bg1eta) < fAccEtaGen && TMath::Abs(bg2eta) < fAccEtaGen
	  && bg1pt > fAccPt && bg2pt > fAccPt 
	  && bm1pt > fAccPt && bm2pt > fAccPt 
	  && TMath::Abs(bm1eta) < fAccEtaRec && TMath::Abs(bm2eta) < fAccEtaRec
	  && bm1gt && bm2gt
	  ) {
	chan = detChan(bm1eta, bm2eta); 
	if (chan == a.index) {
	  ++nreco; // for acceptance
	  ++recoPtA; 
	  if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
	  if (bm > 0) {
	    ++ncand; // for cand efficiency
	  }
	}
      }
    } else if (no) {
      // -- Normalization
      if (TMath::Abs(bg1eta) < fAccEtaGen && TMath::Abs(bg2eta) < fAccEtaGen  && TMath::Abs(bg3eta) < fAccEtaGen
	  && bg1pt > fAccPt && bg2pt > fAccPt && bg3pt > 0.4
	  && bm1pt > fAccPt && bm2pt > fAccPt && bk1pt > 0.5
	  && TMath::Abs(bm1eta) < fAccEtaRec && TMath::Abs(bm2eta) < fAccEtaRec && TMath::Abs(bk1eta) < fAccEtaRec
	  && bm1gt && bm2gt && bk1gt
	  ) {
	chan = detChan(bm1eta, bm2eta); 
	if (chan == a.index) {
	  ++nreco; // for acceptance
	  ++recoPtA; 
	  if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
	  if (bm > 0) {
	    ++ncand; // for cand efficiency
	  }
	}
      }
    } else if (cs) {
      // -- control sample
      if (TMath::Abs(bg1eta) < fAccEtaGen && TMath::Abs(bg2eta) < fAccEtaGen && TMath::Abs(bg3eta) < fAccEtaGen && TMath::Abs(bg4eta) < fAccEtaGen
	  && bg1pt > fAccPt && bg2pt > fAccPt && bg3pt > 0.4 && bg4pt > 0.4
	  && bm1pt > fAccPt && bm2pt > fAccPt && bk1pt > 0.5 && bk2pt > 0.5
	  && TMath::Abs(bm1eta) < fAccEtaRec && TMath::Abs(bm2eta) < fAccEtaRec && TMath::Abs(bk1eta) < fAccEtaRec && TMath::Abs(bk2eta) < fAccEtaRec
	  && bm1gt && bm2gt && bk1gt && bk2gt
	  ) {
	chan = detChan(bm1eta, bm2eta); 
	if (chan == a.index) {
	  ++nreco; // for acceptance
	  ++recoPtA; 
	  if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
	  if (bm > 0) {
	    ++ncand; // for cand efficiency
	  }
	}
      }
    }
  }

  if (recoPtA > 0) {
    a.effPtReco     = static_cast<double>(recoPtB)/static_cast<double>(recoPtA); 
    a.effPtRecoE    = dEff(recoPtB, recoPtA);
  } else {
    a.effPtReco     = 1.; 
    a.effPtRecoE    = 0.;
  }
  a.genAccFileYield = ngen;
  a.genAccYield     = a.genAccFileYield/effFilter; 
  a.recoYield       = nreco; 
  a.candYield       = ncand;

  if (a.genAccYield > 0) {
    a.acc = a.recoYield/a.genAccYield;
    a.accE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genAccYield));
  }  

  if (a.recoYield > 0) {
    a.effCand  = a.candYield/a.recoYield;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.recoYield));
  } 

  cout << "genAccFileYield:  " << a.genAccFileYield << endl;
  cout << "genAccYield:      " << a.genAccYield << endl;
  cout << "recoYield:        " << a.recoYield << endl;
  cout << "candYield:        " << a.candYield << endl;
  cout << "recoPtB:          " << recoPtB << endl;
  cout << "recoPtA:          " << recoPtA << endl;
  cout << "acc:              " << a.acc << endl;
  cout << "effCand:          " << a.effCand << endl;
  cout << "effPtReco:        " << a.effPtReco << endl;
}


// ----------------------------------------------------------------------
void plotClass::init(const char *files, const char *cuts, const char *dir, int mode) {

  gStyle->SetHatchesSpacing(2);

  fDoPrint = true; // create output
  fVerbose = 0; 

  legg = 0; 

  fSize = 0.05; 
  
  fpFunc  = new initFunc(); 

  fMassLo = 4.5; 
  fMassHi = 6.5;

  fNoLo = 5.10;
  fNoHi = 5.40;

  fCsLo = 5.27;
  fCsHi = 5.47;

  fBgLo = 4.9;
  fBgHi = 5.9;

  fSgLo = 5.20;
  fSgHi = 5.45;

  // -- initialize cuts
  cout << "===> Reading cuts from " << Form("anaBmm.%s.cuts", cuts) << endl;
  readCuts(Form("anaBmm.%s.cuts", cuts)); 
  fNchan = fCuts.size(); 

  printCuts(cout); 

  fFont = 42; 
  fMode = mode;  

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (c0 == 0) {
    cout << "TCanvas c0 not found. Creating my own version." << endl;
    c0 = new TCanvas("c0","--c0--",356,0,656,700);
  }
  tl  = new TLatex();
  tl->SetNDC(kTRUE); 
  tl->SetTextSize(fSize);
  tl->SetTextFont(fFont);

  pl  = new TLine();
  pa  = new TArrow();
  box = new TBox();

  f0 = new TF1("f0", f_p1, 0., 7., 2); 

  f1 = new TF1("f1", f_expo, 0., 7.0, 2);

  f2 = new TF1("f2", f_p1aG, 0., 7., 5);
  f2->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f3 = new TF1("f3", f_eaG, 0., 7., 5);
  f3->SetParNames("Area", "Peak", "Sigma", "Offset", "Slope"); 

  f4 = new TF1("f4", f_ea2G, 0., 7., 8);
  f4->SetParNames("Area", "Peak", "Sigma", "Fraction2", "Peak2", "Sigma2", "Offset", "Slope"); 

  fDirectory = dir; 
  fSuffix    = cuts; 
  if (fMode > 0)  fSuffix += Form("-%d", fMode); 
  cout << "==> Dumping output into " << fDirectory << endl;
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".txt";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fOUT.open(fNumbersFileName.c_str(), ios::app);
  
  loadFiles(files);

  fStampString = "CMS, 5fb^{-1}"; 
  if (fDoUseBDT) {
    fStampString = "BDT preliminary"; 
  } else {
    fStampString = "CNC preliminary"; 
  }
  fStampCms = "fill me";

  string sfiles(files);
  if (string::npos != sfiles.find("2011")) {
    fYear = 2011; 
    fStampCms = "#sqrt{s} = 7 TeV";
  } 
  if (string::npos != sfiles.find("2012")) {
    fYear = 2012; 
    fStampCms = "#sqrt{s} = 8 TeV";
  } 

  string hfname  = fDirectory + "/anaBmm." + fSuffix + ".root";

  printCuts(fOUT); 

}


// ----------------------------------------------------------------------
void plotClass::loadFiles(const char *files) {

  cout << "==> Loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files);
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    
    string sbuffer = string(buffer); 

    string::size_type m1 = sbuffer.find("dset="); 
    string stype = sbuffer.substr(5, m1-6); 

    string::size_type m2 = sbuffer.find("lumi="); 
    string sdset = sbuffer.substr(m1+5, m2-m1-6); 
    m1 = m2; 

    m2 = sbuffer.find("file="); 
    string slumi = sbuffer.substr(m1+5, m2-m1-6); 
    string sfile = sbuffer.substr(m2+5); 
    string sname; 
    
    // -- DATA
    TFile *pF(0); 
    if (string::npos != sdset.find("data")) {
      pF = loadFile(sfile); 
      cout << "open Data file "  << sfile  << " as " << sname << " (" << stype << ") with lumi = " << slumi << endl;
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	sname = "SgData"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

       if (string::npos != stype.find("2011") && string::npos != stype.find("sg")) {
	sname = "SgData2011"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }
 
       if (string::npos != stype.find("2012") && string::npos != stype.find("sg")) {
	sname = "SgData2012"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }
 
     if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	sname = "NoData";
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

     if (string::npos != stype.find("2012") && string::npos != stype.find("no")) {
	sname = "NoData2012";
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

     if (string::npos != stype.find("2011") && string::npos != stype.find("no")) {
	sname = "NoData2011";
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	sname = "CsData"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("2012") && string::npos != stype.find("cs")) {
	sname = "CsData2012"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("2011") && string::npos != stype.find("cs")) {
	sname = "CsData2011"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("dstarpi")) {
	sname = "DstarPiData"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Data")); 
      }

      if (string::npos != stype.find("bmt,HT")) {
	sname = "BmtHT"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "HT")); 
      }
    } else {
      string sfilter = sdset; 
      replaceAll(sfilter, "mc,", ""); 
      double effFilter(-1.); 
      effFilter = atof(sfilter.c_str());
      // -- MC
      cout << "open MC file "  << sfile  << " as " << sname << " (" << stype << ") with lumi = " << slumi 
	   << " filter eff = " << effFilter 
	   << endl;
      pF = loadFile(sfile); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	sname = "SgMc"; 
	fF.insert(make_pair(sname, pF)); 	
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("sg")) {
	sname = "SgMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("m57") && string::npos != stype.find("sg")) {
	sname = "SgM57Bs2MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (5.7GeV)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("m53") && string::npos != stype.find("sg")) {
	sname = "SgM53Bs2MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (5.37GeV)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("m51") && string::npos != stype.find("sg")) {
	sname = "SgM51Bs2MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (5.1GeV)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("bd")) {
	sname = "BdMc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.0e-10)); 
	fBFE.insert(make_pair(sname, 0.1)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("bd")) {
	sname = "BdMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.0e-10)); 
	fBFE.insert(make_pair(sname, 0.1)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("bd")) {
	sname = "BdMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.0e-10)); 
	fBFE.insert(make_pair(sname, 0.1)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("bd")) {
	sname = "BdMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.0e-10)); 
	fBFE.insert(make_pair(sname, 0.1)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("bd")) {
	sname = "BdMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.0e-10)); 
	fBFE.insert(make_pair(sname, 0.1)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	sname = "NoMc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.02e-5)); 
	fBFE.insert(make_pair(sname, 0.03322)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	//	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+}")); 
	fName.insert(make_pair(sname, "MC simulation")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("no")) {
	sname = "NoMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("no")) {
	sname = "NoMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("no")) {
	sname = "NoMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("no")) {
	sname = "NoMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("pu") && string::npos != stype.find("no")) {
	sname = "NoMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (PU)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("CMS") && string::npos != stype.find("no")) {
	sname = "NoMcCMS"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (CMS)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	sname = "CsMc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.16e-5)); 
	fBFE.insert(make_pair(sname, 0.24050)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	//	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi")); 
	fName.insert(make_pair(sname, "MC simulation")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("1e33") && string::npos != stype.find("cs")) {
	sname = "CsMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.16e-5)); 
	fBFE.insert(make_pair(sname, 0.24050)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("2e33") && string::npos != stype.find("cs")) {
	sname = "CsMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.16e-5)); 
	fBFE.insert(make_pair(sname, 0.24050)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("3e33") && string::npos != stype.find("cs")) {
	sname = "CsMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.16e-5)); 
	fBFE.insert(make_pair(sname, 0.24050)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("acc") && string::npos != stype.find("cs")) {
	sname = "CsMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.16e-5)); 
	fBFE.insert(make_pair(sname, 0.24050)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("pu") && string::npos != stype.find("cs")) {
	sname = "CsMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.16e-5)); 
	fBFE.insert(make_pair(sname, 0.24050)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (PU)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("default") && string::npos != stype.find("dstarpi")) {
	sname = "DstarPiMc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 7e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow D^{*} #pi")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bs2KK")) {
	sname = "bgBs2KK"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 25.4e-6)); 
	fBFE.insert(make_pair(sname, 0.15)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2KPi")) {
	sname = "bgBs2KPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.0e-6)); 
	fBFE.insert(make_pair(sname, 0.22)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2PiPi")) {
	sname = "bgBs2PiPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 0.73e-6)); 
	fBFE.insert(make_pair(sname, 0.1917808)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2KMuNu")) {
	sname = "bgBs2KMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.42e-4)); 
	fBFE.insert(make_pair(sname, 0.05)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{-}#mu^{+}#nu")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bd2PiMuNu")) {
	sname = "bgBd2PiMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.42e-4)); 
	fBFE.insert(make_pair(sname, 0.05)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{-}#mu^{+}#nu")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2KK")) {
	sname = "bgBd2KK"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 0.13e-6)); 
	fBFE.insert(make_pair(sname, 0.77)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2KPi")) {
	sname = "bgBd2KPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 19.55e-6)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2PiPi")) {
	sname = "bgBd2PiPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.11e-6)); 
	fBFE.insert(make_pair(sname, 0.043052)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2KP")) {
	sname = "bgLb2KP"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.5e-6)); 
	fBFE.insert(make_pair(sname, 0.254545)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2PiP")) {
	sname = "bgLb2PiP"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.5e-6)); 
	fBFE.insert(make_pair(sname, 0.285714)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p #pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2PMuNu")) {
	sname = "bgLb2PMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	//	fBF.insert(make_pair(sname, 1.36e-4)); 
	// 8e-5 to 1.3e-3
	//	fBF.insert(make_pair(sname, 3.0e-4)); 
	//	fBFE.insert(make_pair(sname, 0.33)); 
	fBF.insert(make_pair(sname, 6.5e-4)); 
	fBFE.insert(make_pair(sname, 1.0)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda^{0}_{b} #rightarrow p#mu^{-}#bar{#nu}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bu2PiMuMu")) {
	sname = "bgBu2PiMuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 2.3e-8)); 
	fBFE.insert(make_pair(sname, 0.26)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #pi^{+}#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bu2KMuMu")) {
	sname = "bgBu2KMuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 4.8e-7)); 
	fBFE.insert(make_pair(sname, 0.0833333)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow K^{+}#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bd2Pi0MuMu")) {
	sname = "bgBd2Pi0MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.2e-8)); 
	fBFE.insert(make_pair(sname, 0.5)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{0}#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bd2K0MuMu")) {
	sname = "bgBd2K0MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.8e-7)); 
	fBFE.insert(make_pair(sname, 0.21052)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{0}#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bd2MuMuGamma")) {
	sname = "bgBd2MuMuGamma"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.34e-10)); 
	fBFE.insert(make_pair(sname, 0.2)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} #gamma")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bs2MuMuGamma")) {
	sname = "bgBs2MuMuGamma"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.89e-8)); 
	fBFE.insert(make_pair(sname, 0.2)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0}_{s} #rightarrow #mu^{+}#mu^{-} #gamma")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

    }
  }

  
  // https://docs.google.com/spreadsheet/ccc?key=0AhrdqmQ22qAldHREbHpVOGVMaHhQaU9ULXAydkpZc2c&hl=en#gid=0
  //   fNgen.insert(make_pair("bgLb2KP",    356.e6)); 
  //   fNgen.insert(make_pair("bgLb2PiP",    374.e6)); 
  //   fNgen.insert(make_pair("bgLb2PMuNu", 6910.e6)); 
  
  //   fNgen.insert(make_pair("SgMc",      249.e6)); 
  //   fNgen.insert(make_pair("bgBs2KK",    1392.e6)); 
  //   fNgen.insert(make_pair("bgBs2KPi",    598.e6)); 
  //   fNgen.insert(make_pair("bgBs2PiPi",   148.e6)); 
  //   fNgen.insert(make_pair("bgBs2KMuNu", 6702.e6)); 
  
  //   fNgen.insert(make_pair("BdMc",        40.e6)); 
  //   fNgen.insert(make_pair("bgBd2PiPi",    394.e6)); 
  //   fNgen.insert(make_pair("bgBd2KPi",     996.e6)); 
  //   fNgen.insert(make_pair("bgBd2KK",      200.e6)); 
  //   fNgen.insert(make_pair("bgBd2PiMuNu", 6742.e6)); 
  
  //   fNgen.insert(make_pair("NoMc",15166.e6)); 
  //   fNgen.insert(make_pair("CsMc", 9970.e6)); 
  
  //   fNgen.insert(make_pair("NoMcAcc", 1)); 
  //   fNgen.insert(make_pair("CsMcAcc", 1)); 

}


// ----------------------------------------------------------------------
TFile* plotClass::loadFile(string file) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}


// ----------------------------------------------------------------------
void plotClass::makeAll(int channel) {


}


// ----------------------------------------------------------------------
void plotClass::dumpSamples() {

  fTEX << "% ----------------------------------------------------------------------" << endl;
  string name; 
  double lumi(0.), n(0.), f(0.), bf(0.), bfe(0.), fe(0.); 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    cout << "===> " << imap->first;
    cout << ": " << fF[imap->first]->GetName();
    cout << " -> " << imap->second << endl;

    name = imap->second; 
    lumi = fLumi[imap->first];
    n = fNgen[imap->first];
    bf = fBF[imap->first];
    bfe = bf*fBFE[imap->first];
    fe = fFilterEff[imap->first];
    f = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1);
    replaceAll(name, "#", "\\"); 
    fTEX <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    cout <<  "                       " 
	 << Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    fTEX <<  scientificTex(bf, bfe, Form("%s:bf:%s", fSuffix.c_str(), imap->first.c_str()), -1, 2) << endl;
    fTEX <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    if (fe > 1.e-6) {
      fTEX <<  Form("\\vdef{%s:efilt:%s}   {\\ensuremath{{%5.4f } } }", fSuffix.c_str(), imap->first.c_str(), fe) << endl;
    } else {
      fTEX <<  Form("\\vdef{%s:efilt:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
    if (n>0) {
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
    } else {
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
    if (f>0) {
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), f) << endl;
    } else {
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
  }
  fTEX.flush();

}


// ----------------------------------------------------------------------
void plotClass::bgBlind(TH1 *h, int mode, double lo, double hi) {

  TVirtualFitter::SetMaxIterations(20000);
  
  if (0 == h) { 
    cout << "plotClass::bgBlind(...): No histogram passed! mode = " << mode << endl;
    return;
  }
  
  TF1 *lF1(0), *lF2(0);

  string modifier = (fDoUseBDT?"bdt":"cnc"); 

  double histCount = h->Integral(h->FindBin(fBgLo+0.0001), h->FindBin(fBgHi)); 
  cout << "bgBlind: histCount = " << histCount 
       << " starting at " << fBgLo+0.0001 << " bin(" << h->FindBin(fBgLo+0.0001) << ")"
       << " to " << fBgHi-0.0001 << " bin(" << h->FindBin(fBgHi-0.0001) << ")" 
       << " mode: " << mode 
       << endl;

  fBgHist  = histCount; 
  fBgHistE  = TMath::Sqrt(histCount); 
  double delta = fSgHi-fSgLo; 
  fBgHistExp  = histCount*(delta)/(fBgHi-fBgLo-delta);
  fBgHistLo =  h->Integral(h->FindBin(fBgLo+0.0001), h->FindBin(fSgLo-0.0001)); 
  fBgHistHi =  h->Integral(h->FindBin(fSgHi+0.0001), h->FindBin(fBgHi-0.0001)); 
  if (histCount > 0) {
    fBgHistE = TMath::Sqrt(histCount)/histCount*fBgHist;
  } else {
    fBgHistE  = 0.2; // FIXME?!
    fBsBgExp  = 0.;
    fBsBgExpE = 0.2; 
    fBdBgExp  = 0.;
    fBdBgExpE = 0.2; 
    return;
  }

  double blind = fSgHi - fSgLo; 
  double scaleBlind = (5.45-5.30)/(fBgHi-fBgLo-blind);
  double scaleBs = (fCuts[fChan]->mBsHi-fCuts[fChan]->mBsLo)/(fBgHi-fBgLo-blind);
  double scaleBd = (fCuts[fChan]->mBdHi-fCuts[fChan]->mBdLo)/(fBgHi-fBgLo-blind);

  if (0 == mode) {
    fBsBgExp = fBgHist*scaleBs;
    fBsBgExpE = fBgHistE*scaleBs;
    fBdBgExp = fBgHist*scaleBd;
    fBdBgExpE = fBgHistE*scaleBd;
    h->DrawCopy();
    return;
  } else if (1 == mode) { 
    lF1 = fpFunc->pol0BsBlind(h); 
    lF2 = fpFunc->pol0(h); 
  } else if (2 == mode) {
    lF1 = fpFunc->pol1BsBlind(h); 
    lF2 = fpFunc->pol1(h); 
  } else if (3 == mode) {
    lF1 = fpFunc->expoBsBlind(h); 
    lF2 = fpFunc->expo(h); 
  } else if (4 == mode) {
    fpFunc->fLo = 5.4;
    fpFunc->fHi = 5.9;
    lF1 = fpFunc->pol0BsBlind(h); 
    lF2 = fpFunc->pol0(h); 
  } else if (5 == mode) {
    fpFunc->fLo = 5.4;
    fpFunc->fHi = 5.9;
    lF1 = fpFunc->pol0BsBlind(h); 
    lF2 = fpFunc->pol0(h); 
  }
  
  lF2->SetLineStyle(kDashed);
  h->Fit(lF1, "rl", "", lo, hi); 
  h->DrawCopy();
  lF2->SetLineColor(kBlue);
  lF2->Draw("same");
  lF2->SetParameters(lF1->GetParameters());
  lF2->SetParErrors(lF1->GetParErrors());

  if (!strcmp(gMinuit->fCstatu.Data(), "CONVERGED ")) {
    lF2->Update();
    double integral = lF2->Integral(fCuts[fChan]->mBsLo, fCuts[fChan]->mBsHi, static_cast<const Double_t*>(0), 1.e-15);
    fBsBgExp  = integral/h->GetBinWidth(1); 
    fBdBgExp  = lF2->Integral(fCuts[fChan]->mBdLo, fCuts[fChan]->mBdHi)/h->GetBinWidth(1); 
    fLoBgExp  = lF2->Integral(4.90, 5.20)/h->GetBinWidth(1); 
    fHiBgExp  = lF2->Integral(5.45, 5.90)/h->GetBinWidth(1); 
    fBsBgExpE = fBsBgExp*(fBgHistE/fBgHist);
    fBdBgExpE = fBdBgExp*(fBgHistE/fBgHist);
    cout << "bgBlind: after fit integrals/hist:  lo: " << fLoBgExp << "/" << fBgHistLo
	 << " hi: " << fHiBgExp << "/" << fBgHistHi << endl; 
  } else {
    cout << "+++ Fit did not converge, take flat bg interpretation, fCstatu = ->" << gMinuit->fCstatu.Data() << "<-" << endl;
    fBsBgExp  = (fCuts[fChan]->mBsHi-fCuts[fChan]->mBsLo)/(5.9-4.9-0.25)*fBgHist;
    fBdBgExp  = (fCuts[fChan]->mBdHi-fCuts[fChan]->mBdLo)/(5.9-4.9-0.25)*fBgHist;
    fLoBgExp  = (5.20-4.90)/(5.9-4.9-0.25)*fBgHist;
    fHiBgExp  = (5.90-5.45)/(5.9-4.9-0.25)*fBgHist;
    fBsBgExpE = fBsBgExp*(fBgHistE/fBgHist);
    fBdBgExpE = fBdBgExp*(fBgHistE/fBgHist);
    cout << "+++ fBsBgExp = " << fBsBgExp << " fLoBgExp: " << fLoBgExp << " fHiBgExp: " << fHiBgExp << endl;
  }

  if (4 == mode) {
    string dname = h->GetName();
    string rname =  Form("bslRare_%s", modifier.c_str());
    TH1D *hr = (TH1D*)fHistFile->Get(rname.c_str()); 
    hr = (TH1D*)(hr->Clone(Form("%s-mode4", rname.c_str())));
    setFilledHist(hr, kBlue, kBlue, 3344);
    hr->DrawCopy("same");
    
    // -- build up total bg
    fBsBgExp   = lF2->Integral(fCuts[fChan]->mBsLo, fCuts[fChan]->mBsHi)/h->GetBinWidth(1); 
    fBdBgExp   = lF2->Integral(fCuts[fChan]->mBdLo, fCuts[fChan]->mBsHi)/h->GetBinWidth(1); 
    fLoBgExp   = lF2->Integral(4.90, 5.20)/h->GetBinWidth(1); 
    fHiBgExp   = lF2->Integral(5.45, 5.90)/h->GetBinWidth(1); 
    cout << "flat combinatorial Expectations Lo: " << fLoBgExp << " Bs: " << fBsBgExp << " Hi: " << fHiBgExp << endl;
    
    fBsBgExp  += hr->Integral(hr->FindBin(fCuts[fChan]->mBsLo+fEpsilon), hr->FindBin(fCuts[fChan]->mBsHi-fEpsilon)); 
    fBdBgExp  += hr->Integral(hr->FindBin(fCuts[fChan]->mBdLo+fEpsilon), hr->FindBin(fCuts[fChan]->mBdHi-fEpsilon)); 
    fLoBgExp  += hr->Integral(hr->FindBin(4.90+fEpsilon), hr->FindBin(5.20-fEpsilon)); 
    fHiBgExp  += hr->Integral(hr->FindBin(5.45+fEpsilon), hr->FindBin(5.90-fEpsilon)); 
    cout << "+ rare bg Expectations Lo:          " << fLoBgExp << " Bs: " << fBsBgExp << " Hi: " << fHiBgExp << endl;
  }


  if (5 == mode) {
    if (fBgHistHi < 1) {
      lF2->SetParameter(0, 0.); 
      lF2->SetParError(0, 0.); 
    }

    string dname = h->GetName();
    string rname =  Form("bslRare_%s", modifier.c_str());
    TH1D *hr = (TH1D*)fHistFile->Get(rname.c_str()); 
    cout << "rare bg histogram: " << rname << " at: " << hr << endl;
    hr = (TH1D*)(hr->Clone(Form("%s-mode5", rname.c_str())));

    double rareBgLo    = hr->Integral(hr->FindBin(4.90+fEpsilon), hr->FindBin(5.20-fEpsilon));
    double flatCombLo  = lF2->Integral(4.90, 5.20)/h->GetBinWidth(1);
    if (flatCombLo < 0) {
      flatCombLo = 0.; 
      cout << "+++ negative lF2 integral, take 0 " << endl;
    }
    double scaleRareBg = (fBgHistLo-flatCombLo)/rareBgLo; 

    cout << "+++ rareBgLo: " << rareBgLo << " flatCombLo: " << flatCombLo << " scaleRareBg: " << scaleRareBg << endl;

    // -- Now RESCALE THE RARE SL BACKGROUND!
    hr->Scale(scaleRareBg); 
    setFilledHist(hr, kBlue, kBlue, 3344);
    hr->DrawCopy("same");

    fBsCoBgExp  = lF2->Integral(fCuts[fChan]->mBsLo, fCuts[fChan]->mBsHi)/h->GetBinWidth(1); if (fBsCoBgExp < 0) fBsCoBgExp = 0.;
    fBdCoBgExp  = lF2->Integral(fCuts[fChan]->mBdLo, fCuts[fChan]->mBdHi)/h->GetBinWidth(1); if (fBdCoBgExp < 0) fBdCoBgExp = 0.;
    fLoCoBgExp  = lF2->Integral(4.90, 5.20)/h->GetBinWidth(1); if (fLoCoBgExp < 0) fLoCoBgExp = 0.;
    fHiCoBgExp  = lF2->Integral(5.45, 5.90)/h->GetBinWidth(1); if (fHiCoBgExp < 0) fHiCoBgExp = 0.;

    fBsSlBgExp  = hr->Integral(hr->FindBin(fCuts[fChan]->mBsLo+fEpsilon), hr->FindBin(fCuts[fChan]->mBsHi-fEpsilon)); if (fBsSlBgExp < 0) fBsSlBgExp = 0.;
    fBdSlBgExp  = hr->Integral(hr->FindBin(fCuts[fChan]->mBdLo+fEpsilon), hr->FindBin(fCuts[fChan]->mBdHi-fEpsilon)); if (fBdSlBgExp < 0) fBdSlBgExp = 0.; 
    fLoSlBgExp  = hr->Integral(hr->FindBin(4.90+fEpsilon), hr->FindBin(5.20-fEpsilon)); if (fLoSlBgExp < 0) fLoSlBgExp = 0.;
    fHiSlBgExp  = hr->Integral(hr->FindBin(5.45+fEpsilon), hr->FindBin(5.90-fEpsilon)); if (fHiSlBgExp < 0) fHiSlBgExp = 0.;

    // -- build up total bg
    fBsBgExp = fBsCoBgExp;
    fBdBgExp = fBdCoBgExp;
    fLoBgExp = fLoCoBgExp;
    fHiBgExp = fHiCoBgExp;
    cout << "flat combinatorial Expectations Lo: " << fLoBgExp << " Bs: " << fBsBgExp << " Hi: " << fHiBgExp << endl;
    
    fBsBgExp += fBsSlBgExp; 
    fBdBgExp += fBdSlBgExp; 
    fLoBgExp += fLoSlBgExp; 
    fHiBgExp += fHiSlBgExp; 
    cout << "+ scaled rare bg Expectations Lo:   " << fLoBgExp << " Bs: " << fBsBgExp << " Hi: " << fHiBgExp << endl;

    // FIXME
    fBsBgExpE = 0.5*fBsBgExp;
    fBdBgExpE = 0.5*fBdBgExp;
  }

  delete lF1; 
  delete lF2; 


}


// ----------------------------------------------------------------------
void plotClass::normYield(TH1 *h, int mode, double lo, double hi, double preco) {

//   double pReco = (preco<0? (1 == mode?5.145:5.146): preco);
  double pReco = (preco<0? 5.145: preco);

  h->SetAxisRange(5.0, 5.8);

  TF1 *lF1(0), *lBg(0);

  string name(h->GetName()); 
  double sigma1(0.03);
  if (1 == mode) sigma1 = 0.04; 
  if (string::npos != name.find("hNormC")) sigma1 = 0.02;
  double sigma2(0.1); 
  if (string::npos != name.find("hNormC")) sigma2 = 0.05;
  
  if (0 == mode) { 
    fpFunc->fLo = lo; 
    fpFunc->fHi = hi; 
    //    lF1 = fpFunc->expoErrgauss2c(h, 5.27, 0.03, 0.1, pReco); 
    lF1 = fpFunc->expoErrgauss2(h, 5.28, sigma1, 5.27, sigma2, pReco); 
    lF1->SetNpx(100000);
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  } else {
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoErrGauss(h, 5.27, 0.056, pReco); 
    lF1->SetNpx(100000);
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  }
  h->Fit(lF1, "rem", "", lo, hi); 

  if (0 == mode) {
    cout << "par i = "; 
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(6+i));
      cout << " " << lBg->GetParameter(i);
    }
    cout << endl;
  } else {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  }
  double c  = lF1->GetParameter(0); 
  c = lF1->Integral(5.15, 5.45) - lBg->Integral(5.15, 5.45); 
  double cE = lF1->GetParError(0); 
  double ierr = lF1->IntegralError(5.15, 5.45)/h->GetBinWidth(1); 

  fNoSig = c/h->GetBinWidth(1);
  if (ierr > TMath::Sqrt(fNoSig)) {
    fNoSigE = ierr;
  } else {
    fNoSigE = cE/c*fNoSig;
  }

  cout << "N(Sig) = " << fNoSig << " +/- " << fNoSigE << endl;
  cout << "chi2/dof = " << lF1->GetChisquare() << "/" << lF1->GetNDF() 
       << " prob = " << TMath::Prob(lF1->GetChisquare(), lF1->GetNDF()) 
       << endl;
  
  shrinkPad(0.13, 0.2);
  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muK} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.1, 2.0); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  // -- Overlay BG function
  lBg->SetLineStyle(kDashed);
  lBg->SetLineColor(kRed);
  lBg->SetLineWidth(3);
  lBg->Draw("same");

  tl->SetTextSize(0.07); 
  tl->SetTextColor(kBlack); 
  if (0 == mode) {
    tl->DrawLatex(0.6, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.6, 0.8, "Endcap");   
  } 
  
  stamp(0.20, fStampString, 0.67, fStampCms); 
  if (fDoPrint) {
    
    string pdfname;
    string hname(h->GetName());
    if (string::npos != hname.find("NormC")) {
      pdfname = Form("%s/%s-cncnormC-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode);
      if (fDoUseBDT)  pdfname = Form("%s/bdtnormC-data-chan%d.pdf", fDirectory.c_str(), mode);
    } else {
      pdfname = Form("%s/%s-cncnorm-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode);
      if (fDoUseBDT)  pdfname = Form("%s/bdtnorm-data-chan%d.pdf", fDirectory.c_str(), mode);
    }

    c0->SaveAs(pdfname.c_str());
  }
  

  delete lF1; 
  delete lBg; 

}


// ----------------------------------------------------------------------
// A modified fitting class for the normalization channel.
// Includes the landau peak for Bu2JpsiPi events. The landau has fixed parameters except the amplitude.
// Tested only outside the plotClass framework.
void plotClass::normYield2(TH1 *h, int mode, double lo, double hi, double preco) {

//   double pReco = (preco<0? (1 == mode?5.145:5.146): preco);
  double pReco = (preco<0? 5.145: preco);

  h->SetAxisRange(5.0, 5.8);

  TF1 *lF1(0), *lBg(0);

  TVirtualFitter::SetMaxIterations(20000);

  string name(h->GetName()); 
  double sigma1(0.03);
  if (1 == mode) sigma1 = 0.04; 
  //if (string::npos != name.find("hNormC")) sigma1 = 0.02;
  double sigma2(0.1); 
  //if (string::npos != name.find("hNormC")) sigma2 = 0.05;
  
  if (0 == mode) { 
    fpFunc->fLo = lo; 
    fpFunc->fHi = hi; 
    cout<<sigma1<<" "<<sigma2<<" "<<pReco<<endl;
    //    lF1 = fpFunc->expoErrgauss2c(h, 5.27, 0.03, 0.1, pReco); 
    //lF1 = fpFunc->expoErrgauss2(h, 5.28, sigma1, 5.27, sigma2, pReco); 
    lF1 = fpFunc->expoErrgauss2Landau(h, 5.28, sigma1, 5.27, sigma2, pReco); 
    lF1->SetNpx(100000);
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 

  } else {

    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoErrGaussLandau(h, 5.28, 0.056, pReco); 
    lF1->SetNpx(100000);
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  }
  //  h->Fit(lF1, "rem", "", lo, hi); 
  h->Fit(lF1, "r", "", lo, hi); 

  if (0 == mode) {
    cout << "par i = "; 
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(6+i));
      cout << " " << lBg->GetParameter(i);
    }
    cout << endl;
  } else {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  }

  //get the turnon point for the error func
  if (0 == mode) {
    fNoErrTurnon = lF1->GetParameter(8);
  } else {
    fNoErrTurnon = lF1->GetParameter(5);
  }

  // get the full area of 1 or 2 gaussians
  double binw  = h->GetBinWidth(1);
  double area = 0;
  double cons = lF1->GetParameter(0); 
  double sig1  = lF1->GetParameter(2);
  if(mode==0) {  // barrel
    double frac  = lF1->GetParameter(3);
    double sig2  = lF1->GetParameter(5);
    area  = 2.507 * cons * ( sig1 + frac*sig2);
  } else { // forward
    area  = cons;
  }
  //cout<<" const "<< cons <<" "<<area<<" "<<area/binw<<endl; // this is with full tails 

  // integrals
  double cf  = lF1->Integral(5.15, 5.45); 
  double cb  = lBg->Integral(5.15, 5.45); 
  double cE = lF1->GetParError(0);
  
  // landau
  //double area = 2.507 * cons * (par[2] + par[3]*par[5]); // area of the 2 gaussians
  double area_landau = 0.049 * area;
  double mpvl = 0, sigl=0;
  if(mode==0) {mpvl = lF1->GetParameter(12); sigl = lF1->GetParameter(13);} // bar
  else        {mpvl = lF1->GetParameter(9);  sigl = lF1->GetParameter(10);} // end
  double constl = area_landau/sigl;
  //cout<< mpvl << " " << sigl << " "<<constl<<endl; // landau const = area/sigma

  // area under landau
  TF1 *fl = new TF1("test", fa_landausimp, 4.5, 6, 3);
  fl->SetParameters(mpvl,sigl,constl); 
  double cl  = fl->Integral(5.15, 5.45); // estimate landau contribution 

  double sig = cf - cb - cl; // all - background - landau
  //cout<<" integrals "<< sig <<" "<< cf <<" "<< cb <<" "<< cl <<endl;
  //cout<<" integrals "<<sig/binw<<" "<<cf/binw <<" "<<cb/binw<<" "<<cl/binw<<endl;

  //   double fNoSig = sig/binw;  // convert to counts
  //   double fNoSigE = 0;
  fNoSig = sig/binw;  // convert to counts
  fNoSigE = 0;

  if(true) { // skip
    double ierr = lF1->IntegralError(5.15, 5.45)/binw; // slow
    if (ierr > TMath::Sqrt(fNoSig)) {
      fNoSigE = ierr;
    } else {
      fNoSigE = cE/sig*fNoSig;
    }
  } 

  cout << "N(Sig) = " << fNoSig << " +/- " << fNoSigE << endl;
  cout << "chi2/dof = " << lF1->GetChisquare() << "/" << lF1->GetNDF() 
       << " prob = " << TMath::Prob(lF1->GetChisquare(), lF1->GetNDF()) 
       << endl;


  shrinkPad(0.13, 0.2);
  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muK} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.1, 2.0); 
  h->SetMinimum(0.01); 
  //h->SetLineColor(kBlack); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  //cout<<" after plot histo"<<endl;

  // -- Overlay BG function
  lBg->SetLineStyle(kDashed);
  lBg->SetLineColor(kRed);
  lBg->SetLineWidth(3);
  lBg->DrawCopy("same");

  //cout<<" after plot bg"<<endl;

  // Plot landau
  fl->SetLineStyle(kDotted);
  fl->SetLineColor(kBlue);
  fl->SetLineWidth(3);
  fl->DrawCopy("same");

  //cout<<" after plot landau "<<cl<<" "<<cl/binw<<endl;

  tl->SetNDC(kTRUE); 
  tl->SetTextSize(0.07); 
  tl->SetTextColor(kBlack); 
  if (0 == mode) {
    tl->DrawLatex(0.6, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    cout << "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ printing endcap: " << tl->GetTextSize() << endl;
    tl->DrawLatex(0.6, 0.8, "Endcap");   
  } 
  
  stamp(0.20, fStampString, 0.67, fStampCms); 
  if (fDoPrint) {
    
    string pdfname;
    string hname(h->GetName());
    if (string::npos != hname.find("NormC")) {
      pdfname = Form("%s/%s-cncnormC-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode);
      if (fDoUseBDT)  pdfname = Form("%s/%s-bdtnormC-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode);
    } else {
      pdfname = Form("%s/%s-cncnorm-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode);
      if (fDoUseBDT)  pdfname = Form("%s/%s-bdtnorm-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode);
    }

   c0->SaveAs(pdfname.c_str());
  }

  
  delete lF1; 
  delete lBg; 
  delete fl;

}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// 2 gaussian, 2nd one has a fixed shape and the area related to the area of the 1st gaussian
// fraction: 
// -1 fraction of the 2nd gaussian is left free to be fitted, 
// >=0. the fraction is fixed to the passed value 
// preco:
// -1 the step error function is not used in the fit
// >=0. the step is used with the edge equal to preco
// 
void plotClass::csYield2(TH1 *h, int mode, double lo, double hi, double fraction, double preco) {

  //const fraction = 0.14; // area fraction of the 2nd gauss

  h->SetAxisRange(5.0, 5.8);
 
  TF1 *lF1(0), *lBg(0);

  fpFunc->resetLimits();

  fpFunc->fLo = lo; //5.0;
  fpFunc->fHi = hi; //5.5;
  if (0 == mode) { 
    lF1 = fpFunc->expoErrgauss2f(h, 5.37, 0.040, 5.425, 0.079,fraction, preco);  // 2 gaussians 
    lBg = fpFunc->expoErrGauss(h,5.425,0.079); 
  } else {
    lF1 = fpFunc->expoErrgauss2f(h, 5.37, 0.040, 5.425, 0.079,fraction, preco); 
    lBg = fpFunc->expoErrGauss(h,5.425,0.079); 
  }

  lF1->SetLineColor(kBlue); 
  lF1->Draw();

  //  return;

  h->Fit(lF1, "rm", "", lo, hi); 

  // 
  double ratio  = lF1->GetParameter(3); 
  double area1stGauss = 2.5066 * lF1->GetParameter(0) * lF1->GetParameter(2);
  double area2ndGauss = ratio * area1stGauss;
  //cout<<ratio<<" "<<area1stGauss<<" "<<area2ndGauss<<endl;

  // setup the background function
  lBg->SetParameter(0, area2ndGauss);  // gauss area
  cout << "par " << 0 << ": " << lBg->GetParName(0) << " = " << lBg->GetParameter(0) << endl;
  lBg->SetParameter(1, lF1->GetParameter(4));  // gauss peak
  cout << "par " << 1 << ": " << lBg->GetParName(1) << " = " << lBg->GetParameter(1) << endl;
  lBg->SetParameter(2, lF1->GetParameter(5));  // gauss sigma
  cout << "par " << 2 << ": " << lBg->GetParName(2) << " = " << lBg->GetParameter(2) << endl;
  // test rest
  for (int i = 3; i < lBg->GetNpar(); ++i) {
    lBg->SetParameter(i, lF1->GetParameter(3+i));
    cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
  }

  double c1 = lF1->Integral(5.25, 5.5); 
  double c2 = lBg->Integral(5.25, 5.5); 
  double c = c1 - c2;
  //cout<<c1<<" "<<c2<<" "<<c<<endl;

  double cE = lF1->GetParError(0)/lF1->GetParameter(0); 
  if (0 == mode) {
    cE = (lF1->GetParError(0)/lF1->GetParameter(0))*(lF1->GetParError(0)/lF1->GetParameter(0))
       + (lF1->GetParError(6)/lF1->GetParameter(6))*(lF1->GetParError(6)/lF1->GetParameter(6));
    if(fraction<0.) cE += (lF1->GetParError(3)/lF1->GetParameter(3))*(lF1->GetParError(3)/lF1->GetParameter(3));
    cE = TMath::Sqrt(cE);
  }
  double ierr = lF1->IntegralError(5.25, 5.5)/h->GetBinWidth(1); 
 
  fCsSig = c/h->GetBinWidth(1);
  fCsKstFrac = lF1->GetParameter(3);
  if (ierr > TMath::Sqrt(fCsSig)) {
    cout << "chose integral error" << endl;
    fCsSigE = ierr;
  } else {
    cout << "chose parameter error" << endl;
    fCsSigE = cE*fCsSig;
  }

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  cout << "chi2/dof = " << lF1->GetChisquare() << "/" << lF1->GetNDF() 
       << " prob = " << TMath::Prob(lF1->GetChisquare(), lF1->GetNDF()) 
       << endl;

  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muKK} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.6); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  lF1->SetRange(lo, hi);
  lF1->SetLineColor(kBlack); 
  lF1->SetLineStyle(kSolid);
  lF1->SetLineWidth(3);
  lF1->Draw("same");

  // -- Overlay BG function
  lBg->SetRange(lo, hi);
  lBg->SetLineStyle(kDashed);
  lBg->SetLineColor(kRed);
  lBg->SetLineWidth(3);
  lBg->Draw("same");
  
  // also overlay 2nd Gaussian
  TF1 *lgauss2 = new TF1("gaussian","gaus",lo,hi);
  float gaussmean = area2ndGauss/(lF1->GetParameter(5)*2.5066);
  lgauss2->SetParameter(0, gaussmean ); //const
  lgauss2->SetParameter(1, lF1->GetParameter(4)); //mean 
  lgauss2->SetParameter(2, lF1->GetParameter(5)); //sigma

  lgauss2->SetLineStyle(kDotted);
  lgauss2->SetLineColor(kBlue);
  lgauss2->SetLineWidth(3);
  lgauss2->Draw("same");


  tl->SetTextSize(0.07); 
  tl->SetTextColor(kBlack); 
  if (0 == mode) {
    tl->DrawLatex(0.22, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.22, 0.8, "Endcap");   
  } 

  string hname(h->GetName());
  if (fDoPrint) {
    if (fDoUseBDT) c0->SaveAs(Form("%s/%s-bdtcs-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode));
    else c0->SaveAs(Form("%s/%s-cnccs-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode));
  }

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 
  delete lBg; 

}

// ----------------------------------------------------------------------
void plotClass::csYield(TH1 *h, int mode, double lo, double hi, double preco) {

  h->SetAxisRange(5.0, 5.8);

  TF1 *lF1(0), *lBg(0);

  fpFunc->resetLimits();

  fpFunc->fLo = lo; //5.0;
  fpFunc->fHi = hi; //5.5;
  if (0 == mode) { 
    lF1 = fpFunc->expoGauss(h, 5.37, 0.040); 
    lBg = fpFunc->expo(fpFunc->fLo, fpFunc->fHi); 
  } else {
    lF1 = fpFunc->expoGauss(h, 5.37, 0.06); 
    lBg = fpFunc->expo(fpFunc->fLo, fpFunc->fHi); 
  }

  lF1->SetLineColor(kBlue); 
  lF1->Draw();
  //  return;

  h->Fit(lF1, "rm", "", lo, hi); 


  if (0 == mode) {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  } else {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  }
  

  double c  = lF1->GetParameter(0); 
  c = lF1->Integral(5.25, 5.5) - lBg->Integral(5.25, 5.5); 
  double cE = lF1->GetParError(0)/lF1->GetParameter(0); 
  if (0 == mode) {
    cE = (lF1->GetParError(0)/lF1->GetParameter(0))*(lF1->GetParError(0)/lF1->GetParameter(0))
      + (lF1->GetParError(3)/lF1->GetParameter(3))*(lF1->GetParError(3)/lF1->GetParameter(3));
    cE = TMath::Sqrt(cE);
  }
  double ierr = lF1->IntegralError(5.25, 5.5)/h->GetBinWidth(1); 

  fCsSig = c/h->GetBinWidth(1);
  if (ierr > TMath::Sqrt(fCsSig)) {
    cout << "chose integral error" << endl;
    fCsSigE = ierr;
  } else {
    cout << "chose parameter error" << endl;
    fCsSigE = cE*fCsSig;
  }

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  cout << "chi2/dof = " << lF1->GetChisquare() << "/" << lF1->GetNDF() 
       << " prob = " << TMath::Prob(lF1->GetChisquare(), lF1->GetNDF()) 
       << endl;

  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muKK} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.6); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  // -- Overlay BG function
  lBg->SetLineStyle(kDashed);
  lBg->SetLineColor(kRed);
  lBg->SetLineWidth(3);
  lBg->Draw("same");

  tl->SetTextSize(0.07); 
  tl->SetTextColor(kBlack); 
  if (0 == mode) {
    tl->DrawLatex(0.22, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.22, 0.8, "Endcap");   
  } 

  string hname(h->GetName());
  if (fDoPrint) {
    if (fDoUseBDT) c0->SaveAs(Form("%s/%s-bdtcs-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode));
    else c0->SaveAs(Form("%s/%s-cnccs-data-chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), mode));
  }

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 
  delete lBg; 

}


// ----------------------------------------------------------------------
void plotClass::dpYield(TH1 *h, int mode, double lo, double hi, int bdtc) {

  TF1 *lF1(0), *lBg(0);

  fpFunc->resetLimits();
  fpFunc->fLimit[1] = true;
  fpFunc->fLimitLo[1] = 5.26;
  fpFunc->fLimitHi[1] = 5.30;

  fpFunc->fLimit[2] = true;
  fpFunc->fLimitLo[2] = 0.04;
  fpFunc->fLimitHi[2] = 0.06;

  if (0 == mode) { 
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoGauss(h, 5.279, 0.050); 
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  } else {
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoGauss(h, 5.279, 0.06); 
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  }

  h->Fit(lF1, "rm", "", lo, hi); 


  if (0 == mode) {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  } else {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  }
  
  double intmin = lF1->GetParameter(1)-lF1->GetParameter(2);
  double intmax = lF1->GetParameter(1)+lF1->GetParameter(2);
  double c  = lF1->GetParameter(0); 
  c = lF1->Integral(intmin, intmax) - lBg->Integral(intmin, intmax); 
  double ierr = lF1->IntegralError(intmin, intmax);
  if (ierr/c < 0.001) {
    cout << "ierr too small at " << ierr << endl;
    ierr = lF1->GetParError(0);
    cout << "chosing par error for ierr = " << ierr << endl;
  }

  fDpSig = c/h->GetBinWidth(1);
  fDpSigE = ierr/h->GetBinWidth(1);

  cout << "intmin = " << intmin << " intmax = " << intmax << endl;
  cout << "N(Sig) = " << fDpSig << " +/- " << fDpSigE << endl;
  cout << "chi2/dof = " << lF1->GetChisquare() << "/" << lF1->GetNDF() 
       << " prob = " << TMath::Prob(lF1->GetChisquare(), lF1->GetNDF()) 
       << endl;

  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{D^{*}#pi} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.6); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  // -- Overlay BG function
  lBg->SetLineStyle(kDashed);
  lBg->SetLineColor(kRed);
  lBg->SetLineWidth(3);
  lBg->Draw("same");

  tl->SetTextSize(0.07); 
  tl->SetTextColor(kBlack); 
  if (0 == mode) {
    tl->DrawLatex(0.22, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.22, 0.8, "Endcap");   
  } 

  tl->SetTextSize(0.04); 
  tl->DrawLatex(0.5, 0.85, Form("b > %5.2f", bdtc/100.));   
  tl->DrawLatex(0.5, 0.80, Form("#mu: %5.4f#pm%5.4f", lF1->GetParameter(1), lF1->GetParError(1)));   
  tl->DrawLatex(0.5, 0.75, Form("#sigma: %5.4f#pm%5.4f", lF1->GetParameter(2), lF1->GetParError(2)));   
  tl->DrawLatex(0.5, 0.70, Form("N_{sig}: %5.1f#pm%5.1f", fDpSig, fDpSigE));   
  tl->DrawLatex(0.5, 0.65, Form("#chi^{2}/dof: %3.1f/%i", lF1->GetChisquare(), lF1->GetNDF()));   


  if (fDoPrint) {
    if (fDoUseBDT) c0->SaveAs(Form("%s/bdtdp-data-bdt%d-chan%d.pdf", fDirectory.c_str(), bdtc, mode));
    else c0->SaveAs(Form("%s/dp-data-bdt%d-chan%d.pdf", fDirectory.c_str(), bdtc, mode));
  }

  cout << "N(Sig) = " << fDpSig << " +/- " << fDpSigE << endl;
  
  delete lF1; 

}


// ----------------------------------------------------------------------
void plotClass::printNumbers(numbers &a, ostream &OUT) {
  OUT << "======================================================================" << endl;
  OUT << "numbers for \""  << a.name.c_str() << "\"" << endl;
  OUT << "fitYield        = " << a.fitYield << "+/-" << a.fitYieldE << endl;
  OUT << "fitYieldC       = " << a.fitYieldC << "+/-" << a.fitYieldCE << endl;
  OUT << "bgHist          = " << a.bgObs << endl;
  OUT << "bgHistLo        = " << a.offLo << endl;
  OUT << "bgHistHi        = " << a.offHi << endl;
  OUT << "genAccFileYield = " << a.genAccFileYield << endl;
  OUT << "genAccYield     = " << a.genAccYield << endl;
  OUT << "genFileYield    = " << a.genFileYield << endl;
  OUT << "genYield        = " << a.genYield << endl;
  OUT << "recoYield       = " << a.recoYield << endl;
  OUT << "candYield       = " << a.candYield << endl;
  OUT << "absNoCutsYield  = " << a.absNoCutsYield << endl;
  OUT << "ana0Yield       = " << a.ana0Yield << endl;
  OUT << "anaYield        = " << a.anaYield << endl;
  OUT << "anaMuYield      = " << a.anaMuonYield << endl;
  OUT << "anaTrigYield    = " << a.anaTriggerYield << endl;
  OUT << "anaWmcYield     = " << a.anaWmcYield << endl;
  OUT << "mBsLo           = " << a.mBsLo << endl;
  OUT << "mBsHi           = " << a.mBsHi << endl;
  OUT << "mBdLo           = " << a.mBdLo << endl;
  OUT << "mBdHi           = " << a.mBdHi << endl;
  OUT << "PSS             = " << a.pss << endl;
  OUT << "PDS             = " << a.pds << endl;
  OUT << "PSD             = " << a.psd << endl;
  OUT << "PDD             = " << a.pdd << endl;
  OUT << "bsRare          = " << a.bsRare << endl;
  OUT << "bdRare          = " << a.bdRare << endl;
  OUT << "gen filter      = " << a.effGenFilter << endl;
  OUT << "pt ratio        = " << a.effPtReco << "+/-" << a.effPtRecoE << endl;
  OUT << "acceptance      = " << a.acc << "+/-" << a.accE << endl;
  OUT << "effCand         = " << a.effCand << "+/-" << a.effCandE << endl;
  OUT << "effAna          = " << a.effAna << "+/-" << a.effAnaE << endl; 

  OUT << "effMuidMC       = " << a.effMuidMC << "+/-" << a.effMuidMCE << endl;
  OUT << "effRhoMuidTNPMC = " << a.effRhoMuidTNPMC << "+/-" << a.effRhoMuidTNPMCE << endl;
  OUT << "effMuidTNPMC    = " << a.effMuidTNPMC << "+/-" << a.effMuidTNPMCE << endl;
  OUT << "effRhoMuidTNP   = " << a.effRhoMuidTNP << "+/-" << a.effRhoMuidTNPE << endl;
  OUT << "effMuidTNP      = " << a.effMuidTNP << "+/-" << a.effMuidTNPE << endl;
  OUT << "rhoMuidTNPMC    = " << a.rhoMuidTNPMC << "+/-" << a.rhoMuidTNPMCE << endl;

  OUT << "effTrigMC       = " << a.effTrigMC << "+/-" << a.effTrigMCE << endl;
  OUT << "effRhoTrigTNPMC = " << a.effRhoTrigTNPMC << "+/-" << a.effRhoTrigTNPMCE << endl;
  OUT << "effTrigTNPMC    = " << a.effTrigTNPMC << "+/-" << a.effTrigTNPMCE << endl;
  OUT << "effRhoTrigTNP   = " << a.effRhoTrigTNP << "+/-" << a.effRhoTrigTNPE << endl;
  OUT << "effTrigTNP      = " << a.effTrigTNP << "+/-" << a.effTrigTNPE << endl;
  OUT << "rhoTrigTNPMC    = " << a.rhoTrigTNPMC << "+/-" << a.rhoTrigTNPMCE << endl;

  OUT << "effProd(MC)     = " << a.effProdMC << "+/-" << a.effProdMCE << endl;
  OUT << "effProd(TNP)    = " << a.effProdTNP << "+/-" << a.effProdTNPE << endl;
  OUT << "effProd(TNPMC)  = " << a.effProdTNPMC << "+/-" << a.effProdTNPMCE << endl;
  OUT << "effTot          = " << a.effTot << "+/-" << a.effTotE << endl; 
  OUT << "combGenYield    = " << a.combGenYield << endl; 
  OUT << "prodGenYield    = " << a.prodGenYield << endl; 
  OUT.flush();

}

// ----------------------------------------------------------------------
void plotClass::dumpCutNames(const char *hname) {

  //  ofstream OUT(fNumbersFileName.c_str(), ios::app);
  fTEX << "% ----------------------------------------------------------------------" << endl;
  
  TH1D *h = (TH1D*)gFile->Get(hname); 
  string cut, empty(""), cutstring, cutline, cutvalue; 
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i));
    if (cut == empty) continue;
    string::size_type n1 = cut.find_first_of(":"); 
    string::size_type n2 = cut.find_last_of(":")-1; 
    
    cutstring = cut.substr(0, n1-1); 
    cutline   = cut.substr(n1+3, n2-n1-4);
    replaceAll(cutline, "#", "\\"); 
    cutvalue  = cut.substr(n2+2, cut.size()-n2);

    //    cout << cut << " cutstring==>" << cutstring << "<== cutvalue==>" << cutvalue << " cutline==>" << cutline << "<==" << endl;
    if (string::npos != cutstring.find("JSON")) {
      replaceAll(cutvalue, "_", "\\_"); 
      replaceAll(cutvalue, "cuts/", ""); 
      fTEX <<  Form("\\vdef{%s:%s:cutLine}   {{%s } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:%s:cutValue}  {{%s } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
    } else {
      fTEX <<  Form("\\vdef{%s:%s:cutLine}   {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutline.c_str()) << endl;
      fTEX <<  Form("\\vdef{%s:%s:cutValue}  {\\ensuremath{{%s } } }", fSuffix.c_str(), cutstring.c_str(), cutvalue.c_str()) << endl;
    }

  }
  fTEX.flush();
}



// ----------------------------------------------------------------------
string plotClass::scientificTex(double n, double nE, const char *name, double base, int digits) {

  char line[200]; 
  float a1 = n/TMath::Abs(base); 
  float a2 = nE/TMath::Abs(base); 

  int expo(0); 
  if (base < 0) {
    expo = static_cast<int>(TMath::Log10(n)); 
    if (expo < 0) expo -= 1; 
    a1 = n/TMath::Abs(TMath::Power(10, expo));
    a2 = nE/TMath::Abs(TMath::Power(10, expo));
  } else {
    expo = static_cast<int>(TMath::Log10(base)); 
  }

  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name);
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%i } } }", name, static_cast<int>(n));
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%2.1f \\pm %2.1f)\\times 10^{%i}} } }", name, a1, a2, expo);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%3.2f \\pm %3.2f)\\times 10^{%i}} } }", name, a1, a2, expo);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%4.3f \\pm %4.3f)\\times 10^{%i}} } }", name, a1, a2, expo);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.4f \\pm %5.4f)\\times 10^{%i}} } }", name, a1, a2, expo);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%6.5f \\pm %6.5f)\\times 10^{%i}} } }", name, a1, a2, expo);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%f \\pm %f)\\times 10^{%i}} } }", name, a1, a2, expo);
  }

  return string(line); 
}


// ----------------------------------------------------------------------
string plotClass::formatTex(double n, std::string name, int digits, int sgn) {

  char line[200]; 
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%i } } }", name.c_str(), static_cast<int>(n));
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.1f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.1f } } }", name.c_str(), n);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.2f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.2f } } }", name.c_str(), n);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.3f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.3f } } }", name.c_str(), n);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.4f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.4f } } }", name.c_str(), n);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%6.5f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+6.5f } } }", name.c_str(), n);
  } else if (6 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%7.6f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+7.6f } } }", name.c_str(), n);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+f } } }", name.c_str(), n);
  }

  return string(line); 
}


// ----------------------------------------------------------------------
string plotClass::formatTex(string s, string name) {
  char line[1000]; 
  sprintf(line, "\\vdef{%s}   {%s } ", name.c_str(), s.c_str());
  return string(line); 
}

// ----------------------------------------------------------------------
string plotClass::formatTex(double n, string name, string tag) {
  
  char line[200]; 
  
  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str(), tag.c_str());
    //   } else if ( n > 1.e10) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
    //   } else if ( n > 1.e4) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
  } else if ( n > 100. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%6.0f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.2f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-1) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-3) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
    //   } else if ( n > 1.e-9 ) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
    //   } else if ( n > 1.e-19 ){
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
  } else {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{0.0 } } }", name.c_str(), tag.c_str());
  }
  
  string result(line); 
  return result;
}



// ----------------------------------------------------------------------
void plotClass::makeCanvas(int i) {
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 1000);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}



// ----------------------------------------------------------------------
void plotClass::replaceAll(std::string &s, std::string a, std::string b) {
  
  TString ts(s.c_str()); 
  ts.ReplaceAll(a.c_str(), b.c_str()); 
  s = ts.Data(); 

}


// ----------------------------------------------------------------------
void plotClass::setErrors(TH1D *h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinLowEdge(i) > 5.1 && h->GetBinLowEdge(i+1) <= 5.5) continue;
    if (h->GetBinContent(i) < 1) {
      h->SetBinError(i, 1.); 
    }
  }
}


// ----------------------------------------------------------------------
void plotClass::newLegend(double x1, double y1, double x2, double y2, string title) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}


// ----------------------------------------------------------------------
void plotClass::readCuts(const char *filename) {
  cout << "==> plotClass: Reading " << filename << " for cut settings" << endl;
  vector<string> cutLines; 
  char  buffer[200];
  ifstream is(filename);
  while (is.getline(buffer, 200, '\n')) {
    cutLines.push_back(string(buffer));
  }

  char CutName[100], XmlName[1000];
  float CutValue;
  int dump(1), ok(0);

  cuts *a = 0;
  
  fCuts.clear();

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "index")) {
      ok = 1;
      if (dump) cout << "index:            " << CutValue << endl;
      if (a) fCuts.push_back(a); 
      a = new cuts; 
      a->index = static_cast<int>(CutValue); 
    }
    
    if (!strcmp(CutName, "mBdLo")) {
      a->mBdLo = CutValue; ok = 1;
      if (dump) cout << "mBdLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "bdtpt")) {
      a->bdtpt = CutValue; ok = 1;
      if (dump) cout << "bdtpt:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "bdt")) {
      a->bdt = CutValue; ok = 1;
      if (dump) cout << "bdt:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "bdtMax")) {
      a->bdtMax = CutValue; ok = 1;
      if (dump) cout << "bdtMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBdHi")) {
      a->mBdHi = CutValue; ok = 1;
      if (dump) cout << "mBdHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsLo")) {
      a->mBsLo = CutValue; ok = 1;
      if (dump) cout << "mBsLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsHi")) {
      a->mBsHi = CutValue; ok = 1;
      if (dump) cout << "mBsHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMin")) {
      a->etaMin = CutValue; ok = 1;
      if (dump) cout << "etaMin:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMax")) {
      a->etaMax = CutValue; ok = 1;
      if (dump) cout << "etaMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "pt")) {
      a->pt = CutValue; ok = 1;
      if (dump) cout << "pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1pt")) {
      a->m1pt = CutValue; ok = 1;
      if (dump) cout << "m1pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2pt")) {
      a->m2pt = CutValue; ok = 1;
      if (dump) cout << "m2pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1eta")) {
      a->m1eta = CutValue; ok = 1;
      if (dump) cout << "m1eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2eta")) {
      a->m2eta = CutValue; ok = 1;
      if (dump) cout << "m2eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "iso")) {
      a->iso = CutValue; ok = 1;
      if (dump) cout << "iso:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "chi2dof")) {
      a->chi2dof = CutValue; ok = 1;
      if (dump) cout << "chi2dof:             " << CutValue << endl;
    }

    if (!strcmp(CutName, "alpha")) {
      a->alpha = CutValue; ok = 1;
      if (dump) cout << "alpha:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "fls3d")) {
      a->fls3d = CutValue; ok = 1;
      if (dump) cout << "fls3d:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "docatrk")) {
      a->docatrk = CutValue; ok = 1;
      if (dump) cout << "docatrk:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "closetrk")) {
      a->closetrk = static_cast<int>(CutValue); ok = 1;
      if (dump) cout << "closetrk:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlip")) {
      a->pvlip = CutValue; ok = 1;
      if (dump) cout << "pvlip:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlips")) {
      a->pvlips = CutValue; ok = 1;
      if (dump) cout << "pvlips:                " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlip2")) {
      a->pvlip2 = CutValue; ok = 1;
      if (dump) cout << "pvlip2:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlips2")) {
      a->pvlips2 = CutValue; ok = 1;
      if (dump) cout << "pvlips2:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "maxdoca")) {
      a->maxdoca = CutValue; ok = 1;
      if (dump) cout << "maxdoca:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvip")) {
      a->pvip = CutValue; ok = 1;
      if (dump) cout << "pvip:                    " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvips")) {
      a->pvips = CutValue; ok = 1;
      if (dump) cout << "pvips:                   " << CutValue << endl;
    }

    sscanf(buffer, "%s %s", CutName, XmlName);
    string ctmp = CutName; 
    string sXmlName;
    replaceAll(ctmp, " ", ""); 
    if (!strcmp(ctmp.c_str(), "xml")) {
      a->xmlFile = XmlName; ok = 1; 
      sXmlName = "weights/" + a->xmlFile + "-Events0_BDT.weights.xml"; 
      //      fReaderEvents0.push_back(setupReader(sXmlName, frd)); 
      TMVA::Reader *ar = setupReader(sXmlName, frd); 
      fReaderEvents0[a->index] = ar;
      if (dump) cout << "xml:                   " << sXmlName << endl;
      sXmlName = "weights/" + a->xmlFile + "-Events1_BDT.weights.xml"; 
      //      fReaderEvents1.push_back(setupReader(sXmlName, frd)); 
      ar = setupReader(sXmlName, frd); 
      fReaderEvents1[a->index] = ar;
      if (dump) cout << "xml:                   " << sXmlName << endl;
      sXmlName = "weights/" + a->xmlFile + "-Events2_BDT.weights.xml"; 
      //      fReaderEvents2.push_back(setupReader(sXmlName, frd)); 
      ar = setupReader(sXmlName, frd); 
      fReaderEvents2[a->index] = ar;
      if (dump) cout << "xml:                   " << sXmlName << endl;
    }

    if (!ok) cout << "==> what about " << CutName << endl;
  }

  if (a) fCuts.push_back(a); 

  cout << "==> finished reading cut setting, fCuts.size() =  " << fCuts.size() << endl;
  
}


// ----------------------------------------------------------------------
void plotClass::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  cout << "printCuts ... fCuts.size() = " << fCuts.size() << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    cuts *a = fCuts[i]; 
    OUT << "# -- channel " << a->index << endl;
    OUT << "index    " << a->index << endl;
    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% -- Cuts for channel " << a->index << endl;

    OUT << "xml      " << Form("%s", a->xmlFile.c_str()) << endl;
    OUT << "bdt      " << Form("%4.3f", a->bdt) << endl;
    OUT << "bdtMax   " << Form("%4.3f", a->bdtMax) << endl;
    fTEX <<  Form("\\vdef{%s:bdt:%d}     {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->bdt) << endl;
    fTEX <<  Form("\\vdef{%s:bdtMax:%d}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->bdtMax) << endl;

    OUT << "mBdLo    " << Form("%4.3f", a->mBdLo) << endl;
    OUT << "mBdHi    " << Form("%4.3f", a->mBdHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBdLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBdHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdHi) << endl;

    OUT << "mBsLo    " << Form("%4.3f", a->mBsLo) << endl;
    OUT << "mBsHi    " << Form("%4.3f", a->mBsHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBsLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBsHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsHi) << endl;
 
    OUT << "etaMin   " << Form("%3.1f", a->etaMin) << endl;
    OUT << "etaMax   " << Form("%3.1f", a->etaMax) << endl;
    fTEX <<  Form("\\vdef{%s:etaMin:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMin) << endl;
    fTEX <<  Form("\\vdef{%s:etaMax:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMax) << endl;

    OUT << "pt       " << Form("%3.1f", a->pt) << endl;
    fTEX <<  Form("\\vdef{%s:pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->pt) << endl;
    OUT << "m1pt     " << a->m1pt << endl;
    fTEX <<  Form("\\vdef{%s:m1pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1pt) << endl;
    OUT << "m2pt     " << a->m2pt << endl;
    fTEX <<  Form("\\vdef{%s:m2pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2pt) << endl;
    OUT << "m1eta    " << a->m1eta << endl;
    fTEX <<  Form("\\vdef{%s:m1eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1eta) << endl;
    OUT << "m2eta    " << a->m2eta << endl;
    fTEX <<  Form("\\vdef{%s:m2eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2eta) << endl;

    OUT << "iso      " << a->iso << endl;
    fTEX <<  Form("\\vdef{%s:iso:%d}   {\\ensuremath{{%3.2f } } }", fSuffix.c_str(), a->index, a->iso) << endl;
    OUT << "chi2dof  " << a->chi2dof << endl;
    fTEX <<  Form("\\vdef{%s:chi2dof:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->chi2dof) << endl;
    OUT << "alpha    " << a->alpha << endl;
    fTEX <<  Form("\\vdef{%s:alpha:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->alpha) << endl;
    OUT << "fls3d    " << a->fls3d << endl;
    fTEX <<  Form("\\vdef{%s:fls3d:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->fls3d) << endl;
    OUT << "docatrk  " << a->docatrk << endl;
    fTEX <<  Form("\\vdef{%s:docatrk:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->docatrk) << endl;

    OUT << "closetrk " << a->closetrk << endl;
    fTEX <<  Form("\\vdef{%s:closetrk:%d}   {\\ensuremath{{%d } } }", fSuffix.c_str(), a->index, static_cast<int>(a->closetrk)) << endl;
    OUT << "pvlip    " << a->pvlip << endl;
    fTEX <<  Form("\\vdef{%s:pvlip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlip) << endl;
    OUT << "pvlips   " << a->pvlips << endl;
    fTEX <<  Form("\\vdef{%s:pvlips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlips) << endl;
    OUT << "pvlip2   " << a->pvlip2 << endl;
    fTEX <<  Form("\\vdef{%s:pvlip2:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlip2) << endl;
    OUT << "pvlips2  " << a->pvlips2 << endl;
    fTEX <<  Form("\\vdef{%s:pvlips2:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlips2) << endl;
    OUT << "maxdoca  " << a->maxdoca << endl;
    fTEX <<  Form("\\vdef{%s:maxdoca:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->maxdoca) << endl;
    OUT << "pvip     " << a->pvip << endl;
    fTEX <<  Form("\\vdef{%s:pvip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvip) << endl;
    OUT << "pvips    " << a->pvips << endl;
    fTEX <<  Form("\\vdef{%s:pvips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvips) << endl;

    OUT << "doApplyCowboyVeto  " << fDoApplyCowboyVeto << endl;
    fTEX <<  Form("\\vdef{%s:doApplyCowboyVeto:%d}   {%s }", 
		  fSuffix.c_str(), a->index, fDoApplyCowboyVeto?"yes":"no") << endl;
    OUT << "fDoApplyCowboyVetoAlsoInSignal  " << fDoApplyCowboyVetoAlsoInSignal << endl;
    fTEX <<  Form("\\vdef{%s:fDoApplyCowboyVetoAlsoInSignal:%d}   {%s }", 
		  fSuffix.c_str(), a->index, fDoApplyCowboyVetoAlsoInSignal?"yes":"no") << endl;
    
  }
  OUT.flush();
}



// ----------------------------------------------------------------------
void plotClass::stamp(double x1, string text1, double x2, string text2) {
  cout << "stamp() > " << x1 << " " << text1 << " " << x2 << " " << text2 << endl;
  tl->SetNDC(kTRUE); 
  tl->SetTextSize(fSize); 
  tl->DrawLatex(x1, 0.91, text1.c_str());   
  tl->DrawLatex(x2, 0.91, text2.c_str()); 
  tl->SetTextSize(0.03); 
  tl->SetTextSize(fSize); 
}



// ----------------------------------------------------------------------
void plotClass::drawArrow(double height, int mode, double ylegend) {

  double ylo(0.01); 
  pl->SetLineWidth(static_cast<Width_t>(3.)); 
  
  double d(0.08), y(0.80), x(5.25); 
  
  double yoffset(0.2); 

  if (1 == mode) {
    pl->SetLineColor(kBlue); 
    pl->SetLineColor(kBlue); 
    pl->SetLineStyle(kSolid); 
    pl->DrawLine(fCuts[0]->mBsLo, height, fCuts[0]->mBsHi, height); 
    pl->SetLineWidth(static_cast<Width_t>(2.)); 
    pl->DrawLine(fCuts[0]->mBsLo, height+d, fCuts[0]->mBsLo, height-d); 
    pl->DrawLine(fCuts[0]->mBsHi, height+d, fCuts[0]->mBsHi, height-d); 

    if (1) {
      y = ylegend;
      x = 5.25; 
      d = 0.05;
      pl->SetLineWidth(static_cast<Width_t>(3.)); 
      pl->DrawLine(x, y, x+0.1, y); 
      pl->SetLineWidth(static_cast<Width_t>(2.)); 
      pl->DrawLine(x, y+d, x, y-d); 
      pl->DrawLine(x+0.1, y+d, x+0.1, y-d); 
      tl->SetNDC(kFALSE);
      tl->SetTextSize(0.05); 
      tl->DrawLatex(x+0.15, y-yoffset, "B_{s}^{0} signal window");
    }

  } else if (2 == mode) {
    pl->SetLineColor(kRed); 
    pl->SetLineColor(kRed); 
    pl->SetLineStyle(kDashed); 
    pl->DrawLine(fCuts[0]->mBdLo, height, fCuts[0]->mBdHi, height); 
    pl->SetLineStyle(kSolid); 
    pl->SetLineWidth(static_cast<Width_t>(2.)); 
    pl->DrawLine(fCuts[0]->mBdLo, height+d, fCuts[0]->mBdLo, height-d); 
    pl->DrawLine(fCuts[0]->mBdHi, height+d, fCuts[0]->mBdHi, height-d); 

    if (1) {
      x = 5.25; 
      y = ylegend;
      d = 0.05;
      pl->SetLineWidth(static_cast<Width_t>(3.)); 
      pl->SetLineStyle(kDashed); 
      pl->DrawLine(x, y, x+0.1, y); 
      pl->SetLineStyle(kSolid); 
      pl->SetLineWidth(static_cast<Width_t>(2.)); 
      pl->DrawLine(x, y+d, x, y-d); 
      pl->DrawLine(x+0.1, y+d, x+0.1, y-d); 
      tl->SetNDC(kFALSE);
      tl->SetTextSize(0.05); 
      tl->DrawLatex(x+0.15, y-yoffset, "B^{0} signal window");
    }

  } else if (3 == mode) {
    pa->SetLineColor(kBlack); 
    pa->SetFillColor(kBlack); 
    pa->DrawArrow(4.9, 0.5, 4.9, ylo); 
    pa->DrawArrow(5.2, 0.5, 5.2, ylo); 
    pl->SetLineColor(kBlack); 
    pl->DrawLine(4.9, 0.5, 5.2, 0.5);

    pa->DrawArrow(5.45, 0.5, 5.45, ylo); 
    pa->DrawArrow(5.90, 0.5, 5.90, ylo); 
    pl->SetLineColor(kBlack); 
    pl->DrawLine(5.45, 0.5, 5.90, 0.5);
  }

  
 
}


// ----------------------------------------------------------------------
void plotClass::drawBox(int mode, double hi, double ylo) {

  TBox *b = new TBox; 

  if (1 == mode) {
    b->SetFillColor(kBlue); 
    b->SetLineColor(kBlue); 
    b->SetFillStyle(3004); 
    b->SetFillStyle(3356); 
    b->DrawBox(fCuts[0]->mBsLo, ylo, fCuts[0]->mBsHi, hi); 
  } else if (2 == mode) {
    b->SetLineColor(kRed); 
    b->SetFillColor(kRed); 
    b->SetFillStyle(3005); 
    b->SetFillStyle(3365); 
    b->DrawBox(fCuts[0]->mBdLo, ylo, fCuts[0]->mBdHi, hi); 
  } else if (3 == mode) {
    b->SetLineColor(kBlack); 
    b->SetFillColor(kBlack); 
    b->DrawBox(4.9, ylo, 5.2, hi); 
    b->DrawBox(5.45, ylo, 5.90, hi); 
  }

  
 
}


// ----------------------------------------------------------------------
void plotClass::singleEventPrintout(string suffix, string st, int ievt) {

  fTEX << formatTex(fb.run,      Form("%s:%s%i:run", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.evt,      Form("%s:%s%i:evt", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fChan,     Form("%s:%s%i:chan", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m,        Form("%s:%s%i:m", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.pt,       Form("%s:%s%i:pt", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.phi,      Form("%s:%s%i:phi", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.eta,      Form("%s:%s%i:eta", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << Form("\\vdef{%s:%s%i:channel}   {%s }", suffix.c_str(), st.c_str(), ievt, fChan==0?"barrel":"endcap") << endl;
  fTEX << formatTex((fb.cb?1:0),    Form("%s:%s%i:cowboy", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m1pt,     Form("%s:%s%i:m1pt", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.m2pt,     Form("%s:%s%i:m2pt", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.m1eta,    Form("%s:%s%i:m1eta", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.m2eta,    Form("%s:%s%i:m2eta", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.m1phi,    Form("%s:%s%i:m1phi", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.m2phi,    Form("%s:%s%i:m2phi", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.m1q,       Form("%s:%s%i:m1q", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m2q,       Form("%s:%s%i:m2q", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.iso,      Form("%s:%s%i:iso", suffix.c_str(), st.c_str(), ievt), 3) << endl;
  fTEX << formatTex(fb.alpha,    Form("%s:%s%i:alpha", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.chi2,     Form("%s:%s%i:chi2", suffix.c_str(), st.c_str(), ievt), 2) << endl;
  fTEX << formatTex(fb.dof,      Form("%s:%s%i:dof", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.fls3d,    Form("%s:%s%i:fls3d", suffix.c_str(), st.c_str(), ievt), 2) << endl;
  fTEX << formatTex(fb.fl3d,     Form("%s:%s%i:fl3d", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.fl3dE,    Form("%s:%s%i:fl3dE", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  
  fTEX << formatTex(fb.docatrk,  Form("%s:%s%i:docatrk", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.closetrk, Form("%s:%s%i:closetrk", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.lip,      Form("%s:%s%i:lip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.lipE,     Form("%s:%s%i:lipE", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.tip,      Form("%s:%s%i:tip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.tipE,     Form("%s:%s%i:tipE", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.pvlip,    Form("%s:%s%i:pvlip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.pvlips,   Form("%s:%s%i:pvlips", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.pvip,     Form("%s:%s%i:pvip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.pvips,    Form("%s:%s%i:pvips", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.maxdoca,  Form("%s:%s%i:maxdoca", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fb.pvw8,     Form("%s:%s%i:pvw8", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  fTEX << formatTex(fBDT,        Form("%s:%s%i:bdt", suffix.c_str(), st.c_str(), ievt), 4) << endl;
  
  fTEX << formatTex(fb.m1pix,    Form("%s:%s%i:m1pix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m2pix,    Form("%s:%s%i:m2pix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m1bpix,   Form("%s:%s%i:m1bpix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m2bpix,   Form("%s:%s%i:m2bpix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m1bpixl1, Form("%s:%s%i:m1bpixl1", suffix.c_str(), st.c_str(), ievt), 0) << endl;
  fTEX << formatTex(fb.m2bpixl1, Form("%s:%s%i:m2bpixl1", suffix.c_str(), st.c_str(), ievt), 0) << endl;
}


// ----------------------------------------------------------------------
void plotClass::loopOverTree(TTree *t, std::string mode, int function, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries); 
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0; 
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart; 
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts; 
    } else {
      nend = t->GetEntries();
    }
  }

  nentries = nend - nstart; 

  int step(1000000); 
  if (nentries < 5000000)  step = 500000; 
  if (nentries < 1000000)  step = 100000; 
  if (nentries < 100000)   step = 10000; 
  if (nentries < 10000)    step = 1000; 
  if (nentries < 1000)     step = 100; 
  step = 500000; 
  int imode(0); 
  if (string::npos != mode.find("No")) imode = 10; 
  if (string::npos != mode.find("Cs")) imode = 20; 
  if (string::npos != mode.find("Bs")) imode = 0; 
  if (string::npos != mode.find("Bd")) imode = 1; 
  if (string::npos != mode.find("Bg")) imode = 98; 
  cout << "==> plotClass::loopOverTree> looping in mode " << mode << " -> imode = " << imode 
       << " fDoUseBDT = " << fDoUseBDT
       << " with " << nentries << " entries" 
       << endl;

  string tname; 
  if (98 == imode) {
    tname = fRareName;
  } else {
    tname = mode; 
  }

  TDirectory *dir(0); 
  TTree *small(0); 
  TFile *fLocal(0); 
  if (fSaveSmallTree) {
    dir = gDirectory; 
    if (fSaveLargerTree) {
      fLocal = TFile::Open(Form("%s/larger-%s.root", fDirectory.c_str(), tname.c_str()), "RECREATE"); 
      small = new TTree("t", "t");
    } else {
      fLocal = TFile::Open(Form("%s/small-%s.root", fDirectory.c_str(), tname.c_str()), "RECREATE"); 
      small = new TTree(Form("%s_%s", tname.c_str(), (fDoUseBDT?"bdt":"cnc")), Form("%s_%s", tname.c_str(), (fDoUseBDT?"bdt":"cnc")));
    }
    small->SetDirectory(fLocal);
    small->Branch("run", &fb.run,"run/I");
    small->Branch("evt", &fb.evt,"evt/I");
    small->Branch("ls", &fb.ls,"ls/I");
    small->Branch("m1bdt", &fb.m1rmvabdt ,"m1bdt/D");
    small->Branch("m2bdt", &fb.m2rmvabdt ,"m2bdt/D");
    small->Branch("muid", &fb.gmuid ,"muid/O");
    
    small->Branch("bdt", &fBDT ,"bdt/D");
    if (fSaveLargerTree) {
      small->Branch("hlt", &fb.hlt ,"hlt/O");
      small->Branch("hltm", &fb.hltm ,"hltm/O");
      small->Branch("hltm2", &fb.hltm2 ,"hltm2/O");
      small->Branch("gtqual", &fb.gtqual ,"gtqual/O");
      small->Branch("pt",   &fb.pt ,"hlt/D");
      small->Branch("eta",  &fb.eta ,"eta/D");
      small->Branch("m1pt", &fb.m1pt ,"hlt/D");
      small->Branch("m1q", &fb.m1q ,"m1q/I");
      small->Branch("m2pt", &fb.m2pt ,"hlt/D");
      small->Branch("m2q", &fb.m2q ,"m2q/I");
      small->Branch("pchi2dof", &fb.pchi2dof ,"pchi2dof/D");
      small->Branch("chi2", &fb.chi2 ,"chi2/D");
      small->Branch("dof",  &fb.dof ,"dof/D");
      small->Branch("chi2dof", &fb.chi2dof ,"chi2dof/D");
      small->Branch("fl3d", &fb.fl3d ,"fl3d/D");
      small->Branch("fls3d", &fb.fls3d ,"fls3d/D");
      small->Branch("flsxy", &fb.flsxy ,"flsxy/D");
      small->Branch("pvlip", &fb.pvlip ,"pvlip/D");
      small->Branch("pvlips", &fb.pvlips ,"pvlips/D");
      small->Branch("pvlip2", &fb.pvlip2 ,"pvlip2/D");
      small->Branch("pvlips2", &fb.pvlips2 ,"pvlips2/D");
      small->Branch("pvip", &fb.pvip ,"pvip/D");
      small->Branch("pvips", &fb.pvips ,"pvips/D");
      small->Branch("pvn", &fb.pvn ,"pvn/I");
      small->Branch("pvw8", &fb.pvw8 ,"pvw8/D");
      small->Branch("gtqual", &fb.gtqual ,"gtqual/O");
      small->Branch("q", &fb.q ,"q/I");
      small->Branch("iso", &fb.iso ,"iso/D");
      small->Branch("m1iso", &fb.m1iso ,"m1iso/D");
      small->Branch("m2iso", &fb.m1iso ,"m2iso/D");
      small->Branch("alpha", &fb.alpha ,"alpha/D");
      small->Branch("closetrk", &fb.closetrk ,"closetrk/I");
      small->Branch("docatrk", &fb.docatrk ,"docatrk/D");
      small->Branch("maxdoca", &fb.maxdoca ,"maxdoca/D");
      small->Branch("lip", &fb.lip ,"lip/D");
      small->Branch("tip", &fb.tip ,"tip/D");
    }
    small->Branch("m",     &fb.m,"m/D");
    small->Branch("m3",    &fb.m3,"m3/D");
    small->Branch("m4",    &fb.m4,"m4/D");
    small->Branch("me",    &fb.me,"me/D");
    small->Branch("tau",   &fb.tau ,"tau/D");
    small->Branch("gtau",  &fb.gtau ,"gtau/D");
    small->Branch("m1eta", &fb.m1eta,"m1eta/D");
    small->Branch("m2eta", &fb.m2eta,"m2eta/D");
    small->Branch("eta",   &fb.eta,"eta/D");
  }

  // -- the real loop starts here
  cout << "loopOverTree: nevts = " << nentries << " nstart = " << nstart << " nbegin: " << nbegin << " nend: " << nend << endl;
  bool ptcutOK(false); 
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. Event %8d, run = %6lld evt = %10lld", jentry, 
				       static_cast<Long64_t>(fb.run), static_cast<Long64_t>(fb.evt)) << endl;

    if (fRunMax > fRunMin) {
      if (fb.run < fRunMin) continue; 
      if (fb.run > fRunMax) continue; 
    }
    candAnalysis(imode);
    loopFunction(function, imode);
    ptcutOK = fGoodMuonsPt;
    if (fDoUseBDT) ptcutOK = fGoodBdtPt; 

    if (fSaveSmallTree
	&& fGoodAcceptance 
	&& fGoodQ
	&& fGoodPvAveW8
	&& fGoodTracks 
	&& fGoodTracksPt 
	&& fGoodTracksEta 
	&& ptcutOK
	&& fGoodMuonsEta
	&& fGoodJpsiCuts
	// 	&& (fSaveLargerTree || fGoodMuonsID)
	&& (fSaveLargerTree || fGoodHLT)
	) {
      small->Fill();
    }
  }

  if (fSaveSmallTree) {
    small->Write();
    fLocal->Close();
    dir->cd();
  }
}



// ----------------------------------------------------------------------
TTree* plotClass::getTree(string mode) {
  TTree *t(0);
  cout << "retrieve tree events for mode " << mode << " from file " << fF[mode]->GetName() << endl;
  t = (TTree*)fF[mode]->Get("candAnaMuMu/events"); 
  if (string::npos != mode.find("No")) t = (TTree*)fF[mode]->Get("candAnaBu2JpsiK/events"); 
  if (string::npos != mode.find("Cs")) t = (TTree*)fF[mode]->Get("candAnaBs2JpsiPhi/events"); 
  if (string::npos != mode.find("Sg")) t = (TTree*)fF[mode]->Get("candAnaMuMu/events"); 
  if (string::npos != mode.find("Bd")) t = (TTree*)fF[mode]->Get("candAnaMuMu/events"); 
  return t; 
}

// ----------------------------------------------------------------------
void plotClass::checkAgainstDuplicates(string mode) {

  TTree *t = getTree(mode); 
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("run",1);
  t->SetBranchStatus("evt",1);
  t->SetBranchStatus("m1pt",1);
  t->SetBranchStatus("m2pt",1);

  Long64_t lrun,  levt; 
  double   lm1pt, lm2pt; 


  t->SetBranchAddress("run",&lrun);
  t->SetBranchAddress("evt",&levt);
  t->SetBranchAddress("m1pt",&lm1pt);
  t->SetBranchAddress("m2pt",&lm2pt);

  std::map<string, int> aMap; 
  string ckey; 

  int nentries = Int_t(t->GetEntries());
  cout << "looking at " << nentries << " events" << endl;

  int step(50000); 
  if (nentries < 50000000) step = 5000000; 
  if (nentries < 10000000) step = 1000000; 
  if (nentries < 5000000)  step = 500000; 
  if (nentries < 1000000)  step = 100000; 
  if (nentries < 100000)   step = 10000; 
  if (nentries < 10000)    step = 1000; 
  if (nentries < 1000)     step = 100; 

  for (int jentry = 0; jentry < nentries; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. Event %8d", jentry) << endl;
    
    ckey = Form("%lld:%lld:%5.4f:%5.4f", lrun, levt, lm1pt, lm2pt); 
    aMap.insert(make_pair(ckey, 1));
  }

  
  for (map<string, int>::iterator imap = aMap.begin(); imap != aMap.end(); ++imap) {  
    if (imap->second > 1) {
      cout << "duplicate event: n = " << imap->second << " info: " << imap->first << endl;
    }
  }

}


// ----------------------------------------------------------------------
void plotClass::setupTree(TTree *t, string mode) {

  if (string::npos != mode.find("Mc")) {
    fIsMC = true; 
  } else {
    fIsMC = false; 
  }

  t->SetBranchAddress("pt", &fb.pt);
  t->SetBranchAddress("q", &fb.q);

  t->SetBranchAddress("tau", &fb.tau);
  t->SetBranchAddress("gtau", &fb.gtau);

  t->SetBranchAddress("bdt",&fb.bdt);
  t->SetBranchAddress("bdt",&fBDT);
  t->SetBranchAddress("lip",&fb.lip);
  t->SetBranchAddress("lipE",&fb.lipE);
  t->SetBranchAddress("tip",&fb.tip);
  t->SetBranchAddress("tipE",&fb.tipE);

  t->SetBranchAddress("closetrk",&fb.closetrk);
  t->SetBranchAddress("pvlip",   &fb.pvlip);
  t->SetBranchAddress("pvlips",  &fb.pvlips);
  t->SetBranchAddress("pvlip2",  &fb.pvlip2);
  t->SetBranchAddress("pvlips2", &fb.pvlips2);
  t->SetBranchAddress("maxdoca", &fb.maxdoca);
  t->SetBranchAddress("pvip",    &fb.pvip);
  t->SetBranchAddress("pvips",   &fb.pvips);
  t->SetBranchAddress("pvip3d",  &fb.pvip3d);
  t->SetBranchAddress("pvips3d", &fb.pvips3d);
  t->SetBranchAddress("pvw8",    &fb.pvw8);

  t->SetBranchAddress("m1pix",    &fb.m1pix);
  t->SetBranchAddress("m2pix",    &fb.m2pix);
  t->SetBranchAddress("m1bpix",   &fb.m1bpix);
  t->SetBranchAddress("m2bpix",   &fb.m2bpix);
  t->SetBranchAddress("m1bpixl1", &fb.m1bpixl1);
  t->SetBranchAddress("m2bpixl1", &fb.m2bpixl1);

  t->SetBranchAddress("rr",     &fb.rr);
  t->SetBranchAddress("pvn",    &fb.pvn);
  t->SetBranchAddress("run",    &fb.run);
  t->SetBranchAddress("evt",    &fb.evt);
  t->SetBranchAddress("hlt",    &fb.hlt);
  t->SetBranchAddress("hltm",   &fb.hltm);
  t->SetBranchAddress("hltm2",  &fb.hltm2);
  t->SetBranchAddress("ls",     &fb.ls);
  t->SetBranchAddress("cb",     &fb.cb);
  t->SetBranchAddress("json",   &fb.json);
  t->SetBranchAddress("gmuid",  &fb.gmuid);
  t->SetBranchAddress("gmumvaid", &fb.gmumvaid);
  t->SetBranchAddress("gtqual", &fb.gtqual);
  t->SetBranchAddress("tm",     &fb.tm);
  t->SetBranchAddress("procid", &fb.procid);
  t->SetBranchAddress("m",      &fb.m);
  t->SetBranchAddress("m3",     &fb.m3);
  t->SetBranchAddress("m4",     &fb.m4);
  t->SetBranchAddress("me",     &fb.me);
  t->SetBranchAddress("cm",     &fb.cm);
  t->SetBranchAddress("pt",     &fb.pt);
  t->SetBranchAddress("phi",    &fb.phi);
  t->SetBranchAddress("eta",    &fb.eta);
  t->SetBranchAddress("cosa",   &fb.cosa);
  t->SetBranchAddress("alpha",  &fb.alpha);
  t->SetBranchAddress("iso",    &fb.iso);
  t->SetBranchAddress("chi2",   &fb.chi2);
  t->SetBranchAddress("dof",    &fb.dof);
  t->SetBranchAddress("prob",   &fb.pchi2dof);
  t->SetBranchAddress("chi2dof",&fb.chi2dof);
  t->SetBranchAddress("flsxy",  &fb.flsxy);
  t->SetBranchAddress("fls3d",  &fb.fls3d);
  t->SetBranchAddress("fl3d",   &fb.fl3d);
  t->SetBranchAddress("fl3dE",  &fb.fl3dE);
  t->SetBranchAddress("m1pt",   &fb.m1pt);
  t->SetBranchAddress("m1gt",   &fb.m1gt);
  t->SetBranchAddress("m1eta",  &fb.m1eta);
  t->SetBranchAddress("m1phi",  &fb.m1phi);
  t->SetBranchAddress("m1q",    &fb.m1q);
  t->SetBranchAddress("m2pt",   &fb.m2pt);
  t->SetBranchAddress("m2gt",   &fb.m2gt);
  t->SetBranchAddress("m2eta",  &fb.m2eta);
  t->SetBranchAddress("m2phi",  &fb.m2phi);
  t->SetBranchAddress("m2q",    &fb.m2q);
  t->SetBranchAddress("docatrk",&fb.docatrk);

  t->SetBranchAddress("m1id",   &fb.m1id);
  t->SetBranchAddress("m1rmvaid", &fb.m1rmvaid);
  t->SetBranchAddress("m1trigm",  &fb.m1trigm);
  t->SetBranchAddress("m1rmvabdt", &fb.m1rmvabdt);

  t->SetBranchAddress("m2id",     &fb.m2id);
  t->SetBranchAddress("m2rmvaid", &fb.m2rmvaid);
  t->SetBranchAddress("m2trigm",  &fb.m2trigm);
  t->SetBranchAddress("m2rmvabdt", &fb.m2rmvabdt);


  t->SetBranchAddress("m1iso",     &fb.m1iso);
  t->SetBranchAddress("m2iso",     &fb.m2iso);
  t->SetBranchAddress("closetrks1",&fb.closetrks1);
  t->SetBranchAddress("closetrks2",&fb.closetrks2);
  t->SetBranchAddress("closetrks3",&fb.closetrks3);
  t->SetBranchAddress("othervtx",  &fb.othervtx);
  t->SetBranchAddress("pvdchi2",   &fb.pvdchi2);

  t->SetBranchAddress("g1pt",   &fb.g1pt);
  t->SetBranchAddress("g2pt",   &fb.g2pt);
  t->SetBranchAddress("g1eta",  &fb.g1eta);
  t->SetBranchAddress("g2eta",  &fb.g2eta);
  t->SetBranchAddress("g1id",   &fb.g1id);
  t->SetBranchAddress("g2id",   &fb.g2id);
  if (string::npos != mode.find("No")) {
    if (string::npos != mode.find("Mc")) {
      t->SetBranchAddress("g3pt", &fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
    }
    t->SetBranchAddress("kpt",  &fb.k1pt);
    t->SetBranchAddress("kgt",  &fb.k1gt);
    t->SetBranchAddress("keta", &fb.k1eta);
    t->SetBranchAddress("mpsi", &fb.mpsi);
    t->SetBranchAddress("psipt",&fb.psipt); //FIXME
  }

  if (string::npos != mode.find("Cs")) {
    if (string::npos != mode.find("Mc")) {
      t->SetBranchAddress("g3pt", &fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
      t->SetBranchAddress("g4pt", &fb.g4pt);
      t->SetBranchAddress("g4eta",&fb.g4eta);
    }
    t->SetBranchAddress("psipt",&fb.psipt);   //FIXME
    t->SetBranchAddress("mpsi", &fb.mpsi);
    t->SetBranchAddress("mkk",  &fb.mkk);
    t->SetBranchAddress("dr",   &fb.dr);
    t->SetBranchAddress("k1pt", &fb.k1pt);
    t->SetBranchAddress("k1gt", &fb.k1gt);
    t->SetBranchAddress("k1eta",&fb.k1eta);
    t->SetBranchAddress("k2pt", &fb.k2pt);
    t->SetBranchAddress("k2gt", &fb.k2gt);
    t->SetBranchAddress("k2eta",&fb.k2eta);
  } else {
    fb.mkk = 999.;
    fb.dr = 999.;
  }

  if (string::npos != mode.find("DstarPi")) {
    t->SetBranchAddress("md0",&fb.md0);
    t->SetBranchAddress("dm",&fb.dm);
    t->SetBranchAddress("ptd0",&fb.ptd0);
  }

}


// ----------------------------------------------------------------------
void plotClass::candAnalysis(int mode) {

  cuts *pCuts(0); 
  fChan = detChan(fb.m1eta, fb.m2eta); 
  if (fChan < 0) {
    //    cout << "plotClass::candAnalysis: could not determine channel: " << fb.m1eta << " " << fb.m2eta << endl;
    return;
  }
  pCuts = fCuts[fChan]; 

  bool bp2jpsikp(false), bs2jpsiphi(false); 
  if (10 == mode)  bp2jpsikp = true; 
  if (20 == mode)  bs2jpsiphi = true; 

  // -- reset all
  fBDT = -99.; 
  fGoodHLT = fGoodMuonsID = false;
  fGoodQ = fGoodPvAveW8 = fGoodMaxDoca = fGoodIp = fGoodIpS = fGoodPt = fGoodEta = fGoodAlpha =  fGoodChi2 = fGoodFLS = false;   
  fGoodCloseTrack = fGoodIso = fGoodDocaTrk = fGoodLastCut = fPreselection = false;

  fGoodJpsiCuts = true;

  fGoodAcceptance = true; 
  fGoodBdtPt      = true; 
  fGoodMuonsPt    = true;
  fGoodMuonsEta   = true;
  fGoodTracks     = fb.gtqual;
  fGoodTracksPt   = true;
  fGoodTracksEta  = true;

  fIsCowboy = fb.cb;

  if (fIsMC) {
    if (fb.g1pt < fAccPt) fGoodAcceptance = false; 
    if (fb.g2pt < fAccPt) fGoodAcceptance = false; 
    if (TMath::Abs(fb.g1eta) > 2.5) fGoodAcceptance = false; 
    if (TMath::Abs(fb.g2eta) > 2.5) fGoodAcceptance = false; 
  } else {
    if (!fb.json) {
      return;
    }
  }

  if (fb.m1pt < fAccPt) fGoodAcceptance = false; 
  if (fb.m2pt < fAccPt) fGoodAcceptance = false; 
  if (0 == fb.m1gt)  fGoodAcceptance = false; 
  if (0 == fb.m2gt)  fGoodAcceptance = false; 

  if (fb.m1pt < pCuts->bdtpt) {
    fGoodBdtPt = false; 
  }
  if (fb.m2pt < pCuts->bdtpt) {
    fGoodBdtPt = false; 
  }

  if (fb.m1pt < pCuts->m1pt) {
    fGoodMuonsPt = false; 
  }
  if (fb.m2pt < pCuts->m2pt) {
    fGoodMuonsPt = false; 
  }
  if (TMath::Abs(fb.m1eta) > 2.4) {
    fGoodAcceptance = false; 
  }
  if (TMath::Abs(fb.m2eta) > 2.4) {
    fGoodAcceptance = false; 
  }

  if (TMath::Abs(fb.m1eta) > pCuts->m1eta) {
    fGoodMuonsEta = false; 
  }

  if (TMath::Abs(fb.m2eta) > pCuts->m2eta) {
    fGoodMuonsEta = false; 
  }
  
  if (bp2jpsikp) {
    if (fIsMC) {
      // gen-level cuts for Bu2JpsiKp
      if (fb.g1pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (fb.g2pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (TMath::Abs(fb.g3eta) > 2.5) fGoodAcceptance = false; 
      if (fb.g3pt < 0.4) fGoodAcceptance = false; 
    }
    if (TMath::Abs(fb.k1eta) > 2.4) {
      fGoodAcceptance = false; 
      fGoodTracksEta = false; 
    }
    if (fb.k1pt < 0.5) {
      fGoodAcceptance = false; 
      fGoodTracksPt = false; 
    }
    if (0 == fb.k1gt)  fGoodAcceptance = false; 
  }
  
  if (bs2jpsiphi) {
    if (fIsMC) {
      if (TMath::Abs(fb.g3eta) > 2.5) fGoodAcceptance = false; 
      if (TMath::Abs(fb.g4eta) > 2.5) fGoodAcceptance = false; 
      // gen-level cuts for Bs2JpsiPhi
      if (fb.g1pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (fb.g2pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (fb.g3pt < 0.4) fGoodAcceptance = false; 
      if (fb.g4pt < 0.4) fGoodAcceptance = false; 
    }
    if (TMath::Abs(fb.k1eta) > 2.4) {
      fGoodAcceptance = false; 
      fGoodTracksEta = false; 
    }
    if (TMath::Abs(fb.k2eta) > 2.4) {
      fGoodAcceptance = false; 
      fGoodTracksEta = false; 
    }
    if (fb.k1pt < 0.5) {
      fGoodAcceptance = false; 
      fGoodTracksPt = false; 
    }
    if (fb.k2pt < 0.5) {
      fGoodAcceptance = false; 
      fGoodTracksPt = false; 
    }
    if (0 == fb.k1gt)  fGoodAcceptance = false; 
    if (0 == fb.k2gt)  fGoodAcceptance = false; 

    if (fb.dr   > 0.3) fGoodJpsiCuts = false; 
    if (fb.mkk  < 0.995) fGoodJpsiCuts = false; 
    if (fb.mkk  > 1.045) fGoodJpsiCuts = false; 
  }

  if (bs2jpsiphi || bp2jpsikp) {
    if (fb.mpsi > 3.2) fGoodJpsiCuts = false;
    if (fb.mpsi < 3.0) fGoodJpsiCuts = false;
    if (fb.psipt < 7.0) fGoodJpsiCuts = false;
  } else {
    fGoodJpsiCuts = true; 
  }

  if (fDoUseBDT) {
    if (fGoodAcceptance 
	&& fGoodTracks
	&& fGoodTracksPt
	&& fGoodTracksEta
	&& fGoodBdtPt
	&& fGoodMuonsEta
	&& fGoodJpsiCuts
	) {
      calcBDT(); 
    }
  }

  fGoodMuonsID  = fb.gmuid;

  fW8 = 1.;
  fW8MmuID = fW8Mtrig = fW8DmuID = fW8Dtrig = -1.;
  double w1(-1.), w2(-1.); 
  
  if (fIsMC) {
    PidTable *pT, *pT1, *pT2; 

    // -- Weights with data PidTables
    if (fIsCowboy) {
      pT  = fptCbM; 
      pT1 = fptCbT1;
      pT2 = fptCbT2;
    } else {
      pT  = fptSgM; 
      pT1 = fptSgT1;
      pT2 = fptSgT2;
    }

    double am1eta = TMath::Abs(fb.m1eta);
    double am2eta = TMath::Abs(fb.m2eta);

    w1       = pT->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8DmuID = w1*w2; 
    
    w1       = pT1->effD(fb.m1pt, am1eta, fb.m1phi) * pT2->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT1->effD(fb.m2pt, am2eta, fb.m2phi) * pT2->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8Dtrig = w1*w2; 

    // -- Weights with MC PidTables
    if (fIsCowboy) {
      pT  = fptCbMMC; 
      pT1 = fptCbT1MC;
      pT2 = fptCbT2MC;
    } else {
      pT  = fptSgMMC; 
      pT1 = fptSgT1MC;
      pT2 = fptSgT2MC;
    }

    w1       = pT->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8MmuID = w1*w2; 
    
    w1       = pT1->effD(fb.m1pt, am1eta, fb.m1phi) * pT2->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT1->effD(fb.m2pt, am2eta, fb.m2phi) * pT2->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8Mtrig = w1*w2; 

    if (98 == mode) {
      w1 = w2 = -1.;
      // -- track 1
      if (  13 == fb.g1id) w1 = pT->effD(fb.m1pt, am1eta, 1.);
      if ( -13 == fb.g1id) w1 = pT->effD(fb.m1pt, am1eta, 1.);

      if (  321 == fb.g1id) w1 = fptFakePosKaons->effD(fb.m1pt, am1eta, 1.);
      if ( -321 == fb.g1id) w1 = fptFakeNegKaons->effD(fb.m1pt, am1eta, 1.);

      if (  211 == fb.g1id) w1 = fptFakePosPions->effD(fb.m1pt, am1eta, 1.);
      if ( -211 == fb.g1id) w1 = fptFakeNegPions->effD(fb.m1pt, am1eta, 1.);

      if ( 2212 == fb.g1id) w1 = fptFakePosProtons->effD(fb.m1pt, am1eta, 1.);
      if (-2212 == fb.g1id) w1 = fptFakeNegProtons->effD(fb.m1pt, am1eta, 1.);

      // -- track 2
      if (  13 == fb.g2id) w2 = pT->effD(fb.m2pt, am2eta, 1.);
      if ( -13 == fb.g2id) w2 = pT->effD(fb.m2pt, am2eta, 1.);

      if (  321 == fb.g2id) w2 = fptFakePosKaons->effD(fb.m2pt, am2eta, 1.);
      if ( -321 == fb.g2id) w2 = fptFakeNegKaons->effD(fb.m2pt, am2eta, 1.);

      if (  211 == fb.g2id) w2 = fptFakePosPions->effD(fb.m2pt, am2eta, 1.);
      if ( -211 == fb.g2id) w2 = fptFakeNegPions->effD(fb.m2pt, am2eta, 1.);

      if ( 2212 == fb.g2id) w2 = fptFakePosProtons->effD(fb.m2pt, am2eta, 1.);
      if (-2212 == fb.g2id) w2 = fptFakeNegProtons->effD(fb.m2pt, am2eta, 1.);

      fW8MisId = w1*w2; 
    }

  }
  
  fGoodQ          = (fb.m1q*fb.m2q < 0); 
  fGoodPvAveW8    = (fb.pvw8 > 0.7);
  fGoodMaxDoca    = (TMath::Abs(fb.maxdoca) < pCuts->maxdoca); 
  fGoodIp         = (TMath::Abs(fb.pvip) < pCuts->pvip); 
  fGoodIpS        = (TMath::Abs(fb.pvips) < pCuts->pvips); 

  fGoodLip        = (TMath::Abs(fb.pvlip) < pCuts->pvlip); 
  fGoodLipS       = (TMath::Abs(fb.pvlips) < pCuts->pvlips); 
  
  fGoodPt         = (fb.pt > pCuts->pt);
  fGoodEta        = ((fb.eta > -24.0) && (fb.eta < 24.0)); 
  fGoodAlpha      = (fb.alpha < pCuts->alpha); 
  fGoodChi2       = (fb.chi2/fb.dof < pCuts->chi2dof);
  fGoodFLS        = (fb.fls3d > pCuts->fls3d);
  if (TMath::IsNaN(fb.fls3d)) fGoodFLS = false;
  
  fGoodCloseTrack = (fb.closetrk < pCuts->closetrk); 
  fGoodIso        = (fb.iso > pCuts->iso); 
  fGoodDocaTrk    = (fb.docatrk > pCuts->docatrk);
  fGoodLastCut    = true; 

  fGoodHLT        = fb.hlt && fb.hltm2;

  // -- no trigger matching for rare decays!
  if (98 == mode) fGoodHLT = fb.hlt; 

  fPreselection   = ((fBDT > 0.) && fGoodHLT && fGoodMuonsID ); 

  fAnaCuts.update(); 

}


// ----------------------------------------------------------------------
void plotClass::calcBDT() {
  fBDT = -99.;

  if (!preselection(fb, fChan)) return;

  if (fDoApplyMuonPtCuts) {
    if (1) cout << "applying muon pT cuts?!" << endl;
    if (fb.m1pt < fCuts[fChan]->m1pt) return;
    if (fb.m2pt < fCuts[fChan]->m2pt) return;
  }

  //??  if (5 == mode && 5.2 < mass && mass < 5.45 && fb.iso < 0.7) continue; 
  //  if (rejectInvIso && 5.2 < fb.m && fb.m < 5.45 && fb.iso < 0.7) return;
  //   if (fb.pt > 100) return;
  //   if (fb.pt < 6) return;
  //   if (fb.m1pt < 4) return;
  //   if (fb.m2pt < 4) return;
  //   if (fb.fl3d > 1.5) return;
  //   if (fb.m > 5.9) return;
  //   if (fb.m < 4.9) return;
  
  //   if (!fb.hlt) return;
  //   if (!fb.gmuid) return;
  
  frd.pt = fb.pt; 
  frd.eta = fb.eta; 
  frd.m1eta = fb.m1eta; 
  frd.m2eta = fb.m2eta; 
  frd.m1pt = fb.m1pt; 
  frd.m2pt = fb.m2pt;
  frd.fls3d = fb.fls3d; 
  frd.alpha = fb.alpha; 
  frd.maxdoca = fb.maxdoca;
  frd.pvip = fb.pvip; 
  frd.pvips = fb.pvips; 
  frd.iso = fb.iso; 
  frd.docatrk = fb.docatrk; 
  frd.chi2dof = fb.chi2dof; 
  frd.closetrk = fb.closetrk; 

  frd.m1iso = fb.m1iso; 
  frd.m2iso = fb.m2iso; 

  frd.closetrks1 = fb.closetrks1; 
  frd.closetrks2 = fb.closetrks2; 
  frd.closetrks3 = fb.closetrks3; 

  frd.pvlip2  = fb.pvlip2; 
  frd.pvlips2 = fb.pvlips2; 
  
  frd.m  = fb.m; 
  int remainder = TMath::Abs(fb.evt%3);
  if (0 == remainder) {
    fBDT   = fReaderEvents0[fChan]->EvaluateMVA("BDT"); 
  } else if (1 == remainder) {
    fBDT   = fReaderEvents1[fChan]->EvaluateMVA("BDT"); 
  } else if (2 == remainder) {
    fBDT   = fReaderEvents2[fChan]->EvaluateMVA("BDT"); 
  } else {
    cout << "all hell break loose" << endl;
  }
}




// ----------------------------------------------------------------------
double plotClass::recalcMass(double m1, double m2) {
  
  TLorentzVector p1, p2; 
  p1.SetPtEtaPhiM(fb.k1pt, fb.k1eta, fb.k1phi, m1);
  p2.SetPtEtaPhiM(fb.k2pt, fb.k2eta, fb.k2phi, m2);
  
  TLorentzVector p0 = p1 + p2; 
  
  return p0.M();
  
}


// ----------------------------------------------------------------------
void plotClass::reduceTree(TTree *t) {
  t->SetBranchStatus("*",0);
  t->SetBranchStatus("run", 1); 
  t->SetBranchStatus("evt", 1); 
  t->SetBranchStatus("bdt", 1); 
  t->SetBranchStatus("m", 1); 
  t->SetBranchStatus("m1eta", 1); 
  t->SetBranchStatus("m2eta", 1); 
  t->SetBranchStatus("eta", 1); 
}


// ----------------------------------------------------------------------
TMVA::Reader* plotClass::setupReader(string xmlFile, readerData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  
  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
  //  TString methodName = TString(fBdt) + TString(" method");
  //  TString weightfile = dir + fBdt + "_" + methodNameprefix + TString(".weights.xml");
  TString weightfile = xmlFile;

  // -- read in variables from weight file
  vector<string> allLines; 
  char  buffer[2000];
  cout << "setupReader, open file " << weightfile << endl;
  ifstream is(weightfile); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1); 
  string::size_type m1, m2;
  string stype; 
  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	//	cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
	if (stype == "m1pt") {
	  cout << "  adding m1pt" << endl;
	  reader->AddVariable( "m1pt", &rd.m1pt);
	}
	if (stype == "m2pt") {
	  cout << "  adding m2pt" << endl;
	  reader->AddVariable( "m2pt", &rd.m2pt);
	}
	if (stype == "m1eta") {
	  cout << "  adding m1eta" << endl;
	  reader->AddVariable( "m1eta", &rd.m1eta);
	}
	if (stype == "m2eta") {
	  reader->AddVariable( "m2eta", &rd.m2eta);
	  cout << "  adding m2eta" << endl;
	}
	if (stype == "pt") {
	  cout << "  adding pt" << endl;
	  reader->AddVariable( "pt", &rd.pt);
	}
	if (stype == "eta") {
	  cout << "  adding eta" << endl;
	  reader->AddVariable( "eta", &rd.eta);
	}
	if (stype == "fls3d") {
	  cout << "  adding fls3d" << endl;
	  reader->AddVariable( "fls3d", &rd.fls3d);
	}
	if (stype == "alpha") {
	  cout << "  adding alpha" << endl;
	  reader->AddVariable( "alpha", &rd.alpha);
	}
	if (stype == "maxdoca") {
	  cout << "  adding maxdoca" << endl;
	  reader->AddVariable( "maxdoca", &rd.maxdoca);
	}
	if (stype == "pvip") {
	  cout << "  adding pvip" << endl;
	  reader->AddVariable( "pvip", &rd.pvip);
	}
	if (stype == "pvips") {
	  cout << "  adding pvips" << endl;
	  reader->AddVariable( "pvips", &rd.pvips);
	}
	if (stype == "iso") {
	  cout << "  adding iso" << endl;
	  reader->AddVariable( "iso", &rd.iso);
	}
	if (stype == "docatrk") {
	  cout << "  adding docatrk" << endl;
	  reader->AddVariable( "docatrk", &rd.docatrk);
	}
	if (stype == "closetrk") {
	  cout << "  adding closetrk" << endl;
	  reader->AddVariable( "closetrk", &rd.closetrk);
	}
	if (stype == "chi2dof") {
	  cout << "  adding chi2dof" << endl;
	  reader->AddVariable( "chi2dof", &rd.chi2dof);
	}
	if (stype == "closetrks1") {
	  cout << "  adding closetrks1" << endl;
	  reader->AddVariable( "closetrks1", &rd.closetrks1);
	}
	if (stype == "closetrks2") {
	  cout << "  adding closetrks2" << endl;
	  reader->AddVariable( "closetrks2", &rd.closetrks2);
	}
	if (stype == "closetrks3") {
	  cout << "  adding closetrks3" << endl;
	  reader->AddVariable( "closetrks3", &rd.closetrks3);
	}
	if (stype == "m1iso") {
	  cout << "  adding m1iso" << endl;
	  reader->AddVariable( "m1iso", &rd.m1iso);
	}
	if (stype == "m2iso") {
	  cout << "  adding m2iso" << endl;
	  reader->AddVariable( "m2iso", &rd.m2iso);
	}
	if (stype == "othervtx") {
	  cout << "  adding othervtx" << endl;
	  reader->AddVariable( "othervtx", &rd.othervtx);
	}
	if (stype == "pvdchi2") {
	  cout << "  adding pvdchi2" << endl;
	  reader->AddVariable( "pvdchi2", &rd.pvdchi2);
	}
	if (stype == "pvlip2") {
	  cout << "  adding pvlip2" << endl;
	  reader->AddVariable( "pvlip2", &rd.pvlip2);
	}
	if (stype == "pvlips2") {
	  cout << "  adding pvlips2" << endl;
	  reader->AddVariable( "pvlips2", &rd.pvlips2);
	}
      }
      break;
    }
  }
  
  nvars = -1; 
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add spectators
    if (string::npos != allLines[i].find("Spectators NSpec")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      //      cout << "==> " << stype << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10; 
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2); 
	cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
	if (stype == "m") {
	  cout << "  adding m as spectator" << endl;
	  reader->AddSpectator( "m", &rd.m);  
	}
      }
      break;
    }
  }

  // --- Book the MVA methods
  reader->BookMVA("BDT", weightfile); 
  return reader; 
}

// ----------------------------------------------------------------------
int plotClass::detChan(double m1eta, double m2eta) {
  // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
  if (TMath::Abs(m1eta) < fCuts[0]->etaMax && TMath::Abs(m2eta) < fCuts[0]->etaMax) return 0; 
  if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1; 
  return -1; 
}


// ----------------------------------------------------------------------
double plotClass::getValueByLabel(TH1D *h, string label) {
  string axislabel; 
  if (h) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      axislabel = h->GetXaxis()->GetBinLabel(i);
      if (string::npos != axislabel.find(label)) return h->GetBinContent(i);
    }
  }
  return -999.; 
}


// ----------------------------------------------------------------------=
void plotClass::rmSubString(string &sInput, const string &sub) {
  string::size_type foundpos = sInput.find(sub);
  if (foundpos != string::npos)  sInput.erase(sInput.begin() + foundpos, sInput.begin() + foundpos + sub.length());
}

// ----------------------------------------------------------------------=
void plotClass::rmPath(string &sInput) {
  while(string::size_type foundpos = sInput.find("/")) {
    if (foundpos != string::npos)  sInput.erase(sInput.begin(), sInput.begin() + foundpos);
  }
  sInput.erase(sInput.begin(), sInput.begin() + 1); // delete also leading /
}


// ----------------------------------------------------------------------
double plotClass::quadraticSum(int n, ...) {
  va_list vl;
  va_start(vl, n);
  double a(0.), sum(0.); 
  for (int i = 0; i < n; ++i) {
    a = va_arg(vl, double);
    sum += a*a;
  }

  va_end(vl);
  return TMath::Sqrt(sum); 
}

// ----------------------------------------------------------------------
void plotClass::saveHist(TH1* h, string add2name) {
  if (h->InheritsFrom(TH1::Class())) {
    TH1D *h1 = (TH1D*)(h->Clone(Form("%s", add2name.c_str()))); 
    h1->SetDirectory(fHistFile); 
    h1->Write(); 
  } else if (h->InheritsFrom(TH2::Class())) {
    TH2D *h1 = (TH2D*)(h->Clone(Form("%s", add2name.c_str()))); 
    h1->SetDirectory(fHistFile); 
    h1->Write(); 
  }
}


// ----------------------------------------------------------------------
double plotClass::getMaximum(TH1 *h1, TH1 *h2) {
  if (0 == h1 || 0 == h2) return -1.;
  double max1 = h1->GetMaximum(); 
  double max2 = h2->GetMaximum(); 
  if (max2 > max1) return max2; 
  return max1; 

}
