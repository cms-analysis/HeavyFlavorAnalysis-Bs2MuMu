#include "plotClass.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/util.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/functions.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "../macros/AnalysisDistribution.hh"

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

  TVirtualFitter::SetMaxIterations(50000);

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  // PDG 2010:
  fu  = 0.401;
  fs  = 0.113;
  // -- CMS with PDG input
  fsfu = 0.282;
  fsfuE = 0.037/0.282;
  // -- CMS with LHCb input
  fsfu = 0.267;
  fsfuE = 0.021;
  

  fDoUseBDT = false; 
  fDoPrintSingleEvent = false;

  fDoApplyCowboyVeto = false;   
  fDoApplyCowboyVetoAlsoInSignal = false; 
  fInvertedIso = false; 
  fNormProcessed = false; 

  fStampString = "CMS, 5fb^{-1}"; 
  fStampString = "BDT preliminary"; 
  fStampCms = "#sqrt{s} = 7 TeV";

  fCutsFileName = cuts; 
  init(files, cuts, dir, mode);

  int NBINS = (fMassHi - fMassLo)/0.025;

  int HBINS(15); 
  double HLO(0.), HHI(45.); 

  TH1D *h; 
  TH2D *h2; 

  hMassPiPi = new TH1D("hMassPiPi", "hMassPiPi", 100, 0., 0.6);
  hMassKPi  = new TH1D("hMassKPi",  "hMassKPi", 100, 0.5, 1.0);
  hMassPiK  = new TH1D("hMassPiK",  "hMassPiK", 100, 0.9, 1.4);

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
    
    
    h2 = new TH2D(Form("hBdtMass%d", i), Form("hBdtMass%d", i), NBINS, fMassLo, fMassHi, 100, 0., 1.0);
    fhBdtMass.push_back(h2); 

    h = new TH1D(Form("hMassWithMassCuts%d", i), Form("hMassWithMassCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMassCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithMassCutsManyBins%d", i), Form("fhMassWithMassCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithMassCutsManyBins.push_back(h); 
  
    h = new TH1D(Form("hMassWithAllCuts%d", i), Form("hMassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCuts.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsBlind%d", i), Form("hMassWithAllCutsBlind%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCutsBlind.push_back(h); 
    h = new TH1D(Form("hMassWithAllCutsManyBins%d", i), Form("hMassWithAllCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithAllCutsManyBins.push_back(h); 

    h = new TH1D(Form("hNorm%d", i), Form("hNorm%d", i), 100, 4.9, 5.9);
    fhNorm.push_back(h); 

    h = new TH1D(Form("hDstarPi%d", i), Form("hDstarPi%d", i), 20, 4.9, 5.9);
    fhDstarPi.push_back(h); 

    h = new TH1D(Form("hNormC%d", i), Form("hNormC%d", i), 200, 4.9, 5.9);
    fhNormC.push_back(h); 


    h = new TH1D(Form("hMassWithTriggerCuts%d", i), Form("hMassWithTriggerCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithTriggerCuts.push_back(h); 
    h = new TH1D(Form("hMassWithTriggerCutsManyBins%d", i), Form("hMassWithTriggerCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithTriggerCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassWithMuonCuts%d", i), Form("hMassWithMuonCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMuonCuts.push_back(h); 
    h = new TH1D(Form("hMassWithMuonCutsManyBins%d", i), Form("hMassWithMuonCutsManyBins%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithMuonCutsManyBins.push_back(h); 

    h = new TH1D(Form("fhMassWithAnaCuts%d", i), Form("hMassChan%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAnaCuts.push_back(h); 
    h = new TH1D(Form("fhMassWithAnaCutsManyBins%d", i), Form("hMassChan%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
    fhMassWithAnaCutsManyBins.push_back(h); 

    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h); 
    h = new TH1D(Form("fhMassNoCutsManyBins%d", i), Form("hMassChan%d", i), (fMassHi-fMassLo)*1000, fMassLo, fMassHi);
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

    h = new TH1D(Form("h0TNPTrigger%d", i), Form("hTNPTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0TNPTrigger.push_back(h); 
    h = new TH1D(Form("h1TNPTrigger%d", i), Form("hTNPTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPTrigger.push_back(h); 
    h = new TH1D(Form("h0TNPMuID%d", i), Form("hTNPMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0TNPMuID.push_back(h); 
    h = new TH1D(Form("h1TNPMuID%d", i), Form("hTNPMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPMuID.push_back(h); 

    h = new TH1D(Form("h0TNPMCTrigger%d", i), Form("hTNPMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0TNPMCTrigger.push_back(h); 
    h = new TH1D(Form("h1TNPMCTrigger%d", i), Form("hTNPMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPMCTrigger.push_back(h); 
    h = new TH1D(Form("h0TNPMCMuID%d", i), Form("hTNPMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0TNPMCMuID.push_back(h); 
    h = new TH1D(Form("h1TNPMCMuID%d", i), Form("hTNPMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1TNPMCMuID.push_back(h); 


    h = new TH1D(Form("h0MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0MCTrigger.push_back(h); 
    h = new TH1D(Form("h1MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCTrigger.push_back(h); 
    h = new TH1D(Form("h0MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0MCMuID.push_back(h); 
    h = new TH1D(Form("h1MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCMuID.push_back(h); 
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
}


// ----------------------------------------------------------------------
void plotClass::initNumbers(numbers *a) {

  a->name = "";
  a->genAccYield  = a->genAccFileYield = 0; 
  a->effGenFilter = a->effGenFilterE = 1.;
  a->genFileYield = a->genYield = 0.;
  a->recoYield    = a->muidYield = a->trigYield = a->candYield = a->ana0Yield = a->anaYield = a->anaWmcYield = 0; 
  a->ana0YieldE   = a->anaYieldE = a->anaMuonYieldE = a->anaTriggerYieldE = a->anaWmcYieldE = 0.;
  a->anaMuonYield = a->anaTriggerYield = 0.;
  a->fitYield     = a->fitYieldE = 0.;
  a->fitYieldC    = a->fitYieldCE = 0.;
  a->acc          = a->accE = 0; 
  a->effMuidMC    = a->effMuidMCE = a->effTrigMC = a->effTrigMCE = 0; 
  a->effMuidTNP   = a->effMuidTNPE = a->effTrigTNP = a->effTrigTNPE = 0; 
  a->effMuidTNPMC = a->effMuidTNPMCE = a->effTrigTNPMC = a->effTrigTNPMCE = 0.;
  a->effCand      = a->effCandE = 0; 
  a->effPtReco    = a->effPtRecoE = 0.;
  a->effAna       = a->effAnaE = 0; 
  a->effPtReco    = a->effPtRecoE = 0;
  a->effTot       = a->effTotE = a->aEffProdMC = a->aEffProdMCE = a->effProdMC = a->effProdMCE = a->effProdTNP = a->effProdTNPE = 0.; 
  a->prodGenYield = a->combGenYield;
  // -- this is only relevant for the signal(s)
  a->expSignal = 0.;
  a->pss   = a->pdd  = 1.;
  a->psd   = a->pds  = 1.;
  a->pssE  = a->pddE = 1.;
  a->psdE  = a->pdsE = 1.;
  a->bgObs = a->bgBsExp = a->bgBsExpE = a->bgBdExp = a->bgBdExpE =0; 
  a->bsObs = a->bdObs = 0; 
  a->mBdLo = a->mBdHi = a->mBsLo = a->mBsHi = 0.;

  a->genFileYieldE = a->genYieldE = a->recoYieldE = a->muidYieldE = a->trigYieldE = a->candYieldE = 0;

  a->tauBs = a->tauBsE = a->tauBd = a->tauBdE = 0.; 

  a->offLoRare = a->offLoRareE = a->offHiRare = a->offHiRareE = a->bsRare = a->bsRareE = a->bdRare = a->bdRareE = 0.; 
  a->bsNoScaled = a->bsNoScaledE = a->bdNoScaled = a->bdNoScaledE = 0.;
  a->mBdLo = a->mBdHi = a->mBsLo = a->mBsHi = 0.;

}


// ----------------------------------------------------------------------
void plotClass::loopTree(int mode, int proc) {
  // -- mode:
  // 0  Bs2MuMu MC
  // 1  Bd2MuMu MC
  // 5  Bs2MuMu data
  // 10 Bp2JpsiKp MC
  // 15 Bp2JpsiKp data
  // 20 Bs2JpsiPhi MC
  // 25 Bs2JpsiPhi data
  // 98 preset file/directory, do not change anywhere else 
  // 99 peaking decays

  PidTable *ptT1;
  PidTable *ptT2;
  PidTable *ptM; 

  PidTable *ptT1MC;
  PidTable *ptT2MC;
  PidTable *ptMMC; 
  
  bool bp2jpsikp(false), bs2jpsiphi(false), bd2dstarpi(false), isMC(false); 

  string directory; 
  string fAcc;
  string smode; 
  double effFilter(1.), genFileYield(0.); 

  numbers *aa(0);
  if (0 == mode) {
    isMC = true; 
    directory = "candAnaMuMu"; 
    fF["SgMc"]->cd(directory.c_str());
    fAcc = "SgMcAcc";
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
    smode = "SgMc";
  }  else if (1 == mode) {
    isMC = true; 
    directory = "candAnaMuMu"; 
    fF["BdMc"]->cd(directory.c_str());
    fAcc = "BdMcAcc";
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
    smode = "BdMc";
  } else if (5 == mode) {
    isMC = false;     
    directory = "candAnaMuMu"; 
    fF["SgData"]->cd(directory.c_str());
    smode = "SgData";
  } else if (10 == mode) {
    isMC = true; 
    bp2jpsikp = true; 
    directory = "candAnaBu2JpsiK"; 
    fF["NoMc"]->cd(directory.c_str());
    fAcc = "NoMcAcc";
    effFilter = fFilterEff["NoMc"];
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
    smode = "NoMc";
  } else if (15 == mode) {
    isMC = false;     
    bp2jpsikp = true; 
    directory = "candAnaBu2JpsiK"; 
    fF["NoData"]->cd(directory.c_str());
    smode = "NoData";
  } else if (20 == mode) {
    isMC = true; 
    bs2jpsiphi = true; 
    directory = "candAnaBs2JpsiPhi"; 
    fF["CsMc"]->cd(directory.c_str());
    fAcc = "CsMcAcc";
    effFilter = fFilterEff["CsMc"];
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
    smode = "CsMc";
  } else if (25 == mode) {
    isMC = false;     
    bs2jpsiphi = true; 
    directory = "candAnaBs2JpsiPhi"; 
    fF["CsData"]->cd(directory.c_str());
    smode = "CsData";
  } else if (30 == mode) {
    isMC = false;     
    bd2dstarpi = true; 
    directory = "candAnaBd2DstarPi"; 
    fF["DstarPiMc"]->cd(directory.c_str());
    smode = "DstarPiMc";
  } else if (35 == mode) {
    isMC = false;     
    bd2dstarpi = true; 
    directory = "candAnaBd2DstarPi"; 
    fF["DstarPiData"]->cd(directory.c_str());
    smode = "DstarPiData";
  } else if (98 == mode) {
    directory = "candAnaMuMu"; 
    cout << "mode 98" << endl;
    smode = "SgData";
  } else {
    directory = "candAnaMuMu"; 
    cout << "mode 99" << endl;
    smode = "SgData";
  }

  ptT1 = new PidTable("../macros/pidtables/111210/L1L2Efficiency_VBTF_data_all_histo.dat"); 	
  ptT2 = new PidTable("../macros/pidtables/111210/L3Efficiency_VBTF_data_all_histo.dat"); 	
  ptM  = new PidTable("../macros/pidtables/111210/MuonID_VBTF_data_all_histo.dat"); 

  ptT1MC = new PidTable("../macros/pidtables/111210/L1L2Efficiency_VBTF_datalike_mc_histo.dat"); 	
  ptT2MC = new PidTable("../macros/pidtables/111210/L3Efficiency_VBTF_datalike_mc_histo.dat"); 	
  ptMMC  = new PidTable("../macros/pidtables/111210/MuonID_VBTF_datalike_mc_histo.dat"); 

  //   if (isMC) {
  //     delete ptT1; 
  //     ptT1 = new PidTable("../macros/pidtables/111210/L1L2Efficiency_VBTF_datalike_mc_histo_sg.dat"); 	
  //     delete ptT2; 
  //     ptT2 = new PidTable("../macros/pidtables/111210/L3Efficiency_VBTF_datalike_mc_histo.dat"); 	
  //     delete ptM; 
  //     ptM  = new PidTable("../macros/pidtables/111210/MuonID_VBTF_Track7_datalike_mc_histo_sg.dat"); 
  //   } else {
  //   }

  cout << "--> loopTree with mode " << mode << " fDoUseBDT: " << fDoUseBDT << " inv iso: " << fInvertedIso 
       << " proc = " << proc << "; ";
  gDirectory->pwd();
  if (fDoUseBDT) {
    fCuts[0]->m1pt = 4.0;
    fCuts[0]->m2pt = 4.0; 
    fCuts[1]->m1pt = 4.0;
    fCuts[1]->m2pt = 4.0; 
  }
  cout << "    m1pt> " << fCuts[0]->m1pt << " m2pt> " << fCuts[0]->m2pt << endl;

  // -- reset all histograms
  hMassPiPi->Reset();
  hMassKPi->Reset();
  hMassPiK->Reset();

  for (unsigned int i = 0; i < fNchan; ++i) {
    fhBdtMass[i]->Reset();
    
    fhMuId[i]->Reset();
    fhMuTr[i]->Reset();
    fhMuIdMC[i]->Reset();
    fhMuTrMC[i]->Reset();

    fh0TNPTrigger[i]->Reset();
    fh1TNPTrigger[i]->Reset();
    fh0TNPMuID[i]->Reset();
    fh1TNPMuID[i]->Reset();

    fh0TNPMCTrigger[i]->Reset();
    fh1TNPMCTrigger[i]->Reset();
    fh0TNPMCMuID[i]->Reset();
    fh1TNPMCMuID[i]->Reset();

    fh0MCTrigger[i]->Reset();
    fh1MCTrigger[i]->Reset();
    fh0MCMuID[i]->Reset();
    fh1MCMuID[i]->Reset();

    fhMassAbsNoCuts[i]->Reset();
    fhMassNoCuts[i]->Reset();
    fhMassNoCutsManyBins[i]->Reset();

    fhMassWithAnaCuts[i]->Reset();
    fhMassWithAnaCutsManyBins[i]->Reset();

    fhMassWithMuonCuts[i]->Reset();
    fhMassWithMuonCutsManyBins[i]->Reset();

    fhMassWithTriggerCuts[i]->Reset();
    fhMassWithTriggerCutsManyBins[i]->Reset();

    fhMassWithAllCuts[i]->Reset();
    fhMassWithAllCuts[i]->SetTitle(Form("hMassWithAllCuts%d", i));

    fhMassWithAllCutsBlind[i]->Reset();
    fhMassWithAllCutsManyBins[i]->Reset();

    fhDstarPi[i]->Reset();
    fhNorm[i]->Reset();
    fhNormC[i]->Reset();

    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();
  }

  // -- set up tree
  double mass(0.); 
  TTree *t;
  t = (TTree*)gDirectory->Get("events");
  setupTree(t, smode); 

  TEventList *tlist = new TEventList;
  TEventList *elist = new TEventList;
  double tr1w8(0.), tr2w8(0.), trw8(0.), m1w8(0.), m2w8(0.), mw8(0.0);

  int nentries = Int_t(t->GetEntries());
  int nb(0), ievt(0), bsevt0(0), bdevt0(0), bsevt1(0), bdevt1(0), bgevt0(0), bgevt1(0); 
  cuts *pCuts(0); 

  if (bd2dstarpi) {
    cout << "=============================================================" << endl;
    cout << "BDT cut: " << static_cast<double>(proc)/100. << endl;
    cout << "=============================================================" << endl;
  }

  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    // -- require truth matching when on MC
    if (0 == mode && 0 == fb.tm) continue;
    if (1 == mode && 0 == fb.tm) continue;
    if (10 == mode && 0 == fb.tm) continue;
    if (20 == mode && 0 == fb.tm) continue;

    // -- print any event
    if (fDoPrintSingleEvent) {
      if (fb.run == fPrintSingleRun && fb.evt == fPrintSingleEvt) {
	cout << "====> Printing " << fb.run << "/" << fb.evt << endl;
	string bla("epsSignalEvent"); 
	static int pse(0); 
	singleEventPrintout(bla, bla, pse); 
	++pse;
	return;
      } else {
	continue;
      }
    }



    if (isMC && proc > 0) {
      if (fb.procid != proc) continue;
    }

    // -- channel index
    fChan = -1; 
    pCuts = 0; 
    
    fChan = detChan(fb.m1eta, fb.m2eta); 
    if (fChan > -1) {
      pCuts = fCuts[fChan]; 
    } else {
      //       cout << "event " << jentry << ", fChan = " << fOver << " for eta(B) = " << beta 
      // 	   << " m1eta = " << bm1eta << " m2eta = " << bm2eta
      // 	   << endl;
      continue;
    }

    mass = fb.m; 
    //     if (bs2jpsiphi || bp2jpsikp) {
    //       mass = bcm; 
    //     }     

    fhMassAbsNoCuts[fChan]->Fill(mass);
    // -- require wide mass window ??? FIXME WHY????
    //    this should not be enabled, else one looses the underflow for rare sl decays
    //     if (mass < fMassLo) continue;
    //     if (fMassHi < mass) continue;

    // -- gen-level acceptance cuts
    if (isMC) {
      if (TMath::Abs(fb.g1eta) > 2.5) continue;
      if (TMath::Abs(fb.g2eta) > 2.5) continue;
      if (fb.g1pt < 1.0) continue;
      if (fb.g2pt < 1.0) continue;
      if (bp2jpsikp) {
	// gen-level cuts for Bu2JpsiKp
	if (fb.g1pt < 3.5) continue;
	if (fb.g2pt < 3.5) continue;
	if (TMath::Abs(fb.g3eta) > 2.5) continue;
	if (fb.g3pt < 0.4) continue;
      }
      
      if (bs2jpsiphi) {
	if (TMath::Abs(fb.g3eta) > 2.5) continue;
	if (TMath::Abs(fb.g4eta) > 2.5) continue;
	// gen-level cuts for Bs2JpsiPhi
	if (fb.g1pt < 3.5) continue;
	if (fb.g2pt < 3.5) continue;
	if (fb.g3pt < 0.4) continue;
	if (fb.g4pt < 0.4) continue;
      }
    } else {
      if (!fb.json) {
	//	if (5 == mode) cout << "skipping run: " << brun << " LS: " << bls << endl;
	if (bd2dstarpi) {
	} else {
	  continue;
	}
      }
    }

    // -- immutable cuts: require basic muon and trackQual cuts
    if (false == fb.gtqual)  continue;
    if (TMath::Abs(fb.m1eta) > 2.4) continue;
    if (TMath::Abs(fb.m2eta) > 2.4) continue;

    if (fb.m1pt < 1.0) continue;
    if (fb.m2pt < 1.0) continue;

    if (bp2jpsikp) {
      if (TMath::Abs(fb.k1eta) > 2.4) continue;
      if (fb.k1pt < 0.5) continue;
    }
    if (bs2jpsiphi) {
      if (TMath::Abs(fb.k1eta) > 2.4) continue;
      if (TMath::Abs(fb.k2eta) > 2.4) continue;
      if (fb.k1pt < 0.5) continue;
      if (fb.k2pt < 0.5) continue;
    }

    // -- this is the base, after the raw acceptance cuts
    fhMassNoCuts[fChan]->Fill(mass);
    fhMassNoCutsManyBins[fChan]->Fill(mass); 

    if (fDoUseBDT) {
      // -- skip inverted isolation events
      //      if (5 == mode && 5.2 < mass && mass < 5.45) continue; 
      
      //       if(5 == mode || 15 == mode || 25 == mode) {
      // 	if (false == fb.hlt) continue;
      // 	if (false == fb.gmuid) continue;
      //       }

      if (bs2jpsiphi && fb.dr > 0.25) continue;
      if (bs2jpsiphi && fb.mkk < 0.995) continue;
      if (bs2jpsiphi && fb.mkk > 1.045) continue;
      if (bd2dstarpi && fb.md0 < 1.82) continue;
      if (bd2dstarpi && fb.md0 > 1.92) continue;
      if (bd2dstarpi && fb.dm  < 0.1445) continue;
      if (bd2dstarpi && fb.dm  > 0.1470) continue;
      if (bd2dstarpi && fb.ptd0 < 15) continue;
      if (bd2dstarpi && fb.pt   < 15) continue;

      if (bs2jpsiphi || bp2jpsikp) {
	if (fb.mpsi > 3.2) continue;
	if (fb.mpsi < 3.0) continue;
	// -- cowboy veto 
	if (fDoApplyCowboyVeto && fb.cb) continue;
	if (fb.psipt < 7) continue;
      } 

      if (5 == mode) {
	calcBDT(true); 
      } else {
	calcBDT(false);
      }

      // -- moved small tree filling here to have more statistics
      if (fBDT > 0 && fb.hlt && fb.gmuid) {
	elist->Enter(jentry); 
      }
      if (bd2dstarpi) {
	if (fBDT < static_cast<double>(proc)/100.) {
	  continue;
	}
      } else if (fBDT < pCuts->bdt) {
	continue;
      }
      //      cout << "keeping event: " << fBDT << endl;
    } else {

      // -- analysis cuts
      if (fb.m1q*fb.m2q > 0) continue;
      if (fb.m1pt < pCuts->m1pt) continue; 
      if (fb.m2pt < pCuts->m2pt) continue; 
      
      if (fb.fl3d > 1.5) continue;
      if (fb.pvw8 < 0.6) continue;
      
      if (fb.pt < pCuts->pt) continue; 
      if (fb.iso < pCuts->iso) continue; 
      if (fb.chi2/fb.dof > pCuts->chi2dof) continue;
      if (TMath::IsNaN(fb.fls3d)) continue;
      if (fb.fls3d < pCuts->fls3d) continue;
      if (fb.alpha > pCuts->alpha) continue;
      if (fb.docatrk < pCuts->docatrk) continue;
      
      if (bs2jpsiphi && fb.dr > 0.25) continue;
      if (bs2jpsiphi && fb.mkk < 0.995) continue;
      if (bs2jpsiphi && fb.mkk > 1.045) continue;
      
      // -- new cuts
      if (fb.closetrk >= pCuts->closetrk) continue;
      if (TMath::Abs(fb.pvlip) > pCuts->pvlip) continue;
      if (TMath::Abs(fb.pvlips) > pCuts->pvlips) continue;
      if (TMath::Abs(fb.pvlip2) < pCuts->pvlip2) continue;
      if (TMath::Abs(fb.pvlips2) < pCuts->pvlips2) continue;
      if (fb.maxdoca > pCuts->maxdoca) continue;
      if (fb.pvips > pCuts->pvips) continue;
      if (fb.pvip > pCuts->pvip) continue;
      
      if (bs2jpsiphi || bp2jpsikp) {
	if (fb.mpsi > 3.2) continue;
	if (fb.mpsi < 3.0) continue;
	// -- cowboy veto 
	if (fDoApplyCowboyVeto && fb.cb) continue;
	if (fb.psipt < 7) continue;
      } 
      
      if (fDoApplyCowboyVetoAlsoInSignal && fb.cb) continue;
    }


    fhMassWithAnaCuts[fChan]->Fill(mass); 
    fhMassWithAnaCutsManyBins[fChan]->Fill(mass); 

    // ----------------------------
    // -- Intermezzo with PidTables
    // ----------------------------

    // -- muon ID: Data PidTables
    m1w8 = ptM->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
    m2w8 = ptM->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
    mw8  = m1w8*m2w8; 
    fhMuId[fChan]->Fill(mw8); 
    fh0TNPMuID[fChan]->Fill(fb.pt, mw8); 
    fh1TNPMuID[fChan]->Fill(fb.pt); 

    // -- MC for comparison 
    m1w8 = ptMMC->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
    m2w8 = ptMMC->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
    mw8  = m1w8*m2w8; 
    if (mw8 > 0.) {
      fhMuIdMC[fChan]->Fill(mw8, 1./mw8); 
    }
    fh0TNPMCMuID[fChan]->Fill(fb.pt, mw8); 
    fh1TNPMCMuID[fChan]->Fill(fb.pt); 


    // -- muon trigger: Data PidTables
    tr1w8 = ptT1->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.)*ptT2->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
    tr2w8 = ptT1->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.)*ptT2->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (bs2jpsiphi || bp2jpsikp) {
      if (fb.rr >= 3) {
	if (TMath::Abs(fb.m1eta) > 2.2) trw8 = 0; 
	if (TMath::Abs(fb.m2eta) > 2.2) trw8 = 0; 
      }
    }
    fhMuTr[fChan]->Fill(trw8); 
    fh0TNPTrigger[fChan]->Fill(fb.pt, trw8); 
    fh1TNPTrigger[fChan]->Fill(fb.pt); 

    // -- MC for comparison 
    tr1w8 = ptT1MC->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.)*ptT2->effD(fb.m1pt, TMath::Abs(fb.m1eta), 0.);
    tr2w8 = ptT1MC->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.)*ptT2->effD(fb.m2pt, TMath::Abs(fb.m2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (bs2jpsiphi || bp2jpsikp) {
      if (fb.rr >= 3) {
	if (TMath::Abs(fb.m1eta) > 2.2) trw8 = 0; 
	if (TMath::Abs(fb.m2eta) > 2.2) trw8 = 0; 
      }
    }
    fhMuTrMC[fChan]->Fill(trw8); 
    fh0TNPMCTrigger[fChan]->Fill(fb.pt, trw8); 
    fh1TNPMCTrigger[fChan]->Fill(fb.pt); 


    // ---------------------------
    // -- now to the MC simulation
    // ---------------------------

    // -- MUON ID
    if (false == fb.gmuid) continue;
    fhMassWithMuonCuts[fChan]->Fill(mass); 
    fhMassWithMuonCutsManyBins[fChan]->Fill(mass); 

    // -- TRIGGER
    fh0MCTrigger[fChan]->Fill(fb.pt); 
    if (false == fb.hlt) continue;
    fh1MCTrigger[fChan]->Fill(fb.pt); 
    fhMassWithTriggerCuts[fChan]->Fill(mass); 
    fhMassWithTriggerCutsManyBins[fChan]->Fill(mass); 

    fhMassWithAllCuts[fChan]->Fill(mass); 
    if (5 == mode && !(5.2 < mass && mass < 5.45)) {
      fhMassWithAllCutsBlind[fChan]->Fill(mass); 
    }

    fhMassWithAllCutsManyBins[fChan]->Fill(mass); 

    fhNorm[fChan]->Fill(mass);
    fhDstarPi[fChan]->Fill(mass);
    fhNormC[fChan]->Fill(fb.cm);
    //    fhNorm[fChan]->Fill(mass);
    
    if (!fDoUseBDT) elist->Enter(jentry); 

    if (0 == mode && mass < pCuts->mBsLo) continue;
    if (0 == mode && mass > pCuts->mBsHi) continue;
    if (1 == mode && mass < pCuts->mBdLo) continue;
    if (1 == mode && mass > pCuts->mBdHi) continue;
    if (10 == mode && mass < fNoLo) continue;
    if (10 == mode && mass > fNoHi) continue;
    if (20 == mode && mass < fCsLo) continue;
    if (20 == mode && mass > fCsHi) continue;

    fhMassWithMassCuts[fChan]->Fill(mass);
    fhMassWithMassCutsManyBins[fChan]->Fill(mass); 


    fhBdtMass[fChan]->Fill(mass, fBDT); 

    if (fDoPrint && 5 == mode && mass > 4.9 && mass < 5.9) {
      tlist->Enter(jentry); 
      cout << Form("%d m = %4.3f pT = %3.1f eta = %3.2f", fChan, mass, fb.pt, fb.eta)
	//	   <<	" r = " << brun << "/" << bevt
	   << Form(" mpt = %3.1f,%3.1f", fb.m1pt, fb.m2pt)
	   << Form(" meta = %3.2f,%3.2f", TMath::Abs(fb.m1eta), TMath::Abs(fb.m2eta))
	   << Form(" a = %4.3f iso = %3.2f chi2 = %3.1f fls3d = %3.1f, fl/E=%3.1f/%3.2f", 
		   fb.alpha, fb.iso, fb.chi2/fb.dof, fb.fls3d, fb.fl3d, fb.fl3dE)
	   << Form(" pv1: %4.3f/%3.1f", fb.pvlip, fb.pvlips) 
	   << Form(" d/s: %5.4f/%3.2f md: %4.3f d: %4.3f", fb.pvip, fb.pvips, fb.maxdoca, fb.docatrk) 
	   << endl;
      fOUT << Form("%d m = %4.3f pT = %3.1f eta = %3.2f", fChan, mass, fb.pt, fb.eta)
	//	   <<	" run = " << brun << " event = " << bevt
	   << Form(" mpt = %3.1f,%3.1f", fb.m1pt, fb.m2pt)
	   << Form(" meta = %3.2f,%3.2f", TMath::Abs(fb.m1eta), TMath::Abs(fb.m2eta))
	   << Form(" a = %4.3f iso = %3.2f chi2 = %3.1f fls3d = %3.1f, fl/E=%3.1f/%3.2f", 
		   fb.alpha, fb.iso, fb.chi2/fb.dof, fb.fls3d, fb.fl3d, fb.fl3dE)
	   << Form(" pv1: %4.3f/%3.1f", fb.pvlip, fb.pvlips) 
	   << Form(" d/s: %5.4f/%3.2f md: %4.3f d: %4.3f", fb.pvip, fb.pvips, fb.maxdoca, fb.docatrk) 
	   << endl;


      string st("SgEvt"); 

      if (pCuts->mBsLo < mass  && mass < pCuts->mBsHi) {
	st = "BsSgEvt"; 
	if (0 == fChan) {
	  ievt = bsevt0;
	  ++bsevt0; 
	} else {
	  ievt = bsevt1;
	  ++bsevt1; 
	}
      } else if (pCuts->mBdLo < mass  && mass < pCuts->mBdHi) {
	st = "BdSgEvt"; 
	if (0 == fChan) {
	  ievt = bdevt0;
	  ++bdevt0; 
	} else {
	  ievt = bdevt1;
	  ++bdevt1; 
	}
      } else {
	if (0 == fChan) {
	  ievt = bgevt0;
	  ++bgevt0; 
	} else {
	  ievt = bgevt1;
	  ++bgevt1; 
	}
      }

      if (0 == fChan) st += "0:";
      if (1 == fChan) st += "1:";

      string suffix(fSuffix); 
      if (fDoUseBDT) suffix += "Bdt"; 

      singleEventPrintout(suffix, st, ievt);

    }
    
  }


  // -- dump output
  t->SetEventList(elist);
  string tname; 
  if (99 == mode) {
    tname = fRareName;
  } else {
    tname = smode; 
  }
  TFile *fLocal = TFile::Open(Form("small-%s.root", tname.c_str()), "RECREATE"); 
  reduceTree(t);

  TTree *small = t->CopyTree(""); 
  cout << "--> Writing small tree with name " << Form("%s_%s", tname.c_str(), (fDoUseBDT?"bdt":"cnc")) << endl;
  small->SetName(Form("%s_%s", tname.c_str(), (fDoUseBDT?"bdt":"cnc")));
  fLocal->Write();
  fLocal->Close();
  t->SetBranchStatus("*",1);

  cout << gDirectory->GetName() << ": " << fhMassWithAllCuts[0]->GetSumOfWeights() 
       << " (" << fhMassWithAllCuts[0]->Integral(0, fhMassWithAllCuts[0]->GetNbinsX())  << ") " 
       << "  " << fhMassWithAllCuts[1]->GetSumOfWeights()
       << " (" << fhMassWithAllCuts[1]->Integral(0, fhMassWithAllCuts[1]->GetNbinsX())  << ") " 
       << endl;

  if (99 == mode) {
    // -- Cache the pwd...
    TDirectory *pD = gDirectory; 
    fHistFile->cd();
    for (unsigned int i = 0; i < fNchan; ++i) {
      string modifier = (fDoUseBDT?"bdt":"cnc"); 
      fhMassWithAllCuts[i]->SetName(Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), mode, i)); 
      fhMassWithAllCuts[i]->SetTitle(Form("hMassWithAllCuts_%s_%d_chan%d %s", modifier.c_str(), mode, i, pD->GetFile()->GetName())); 
      fhMassWithAllCuts[i]->Write();
      //      fhMassWithAllCuts[i]->SetDirectory(fHistFile);
    }
    pD->cd();

    return;
  }


  //   c0->Clear();
  //   c0->Divide(1,2);
  for (unsigned int i = 0; i < fNchan; ++i) {
    fChan = i; 
    //     c0->cd(i+1);
    pCuts = fCuts[i]; 
    aa = 0; 
    if (0 == mode || 5 == mode) {
      aa = fNumbersBs[i];
      aa->mBdHi = pCuts->mBdHi;
      aa->mBdLo = pCuts->mBdLo;
      aa->mBsHi = pCuts->mBsHi;
      aa->mBsLo = pCuts->mBsLo;
    }
    // -- Bd2MuMu only for the MC numbers!
    if (1 == mode) {
      aa = fNumbersBd[i];
      aa->mBdHi = pCuts->mBdHi;
      aa->mBdLo = pCuts->mBdLo;
      aa->mBsHi = pCuts->mBsHi;
      aa->mBsLo = pCuts->mBsLo;
    }

    if (10 == mode || 15 == mode) {
      aa = fNumbersNo[i];
    }

    if (20 == mode || 25 == mode) {
      aa = fNumbersCs[i];
    }
    
    if (98 == mode) {
      aa = fNumbersBla[i]; 
    }

    if (30 == mode || 35 == mode) {
      aa = fNumbersBla[i]; 
    }

    if (0 == aa) { 
      cout << "++++++++++++????????????++++++++++ aa unset?????????????" << endl;
      cout << "mode = " << mode << endl;
      cout << "i    = " << i << endl;
      cout << "++++++++++++????????????++++++++++ aa unset?????????????" << endl;

    }


    // -- signal cross feed 
    double tot   = fhMassWithAllCutsManyBins[i]->GetSumOfWeights();
    double bd    = fhMassWithAllCutsManyBins[i]->Integral(fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBdLo), 
						       fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBdHi));
    double bs    = fhMassWithAllCutsManyBins[i]->Integral(fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBsLo), 
						       fhMassWithAllCutsManyBins[i]->FindBin(pCuts->mBsHi));
    
    fhMuId[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtmuid-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/muid-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }
    fhMuTr[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtmutr-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else  c0->SaveAs(Form("%s/mutr-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    fhMassAbsNoCuts[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtanc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/anc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    fhMassNoCuts[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtnoc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/noc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    fhMassWithAllCutsManyBins[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtwac-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/wac-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    fhMassWithMuonCutsManyBins[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtmuc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/muc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    fhMassWithTriggerCutsManyBins[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdttrc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/trc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    fhMassWithMassCutsManyBins[i]->Draw();
    if (fDoPrint) {
      if (fDoUseBDT) c0->SaveAs(Form("%s/bdtwmc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
      else c0->SaveAs(Form("%s/wmc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));
    }

    // -- Efficiency and acceptance
    if (isMC) {
      accEffFromEffTree(fAcc, directory, *aa, *fCuts[i], proc);
      double a = fhMassNoCuts[i]->GetSumOfWeights(); 
      double b = fhMassWithAnaCuts[i]->GetSumOfWeights();
      double c = fhMassWithMuonCuts[i]->GetSumOfWeights();
      double d = fhMassWithTriggerCuts[i]->GetSumOfWeights();
      double e = fhMassWithAllCuts[i]->GetSumOfWeights();
      //       double e2 = fhMassWithAllCuts[i]->Integral();
      //       cout << "=================> " << e << " .. " << e2 << endl;
      double f = fhMassWithMassCuts[i]->GetSumOfWeights();
      //WRONG if effTree is from different file than events
      //aa->effCand          = a/aa->recoYield;
      //aa->effCandE         = dEff(static_cast<int>(a), static_cast<int>(aa->recoYield));
      aa->ana0Yield        = a;
      aa->ana0YieldE       = TMath::Sqrt(a);
      aa->anaYield         = b; 
      aa->anaYieldE        = TMath::Sqrt(b); 
      aa->anaMuonYield     = c; 
      aa->anaMuonYieldE    = TMath::Sqrt(c); 
      aa->anaTriggerYield  = d; 
      aa->anaTriggerYieldE = TMath::Sqrt(d); 
      aa->anaWmcYield      = f; 
      aa->anaWmcYieldE     = TMath::Sqrt(f);
      aa->effAna           = b/a*aa->effPtReco;
      aa->effAnaE          = dEff(static_cast<int>(b), static_cast<int>(a)); // FIXME add error from effPtReco
      aa->effMuidMC        = c/b;
      aa->effMuidMCE       = dEff(static_cast<int>(c), static_cast<int>(b));
      aa->effMuidTNP       = fhMuId[i]->GetMean();
      aa->effMuidTNPE      = fhMuId[i]->GetMeanError();
      aa->effMuidTNPMC     = fhMuIdMC[i]->GetMean();
      aa->effMuidTNPMCE    = fhMuIdMC[i]->GetMeanError();
      aa->effTrigMC        = d/c;
      aa->effTrigMCE       = dEff(static_cast<int>(d), static_cast<int>(c));
      aa->effTrigTNP       = fhMuTr[i]->GetMean();
      aa->effTrigTNPE      = fhMuTr[i]->GetMeanError();
      aa->effTrigTNPMC     = fhMuTrMC[i]->GetMean();
      aa->effTrigTNPMCE    = fhMuTrMC[i]->GetMeanError();
      aa->genFileYield     = genFileYield;
      aa->effGenFilter     = effFilter; 
      aa->genYield         = aa->genFileYield/aa->effGenFilter;
      aa->effTot           = e/aa->genYield;
      aa->effTotE          = dEff(static_cast<int>(e), static_cast<int>(aa->genYield));
      aa->effProdMC        = aa->effCand * aa->effAna * aa->effMuidMC * aa->effTrigMC;
      aa->effProdMCE       = 0.;
      aa->effProdTNP       = aa->effCand * aa->effAna * aa->effMuidTNP * aa->effTrigTNP;
      aa->effProdTNPE      = 0.;

      aa->combGenYield     = e/(aa->acc * aa->effProdMC);
      aa->prodGenYield     = e/(aa->effTot); 
    }

    if (fDoPrint) {
      cout << "bs:  " << bs << endl;
      cout << "bd:  " << bd << endl;
      cout << "tot: " << tot << endl;
      cout << "wmc: " << aa->anaWmcYield << endl;
    }

    // -- compute PSD, PDS, etc
    if (0 == mode) {
      aa->pss = bs/tot;
      aa->pssE = dEff(static_cast<int>(bs), static_cast<int>(tot)); 
      aa->pds = bd/tot;
      aa->pdsE = dEff(static_cast<int>(bd), static_cast<int>(tot)); 
    } 
    
    if (1 == mode) {
      aa->pdd = bd/tot;
      aa->pddE = dEff(static_cast<int>(bd), static_cast<int>(tot)); 
      aa->psd = bs/tot;
      aa->psdE = dEff(static_cast<int>(bs), static_cast<int>(tot)); 
    } 

    if (0 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC SIGNAL, channel " << i << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtsig-mc-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/sig-mc-chan%d.pdf", fDirectory.c_str(), i));
      }
      cout << "----> "; gDirectory->pwd(); 
    } else if (5 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA SIGNAL, channel " << i  << endl;
      bgBlind(fhMassWithAllCutsBlind[i], 1, fBgLo, fBgHi);
      cout << "fBgHist = " << fBgHist << "+/-" << fBgHistE 
	   << " lo: " << fBgHistLo << " hi: " << fBgHistHi 
	   << endl;
      aa->bgObs = fBgHist;
      aa->offLo = fBgHistLo;
      aa->offHi = fBgHistHi;
      double blind = fSgHi - fSgLo; 
      double scaleBs = (aa->mBsHi-aa->mBsLo)/(fBgHi-fBgLo-blind);
      double scaleBd = (aa->mBdHi-aa->mBdLo)/(fBgHi-fBgLo-blind);
      aa->tauBs    = scaleBs; 
      aa->tauBsE   = 0.04*scaleBs; 
      aa->tauBd    = scaleBd; 
      aa->tauBdE   = 0.04*scaleBd; 
      cout << "CCCCCCCCCCCC  " << scaleBs << " " << aa->tauBs << "+/-" << aa->tauBsE << endl;
      double combBg  = aa->bgObs - aa->offLoRare - aa->offHiRare;
      double combBgE = 0.2*0.2*(aa->offLoRare + aa->offHiRare)*(aa->offLoRare + aa->offHiRare)  + aa->bgObs;
      combBgE        = TMath::Sqrt(combBgE);
      aa->bgBsExp    = scaleBs*combBg;
      aa->bgBsExpE   = scaleBs*combBgE;
      aa->bgBdExp    = scaleBd*combBg;
      aa->bgBdExpE   = scaleBd*combBgE;

      double cnt = fhMassWithAllCuts[i]->Integral(fhMassWithAllCuts[i]->FindBin(aa->mBsLo), 
						  fhMassWithAllCuts[i]->FindBin(aa->mBsHi));
      aa->bsObs = cnt;

      cnt = fhMassWithAllCuts[i]->Integral(fhMassWithAllCuts[i]->FindBin(aa->mBdLo+0.001), 
					   fhMassWithAllCuts[i]->FindBin(aa->mBdHi-0.001));
      aa->bdObs = cnt;

      // -- blinded version
      TH1D *h = fhMassWithAllCutsBlind[i]; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#mu} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.3); 
      h->SetMinimum(0.01); 
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->SetAxisRange(4.9, 5.9, "X"); 
      //      h->SetMaximum(.2);
      h->SetMaximum(h->GetBinContent(h->GetMaximumBin())+2); 

      double ymax0 = (fhMassWithAllCutsBlind[0]->GetBinContent(fhMassWithAllCutsBlind[0]->GetMaximumBin()) >
		      fhMassWithAllCutsBlind[1]->GetBinContent(fhMassWithAllCutsBlind[1]->GetMaximumBin()) ? 
		      fhMassWithAllCutsBlind[0]->GetBinContent(fhMassWithAllCutsBlind[0]->GetMaximumBin())+2:
		      fhMassWithAllCutsBlind[1]->GetBinContent(fhMassWithAllCutsBlind[1]->GetMaximumBin())+2);
      h->SetMaximum(ymax0); 
      h->Draw();
      //       drawArrow(0.5, 2); 
      //       drawArrow(0.5, 1); 
      drawBox(2, 0.5); 
      drawBox(1, 0.5); 
      
      TH1D *dummy1 = new TH1D("dummy1", "", 10, 0., 10.); setFilledHist(dummy1, kBlue, kBlue, 3356); 
      TH1D *dummy2 = new TH1D("dummy2", "", 10, 0., 10.); setFilledHist(dummy2, kRed, kRed, 3365); 
      
      newLegend(0.4, 0.65, 0.8, 0.8); 
      legg->SetTextSize(0.04);  

      legg->AddEntry(dummy1, "B_{s}^{0} signal window", "f"); 
      legg->AddEntry(dummy2, "B^{0} signal window", "f"); 
      if (0 == i) {
	legg->SetHeader("Barrel");   
      } else {
	legg->SetHeader("Endcap");   
      }      

      legg->Draw();


      stamp(0.18, fStampString, 0.67, fStampCms); 
      if (fDoPrint)  {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtsig-data-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/sig-data-chan%d.pdf", fDirectory.c_str(), i));
      }
      // -- unblinded version
      h = fhMassWithAllCuts[i]; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#mu} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.3); 
      h->SetMinimum(0.01); 
      ymax0 = (fhMassWithAllCuts[0]->GetBinContent(fhMassWithAllCuts[0]->GetMaximumBin()) >
	       fhMassWithAllCuts[1]->GetBinContent(fhMassWithAllCuts[1]->GetMaximumBin()) ? 
	       fhMassWithAllCuts[0]->GetBinContent(fhMassWithAllCuts[0]->GetMaximumBin())+2:
	       fhMassWithAllCuts[1]->GetBinContent(fhMassWithAllCuts[1]->GetMaximumBin())+2);
      h->SetMaximum(ymax0); 
      h->SetAxisRange(4.9, 5.9, "X"); 
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->Draw();
      double yleg = ymax0-1.1;
      
      drawArrow(0.6, 1, yleg); 
      drawArrow(0.4, 2, yleg-0.35); 

      tl->SetNDC(kTRUE); 
      tl->SetTextSize(0.07); 
      tl->SetTextColor(kBlack); 
      if (0 == i) {
	tl->DrawLatex(0.6, 0.8, "Barrel");   
      } 
      
      if (1 == i) {
	tl->DrawLatex(0.6, 0.8, "Endcap");   
      } 
      
      stamp(0.18, fStampString, 0.67, fStampCms); 
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtunblinded-sig-data-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/unblinded-sig-data-chan%d.pdf", fDirectory.c_str(), i));
      }

      delete dummy1; 
      delete dummy2; 
      
      fhBdtMass[i]->Draw("colz");
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdt-mass-sig-data-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/bdt-mass-sig-data-chan%d.pdf", fDirectory.c_str(), i));
      }


    } else if (10 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC NORMALIZATION, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtnorm-mc-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/norm-mc-chan%d.pdf", fDirectory.c_str(), i));
      }
    } else if (15 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA NORMALIZATION, channel " << i  << endl;
      //      TH1D *h = fhMassWithAllCutsManyBins[i]; 
      //      TH1D *h = fhMassWithAllCuts[i]; 
      //      normYield(h, mode, 5.10, 5.5);
      TH1D *h = fhNorm[i];
      //      normYield(h, 0, 4.9, 5.9, 5.145);

      //normYield(h, fChan, 5.0);
      normYield2(h, fChan, 5.0);   // Includes the landau for Bu2JpsiPi

      aa->fitYield  = fNoSig; 
      aa->fitYieldE = fNoSigE; 

      h = fhNormC[i];
      //      normYield(h, 0, 5.0, 5.9, 5.145);
      aa->fitYieldC  = fNoSig; 
      aa->fitYieldCE = fNoSigE; 
    } else if (20 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC CONTROL SAMPLE, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtcs-mc-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/cs-mc-chan%d.pdf", fDirectory.c_str(), i));
      }
    } else if (25 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA CONTROL SAMPLE, channel " << i  << endl;
      //      TH1D *h = fhMassWithAllCuts[i]; 
      //      csYield(h, mode, 5.20, 5.6);
      TH1D *h = fhNorm[i];
      //      csYield(h, i, 5.1, 6.0);
      csYield(h, fChan, 5.1);  // standard control fit
      //const double fraction = 0.14; // fix the ratio of the 2 gaussians, -1. will let it float
      //csYield2(h, fChan, 5.1,5.6,fraction,-1.); // new control fit
      aa->fitYield  = fCsSig; 
      aa->fitYieldE = fCsSigE; 
    } else if (30 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: Bd2DstarPi SAMPLE, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      TH1D *h = fhDstarPi[i];
      h->Draw();
      tl->SetTextSize(0.07); 
      tl->SetTextColor(kBlack); 
      if (0 == mode) {
	tl->DrawLatex(0.22, 0.8, "Barrel");   
      } 
      
      if (1 == mode) {
	tl->DrawLatex(0.22, 0.8, "Endcap");   
      } 
      
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtdp-mc-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/dp-mc-chan%d.pdf", fDirectory.c_str(), i));
      }
    } else if (35 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: Bd2DstarPi SAMPLE, channel " << i  << endl;
      TH1D *h = fhDstarPi[i];
      //      dpYield(h, i, 5.05, 6.0, proc);
      dpYield(h, i, 4.9, 5.65, proc);
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
    } else if (98 == mode) {
      double a = fhMassNoCuts[i]->GetSumOfWeights(); 
      double b = fhMassWithAnaCuts[i]->GetSumOfWeights();
      double c = fhMassWithMuonCuts[i]->GetSumOfWeights();
      double d = fhMassWithTriggerCuts[i]->GetSumOfWeights();
      double e = fhMassWithAllCuts[i]->GetSumOfWeights();
      double f = fhMassWithMassCuts[i]->GetSumOfWeights();
      aa->ana0Yield        = a;
      aa->ana0YieldE       = TMath::Sqrt(a);
      aa->anaYield         = b; 
      aa->anaYieldE        = TMath::Sqrt(b); 
      aa->anaMuonYield     = c; 
      aa->anaMuonYieldE    = TMath::Sqrt(c); 
      aa->anaTriggerYield  = d; 
      aa->anaTriggerYieldE = TMath::Sqrt(d); 
      aa->anaWmcYield      = f; 
      aa->anaWmcYieldE     = TMath::Sqrt(f);
      aa->effAna           = b/a;
      aa->effAnaE          = dEff(static_cast<int>(b), static_cast<int>(a));
      aa->effMuidMC        = c/b;
      aa->effMuidMCE       = dEff(static_cast<int>(c), static_cast<int>(b));
      aa->effTrigMC        = d/c;
      aa->effTrigMCE       = dEff(static_cast<int>(d), static_cast<int>(c));
      aa->effTot           = e/aa->genYield;
      aa->effTotE          = dEff(static_cast<int>(e), static_cast<int>(aa->genYield));
      continue;
    }

    if (0 == mode || 5 == mode || 10 == mode ||15 == mode || 20 == mode ||25 == mode) {
      printNumbers(*aa, cout); 
    }
    printNumbers(*aa, fOUT); 

    // -- Cache the pwd...
    TDirectory *pD = gDirectory; 
    cout << "fHistFile: " << fHistFile  << endl;
    fHistFile->cd();
    string modifier = (fDoUseBDT?"bdt":"cnc"); 
    fhBdtMass[i]->SetName(Form("hBdtMass%d_chan%d", mode, i)); fhBdtMass[i]->Write();

    fhMassWithMassCutsManyBins[i]->SetName(Form("hMassWithMassCutsManyBins%d_chan%d", mode, i)); fhMassWithMassCutsManyBins[i]->Write();
    fhMassWithMassCuts[i]->SetName(Form("hMassWithMassCuts%d_chan%d", mode, i)); fhMassWithMassCuts[i]->Write();
    fhMassWithAllCutsManyBins[i]->SetName(Form("hMassWithAllCutsManyBins_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhMassWithAllCutsManyBins[i]->SetTitle(Form("hMassWithAllCutsManyBins_%s_%d_chan%d %s", modifier.c_str(), mode, i, pD->GetName())); 
    fhMassWithAllCutsManyBins[i]->Write();
    //     fhMassWithAllCutsManyBins[i]->SetDirectory(fHistFile);
    fhMassWithAllCuts[i]->SetName(Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhMassWithAllCuts[i]->Write();
    //     fhMassWithAllCuts[i]->SetDirectory(fHistFile);
    fhNorm[i]->SetName(Form("hNorm_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhNorm[i]->Write();
    //    fhNorm[i]->SetDirectory(fHistFile);
    fhNormC[i]->SetName(Form("hNormC_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhNormC[i]->Write();
    //    fhNormC[i]->SetDirectory(fHistFile);
    fhDstarPi[i]->SetName(Form("hDstarPi_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhDstarPi[i]->Write();
    //    fhDstarPi[i]->SetDirectory(fHistFile);

    if (5 == mode) {
      t->SetEventList(tlist);
      reduceTree(t);
      TTree *small = t->CopyTree(""); 
      t->SetName(Form("%s", (fDoUseBDT?"bdt":"cnc")));
      small->Write(); 
      t->SetBranchStatus("*",1);
      //      small->SetDirectory(fHistFile); 
    }
    //    fhMassNoCuts[i]->SetName(Form("hMassNoCuts%d_chan%d", mode, i)); fhMassNoCuts[i]->Write();
    //    fhMassAbsNoCuts[i]->SetName(Form("hMassAbsNoCuts%d_chan%d", mode, i)); fhMassAbsNoCuts[i]->Write();
    //    fhMuTr[i]->SetName(Form("hMuTr%d_chan%d", mode, i)); fhMuTr[i]->Write();
    //    fhMuId[i]->SetName(Form("hMuId%d_chan%d", mode, i)); fhMuId[i]->Write();
    // -- and get back to it
    pD->cd();

  }


  delete ptT1;
  delete ptT2;
  delete ptM;
  
  delete ptT1MC; 
  delete ptT2MC; 
  delete ptMMC; 

  return;
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

  int nb(0); 
  int ngen(0), ngenlevel(0); 
  int nentries = Int_t(t->GetEntries());
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);

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
  int nb(0); 
  int ngen(0), nchangen(0), nreco(0), nchan(0), nmuid(0), nhlt(0), ncand(0); 
  int chan(-1); 
  int recoPtA(0), recoPtB(0); 
  cout << "channel = " << a.index << endl;
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    if (bidx < 0) continue;
    ++ngen;
    if (proc > 0 && bprocid != proc) continue;
    if (sg) {
      // -- Signal
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1.) {
	    if (bm1pt > 1. && bm2pt > 1.
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4
		&& bm1gt && bm2gt
		) {
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt 
		    ) {
		  ++nchan; 
		  if (bm1id && bm2id) {
		    ++nmuid;
		    if (bhlt) {
		      ++nhlt;
		      if (bm > 0) {
			++ncand;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } else if (no) {
      // -- Normalization
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4) {
	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.5
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4
		&& bm1gt && bm2gt && bk1gt
		) {
	      ++recoPtA; 
	      if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt
		    ) {
		  ++nchan; 
		  if (bm > 0) {
		    ++ncand;
		  }
		  if (bm1id && bm2id) {
		    ++nmuid;
		    if (bhlt) {
		      ++nhlt;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } else if (cs) {
      // -- Control sample
      chan = detChan(bg1eta, bg2eta); 
      if (chan == a.index) {
	++nchangen;
	if (TMath::Abs(bg1eta) < 2.5 && TMath::Abs(bg2eta) < 2.5 && TMath::Abs(bg3eta) < 2.5 && TMath::Abs(bg4eta) < 2.5) {
	  if (bg1pt > 1. && bg2pt > 1. && bg3pt > 0.4 && bg4pt > 0.4) {
	    if (bm1pt > 1. && bm2pt > 1. && bk1pt > 0.5 && bk2pt > 0.5
		&& TMath::Abs(bm1eta) < 2.4 && TMath::Abs(bm2eta) < 2.4 && TMath::Abs(bk1eta) < 2.4 && TMath::Abs(bk2eta) < 2.4
		&& bm1gt && bm2gt && bk1gt && bk2gt
		) {
	      ++recoPtA; 
	      if (bm1pt > 3.5 && bm2pt > 3.5) ++recoPtB;
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt) {
		  ++nchan; 
		  if (bm > 0) {
		    ++ncand;
		  }
		  if (bm1id && bm2id) {
		    ++nmuid;
		    if (bhlt) {
		      ++nhlt;
		    }
		  }
		}
	      }
	    }
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
  a.recoYield     = nreco; // reco'ed in chan, basic global reconstruction cuts 
  a.muidYield     = nmuid;
  a.trigYield     = nhlt;
  a.candYield     = ncand;

  if (a.genAccYield > 0) {
    a.acc = a.recoYield/a.genAccYield;
    a.accE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genAccYield));
  }  

  if (a.trigYield > 0) {
    a.effCand  = a.candYield/a.trigYield;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.trigYield));
    a.effCand  = a.candYield/nchan;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(nchan));
    a.effCand  = 0.98; // estimate
  } 

  cout << "NGEN      = " << ngen << endl;
  cout << "NGENCHAN  = " << nchangen << endl;
  cout << "NRECO     = " << nreco << endl;
  cout << "NRECOCHAN = " << nchan << endl;
  cout << "NMUID     = " << nmuid << endl;
  cout << "NHLT      = " << nhlt << endl;
  cout << "NCAND     = " << ncand << endl;


//   a.genChanYield  = nchangen; 
//   a.chanYield     = nchan; // reco'ed in chan, with channel-dependent (pT) cuts

//   if (a.genAccYield > 0) {
//     a.cFrac  = a.genChanYield/a.genAccYield;
//     a.cFracE = dEff(static_cast<int>(a.genChanYield), static_cast<int>(a.genAccYield));
//   }

//   if (a.genChanYield > 0) {
//     a.accChan = a.recoYield/a.genChanYield;
//     a.accChanE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genChanYield));
//   }
  
//   if (a.recoYield > 0) {
//     a.effChan = a.chanYield/a.recoYield;
//     a.effChanE = dEff(static_cast<int>(a.chanYield), static_cast<int>(a.recoYield));
//   }

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

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	sname = "CsData"; 
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
      if (string::npos != stype.find("bmt,Photon")) {
	sname = "BmtPhoton"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Photon")); 
      }
      if (string::npos != stype.find("bmt,Jet")) {
	sname = "BmtJet"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "Jet")); 
      }
      cout << "open Data file "  << sfile  << " as " << sname << " (" << stype << ") with lumi = " << slumi << endl;
    } else {
      string sfilter = sdset; 
      replaceAll(sfilter, "mc,", ""); 
      double effFilter = atof(sfilter.c_str());
      // -- MC
      pF = loadFile(sfile); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	sname = "SgMc"; 
	fF.insert(make_pair(sname, pF)); 	
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("sg")) {
	sname = "SgMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("sg")) {
	sname = "SgMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("sg")) {
	sname = "SgMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("pu") && string::npos != stype.find("sg")) {
	sname = "SgMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (PU)")); 
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
      if (string::npos != stype.find("Bx2MuMu") && string::npos != stype.find("sg")) {
	sname = "SgBx2MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (5.7GeV)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("By2MuMu") && string::npos != stype.find("sg")) {
	sname = "SgBy2MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (5.1GeV)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("Bs2MuMu") && string::npos != stype.find("sg")) {
	sname = "SgBs2MuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.2e-9)); 
	fBFE.insert(make_pair(sname, 0.06)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (5.37GeV)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("bd")) {
	sname = "BdMc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.0e-10)); 
	fBFE.insert(make_pair(sname, 0.1)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
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
	fBF.insert(make_pair(sname, 6.0e-5)); 
	fBFE.insert(make_pair(sname, 0.03)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	//	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (MC)")); 
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
	fBF.insert(make_pair(sname, 3.8e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	//	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (MC)")); 
	fName.insert(make_pair(sname, "MC simulation")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("1e33") && string::npos != stype.find("cs")) {
	sname = "CsMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.8e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("2e33") && string::npos != stype.find("cs")) {
	sname = "CsMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.8e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("3e33") && string::npos != stype.find("cs")) {
	sname = "CsMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.8e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("acc") && string::npos != stype.find("cs")) {
	sname = "CsMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.8e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("pu") && string::npos != stype.find("cs")) {
	sname = "CsMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.8e-5)); 
	fBFE.insert(make_pair(sname, 0.32)); 
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
	fBF.insert(make_pair(sname, 23.e-6)); 
	fBFE.insert(make_pair(sname, 0.15)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2KPi")) {
	sname = "bgBs2KPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.4e-6)); 
	fBFE.insert(make_pair(sname, 0.22)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2PiPi")) {
	sname = "bgBs2PiPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 0.95e-6)); 
	fBFE.insert(make_pair(sname, 0.9)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2KMuNu")) {
	sname = "bgBs2KMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.36e-4)); 
	fBFE.insert(make_pair(sname, 0.05)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{-}#mu^{+}#nu")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bd2PiMuNu")) {
	sname = "bgBd2PiMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.36e-4)); 
	fBFE.insert(make_pair(sname, 0.05)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{-}#mu^{+}#nu")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2KK")) {
	sname = "bgBd2KK"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 0.11e-6)); 
	fBFE.insert(make_pair(sname, 0.73)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2KPi")) {
	sname = "bgBd2KPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 19.4e-6)); 
	fBFE.insert(make_pair(sname, 0.05)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2PiPi")) {
	sname = "bgBd2PiPi"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.1e-6)); 
	fBFE.insert(make_pair(sname, 0.04)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2KP")) {
	sname = "bgLb2KP"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.6e-6)); 
	fBFE.insert(make_pair(sname, 0.3)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2PiP")) {
	sname = "bgLb2PiP"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 3.5e-6)); 
	fBFE.insert(make_pair(sname, 0.31)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p #pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2PMuNu")) {
	sname = "bgLb2PMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.36e-4)); 
	fBFE.insert(make_pair(sname, 0.05)); 
	fProdR.insert(make_pair(sname, fsfu)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda^{0}_{b} #rightarrow p#mu^{-}#bar{#nu}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bu2PiMuMu")) {
	sname = "bgBu2PiMuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 4.9e-8)); 
	fBFE.insert(make_pair(sname, 0.2)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow #pi^{+}#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bu2KMuMu")) {
	sname = "bgBu2KMuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 5.0e-7)); 
	fBFE.insert(make_pair(sname, 0.2)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow K^{+}#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bs2PhiMuMu")) {
	sname = "bgBs2PhiMuMu"; 
	fF.insert(make_pair(sname, pF)); 
	fBF.insert(make_pair(sname, 1.4e-6)); 
	fBFE.insert(make_pair(sname, 0.4)); 
	fProdR.insert(make_pair(sname, 1.0)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0}_{s} #rightarrow #phi#mu^{+}#mu^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      cout << "open MC file "  << sfile  << " as " << sname << " (" << stype << ") with lumi = " << slumi << endl;
    }
  }

  
  // https://docs.google.com/spreadsheet/ccc?key=0AhrdqmQ22qAldHREbHpVOGVMaHhQaU9ULXAydkpZc2c&hl=en#gid=0
  fNgen.insert(make_pair("bgLb2KP",    356.e6)); 
  fNgen.insert(make_pair("bgLb2PiP",    374.e6)); 
  fNgen.insert(make_pair("bgLb2PMuNu", 6910.e6)); 

  fNgen.insert(make_pair("SgMc",      249.e6)); 
  fNgen.insert(make_pair("bgBs2KK",    1392.e6)); 
  fNgen.insert(make_pair("bgBs2KPi",    598.e6)); 
  fNgen.insert(make_pair("bgBs2PiPi",   148.e6)); 
  fNgen.insert(make_pair("bgBs2KMuNu", 6702.e6)); 

  fNgen.insert(make_pair("BdMc",        40.e6)); 
  fNgen.insert(make_pair("bgBd2PiPi",    394.e6)); 
  fNgen.insert(make_pair("bgBd2KPi",     996.e6)); 
  fNgen.insert(make_pair("bgBd2KK",      200.e6)); 
  fNgen.insert(make_pair("bgBd2PiMuNu", 6742.e6)); 

  fNgen.insert(make_pair("NoMc",15166.e6)); 
  fNgen.insert(make_pair("CsMc", 9970.e6)); 

  fNgen.insert(make_pair("NoMcAcc", 1)); 
  fNgen.insert(make_pair("CsMcAcc", 1)); 

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
  double lumi(0), n(0), f(0); 
  for (map<string, string>::iterator imap = fName.begin(); imap != fName.end(); ++imap) {  
    cout << "===> " << imap->first;
    cout << ": " << fF[imap->first]->GetName();
    cout << " -> " << imap->second << endl;
    name = imap->second; 
    lumi = fLumi[imap->first];
    n = fNgen[imap->first];
    f = ((TH1D*)fF[imap->first]->Get("monEvents"))->GetBinContent(1);
    replaceAll(name, "#", "\\"); 
    //    cout <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    fTEX <<  Form("\\vdef{%s:sampleName:%s}   {\\ensuremath{{%s } } }", fSuffix.c_str(), imap->first.c_str(), name.c_str()) << endl;
    //    cout <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    fTEX <<  Form("\\vdef{%s:lumi:%s}   {\\ensuremath{{%4.1f } } }", fSuffix.c_str(), imap->first.c_str(), lumi) << endl;
    if (n>0) {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f } } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
    } else {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str() ) << endl;
      fTEX <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
    if (f>0) {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{%4.0f } } }", fSuffix.c_str(), imap->first.c_str(), n) << endl;
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{%4.0f} } }", fSuffix.c_str(), imap->first.c_str(), f) << endl;
    } else {
      //      cout <<  Form("\\vdef{%s:ngen:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str() ) << endl;
      fTEX <<  Form("\\vdef{%s:nfilt:%s}   {\\ensuremath{{n/a } } }", fSuffix.c_str(), imap->first.c_str()) << endl;
    }
  }
  fTEX.flush();

}


// ----------------------------------------------------------------------
void plotClass::bgBlind(TH1 *h, int mode, double lo, double hi) {
  
  if (0 == h) { 
    cout << "plotClass::bgBlind(...): No histogram passed! mode = " << mode << endl;
    return;
  }
  
  TF1 *lF1(0);

  double histCount = h->Integral(h->FindBin(fBgLo+0.0001), h->FindBin(fBgHi-0.0001)); 
  cout << "bgBlind: histCount = " << histCount 
       << " starting at " << h->FindBin(fBgLo+0.0001) << " to " << h->FindBin(fBgHi-0.0001) << endl;
  fBgHist  = histCount; 
  fBgHistE  = TMath::Sqrt(histCount); 
  double delta = fSgHi-fSgLo; 
  fBgHistExp  = histCount*(delta)/(fBgHi-fBgLo-delta);
  fBgHistLo =  h->Integral(h->FindBin(fBgLo+0.0001), h->FindBin(fSgLo-0.0001)); 
  fBgHistHi =  h->Integral(h->FindBin(fSgHi+0.0001), h->FindBin(fBgHi-0.0001)); 
  if (histCount > 0) {
    fBgHistE = TMath::Sqrt(histCount)/histCount*fBgHist;
  } else {
    fBgHistE = 0.2; // FIXME?!
    fBgExp = 0.;
    fBgExpE = 0.2; 
    return;
  }

  if (0 == mode) {
    fBgExp = fBgHist;
    fBgExpE = fBgHistE;
    return;
  }
  else if (1 == mode) { 
    //    lF1 = fpFunc->pol0(h); 
    lF1 = fpFunc->pol0BsBlind(h); 
  } else if (2 == mode) {
    lF1 = fpFunc->pol1BsBlind(h); 
  }
  
  //  setErrors(h);
  h->Fit(lF1, "lr", "", lo, hi); 

  double c  = h->GetFunction("f1")->GetParameter(0); 
  double cE = h->GetFunction("f1")->GetParError(0); 
  fBgExp  = c * (fSgHi - fSgLo)/h->GetBinWidth(1);
  fBgExpE = cE/c*fBgExp; 
  cout << "bgBlind: c = " << c << " sig width = " << (fSgHi - fSgLo) << endl;
}


// ----------------------------------------------------------------------
void plotClass::normYield(TH1 *h, int mode, double lo, double hi, double preco) {

//   double pReco = (preco<0? (1 == mode?5.145:5.146): preco);
  double pReco = (preco<0? 5.145: preco);

  TF1 *lF1(0), *lBg(0);

  string name(h->GetName()); 
  double sigma1(0.03);
  if (1 == fChan) sigma1 = 0.04; 
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
  if (0 == fChan) {
    tl->DrawLatex(0.6, 0.8, "Barrel");   
  } 

  if (1 == fChan) {
    tl->DrawLatex(0.6, 0.8, "Endcap");   
  } 
  
  stamp(0.20, fStampString, 0.67, fStampCms); 
  if (fDoPrint) {
    
    string pdfname;
    string hname(h->GetName());
    if (string::npos != hname.find("NormC")) {
      pdfname = Form("%s/normC-data-chan%d.pdf", fDirectory.c_str(), fChan);
      if (fDoUseBDT)  pdfname = Form("%s/bdtnormC-data-chan%d.pdf", fDirectory.c_str(), fChan);
    } else {
      pdfname = Form("%s/norm-data-chan%d.pdf", fDirectory.c_str(), fChan);
      if (fDoUseBDT)  pdfname = Form("%s/bdtnorm-data-chan%d.pdf", fDirectory.c_str(), fChan);
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
    lF1 = fpFunc->expoErrGaussLandau(h, 5.27, 0.056, pReco); 
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

  tl->SetTextSize(0.07); 
  tl->SetTextColor(kBlack); 
  if (0 == fChan) {
    tl->DrawLatex(0.6, 0.8, "Barrel");   
  } 

  if (1 == fChan) {
    tl->DrawLatex(0.6, 0.8, "Endcap");   
  } 
  
  stamp(0.20, fStampString, 0.67, fStampCms); 
  if (fDoPrint) {
    
    string pdfname;
    string hname(h->GetName());
    if (string::npos != hname.find("NormC")) {
      pdfname = Form("%s/normC-data-chan%d.pdf", fDirectory.c_str(), fChan);
      if (fDoUseBDT)  pdfname = Form("%s/bdtnormC-data-chan%d.pdf", fDirectory.c_str(), fChan);
    } else {
      pdfname = Form("%s/norm-data-chan%d.pdf", fDirectory.c_str(), fChan);
      if (fDoUseBDT)  pdfname = Form("%s/bdtnorm-data-chan%d.pdf", fDirectory.c_str(), fChan);
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

  zone(1,2);
  h->DrawCopy();
  c0->cd(2);
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

  //lF1->SetLineColor(kBlue); 
  //lF1->SetLineStyle(kDotted);
  //lF1->SetLineWidth(3);
  //lF1->Draw("same");

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


  if (fDoPrint) {
    if (fDoUseBDT) c0->SaveAs(Form("%s/bdtcs-data-chan%d.pdf", fDirectory.c_str(), mode));
    else c0->SaveAs(Form("%s/cs-data-chan%d.pdf", fDirectory.c_str(), mode));
  }

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 
  delete lBg; 

}

// ----------------------------------------------------------------------
void plotClass::csYield(TH1 *h, int mode, double lo, double hi, double preco) {

  double pReco = (preco<0? 5.2: preco);

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

  zone(1,2);
  h->DrawCopy();
  c0->cd(2);
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


  if (fDoPrint) {
    if (fDoUseBDT) c0->SaveAs(Form("%s/bdtcs-data-chan%d.pdf", fDirectory.c_str(), mode));
    else c0->SaveAs(Form("%s/cs-data-chan%d.pdf", fDirectory.c_str(), mode));
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
  OUT << "muidYield       = " << a.muidYield << endl;
  OUT << "trigYield       = " << a.trigYield << endl;
  OUT << "candYield       = " << a.candYield << endl;
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
  OUT << "effMuidTNP      = " << a.effMuidTNP << "+/-" << a.effMuidTNPE << endl;
  OUT << "effMuidTNPMC    = " << a.effMuidTNPMC << "+/-" << a.effMuidTNPMCE << endl;
  OUT << "effTrigMC       = " << a.effTrigMC << "+/-" << a.effTrigMCE << endl;
  OUT << "effTrigTNP      = " << a.effTrigTNP << "+/-" << a.effTrigTNPE << endl;
  OUT << "effTrigTNPMC    = " << a.effTrigTNPMC << "+/-" << a.effTrigTNPMCE << endl;
  OUT << "effProd(MC)     = " << a.effProdMC << endl;
  OUT << "effProd(TNP)    = " << a.effProdTNP << endl;
  OUT << "effProd(MC)A    = " << a.effProdMC*a.acc << endl;
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
  double value; 
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    cut = string(h->GetXaxis()->GetBinLabel(i));
    value = h->GetBinContent(i);
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
string plotClass::scientificTex(double n, double nE, std::string name, double base, int digits) {

  char line[200]; 
  double a1 = n/base; 
  double a2 = nE/base; 
  int  expo = TMath::Log10(base); 

  if ( isnan(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%i } } }", name.c_str(), static_cast<int>(n));
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%2.1f \\pm %2.1f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%3.2f \\pm %3.2f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%4.3f \\pm %4.3f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.4f \\pm %5.4f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%6.5f \\pm %6.5f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%f \\pm %f)\\times 10^{%i}} } }", name.c_str(), a1, a2, expo);
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

    if (!strcmp(CutName, "bdt")) {
      a->bdt = CutValue; ok = 1;
      if (dump) cout << "bdt:              " << CutValue << endl;
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
      a->xmlFile = XmlName;
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
    fTEX <<  Form("\\vdef{%s:bdt:%d}     {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->bdt) << endl;

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
  tl->SetTextSize(fSize); 
  tl->DrawLatex(x1, 0.91, text1.c_str());   
  tl->DrawLatex(x2, 0.91, text2.c_str()); 
  tl->SetTextSize(0.03); 
  tl->SetTextSize(fSize); 
}



// ----------------------------------------------------------------------
void plotClass::drawArrow(double height, int mode, double ylegend) {

  double ylo(0.01); 
  pl->SetLineWidth(3.); 
  
  double d(0.08), y(0.80), x(5.25); 
  
  if (1 == mode) {
    pl->SetLineColor(kBlue); 
    pl->SetLineColor(kBlue); 
    pl->SetLineStyle(kSolid); 
    pl->DrawLine(fCuts[0]->mBsLo, height, fCuts[0]->mBsHi, height); 
    pl->SetLineWidth(2.); 
    pl->DrawLine(fCuts[0]->mBsLo, height+d, fCuts[0]->mBsLo, height-d); 
    pl->DrawLine(fCuts[0]->mBsHi, height+d, fCuts[0]->mBsHi, height-d); 

    if (1) {
      y = ylegend;
      x = 5.25; 
      d = 0.05;
      pl->SetLineWidth(3.); 
      pl->DrawLine(x, y, x+0.1, y); 
      pl->SetLineWidth(2.); 
      pl->DrawLine(x, y+d, x, y-d); 
      pl->DrawLine(x+0.1, y+d, x+0.1, y-d); 
      tl->SetNDC(kFALSE);
      tl->SetTextSize(0.05); 
      tl->DrawLatex(x+0.15, y-0.05, "B_{s}^{0} signal window");
    }

  } else if (2 == mode) {
    pl->SetLineColor(kRed); 
    pl->SetLineColor(kRed); 
    pl->SetLineStyle(kDashed); 
    pl->DrawLine(fCuts[0]->mBdLo, height, fCuts[0]->mBdHi, height); 
    pl->SetLineStyle(kSolid); 
    pl->SetLineWidth(2.); 
    pl->DrawLine(fCuts[0]->mBdLo, height+d, fCuts[0]->mBdLo, height-d); 
    pl->DrawLine(fCuts[0]->mBdHi, height+d, fCuts[0]->mBdHi, height-d); 

    if (1) {
      x = 5.25; 
      y = ylegend;
      d = 0.05;
      pl->SetLineWidth(3.); 
      pl->SetLineStyle(kDashed); 
      pl->DrawLine(x, y, x+0.1, y); 
      pl->SetLineStyle(kSolid); 
      pl->SetLineWidth(2.); 
      pl->DrawLine(x, y+d, x, y-d); 
      pl->DrawLine(x+0.1, y+d, x+0.1, y-d); 
      tl->SetNDC(kFALSE);
      tl->SetTextSize(0.05); 
      tl->DrawLatex(x+0.15, y-0.05, "B^{0} signal window");
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
void plotClass::drawBox(int mode, double hi, int ylo) {

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
void plotClass::loopOverTree(TTree *t, std::string mode, int function, int nevts) {
  int nentries = Int_t(t->GetEntries());
  if (nevts > 0) nentries = nevts;
  //  nentries = 100000; 

  hMassPiPi->Reset(); hMassPiPi->SetTitle(Form("mPiPi %s", mode.c_str()));
  hMassKPi->Reset();  hMassKPi->SetTitle(Form("mKPi %s", mode.c_str()));
  hMassPiK->Reset();  hMassPiK->SetTitle(Form("mPiK %s", mode.c_str()));

  int nb(0);
  int step(50000); 
  if (nentries < 50000000) step = 5000000; 
  if (nentries < 10000000) step = 1000000; 
  if (nentries < 5000000)  step = 500000; 
  if (nentries < 1000000)  step = 100000; 
  if (nentries < 100000)   step = 10000; 
  if (nentries < 10000)    step = 1000; 
  if (nentries < 1000)     step = 100; 
  int imode(0); 
  if (string::npos != mode.find("No")) imode = 10; 
  if (string::npos != mode.find("Cs")) imode = 20; 
  cout << "==> plotClass::loopOverTree> looping in mode " << mode << " -> imode = " << imode 
       << " with " << nentries << " entries" 
       << endl;

  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. Event %8d", jentry) << endl;
    candAnalysis(imode);
    loopFunction(function);
  }

  //   hMassPiPi->Write();
  //   hMassKPi->Write();
  //   hMassPiK->Write();
  
}

// ----------------------------------------------------------------------
TTree* plotClass::getTree(string mode) {
  TTree *t(0);
  cout << "retrieve tree events for mode " << mode << " from file " << fF[mode]->GetName() << endl;
  if (string::npos != mode.find("No")) t = (TTree*)fF[mode]->Get("candAnaBu2JpsiK/events"); 
  if (string::npos != mode.find("Cs")) t = (TTree*)fF[mode]->Get("candAnaBs2JpsiPhi/events"); 
  if (string::npos != mode.find("Sg")) t = (TTree*)fF[mode]->Get("candAnaMuMu/events"); 
  return t; 
}

// ----------------------------------------------------------------------
void plotClass::setupTree(TTree *t, string mode) {
  t->SetBranchAddress("pt", &fb.pt);

  t->SetBranchAddress("bdt",&fb.bdt);
  t->SetBranchAddress("bdt",&fBDT);
  t->SetBranchAddress("lip",&fb.lip);
  t->SetBranchAddress("lipE",&fb.lipE);
  t->SetBranchAddress("tip",&fb.tip);
  t->SetBranchAddress("tipE",&fb.tipE);

  t->SetBranchAddress("closetrk",&fb.closetrk);
  t->SetBranchAddress("pvlip",&fb.pvlip);
  t->SetBranchAddress("pvlips",&fb.pvlips);
  t->SetBranchAddress("maxdoca",&fb.maxdoca);
  t->SetBranchAddress("pvip",&fb.pvip);
  t->SetBranchAddress("pvips",&fb.pvips);
  t->SetBranchAddress("pvw8",&fb.pvw8);

  t->SetBranchAddress("m1pix",&fb.m1pix);
  t->SetBranchAddress("m2pix",&fb.m2pix);
  t->SetBranchAddress("m1bpix",&fb.m1bpix);
  t->SetBranchAddress("m2bpix",&fb.m2bpix);
  t->SetBranchAddress("m1bpixl1",&fb.m1bpixl1);
  t->SetBranchAddress("m2bpixl1",&fb.m2bpixl1);

  t->SetBranchAddress("rr",&fb.rr);
  t->SetBranchAddress("pvn",&fb.pvn);
  t->SetBranchAddress("run",&fb.run);
  t->SetBranchAddress("evt",&fb.evt);
  t->SetBranchAddress("hlt",&fb.hlt);
  t->SetBranchAddress("ls",&fb.ls);
  t->SetBranchAddress("cb",&fb.cb);
  t->SetBranchAddress("json",&fb.json);
  t->SetBranchAddress("gmuid",&fb.gmuid);
  t->SetBranchAddress("gtqual",&fb.gtqual);
  t->SetBranchAddress("w8mu",&fb.w8mu);
  t->SetBranchAddress("w8tr",&fb.w8tr);
  t->SetBranchAddress("tm",&fb.tm);
  t->SetBranchAddress("procid",&fb.procid);
  t->SetBranchAddress("m",&fb.m);
  t->SetBranchAddress("cm",&fb.cm);
  t->SetBranchAddress("pt",&fb.pt);
  t->SetBranchAddress("phi",&fb.phi);
  t->SetBranchAddress("eta",&fb.eta);
  t->SetBranchAddress("cosa",&fb.cosa);
  t->SetBranchAddress("alpha",&fb.alpha);
  t->SetBranchAddress("iso",&fb.iso);
  t->SetBranchAddress("chi2",&fb.chi2);
  t->SetBranchAddress("dof",&fb.dof);
  t->SetBranchAddress("prob",&fb.pchi2dof);
  t->SetBranchAddress("flsxy",&fb.flsxy);
  t->SetBranchAddress("fls3d",&fb.fls3d);
  t->SetBranchAddress("fl3d",&fb.fl3d);
  t->SetBranchAddress("fl3dE",&fb.fl3dE);
  t->SetBranchAddress("m1pt",&fb.m1pt);
  t->SetBranchAddress("m1eta",&fb.m1eta);
  t->SetBranchAddress("m1phi",&fb.m1phi);
  t->SetBranchAddress("m1q",&fb.m1q);
  t->SetBranchAddress("m2pt",&fb.m2pt);
  t->SetBranchAddress("m2eta",&fb.m2eta);
  t->SetBranchAddress("m2phi",&fb.m2phi);
  t->SetBranchAddress("m2q",&fb.m2q);
  t->SetBranchAddress("docatrk",&fb.docatrk);

  t->SetBranchAddress("g1pt",&fb.g1pt);
  t->SetBranchAddress("g2pt",&fb.g2pt);
  t->SetBranchAddress("g1eta",&fb.g1eta);
  t->SetBranchAddress("g2eta",&fb.g2eta);
  if (string::npos != mode.find("No")) {
    if (string::npos != mode.find("Mc")) {
      t->SetBranchAddress("g3pt",&fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
    }
    t->SetBranchAddress("kpt",&fb.k1pt);
    t->SetBranchAddress("keta",&fb.k1eta);
    t->SetBranchAddress("mpsi",&fb.mpsi);
    t->SetBranchAddress("psipt",&fb.psipt); //FIXME
  }

  if (string::npos != mode.find("Cs")) {
    if (string::npos != mode.find("Mc")) {
      t->SetBranchAddress("g3pt",&fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
      t->SetBranchAddress("g4pt",&fb.g4pt);
      t->SetBranchAddress("g4eta",&fb.g4eta);
    }
    t->SetBranchAddress("psipt",&fb.psipt);   //FIXME
    t->SetBranchAddress("mpsi",&fb.mpsi);
    t->SetBranchAddress("mkk",&fb.mkk);
    t->SetBranchAddress("dr",&fb.dr);
    t->SetBranchAddress("k1pt",&fb.k1pt);
    t->SetBranchAddress("k1eta",&fb.k1eta);
    t->SetBranchAddress("k2pt",&fb.k2pt);
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

  bool bp2jpsikp(false), bs2jpsiphi(false); 
  if (10 == mode)  bp2jpsikp = true; 
  if (20 == mode)  bs2jpsiphi = true; 

  // -- reset all
  fBDT = -99.; 
  fGoodHLT = fGoodMuonsID = fGoodMuonsPt = fGoodMuonsEta = fGoodTracks = fGoodTracksPt = fGoodTracksEta = false;
  fGoodQ = fGoodPvAveW8 = fGoodMaxDoca = fGoodIp = fGoodIpS = fGoodPt = fGoodEta = fGoodAlpha =  fGoodChi2 = fGoodFLS = false;   
  fGoodCloseTrack = fGoodIso = fGoodDocaTrk = fGoodLastCut = fPreselection = false;

  if (fIsMC) {
    if (fb.g1pt < 3.5) return;
    if (fb.g2pt < 3.5) return;
    if (TMath::Abs(fb.g1eta) > 2.5) return;
    if (TMath::Abs(fb.g2eta) > 2.5) return;
    if (!fb.json) {
      return;
    }
  }

  if (fIsMC && bp2jpsikp) {
    // gen-level cuts for Bu2JpsiKp
    if (fb.g1pt < 3.5) return;
    if (fb.g2pt < 3.5) return;
    if (TMath::Abs(fb.g3eta) > 2.5) return;
    if (fb.g3pt < 0.4) return;
  }
  
  if (fIsMC && bs2jpsiphi) {
    if (TMath::Abs(fb.g3eta) > 2.5) return;
    if (TMath::Abs(fb.g4eta) > 2.5) return;
    // gen-level cuts for Bs2JpsiPhi
    if (fb.g1pt < 3.5) return;
    if (fb.g2pt < 3.5) return;
    if (fb.g3pt < 0.4) return;
    if (fb.g4pt < 0.4) return;
  }

  if (bp2jpsikp) {
    if (TMath::Abs(fb.k1eta) > 2.4) return;
    if (fb.k1pt < 0.5) return;
  }
  if (bs2jpsiphi) {
    if (TMath::Abs(fb.k1eta) > 2.4) return;
    if (TMath::Abs(fb.k2eta) > 2.4) return;
    if (fb.k1pt < 0.5) return;
    if (fb.k2pt < 0.5) return;
    if (fb.dr >0.3) return;
    if (fb.mkk < 0.995) return;
    if (fb.mkk > 1.045) return;
  }

  if (bs2jpsiphi || bp2jpsikp) {
    if (fb.mpsi > 3.2) return;
    if (fb.mpsi < 3.0) return;
    // -- cowboy veto 
    if (fDoApplyCowboyVeto && fb.cb) return;
    if (fb.psipt < 7) return;
  } 

  hMassPiPi->Fill(recalcMass(MPION, MPION)); 
  hMassKPi->Fill(recalcMass(MKAON, MPION)); 
  hMassPiK->Fill(recalcMass(MKAON, MKAON)); 

  if (fChan > -1) {
    pCuts = fCuts[fChan]; 
  } else {
    //       cout << "event " << jentry << ", fChan = " << fOver << " for eta(B) = " << beta 
    // 	   << " m1eta = " << bm1eta << " m2eta = " << bm2eta
    // 	   << endl;
    return;
  }

  if (fDoUseBDT) {
    if (TMath::IsNaN(fb.fls3d)) return; 
    calcBDT(true); 
  }

  fGoodMuonsID    = fb.gmuid;
  fGoodMuonsPt    = ((fb.m1pt > pCuts->m1pt)  && (fb.m2pt > pCuts->m2pt));
  fGoodMuonsEta   = ((fb.m1eta > -2.4) && (fb.m1eta < 2.4) && (fb.m2eta > -2.4) && (fb.m2eta < 2.4));
  fGoodTracks     = fb.gtqual;
  fGoodTracksPt   = ((fb.m1pt > 3.0) && (fb.m1pt < 100.0) && (fb.m2pt > 3.0) && (fb.m2pt < 100.));
  fGoodTracksEta  = ((fb.m1eta > -2.4) && (fb.m1eta < 2.4) && (fb.m2eta > -2.4) && (fb.m2eta < 2.4));
  
  fGoodQ          = (fb.m1q*fb.m2q < 0); 
  fGoodPvAveW8    = (fb.pvw8 > 0.7);
  fGoodMaxDoca    = (TMath::Abs(fb.maxdoca) < pCuts->maxdoca); 
  fGoodIp         = (TMath::Abs(fb.pvip) < pCuts->pvip); 
  fGoodIpS        = (TMath::Abs(fb.pvips) < pCuts->pvips); 
  
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

  // -- different preselection cuts used for the plotting of various (sbs) distributions
  //  fGoodHLT        = fb.hlt && fGoodMuonsID && (fb.fls3d > 4) && (fb.chi2/fb.dof < 5) && (fBDT > -1.);

  // -- The following is applied for the data tree and MUST be applied offline as well
  //  fPreselection = fPreselection && fGoodPvLip && fGoodPvLipS && fGoodQ; 
  //  fPreselection = fPreselection && (fCandPt > 5) && (fCandA < 0.2) && (fCandFLS3d > 5) && (fCandChi2/fCandDof < 5); 

  fGoodHLT        = fb.hlt;
  fPreselection   = ((fBDT > 0.) && fb.hlt && fGoodMuonsID ); 

  fAnaCuts.update(); 

}


// ----------------------------------------------------------------------
void plotClass::calcBDT(bool rejectInvIso) {
  fBDT = -99.;

  if (!preselection(fb, fChan, rejectInvIso)) return;

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
  frd.chi2dof = fb.chi2/fb.dof; 
  frd.closetrk = fb.closetrk; 
  
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
	if (stype == "chi2/dof") {
	  cout << "  adding chi2/dof" << endl;
	  reader->AddVariable( "chi2dof := chi2/dof", &rd.chi2dof);
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
