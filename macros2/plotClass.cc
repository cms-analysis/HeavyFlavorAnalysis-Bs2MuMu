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
plotClass::plotClass(const char *files, const char *cuts, const char *dir, int mode) { 

  fFiles = files; 

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  init(files, cuts, dir, mode);

  int NBINS = (fMassHi - fMassLo)/0.025;

  int HBINS(15); 
  double HLO(0.), HHI(45.); 

  TH1D *h; 
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
//   fHistFile->Write();
//   fHistFile->Close();
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
int plotClass::detChan(double m1eta, double m2eta) {
  // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
  if (TMath::Abs(m1eta) < fCuts[0]->etaMax && TMath::Abs(m2eta) < fCuts[0]->etaMax) return 0; 
  if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1; 
  return -1; 
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
  
  bool bp2jpsikp(false), bs2jpsiphi(false), isMC(false); 

  string directory; 
  string fAcc;
  
  double effFilter(1.), genFileYield(0.); 

  numbers *aa(0);
  if (0 == mode) {
    isMC = true; 
    directory = "candAnaMuMu"; 
    fF["SgMc"]->cd(directory.c_str());
    fAcc = "SgMcAcc";
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
  }  else if (1 == mode) {
    isMC = true; 
    directory = "candAnaMuMu"; 
    fF["BdMc"]->cd(directory.c_str());
    fAcc = "BdMcAcc";
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
  } else if (5 == mode) {
    isMC = false;     
    directory = "candAnaMuMu"; 
    fF["SgData"]->cd(directory.c_str());
  } else if (10 == mode) {
    isMC = true; 
    bp2jpsikp = true; 
    directory = "candAnaBu2JpsiK"; 
    fF["NoMc"]->cd(directory.c_str());
    fAcc = "NoMcAcc";
    effFilter = fFilterEff["NoMc"];
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
  } else if (15 == mode) {
    isMC = false;     
    bp2jpsikp = true; 
    directory = "candAnaBu2JpsiK"; 
    fF["NoData"]->cd(directory.c_str());
  } else if (20 == mode) {
    isMC = true; 
    bs2jpsiphi = true; 
    directory = "candAnaBs2JpsiPhi"; 
    fF["CsMc"]->cd(directory.c_str());
    fAcc = "CsMcAcc";
    effFilter = fFilterEff["CsMc"];
    genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();
  } else if (25 == mode) {
    isMC = false;     
    bs2jpsiphi = true; 
    directory = "candAnaBs2JpsiPhi"; 
    fF["CsData"]->cd(directory.c_str());
  } else if (98 == mode) {
    directory = "candAnaMuMu"; 
    cout << "mode 98" << endl;
  } else {
    directory = "candAnaMuMu"; 
    cout << "mode 99" << endl;
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

  cout << "--> loopTree with mode " << mode << " proc = " << proc << " on file ";
  gDirectory->pwd();

  // -- reset all histograms
  for (unsigned int i = 0; i < fNchan; ++i) {
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

    fhNorm[i]->Reset();
    fhNormC[i]->Reset();

    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();
  }

  // -- set up tree
  double mass(0.); 
  TTree *t;
  t = (TTree*)gDirectory->Get("events");
  TEventList *tlist = new TEventList;
  int brr, brun, bevt, bls, btm, bq1, bq2, bprocid; 
  double bg1pt, bg2pt, bg1eta, bg2eta;
  double bbdt, bbdt2; 
  double bm, bcm, bpt, beta, bphi, bcosa, balpha, biso, bchi2, bdof, bdocatrk, bfls3d, bfl3dE, bfl3d;
  double bm1pt, bm1eta, bm2pt, bm2eta, bm1phi, bm2phi;
  double bk1pt, bk1eta, bk2pt, bk2eta; 
  double bg3pt, bg3eta, bg4pt, bg4eta; 
  double bptpsi, bmpsi, bmkk, bdr;
  double bw8mu, bw8tr;
  bool bhlt, bgmuid, bgtqual, bjson, bcb;
  double tr1w8(0.), tr2w8(0.), trw8(0.), m1w8(0.), m2w8(0.), mw8(0.0);

  double blip, blipE, btip, btipE; 
  int bm1pix, bm2pix, bm1bpix, bm2bpix, bm1bpixl1, bm2bpixl1;

  int bclosetrk; 
  double bpvlip, bpvlips, bpvlip2, bpvlips2, bmaxdoca, bpvip, bpvips, bpvw8; 

  t->SetBranchAddress("bdt",&bbdt);
  t->SetBranchAddress("bdt2",&bbdt2);
  t->SetBranchAddress("lip",&blip);
  t->SetBranchAddress("lipE",&blipE);
  t->SetBranchAddress("tip",&btip);
  t->SetBranchAddress("tipE",&btipE);

  t->SetBranchAddress("closetrk",&bclosetrk);
  t->SetBranchAddress("pvlip",&bpvlip);
  t->SetBranchAddress("pvlips",&bpvlips);
  t->SetBranchAddress("pvlip2",&bpvlip2);
  t->SetBranchAddress("pvlips2",&bpvlips2);
  t->SetBranchAddress("maxdoca",&bmaxdoca);
  t->SetBranchAddress("pvip",&bpvip);
  t->SetBranchAddress("pvips",&bpvips);
  t->SetBranchAddress("pvw8",&bpvw8);

  t->SetBranchAddress("m1pix",&bm1pix);
  t->SetBranchAddress("m2pix",&bm2pix);
  t->SetBranchAddress("m1bpix",&bm1bpix);
  t->SetBranchAddress("m2bpix",&bm2bpix);
  t->SetBranchAddress("m1bpixl1",&bm1bpixl1);
  t->SetBranchAddress("m2bpixl1",&bm2bpixl1);

  t->SetBranchAddress("rr",&brr);
  t->SetBranchAddress("run",&brun);
  t->SetBranchAddress("evt",&bevt);
  t->SetBranchAddress("hlt",&bhlt);
  t->SetBranchAddress("ls",&bls);
  t->SetBranchAddress("cb",&bcb);
  t->SetBranchAddress("json",&bjson);
  t->SetBranchAddress("gmuid",&bgmuid);
  t->SetBranchAddress("gtqual",&bgtqual);
  t->SetBranchAddress("w8mu",&bw8mu);
  t->SetBranchAddress("w8tr",&bw8tr);
  t->SetBranchAddress("tm",&btm);
  t->SetBranchAddress("procid",&bprocid);
  t->SetBranchAddress("m",&bm);
  t->SetBranchAddress("cm",&bcm);
  t->SetBranchAddress("pt",&bpt);
  t->SetBranchAddress("phi",&bphi);
  t->SetBranchAddress("eta",&beta);
  t->SetBranchAddress("cosa",&bcosa);
  t->SetBranchAddress("alpha",&balpha);
  t->SetBranchAddress("iso",&biso);
  t->SetBranchAddress("chi2",&bchi2);
  t->SetBranchAddress("dof",&bdof);
  t->SetBranchAddress("fls3d",&bfls3d);
  t->SetBranchAddress("fl3d",&bfl3d);
  t->SetBranchAddress("fl3dE",&bfl3dE);
  t->SetBranchAddress("m1pt",&bm1pt);
  t->SetBranchAddress("m1eta",&bm1eta);
  t->SetBranchAddress("m1phi",&bm1phi);
  t->SetBranchAddress("m1q",&bq1);
  t->SetBranchAddress("m2pt",&bm2pt);
  t->SetBranchAddress("m2eta",&bm2eta);
  t->SetBranchAddress("m2phi",&bm2phi);
  t->SetBranchAddress("m2q",&bq2);
  t->SetBranchAddress("docatrk",&bdocatrk);

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);
  if (bp2jpsikp) {
    if (isMC) {
      t->SetBranchAddress("g3pt",&bg3pt);
      t->SetBranchAddress("g3eta",&bg3eta);
    }
    t->SetBranchAddress("kpt",&bk1pt);
    t->SetBranchAddress("keta",&bk1eta);
    t->SetBranchAddress("mpsi",&bmpsi);
    t->SetBranchAddress("psipt",&bptpsi); //FIXME
  }

  if (bs2jpsiphi) {
    if (isMC) {
      t->SetBranchAddress("g3pt",&bg3pt);
      t->SetBranchAddress("g3eta",&bg3eta);
      t->SetBranchAddress("g4pt",&bg4pt);
      t->SetBranchAddress("g4eta",&bg4eta);
    }
    t->SetBranchAddress("psipt",&bptpsi);   //FIXME
    t->SetBranchAddress("mpsi",&bmpsi);
    t->SetBranchAddress("mkk",&bmkk);
    t->SetBranchAddress("dr",&bdr);
    t->SetBranchAddress("k1pt",&bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k2pt",&bk2pt);
    t->SetBranchAddress("k2eta",&bk2eta);
  } else {
    bmkk = 999.;
    bdr = 999.;
  }

  int nentries = Int_t(t->GetEntries());
  int nb(0), ievt(0), bsevt(0), bdevt(0), bgevt(0); 
  cuts *pCuts(0); 
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);
    // -- require truth matching when on MC
    if (0 == mode && 0 == btm) continue;
    if (1 == mode && 0 == btm) continue;
    if (10 == mode && 0 == btm) continue;
    if (20 == mode && 0 == btm) continue;

    if (isMC && proc > 0) {
      if (bprocid != proc) continue;
    }

    // -- channel index
    fChan = -1; 
    pCuts = 0; 
    
    fChan = detChan(bm1eta, bm2eta); 
    if (fChan > -1) {
      pCuts = fCuts[fChan]; 
    } else {
      //       cout << "event " << jentry << ", fChan = " << fOver << " for eta(B) = " << beta 
      // 	   << " m1eta = " << bm1eta << " m2eta = " << bm2eta
      // 	   << endl;
      continue;
    }

    mass = bm; 
//     if (bs2jpsiphi || bp2jpsikp) {
//       mass = bcm; 
//     }     

    fhMassAbsNoCuts[fChan]->Fill(mass);
    // -- require wide mass window ??? FIXME WHY????
    //     if (mass < fMassLo) continue;
    //     if (fMassHi < mass) continue;

    // -- gen-level acceptance cuts
    if (isMC) {
      if (TMath::Abs(bg1eta) > 2.5) continue;
      if (TMath::Abs(bg2eta) > 2.5) continue;
      if (bg1pt < 1.0) continue;
      if (bg2pt < 1.0) continue;
      if (bp2jpsikp) {
	// gen-level cuts for Bu2JpsiKp
	if (bg1pt < 3.5) continue;
	if (bg2pt < 3.5) continue;
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (bg3pt < 0.4) continue;
      }
      
      if (bs2jpsiphi) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (TMath::Abs(bg4eta) > 2.5) continue;
	// gen-level cuts for Bs2JpsiPhi
	if (bg1pt < 3.5) continue;
	if (bg2pt < 3.5) continue;
	if (bg3pt < 0.4) continue;
	if (bg4pt < 0.4) continue;
      }
    } else {
      if (!bjson) {
	//	if (5 == mode) cout << "skipping run: " << brun << " LS: " << bls << endl;
	continue;
      }
    }

    // -- immutable cuts: require basic muon and trackQual cuts
    if (false == bgtqual)  continue;
    if (TMath::Abs(bm1eta) > 2.4) continue;
    if (TMath::Abs(bm2eta) > 2.4) continue;

    if (bm1pt < 1.0) continue;
    if (bm2pt < 1.0) continue;

    if (bp2jpsikp) {
      if (TMath::Abs(bk1eta) > 2.4) continue;
      if (bk1pt < 0.5) continue;
    }
    if (bs2jpsiphi) {
      if (TMath::Abs(bk1eta) > 2.4) continue;
      if (TMath::Abs(bk2eta) > 2.4) continue;
      if (bk1pt < 0.5) continue;
      if (bk2pt < 0.5) continue;
    }

    // -- this is the base, after the raw acceptance cuts
    fhMassNoCuts[fChan]->Fill(mass);
    fhMassNoCutsManyBins[fChan]->Fill(mass); 

    if (fDoUseBDT) {
      // Gemma's cuts:
      // //for MC
      //  mycuts="(m1pt>4.5) && (m2pt>4) && (pt>5) && (alpha<0.2) && (fls3d>5)  && (json) && (gmupt) && (gmueta) 
      //  && gteta && gtpt && gtqual  && (m1q*m2q<0) && (m>4.9) && (m<5.9) && (tm>0) && fabs(m1eta)<2.4 && fabs(m2eta)<2.4
      //  && fabs(pvlip)<0.015 && fabs(pvlips)<3 && fl3d<2";

      // //for Data   (should be the same but hlt and gmuid required in addition)
      //  mycutb="hlt  && (m1pt>4.5) && (m2pt>4) && (pt>5) && (alpha<0.2) && (fls3d>5) && json  && gmuid
      //  && gmupt && gmueta && gteta && gtpt && gtqual  && (m1q*m2q<0) 
      //  && ( ( (m>4.9) && (m<=5.2) )|| ((m>=5.45) && (m<5.9))) && fabs(m1eta)<2.4 && fabs(m2eta)<2.4
      //  && fabs(pvlip)<0.015 && fabs(pvlips)<3 && fl3d<2";

      if (bm1pt < 4.5) continue; 
      if (bm2pt < 4.0) continue; 
      if (bq1*bq2 > 0) continue;
      if (bfl3d > 1.5) continue;
      if (bfls3d < 5) continue;
      if (balpha > 0.2) continue;
      if (bpvlip > 0.015) continue;
      if (bpvlips > 3) continue;
      if (!bjson) continue;

      // -- skip inverted isolation events
      if (5 == mode && 5.2 < mass && mass < 5.45 && biso < 0.7) continue; 
      
      if(5 == mode || 15 == mode || 25 == mode) {
	if (false == bhlt) continue;
	if (false == bgmuid) continue;
      }

      if (TMath::IsNaN(bfls3d)) continue;
      if (bs2jpsiphi && bdr >0.3) continue;
      if (bs2jpsiphi && bmkk < 0.995) continue;
      if (bs2jpsiphi && bmkk > 1.045) continue;

      if (bs2jpsiphi || bp2jpsikp) {
	if (bmpsi > 3.2) continue;
	if (bmpsi < 3.0) continue;
	// -- cowboy veto 
	if (fDoApplyCowboyVeto && bcb) continue;
	if (bptpsi < 7) continue;
      } 

      if (0 == fChan && bbdt < pCuts->bdt) continue;
      if (1 == fChan && bbdt2 < pCuts->bdt) continue;
    } else {
    
      // -- analysis cuts
      if (bq1*bq2 > 0) continue;
      if (bm1pt < pCuts->m1pt) continue; 
      if (bm2pt < pCuts->m2pt) continue; 
      
      if (bfl3d > 1.5) continue;
      if (bpvw8 < 0.6) continue;
      
      if (bpt < pCuts->pt) continue; 
      if (biso < pCuts->iso) continue; 
      if (bchi2/bdof > pCuts->chi2dof) continue;
      if (TMath::IsNaN(bfls3d)) continue;
      if (bfls3d < pCuts->fls3d) continue;
      if (balpha > pCuts->alpha) continue;
      if (bdocatrk < pCuts->docatrk) continue;
      
      if (bs2jpsiphi && bdr >0.3) continue;
      if (bs2jpsiphi && bmkk < 0.995) continue;
      if (bs2jpsiphi && bmkk > 1.045) continue;
      
      // -- new cuts
      if (bclosetrk >= pCuts->closetrk) continue;
      if (TMath::Abs(bpvlip) > pCuts->pvlip) continue;
      if (TMath::Abs(bpvlips) > pCuts->pvlips) continue;
      if (TMath::Abs(bpvlip2) < pCuts->pvlip2) continue;
      if (TMath::Abs(bpvlips2) < pCuts->pvlips2) continue;
      if (bmaxdoca > pCuts->maxdoca) continue;
      if (bpvips > pCuts->pvips) continue;
      if (bpvip > pCuts->pvip) continue;
      
      if (bs2jpsiphi || bp2jpsikp) {
	if (bmpsi > 3.2) continue;
	if (bmpsi < 3.0) continue;
	// -- cowboy veto 
	if (fDoApplyCowboyVeto && bcb) continue;
	if (bptpsi < 7) continue;
      } 
      
      if (fDoApplyCowboyVetoAlsoInSignal && bcb) continue;
    }


    fhMassWithAnaCuts[fChan]->Fill(mass); 
    fhMassWithAnaCutsManyBins[fChan]->Fill(mass); 

    // ----------------------------
    // -- Intermezzo with PidTables
    // ----------------------------

    // -- muon ID: Data PidTables
    m1w8 = ptM->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    m2w8 = ptM->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    mw8  = m1w8*m2w8; 
    fhMuId[fChan]->Fill(mw8); 
    fh0TNPMuID[fChan]->Fill(bpt, mw8); 
    fh1TNPMuID[fChan]->Fill(bpt); 

    // -- MC for comparison 
    m1w8 = ptMMC->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    m2w8 = ptMMC->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    mw8  = m1w8*m2w8; 
    if (mw8 > 0.) {
      fhMuIdMC[fChan]->Fill(mw8, 1./mw8); 
    }
    fh0TNPMCMuID[fChan]->Fill(bpt, mw8); 
    fh1TNPMCMuID[fChan]->Fill(bpt); 


    // -- muon trigger: Data PidTables
    tr1w8 = ptT1->effD(bm1pt, TMath::Abs(bm1eta), 0.)*ptT2->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    tr2w8 = ptT1->effD(bm2pt, TMath::Abs(bm2eta), 0.)*ptT2->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (bs2jpsiphi || bp2jpsikp) {
      if (brr >= 3) {
	if (TMath::Abs(bm1eta) > 2.2) trw8 = 0; 
	if (TMath::Abs(bm2eta) > 2.2) trw8 = 0; 
      }
    }
    fhMuTr[fChan]->Fill(trw8); 
    fh0TNPTrigger[fChan]->Fill(bpt, trw8); 
    fh1TNPTrigger[fChan]->Fill(bpt); 

    // -- MC for comparison 
    tr1w8 = ptT1MC->effD(bm1pt, TMath::Abs(bm1eta), 0.)*ptT2->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    tr2w8 = ptT1MC->effD(bm2pt, TMath::Abs(bm2eta), 0.)*ptT2->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (bs2jpsiphi || bp2jpsikp) {
      if (brr >= 3) {
	if (TMath::Abs(bm1eta) > 2.2) trw8 = 0; 
	if (TMath::Abs(bm2eta) > 2.2) trw8 = 0; 
      }
    }
    fhMuTrMC[fChan]->Fill(trw8); 
    fh0TNPMCTrigger[fChan]->Fill(bpt, trw8); 
    fh1TNPMCTrigger[fChan]->Fill(bpt); 


    // ---------------------------
    // -- now to the MC simulation
    // ---------------------------

    // -- MUON ID
    if (false == bgmuid) continue;
    fhMassWithMuonCuts[fChan]->Fill(mass); 
    fhMassWithMuonCutsManyBins[fChan]->Fill(mass); 

    // -- TRIGGER
    fh0MCTrigger[fChan]->Fill(bpt); 
    if (false == bhlt) continue;
    fh1MCTrigger[fChan]->Fill(bpt); 
    fhMassWithTriggerCuts[fChan]->Fill(mass); 
    fhMassWithTriggerCutsManyBins[fChan]->Fill(mass); 

    fhMassWithAllCuts[fChan]->Fill(mass); 
    if (5 == mode && !(5.2 < mass && mass < 5.45)) {
      fhMassWithAllCutsBlind[fChan]->Fill(mass); 
    }

    fhMassWithAllCutsManyBins[fChan]->Fill(mass); 

    fhNorm[fChan]->Fill(mass);
    fhNormC[fChan]->Fill(bcm);
    //    fhNorm[fChan]->Fill(mass);
    
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

    if (fDoPrint && 5 == mode && mass > 4.9 && mass < 5.9) {
      tlist->Enter(jentry); 
      cout << Form("%d m = %4.3f pT = %3.1f eta = %3.2f", fChan, mass, bpt, beta)
	//	   <<	" r = " << brun << "/" << bevt
	   << Form(" mpt = %3.1f,%3.1f", bm1pt, bm2pt)
	   << Form(" meta = %3.2f,%3.2f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %3.2f chi2 = %3.1f fls3d = %3.1f, fl/E=%3.1f/%3.2f", 
		   TMath::ACos(bcosa), biso, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
	   << Form(" pv1: %4.3f/%3.1f", bpvlip, bpvlips) 
	   << Form(" d/s: %5.4f/%3.2f md: %4.3f d: %4.3f", bpvip, bpvips, bmaxdoca, bdocatrk) 
	   << endl;
      fOUT << Form("%d m = %4.3f pT = %3.1f eta = %3.2f", fChan, mass, bpt, beta)
	//	   <<	" run = " << brun << " event = " << bevt
	   << Form(" mpt = %3.1f,%3.1f", bm1pt, bm2pt)
	   << Form(" meta = %3.2f,%3.2f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %3.2f chi2 = %3.1f fls3d = %3.1f, fl/E=%3.1f/%3.2f", 
		   TMath::ACos(bcosa), biso, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
	   << Form(" pv1: %4.3f/%3.1f", bpvlip, bpvlips) 
	   << Form(" d/s: %5.4f/%3.2f md: %4.3f d: %4.3f", bpvip, bpvips, bmaxdoca, bdocatrk) 
	   << endl;


      string st("SgEvt"); 

      if (pCuts->mBsLo < mass  && mass < pCuts->mBsHi) {
	st = "BsSgEvt"; 
	ievt = bsevt;
	++bsevt; 
      } else if (pCuts->mBdLo < mass  && mass < pCuts->mBdHi) {
	st = "BdSgEvt"; 
	ievt = bdevt;
	++bdevt; 
      } else {
	ievt = bgevt;
	++bgevt; 
      }

      if (0 == fChan) st += "0";
      if (1 == fChan) st += "1";

      string suffix(fSuffix); 
      if (fDoUseBDT) suffix += "Bdt"; 

      fTEX << formatTex(brun,      Form("%s:%s%i:run", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bevt,      Form("%s:%s%i:evt", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(fChan,     Form("%s:%s%i:chan", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm,        Form("%s:%s%i:m", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bpt,       Form("%s:%s%i:pt", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bphi,      Form("%s:%s%i:phi", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(beta,      Form("%s:%s%i:eta", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << Form("\\vdef{%s:%s%i:channel}   {%s }", suffix.c_str(), st.c_str(), ievt, fChan==0?"barrel":"endcap") << endl;
      fTEX << formatTex((bcb?1:0),    Form("%s:%s%i:cowboy", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm1pt,     Form("%s:%s%i:m1pt", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm2pt,     Form("%s:%s%i:m2pt", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm1eta,    Form("%s:%s%i:m1eta", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm2eta,    Form("%s:%s%i:m2eta", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm1phi,    Form("%s:%s%i:m1phi", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm2phi,    Form("%s:%s%i:m2phi", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bq1,       Form("%s:%s%i:m1q", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bq2,       Form("%s:%s%i:m2q", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(biso,      Form("%s:%s%i:iso", suffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(balpha,    Form("%s:%s%i:alpha", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bchi2,     Form("%s:%s%i:chi2", suffix.c_str(), st.c_str(), ievt), 2) << endl;
      fTEX << formatTex(bdof,      Form("%s:%s%i:dof", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bfls3d,    Form("%s:%s%i:fls3d", suffix.c_str(), st.c_str(), ievt), 2) << endl;
      fTEX << formatTex(bfl3d,     Form("%s:%s%i:fl3d", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bfl3dE,    Form("%s:%s%i:fl3dE", suffix.c_str(), st.c_str(), ievt), 4) << endl;

      fTEX << formatTex(bdocatrk,  Form("%s:%s%i:docatrk", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(blip,      Form("%s:%s%i:lip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(blipE,     Form("%s:%s%i:lipE", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(btip,      Form("%s:%s%i:tip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(btipE,     Form("%s:%s%i:tipE", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bpvlip,    Form("%s:%s%i:pvlip", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bpvlips,   Form("%s:%s%i:pvlips", suffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bmaxdoca,  Form("%s:%s%i:maxdoca", suffix.c_str(), st.c_str(), ievt), 4) << endl;

      fTEX << formatTex(bm1pix,    Form("%s:%s%i:m1pix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm2pix,    Form("%s:%s%i:m2pix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm1bpix,   Form("%s:%s%i:m1bpix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm2bpix,   Form("%s:%s%i:m2bpix", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm1bpixl1, Form("%s:%s%i:m1bpixl1", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm2bpixl1, Form("%s:%s%i:m2bpixl1", suffix.c_str(), st.c_str(), ievt), 0) << endl;
      
    }
    
  }

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
    }
    pD->cd();

    return;
  }

  //   c0->Clear();
  //   c0->Divide(1,2);
  for (unsigned int i = 0; i < fNchan; ++i) {
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
      h->Draw();
      //       drawArrow(0.5, 2); 
      //       drawArrow(0.5, 1); 
      drawBox(2, 0.5); 
      drawBox(1, 0.5); 
      
      TH1D *dummy1 = new TH1D("dummy1", "", 10, 0., 10.); setFilledHist(dummy1, kBlue, kBlue, 3356); 
      TH1D *dummy2 = new TH1D("dummy2", "", 10, 0., 10.); setFilledHist(dummy2, kRed, kRed, 3365); 
      
      newLegend(0.4, 0.65, 0.8, 0.8); 
      legg->SetTextSize(0.045);  
      legg->AddEntry(dummy1, "B_{s}^{0} signal window", "f"); 
      legg->AddEntry(dummy2, "B^{0} signal window", "f"); 
      legg->Draw();
      
      //      stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
      if (fDoPrint)  {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtsig-data-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/sig-data-chan%d.pdf", fDirectory.c_str(), i));
      }
      // -- unblinded version
      h = fhMassWithAllCuts[i]; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#mu} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.3); 
      h->SetMinimum(0.01); 
      h->SetNdivisions(003, "Y");
      h->SetAxisRange(4.9, 5.9, "X"); 
      //      h->SetMaximum(2.2);
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->Draw();
      drawArrow(0.6, 1); 
      drawArrow(0.4, 2); 

      tl->SetNDC(kTRUE); 
      tl->SetTextSize(0.07); 
      if (0 == i) {
	tl->DrawLatex(0.6, 0.8, "Barrel");   
      } 
      
      if (1 == i) {
	tl->DrawLatex(0.6, 0.8, "Endcap");   
      } 
      

//       drawBox(2, 0.5); 
//       drawBox(1, 0.5); 

//       newLegend(0.4, 0.65, 0.8, 0.8); 
//       legg->SetTextSize(0.045);  
//       legg->AddEntry(dummy1, "B_{s}^{0} signal window", "f"); 
//       legg->AddEntry(dummy2, "B^{0} signal window", "f"); 
//       legg->Draw();


//      stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
      if (fDoPrint) {
	if (fDoUseBDT) c0->SaveAs(Form("%s/bdtunblinded-sig-data-chan%d.pdf", fDirectory.c_str(), i));
	else c0->SaveAs(Form("%s/unblinded-sig-data-chan%d.pdf", fDirectory.c_str(), i));
      }

      delete dummy1; 
      delete dummy2; 

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
      normYield(h, i, 4.8, 6.0);
      aa->fitYield  = fNoSig; 
      aa->fitYieldE = fNoSigE; 
      h = fhNormC[i];
      normYield(h, 0, 5.16, 5.6, 4.15);
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
      csYield(h, i, 5.1, 6.0);
      aa->fitYield  = fCsSig; 
      aa->fitYieldE = fCsSigE; 
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


    printNumbers(*aa, cout); 
    printNumbers(*aa, fOUT); 
    
    // -- Cache the pwd...
    TDirectory *pD = gDirectory; 
    fHistFile->cd();
    string modifier = (fDoUseBDT?"bdt":"cnc"); 
    //    fhMassWithMassCutsManyBins[i]->SetName(Form("hMassWithMassCutsManyBins%d_chan%d", mode, i)); fhMassWithMassCutsManyBins[i]->Write();
    //    fhMassWithMassCuts[i]->SetName(Form("hMassWithMassCuts%d_chan%d", mode, i)); fhMassWithMassCuts[i]->Write();
    fhMassWithAllCutsManyBins[i]->SetName(Form("hMassWithAllCutsManyBins_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhMassWithAllCutsManyBins[i]->SetTitle(Form("hMassWithAllCutsManyBins_%s_%d_chan%d %s", modifier.c_str(), mode, i, pD->GetName())); 
    fhMassWithAllCutsManyBins[i]->Write();
    fhMassWithAllCuts[i]->SetName(Form("hMassWithAllCuts_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhMassWithAllCuts[i]->Write();
    fhNorm[i]->SetName(Form("hNorm_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhNorm[i]->Write();
    fhNormC[i]->SetName(Form("hNormC_%s_%d_chan%d", modifier.c_str(), mode, i)); 
    fhNormC[i]->Write();

    if (5 == mode) {
      t->SetEventList(tlist);
      TTree *small = t->CopyTree(""); 
      t->SetName(Form("%s", (fDoUseBDT?"bdt":"cnc")));
      small->Write(); 
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
    a.effCandE = 0.01; 
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
  cout << "--> Dumping output into " << fDirectory << endl;
  fNumbersFileName = fDirectory + "/anaBmm." + fSuffix + ".txt";
  // FIXME
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fOUT.open(fNumbersFileName.c_str(), ios::app);
  
  loadFiles(files);
  string hfname  = fDirectory + "/anaBmm." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "RECREATE");

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
    } else {
      string sfilter = sdset; 
      replaceAll(sfilter, "mc,", ""); 
      double effFilter = atof(sfilter.c_str());
      // -- MC
      pF = loadFile(sfile); 
      if (string::npos != stype.find("default") && string::npos != stype.find("sg")) {
	sname = "SgMc"; 
	fF.insert(make_pair(sname, pF)); 	
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("sg")) {
	sname = "SgMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("sg")) {
	sname = "SgMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("sg")) {
	sname = "SgMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("pu") && string::npos != stype.find("sg")) {
	sname = "SgMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (PU)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("sg")) {
	sname = "SgMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #mu^{+}#mu^{-} (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("bd")) {
	sname = "BdMc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("bd")) {
	sname = "BdMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("bd")) {
	sname = "BdMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("bd")) {
	sname = "BdMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("bd")) {
	sname = "BdMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #mu^{+}#mu^{-} (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }

      if (string::npos != stype.find("default") && string::npos != stype.find("no")) {
	sname = "NoMc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("2e33") && string::npos != stype.find("no")) {
	sname = "NoMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("1e33") && string::npos != stype.find("no")) {
	sname = "NoMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("3e33") && string::npos != stype.find("no")) {
	sname = "NoMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("acc") && string::npos != stype.find("no")) {
	sname = "NoMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }
      if (string::npos != stype.find("pu") && string::npos != stype.find("no")) {
	sname = "NoMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{+} #rightarrow J/#psi K^{+} (PU)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("default") && string::npos != stype.find("cs")) {
	sname = "CsMc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (MC)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("1e33") && string::npos != stype.find("cs")) {
	sname = "CsMc1e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (1e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("2e33") && string::npos != stype.find("cs")) {
	sname = "CsMc2e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (2e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("3e33") && string::npos != stype.find("cs")) {
	sname = "CsMc3e33"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (3e33)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("acc") && string::npos != stype.find("cs")) {
	sname = "CsMcAcc"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (acc)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("pu") && string::npos != stype.find("cs")) {
	sname = "CsMcPU"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow J/#psi #phi (PU)")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bs2KK")) {
	sname = "bgBs2KK"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2KPi")) {
	sname = "bgBs2KPi"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2PiPi")) {
	sname = "bgBs2PiPi"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow #pi^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bs2KMuNu")) {
	sname = "bgBs2KMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B_{s}^{0} #rightarrow K^{-}#mu^{+}#nu")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	

      if (string::npos != stype.find("bg,Bd2PiMuNu")) {
	sname = "bgBd2PiMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{-}#mu^{+}#nu")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2KK")) {
	sname = "bgBd2KK"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2KPi")) {
	sname = "bgBd2KPi"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow K^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Bd2PiPi")) {
	sname = "bgBd2PiPi"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "B^{0} #rightarrow #pi^{+}#pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2KP")) {
	sname = "bgLb2KP"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p K^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2PiP")) {
	sname = "bgLb2PiP"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda_{b}^{0} #rightarrow p #pi^{-}")); 
	fFilterEff.insert(make_pair(sname, effFilter)); 
      }	
      if (string::npos != stype.find("bg,Lb2PMuNu")) {
	sname = "bgLb2PMuNu"; 
	fF.insert(make_pair(sname, pF)); 
	fLumi.insert(make_pair(sname, atof(slumi.c_str()))); 
	fName.insert(make_pair(sname, "#Lambda^{0}_{b} #rightarrow p#mu^{-}#bar{#nu}")); 
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

  double pReco = (preco<0? (1 == mode?5.145:5.146): preco);

  TF1 *lF1(0), *lBg(0);
  
  if (0 == mode) { 
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoErrgauss2c(h, 5.27, 0.034, pReco); 
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  } else {
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoErrGauss(h, 5.27, 0.056, pReco); 
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  }
  h->Fit(lF1, "rm", "", lo, hi); 

  if (0 == mode) {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(5+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  } else {
    for (int i = 0; i < lBg->GetNpar(); ++i) {
      lBg->SetParameter(i, lF1->GetParameter(3+i));
      cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
    }
  }
  double c  = lF1->GetParameter(0); 
  c = lF1->Integral(5.15, 5.4) - lBg->Integral(5.15, 5.4); 
  double cE = lF1->GetParError(0); 
  double ierr = lF1->IntegralError(5.15, 5.4)/h->GetBinWidth(1); 

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



  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muK} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.6); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  // -- Overlay BG function
  lBg->SetLineStyle(kDashed);
  lBg->SetLineColor(kRed);
  lBg->SetLineWidth(3);
  lBg->Draw("same");

//   tl->SetTextSize(0.07); 
//   if (0 == mode) {
//     tl->DrawLatex(0.6, 0.8, "Barrel");   
//   } 

//   if (1 == mode) {
//     tl->DrawLatex(0.6, 0.8, "Endcap");   
//   } 

  //  stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint) {
    
    string pdfname;
    string hname(h->GetName());
    if (string::npos != hname.find("NormC")) {
      pdfname = Form("%s/normC-data-chan%d.pdf", fDirectory.c_str(), mode);
    } else {
      pdfname = Form("%s/norm-data-chan%d.pdf", fDirectory.c_str(), mode);
    }

    if (fDoUseBDT)  pdfname = Form("%s/bdtnorm-data-chan%d.pdf", fDirectory.c_str(), mode);
    c0->SaveAs(pdfname.c_str());
  }
  

  //   double c  = h->GetFunction("f1")->GetParameter(0); 
  //   double cE = h->GetFunction("f1")->GetParError(0); 
  
  delete lF1; 

}


// ----------------------------------------------------------------------
void plotClass::csYield(TH1 *h, int mode, double lo, double hi, double preco) {

  double pReco = (preco<0? 5.2: preco);

  TF1 *lF1(0), *lBg(0);

  if (0 == mode) { 
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoGauss(h, 5.37, 0.030); 
    lBg = fpFunc->expoErr(fpFunc->fLo, fpFunc->fHi); 
  } else {
    fpFunc->fLo = lo; //5.0;
    fpFunc->fHi = hi; //5.5;
    lF1 = fpFunc->expoGauss(h, 5.37, 0.056); 
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
  if (0 == mode) {
    tl->DrawLatex(0.22, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.22, 0.8, "Endcap");   
  } 


  //  stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint) {
    if (fDoUseBDT) c0->SaveAs(Form("%s/bdtcs-data-chan%d.pdf", fDirectory.c_str(), mode));
    else c0->SaveAs(Form("%s/cs-data-chan%d.pdf", fDirectory.c_str(), mode));
  }

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 

}


// ----------------------------------------------------------------------
void plotClass::printNumbers(numbers &a, ostream &OUT) {
  OUT << "======================================================================" << endl;
  OUT << "numbers for \""  << a.name.c_str() << "\"" << endl;
  OUT << "fitYield        = " << a.fitYield << "+/-" << a.fitYieldE << endl;
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
void plotClass::newLegend(double x1, double y1, double x2, double y2) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2);
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

  char CutName[100];
  float CutValue;
  int dump(0), ok(0);

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

  }

  if (a) fCuts.push_back(a); 

  if (!ok) cout << "==> what about " << CutName << endl;
  

}


// ----------------------------------------------------------------------
void plotClass::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    cuts *a = fCuts[i]; 
    OUT << "# -- channel " << a->index << endl;
    OUT << "index    " << a->index << endl;
    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% -- Cuts for channel " << a->index << endl;

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
  //  tl->DrawLatex(x2, 0.85, text2.c_str()); 
}



// ----------------------------------------------------------------------
void plotClass::drawArrow(double height, int mode, int color) {

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
      y = 1.75;
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
      y = 1.55;
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
