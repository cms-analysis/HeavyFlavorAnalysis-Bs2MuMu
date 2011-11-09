#include "plotReducedTree.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std; 
using std::string; 

ClassImp(plotReducedTree)


// ----------------------------------------------------------------------
plotReducedTree::plotReducedTree(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  int NBINS = (fMassHi - fMassLo)/0.025;

  int HBINS(15); 
  double HLO(0.), HHI(45.); 

  TH1D *h; 
  numbers *a(0); 
  for (unsigned int i = 0; i < fNchan; ++i) {
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

    h = new TH1D(Form("hNorm%d", i), Form("hNorm%d", i), 60, 5.0, 5.6);
    fhNorm.push_back(h); 


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

    h = new TH1D(Form("h0PidTrigger%d", i), Form("hPidTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0PidTrigger.push_back(h); 
    h = new TH1D(Form("h1PidTrigger%d", i), Form("hPidTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1PidTrigger.push_back(h); 
    h = new TH1D(Form("h0PidMuID%d", i), Form("hPidMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0PidMuID.push_back(h); 
    h = new TH1D(Form("h1PidMuID%d", i), Form("hPidMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1PidMuID.push_back(h); 

    h = new TH1D(Form("h0PidMCTrigger%d", i), Form("hPidMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0PidMCTrigger.push_back(h); 
    h = new TH1D(Form("h1PidMCTrigger%d", i), Form("hPidMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1PidMCTrigger.push_back(h); 
    h = new TH1D(Form("h0PidMCMuID%d", i), Form("hPidMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0PidMCMuID.push_back(h); 
    h = new TH1D(Form("h1PidMCMuID%d", i), Form("hPidMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1PidMCMuID.push_back(h); 


    h = new TH1D(Form("h0MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI); h->Sumw2();
    fh0MCTrigger.push_back(h); 
    h = new TH1D(Form("h1MCTrigger%d", i), Form("hMCTrigger%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCTrigger.push_back(h); 
    h = new TH1D(Form("h0MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh0MCMuID.push_back(h); 
    h = new TH1D(Form("h1MCMuID%d", i), Form("hMCMuID%d", i), HBINS, HLO, HHI);  h->Sumw2();
    fh1MCMuID.push_back(h); 
  }


}

// ----------------------------------------------------------------------
plotReducedTree::~plotReducedTree() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedTree::makeAll(int channels) {
  


}




// ----------------------------------------------------------------------
void plotReducedTree::loopTree(int mode, int proc) {
  // -- mode:
  // 0  Bs2MuMu MC
  // 1  Bd2MuMu MC
  // 5  Bs2MuMu data
  // 10 Bp2JpsiKp MC
  // 11 Bp2JpsiKp data
  // 20 Bs2JpsiPhi MC
  // 21 Bs2JpsiPhi data

  PidTable *ptT1;
  PidTable *ptT2;
  PidTable *ptM; 

  PidTable *ptT1MC;
  PidTable *ptT2MC;
  PidTable *ptMMC; 
  
  bool bp2jpsikp(false), bs2jpsiphi(false), isMC(false); 

  string fAcc;

  numbers *aa(0);
  if (0 == mode) {
    isMC = true; 
    fF["SgMc"]->cd();
    //    fpMc[fSgMc]->cd(); 
    fAcc = "SgMcAcc";
  }  else if (1 == mode) {
    isMC = true; 
    fF["BdMc"]->cd();
    //    fpMc[fBdMc]->cd(); 
    fAcc = "BdMcAcc";
  } else if (5 == mode) {
    isMC = false;     
    fF["SgData"]->cd();
    //    fpData[fSgData]->cd(); 
  } else if (10 == mode) {
    isMC = true; 
    bp2jpsikp = true; 
    fF["NoMc"]->cd();
    //    fpMc[fNoMc]->cd(); 
    fAcc = "NoMcAcc";
  } else if (11 == mode) {
    isMC = false;     
    bp2jpsikp = true; 
    fF["NoData"]->cd();
    //    fpData[fNoData]->cd(); 
  } else if (20 == mode) {
    isMC = true; 
    bs2jpsiphi = true; 
    fF["CsMc"]->cd();
    //    fpMc[fCsMc]->cd(); 
    fAcc = "CsMcAcc";
  } else if (21 == mode) {
    isMC = false;     
    bs2jpsiphi = true; 
    fF["CsData"]->cd();
    //    fpData[fCsData]->cd(); 
  } else if (98 == mode) {
    cout << "mode 98" << endl;
  } else {
    cout << "mode 99" << endl;
  }

  ptT1 = new PidTable("pidtables/110606/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 	
  ptT2 = new PidTable("pidtables/110606/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 	
  ptM  = new PidTable("pidtables/110606/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_data_all.dat"); 

  ptT1MC = new PidTable("pidtables/110625/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 	
  ptT2MC = new PidTable("pidtables/110625/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 	
  ptMMC  = new PidTable("pidtables/110625/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 

  //   if (isMC) {
  //     ptT1 = new PidTable("pidtables/110625/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 	
  //     ptT2 = new PidTable("pidtables/110625/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 	
  //     ptM  = new PidTable("pidtables/110625/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_mc_MC.dat"); 
  //     //    ptM  = new PidTable("pidtables/110701/H2D_TagPt_MuonIDEfficiency_GlbTM_TagHighPt_mc.dat"); 
  //     //     ptT1 = new PidTable("pidtables/110625/H2D_L1L2Efficiency_GlbTM_ProbeTrackMatched_mc_MCTRUTH.dat"); 	
  //     //     ptT2 = new PidTable("pidtables/110625/H2D_L3Efficiency_GlbTM_ProbeTrackMatched_mc_MCTRUTH.dat"); 	
  //     //     ptM  = new PidTable("pidtables/110625/H2D_MuonIDEfficiency_GlbTM_ProbeTrackMatched_mc_MCTRUTH.dat"); 
  //   } else {
  //   }

  cout << "--> loopTree with mode " << mode << " proc = " << proc << " on file ";
  gFile->pwd();

  // -- reset all histograms
  for (unsigned int i = 0; i < fNchan; ++i) {
    fhMuId[i]->Reset();
    fhMuTr[i]->Reset();
    fhMuIdMC[i]->Reset();
    fhMuTrMC[i]->Reset();

    fh0PidTrigger[i]->Reset();
    fh1PidTrigger[i]->Reset();
    fh0PidMuID[i]->Reset();
    fh1PidMuID[i]->Reset();

    fh0PidMCTrigger[i]->Reset();
    fh1PidMCTrigger[i]->Reset();
    fh0PidMCMuID[i]->Reset();
    fh1PidMCMuID[i]->Reset();

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
    fhMassWithAllCutsBlind[i]->Reset();
    fhMassWithAllCutsManyBins[i]->Reset();

    fhNorm[i]->Reset();

    fhMassWithMassCutsManyBins[i]->Reset();
    fhMassWithMassCuts[i]->Reset();
  }

  // -- set up tree
  double mass(0.); 
  TTree *t;
  t = (TTree*)gFile->Get("events");
  int brun, bevt, bls, btm, bq1, bq2, bprocid; 
  double bg1pt, bg2pt, bg1eta, bg2eta;
  double bm, bcm, bpt, beta, bphi, bcosa, balpha, biso, bchi2, bdof, bdocatrk, bfls3d, bfl3dE, bfl3d;
  double bm1pt, bm1eta, bm2pt, bm2eta, bm1phi, bm2phi;
  double bk1pt, bk1eta, bk2pt, bk2eta; 
  double bg3pt, bg3eta, bg4pt, bg4eta; 
  double bmpsi, bmkk, bdr;
  double bw8mu, bw8tr;
  bool bhlt, bgmuid, bgtqual, bjson;
  double tr1w8(0.), tr2w8(0.), trw8(0.), m1w8(0.), m2w8(0.), mw8(0.0);

  double blip, blipE, btip, btipE; 
  int bm1pix, bm2pix, bm1bpix, bm2bpix, bm1bpixl1, bm2bpixl1;

  t->SetBranchAddress("lip",&blip);
  t->SetBranchAddress("lipE",&blipE);
  t->SetBranchAddress("tip",&btip);
  t->SetBranchAddress("tipE",&btipE);

  t->SetBranchAddress("m1pix",&bm1pix);
  t->SetBranchAddress("m2pix",&bm2pix);
  t->SetBranchAddress("m1bpix",&bm1bpix);
  t->SetBranchAddress("m2bpix",&bm2bpix);
  t->SetBranchAddress("m1bpixl1",&bm1bpixl1);
  t->SetBranchAddress("m2bpixl1",&bm2bpixl1);

  t->SetBranchAddress("run",&brun);
  t->SetBranchAddress("evt",&bevt);
  t->SetBranchAddress("hlt",&bhlt);
  t->SetBranchAddress("ls",&bls);
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
  }

  if (bs2jpsiphi) {
    if (isMC) {
      t->SetBranchAddress("g3pt",&bg3pt);
      t->SetBranchAddress("g3eta",&bg3eta);
      t->SetBranchAddress("g4pt",&bg4pt);
      t->SetBranchAddress("g4eta",&bg4eta);
    }
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
  bool cowboy(false); 
  double dphi; 
  TLorentzVector vm1, vm2, vpsi; 
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
    //     if (10 == mode || 11 == mode) {
    //       mass = bcm; 
    //     } else if (20 == mode || 21 == mode) {
    //       mass = bcm; 
    //     }
    

    fhMassAbsNoCuts[fChan]->Fill(mass);
    // -- require wide mass window
    if (mass < fMassLo) continue;
    if (fMassHi < mass) continue;

    // -- gen-level acceptance cuts
    if (isMC) {
      if (TMath::Abs(bg1eta) > 2.5) continue;
      if (TMath::Abs(bg2eta) > 2.5) continue;
      if (bg1pt < 1.0) continue;
      if (bg2pt < 1.0) continue;
      if (bp2jpsikp) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (bg3pt < 0.4) continue;
      }
      
      if (bs2jpsiphi) {
	if (TMath::Abs(bg3eta) > 2.5) continue;
	if (TMath::Abs(bg4eta) > 2.5) continue;
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
    if (false == bgtqual) continue;
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
    
    // -- analysis cuts
    if (bq1*bq2 > 0) continue;
    if (bm1pt < pCuts->m1pt) continue; 
    if (bm2pt < pCuts->m2pt) continue; 

    if (bfl3d > 2) continue;

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

    cowboy = false; 
    vm1.SetPtEtaPhiM(bm1pt, bm1eta, bm1phi, MMUON);
    vm2.SetPtEtaPhiM(bm2pt, bm2eta, bm2phi, MMUON);
    dphi = vm1.DeltaPhi(vm2); 
    cowboy = bq1*dphi > 0; 
    if (bs2jpsiphi || bp2jpsikp) {
      vpsi = vm1 + vm2; 
      if (bmpsi > 3.2) continue;
      if (bmpsi < 3.0) continue;
      // -- cowboy veto 
      if (cowboy) {
 	continue;
      }
      if (vpsi.Perp() < 7) {
	continue;
      }
    } 

    fhMassWithAnaCuts[fChan]->Fill(mass); 
    fhMassWithAnaCutsManyBins[fChan]->Fill(mass); 

    // -- MUON ID
    if (false == bgmuid) continue;

    // -- Data PidTables
    m1w8 = ptM->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    m2w8 = ptM->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    mw8  = m1w8*m2w8; 

    if (mw8 > 0.) {
      fhMuId[fChan]->Fill(mw8, 1./mw8); 
    } else {
      cout << "mw8 = " << mw8 << " for muons = " << bm1pt << " " << bm1eta << " " << bm2pt << " " << bm2eta << endl;
    }

    fh0PidMuID[fChan]->Fill(bpt, mw8); 
    fh1PidMuID[fChan]->Fill(bpt); 

    fhMassWithMuonCuts[fChan]->Fill(mass); 
    fhMassWithMuonCutsManyBins[fChan]->Fill(mass); 

    // -- MC for comparison 
    m1w8 = ptMMC->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    m2w8 = ptMMC->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    mw8  = tr1w8*tr2w8; 
    if (mw8 > 0.) {
      fhMuIdMC[fChan]->Fill(mw8, 1./mw8); 
    }
    fh0PidMCMuID[fChan]->Fill(bpt, mw8); 
    fh1PidMCMuID[fChan]->Fill(bpt); 

    // -- TRIGGER
    fh1MCTrigger[fChan]->Fill(bpt); 
    if (false == bhlt) continue;
    fhMassWithTriggerCuts[fChan]->Fill(mass); 
    fhMassWithTriggerCutsManyBins[fChan]->Fill(mass); 

    // -- Data PidTables
    tr1w8 = ptT1->effD(bm1pt, TMath::Abs(bm1eta), 0.)*ptT2->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    tr2w8 = ptT1->effD(bm2pt, TMath::Abs(bm2eta), 0.)*ptT2->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (trw8 > 0.) {
      fhMuTr[fChan]->Fill(trw8, 1./trw8); 
    }

    fh0PidTrigger[fChan]->Fill(bpt, trw8); 
    fh1PidTrigger[fChan]->Fill(bpt); 

    // -- MC for comparison 
    tr1w8 = ptT1MC->effD(bm1pt, TMath::Abs(bm1eta), 0.)*ptT2->effD(bm1pt, TMath::Abs(bm1eta), 0.);
    tr2w8 = ptT1MC->effD(bm2pt, TMath::Abs(bm2eta), 0.)*ptT2->effD(bm2pt, TMath::Abs(bm2eta), 0.);
    trw8  = tr1w8*tr2w8; 
    if (trw8 > 0.) {
      fhMuTrMC[fChan]->Fill(trw8, 1./trw8); 
    }
    fh0PidMCTrigger[fChan]->Fill(bpt, trw8); 
    fh1PidMCTrigger[fChan]->Fill(bpt); 

    fh0MCTrigger[fChan]->Fill(bpt); 

    fhMassWithAllCuts[fChan]->Fill(mass); 
    if (5 == mode && !(5.2 < mass && mass < 5.45)) {
      fhMassWithAllCutsBlind[fChan]->Fill(mass); 
    }

    fhMassWithAllCutsManyBins[fChan]->Fill(mass); 

    fhNorm[fChan]->Fill(mass);
    
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
      cout << Form("m = %4.3f pT = %4.3f eta = %4.3f", mass, bpt, beta)
	   <<	" r = " << brun << "/" << bevt
	   << " chan = " << fChan 
	   << Form(" mpt = %4.3f,%4.3f", bm1pt, bm2pt)
	   << Form(" meta = %4.3f,%4.3f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %4.3f chi2 = %4.3f fls3d = %4.3f, fl/E=%4.3f/%4.3f", 
		   TMath::ACos(bcosa), biso, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
	   << endl;
      fOUT << Form("m = %4.3f pT = %4.3f eta = %4.3f", mass, bpt, beta)
	//	   <<	" run = " << brun << " event = " << bevt
	   << " chan = " << fChan 
	   << Form(" mpt = %4.3f,%4.3f", bm1pt, bm2pt)
	   << Form(" meta = %4.3f,%4.3f", TMath::Abs(bm1eta), TMath::Abs(bm2eta))
	   << Form(" a = %4.3f iso = %4.3f chi2 = %4.3f fls3d = %4.3f, fl/E=%4.3f/%4.3f", 
		   TMath::ACos(bcosa), biso, bchi2/bdof, bfls3d, bfl3d, bfl3dE)
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

      fTEX << formatTex(brun,      Form("%s:%s%i:run", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bevt,      Form("%s:%s%i:evt", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(fChan,     Form("%s:%s%i:chan", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm,        Form("%s:%s%i:m", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bpt,       Form("%s:%s%i:pt", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bphi,      Form("%s:%s%i:phi", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(beta,      Form("%s:%s%i:eta", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << Form("\\vdef{%s:%s%i:channel}   {%s }", fSuffix.c_str(), st.c_str(), ievt, fChan==0?"barrel":"endcap") << endl;
      fTEX << formatTex((cowboy?1:0),    Form("%s:%s%i:cowboy", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm1pt,     Form("%s:%s%i:m1pt", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm2pt,     Form("%s:%s%i:m2pt", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm1eta,    Form("%s:%s%i:m1eta", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm2eta,    Form("%s:%s%i:m2eta", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm1phi,    Form("%s:%s%i:m1phi", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bm2phi,    Form("%s:%s%i:m2phi", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(bq1,       Form("%s:%s%i:m1q", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bq2,       Form("%s:%s%i:m2q", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(biso,      Form("%s:%s%i:iso", fSuffix.c_str(), st.c_str(), ievt), 3) << endl;
      fTEX << formatTex(balpha,    Form("%s:%s%i:alpha", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bchi2,     Form("%s:%s%i:chi2", fSuffix.c_str(), st.c_str(), ievt), 2) << endl;
      fTEX << formatTex(bdof,      Form("%s:%s%i:dof", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bfls3d,    Form("%s:%s%i:fls3d", fSuffix.c_str(), st.c_str(), ievt), 2) << endl;
      fTEX << formatTex(bfl3d,     Form("%s:%s%i:fl3d", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(bfl3dE,    Form("%s:%s%i:fl3dE", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;

      fTEX << formatTex(bdocatrk,  Form("%s:%s%i:docatrk", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(blip,      Form("%s:%s%i:lip", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(blipE,     Form("%s:%s%i:lipE", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(btip,      Form("%s:%s%i:tip", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;
      fTEX << formatTex(btipE,     Form("%s:%s%i:tipE", fSuffix.c_str(), st.c_str(), ievt), 4) << endl;

      fTEX << formatTex(bm1pix,    Form("%s:%s%i:m1pix", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm2pix,    Form("%s:%s%i:m2pix", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm1bpix,   Form("%s:%s%i:m1bpix", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm2bpix,   Form("%s:%s%i:m2bpix", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm1bpixl1, Form("%s:%s%i:m1bpixl1", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      fTEX << formatTex(bm2bpixl1, Form("%s:%s%i:m2bpixl1", fSuffix.c_str(), st.c_str(), ievt), 0) << endl;
      
    }
    
  }

  cout << gFile->GetName() << ": " << fhMassWithAllCuts[0]->GetSumOfWeights() 
       << "  " << fhMassWithAllCuts[1]->GetSumOfWeights()
       << endl;
  if (98 == mode) {
    if (fhMassWithAllCuts[0]->GetSumOfWeights() > 0) {
      fhMassNoCuts[0]->Scale(fhMassWithAllCuts[0]->GetSumOfWeights()/fhMassNoCuts[0]->GetSumOfWeights());
    } else {
      fhMassNoCuts[0]->Scale(2.3/fhMassNoCuts[0]->GetSumOfWeights());
    }
    if (fhMassWithAllCuts[1]->GetSumOfWeights() > 0) {
      fhMassNoCuts[1]->Scale(fhMassWithAllCuts[1]->GetSumOfWeights()/fhMassNoCuts[1]->GetSumOfWeights());
    } else {
      fhMassNoCuts[1]->Scale(2.3/fhMassNoCuts[1]->GetSumOfWeights());
    }

    if (fhMassWithAllCutsManyBins[0]->GetSumOfWeights() > 0) {
      fhMassNoCutsManyBins[0]->Scale(fhMassWithAllCutsManyBins[0]->GetSumOfWeights()/fhMassNoCutsManyBins[0]->GetSumOfWeights());
    } else {
      fhMassNoCutsManyBins[0]->Scale(2.3/fhMassNoCutsManyBins[0]->GetSumOfWeights());
    }
    if (fhMassWithAllCutsManyBins[1]->GetSumOfWeights() > 0) {
      fhMassNoCutsManyBins[1]->Scale(fhMassWithAllCutsManyBins[1]->GetSumOfWeights()/fhMassNoCutsManyBins[1]->GetSumOfWeights());
    } else {
      fhMassNoCutsManyBins[1]->Scale(2.3/fhMassNoCutsManyBins[1]->GetSumOfWeights());
    }

//     c0->Clear(); 
//     c0->Divide(1,2);
//     c0->cd(1);
//     fhMassNoCuts[0]->Draw();
//     c0->cd(2);
//     fhMassNoCuts[1]->Draw();
    
    return;
  }

  if (99 == mode) return;

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

    if (10 == mode || 11 == mode) {
      aa = fNumbersNo[i];
    }

    if (20 == mode || 21 == mode) {
      aa = fNumbersCs[i];
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
    if (fDoPrint)    c0->SaveAs(Form("%s/muid-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMuTr[i]->Draw();
    if (fDoPrint)   c0->SaveAs(Form("%s/mutr-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassAbsNoCuts[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/anc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassNoCuts[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/noc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithAllCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/wac-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithMuonCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/muc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithTriggerCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/trc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    fhMassWithMassCutsManyBins[i]->Draw();
    if (fDoPrint)    c0->SaveAs(Form("%s/wmc-mode-%d-chan%d.pdf", fDirectory.c_str(), mode, i));

    // -- Efficiency and acceptance
    if (isMC) {
      accEffFromEffTree(fAcc, *aa, *fCuts[i], proc);
      //      accEffFromEffTree(gFile, *aa, *fCuts[i], 0, 1, proc);
      //      double a = fhMassAbsNoCuts[i]->GetSumOfWeights(); 
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
      aa->effAna           = b/a;
      aa->effAnaE          = dEff(static_cast<int>(b), static_cast<int>(a));
      aa->effMuidMC        = c/b;
      aa->effMuidMCE       = dEff(static_cast<int>(c), static_cast<int>(b));
      aa->effMuidPid       = fhMuId[i]->GetMean();
      aa->effMuidPidE      = fhMuId[i]->GetMeanError();
      aa->effMuidPidMC     = fhMuIdMC[i]->GetMean();
      aa->effMuidPidMCE    = fhMuIdMC[i]->GetMeanError();
      aa->effTrigMC        = d/c;
      aa->effTrigMCE       = dEff(static_cast<int>(d), static_cast<int>(c));
      aa->effTrigPid       = fhMuTr[i]->GetMean();
      aa->effTrigPidE      = fhMuTr[i]->GetMeanError();
      aa->effTrigPidMC     = fhMuTrMC[i]->GetMean();
      aa->effTrigPidMCE    = fhMuTrMC[i]->GetMeanError();
      aa->effTot           = e/(aa->genYield);
      aa->effTotE          = dEff(static_cast<int>(e), static_cast<int>(aa->genYield));
      aa->effTotChan       = e/(aa->genChanYield);
      aa->effTotChanE      = dEff(static_cast<int>(e), static_cast<int>(aa->genChanYield));
      aa->effProdMC        = aa->effCand * aa->effAna * aa->effMuidMC * aa->effTrigMC;
      aa->effProdMCE       = 0.;
      aa->effProdPid       = aa->effCand * aa->effAna * aa->effMuidPid * aa->effTrigPid;
      aa->effProdPidE      = 0.;
      aa->aEffProdMC       = aa->effProdMC * aa->accChan * aa->cFrac;
      aa->aEffProdMCE      = 0.;

      aa->combGenYield     = e/(aa->acc * aa->effProdMC);
      aa->chanGenYield     = e/(aa->accChan * aa->effProdMC); 
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
      if (fDoPrint)      c0->SaveAs(Form("%s/sig-mc-chan%d.pdf", fDirectory.c_str(), i));
      cout << "----> "; gFile->pwd(); 
    } else if (5 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA SIGNAL, channel " << i  << endl;
      bgBlind(fhMassWithAllCutsBlind[i], 1, fBgLo, fBgHi);
      cout << "fBgHist = " << fBgHist << "+/-" << fBgHistE << endl;
      aa->bgObs = fBgHist;
      double blind = 5.45 - 5.20; 
      double scaleBs = (aa->mBsHi-aa->mBsLo)/(fBgHi-fBgLo-blind);
      double scaleBd = (aa->mBdHi-aa->mBdLo)/(fBgHi-fBgLo-blind);
      aa->tauBs    = scaleBs; 
      aa->tauBsE   = 0.04*scaleBs; 
      aa->tauBd    = scaleBd; 
      aa->tauBdE   = 0.04*scaleBd; 
      cout << "CCCCCCCCCCCC  " << scaleBs << " " << aa->tauBs << "+/-" << aa->tauBsE << endl;
      aa->bgBsExp  = scaleBs*aa->bgObs;
      aa->bgBsExpE = scaleBs*TMath::Sqrt(aa->bgObs);
      aa->bgBdExp  = scaleBd*aa->bgObs;
      aa->bgBdExpE = scaleBd*TMath::Sqrt(aa->bgObs);

      double cnt = fhMassWithAllCuts[i]->Integral(fhMassWithAllCuts[i]->FindBin(aa->mBsLo), 
						  fhMassWithAllCuts[i]->FindBin(aa->mBsHi));
      aa->bsObs = cnt;

      cnt = fhMassWithAllCuts[i]->Integral(fhMassWithAllCuts[i]->FindBin(aa->mBdLo), 
					   fhMassWithAllCuts[i]->FindBin(aa->mBdHi));
      aa->bdObs = cnt;

      // -- blinded version
      TH1D *h = fhMassWithAllCutsBlind[i]; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#mu} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.3); 
      h->SetMinimum(0.01); 
      gStyle->SetOptStat(0); 
      gStyle->SetOptTitle(0); 
      h->SetAxisRange(4.9, 5.9, "X"); 
      h->SetMaximum(2.2);
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
      
      stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
      if (fDoPrint)  c0->SaveAs(Form("%s/sig-data-chan%d.pdf", fDirectory.c_str(), i));
      // -- unblinded version
      h = fhMassWithAllCuts[i]; 
      setHist(h, kBlack, 20, 1.); 
      setTitles(h, "m_{#mu#mu} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.3); 
      h->SetMinimum(0.01); 
      h->SetNdivisions(003, "Y");
      h->SetAxisRange(4.9, 5.9, "X"); 
      h->SetMaximum(2.2);
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


      stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
      if (fDoPrint)  c0->SaveAs(Form("%s/unblinded-sig-data-chan%d.pdf", fDirectory.c_str(), i));
    } else if (10 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC NORMALIZATION, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) c0->SaveAs(Form("%s/norm-mc-chan%d.pdf", fDirectory.c_str(), i));
    } else if (11 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA NORMALIZATION, channel " << i  << endl;
      //      TH1D *h = fhMassWithAllCutsManyBins[i]; 
      //      TH1D *h = fhMassWithAllCuts[i]; 
      //      normYield(h, mode, 5.10, 5.5);
      TH1D *h = fhNorm[i];
      normYield(h, i, 5.0, 5.6);
      aa->fitYield  = fNoSig; 
      aa->fitYieldE = fNoSigE; 
    } else if (20 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: MC CONTROL SAMPLE, channel " << i  << endl;
      fhMassNoCuts[i]->Draw();
      fhMassWithAnaCuts[i]->Draw("same");
      tl->DrawLatex(0.2, 0.77, Form("w/o cuts: %5.1f", fhMassNoCuts[i]->GetSumOfWeights())); 
      tl->DrawLatex(0.2, 0.7, Form("w/   cuts: %5.1f", fhMassWithAnaCuts[i]->GetSumOfWeights())); 
      if (fDoPrint) c0->SaveAs(Form("%s/cs-mc-chan%d.pdf", fDirectory.c_str(), i));
    } else if (21 == mode) {
      cout << "----------------------------------------------------------------------" << endl;
      cout << "==> loopTree: DATA CONTROL SAMPLE, channel " << i  << endl;
      //      TH1D *h = fhMassWithAllCuts[i]; 
      //      csYield(h, mode, 5.20, 5.6);
      TH1D *h = fhNorm[i];
      csYield(h, i, 5.0, 5.6);
      aa->fitYield  = fCsSig; 
      aa->fitYieldE = fCsSigE; 
    } 

    printNumbers(*aa, cout); 
    printNumbers(*aa, fOUT); 
    
    // -- Cache the pwd...
    TDirectory *pD = gFile; 
    fHistFile->cd();
    cout << "==> " << i << "  " << mode << " hMassWithMassCuts[i] = " << fhMassWithMassCuts[i] << endl;
    fhMassWithMassCutsManyBins[i]->SetName(Form("hMassWithMassCutsManyBins%d_chan%d", mode, i)); fhMassWithMassCutsManyBins[i]->Write();
    fhMassWithAllCutsManyBins[i]->SetName(Form("hMassWithAllCutsManyBins%d_chan%d", mode, i)); fhMassWithAllCutsManyBins[i]->Write();
    fhMassWithMassCuts[i]->SetName(Form("hMassWithMassCuts%d_chan%d", mode, i)); fhMassWithMassCuts[i]->Write();
    fhMassWithAllCuts[i]->SetName(Form("hMassWithAllCuts%d_chan%d", mode, i)); fhMassWithAllCuts[i]->Write();
    fhNorm[i]->SetName(Form("hNorm%d_chan%d", mode, i)); fhNorm[i]->Write();
    fhMassNoCuts[i]->SetName(Form("hMassNoCuts%d_chan%d", mode, i)); fhMassNoCuts[i]->Write();
    fhMassAbsNoCuts[i]->SetName(Form("hMassAbsNoCuts%d_chan%d", mode, i)); fhMassAbsNoCuts[i]->Write();
    fhMuTr[i]->SetName(Form("hMuTr%d_chan%d", mode, i)); fhMuTr[i]->Write();
    fhMuId[i]->SetName(Form("hMuId%d_chan%d", mode, i)); fhMuId[i]->Write();
    // -- and get back to it
    pD->cd();
  }

  delete ptT1;
  delete ptT2;
  delete ptM;

  return;
}


// ----------------------------------------------------------------------
void plotReducedTree::initNumbers(numbers *a) {

  a->name = "";
  a->effGenFilter = a->effGenFilterE = 1.;
  a->fitYield = a->fitYieldE = 0.;
  a->genFileYield = a->genYield = a->recoYield = a->muidYield = a->trigYield = a->candYield = a->ana0Yield = a->anaYield = a->anaWmcYield = 0; 
  a->acc = a->accE = 0; 
  a->effMuidMC =  a->effMuidMCE = a->effTrigMC = a->effTrigMCE = 0; 
  a->effMuidPid = a->effMuidPidE = a->effTrigPid = a->effTrigPidE = 0; 
  a->effCand = a->effCandE = 0; 
  a->effAna = a->effAnaE = 0; 
  // -- this is only relevant for the signal(s)
  a->pss   = a->pdd = 1.;
  a->psd   = a->pds = 1.;
  a->bgObs = a->bgBsExp = a->bgBsExpE = a->bgBdExp = a->bgBdExpE =0; 
  a->bsObs = a->bdObs = 0; 
  a->mBdLo = a->mBdHi = a->mBsLo = a->mBsHi = 0.;
}


// ----------------------------------------------------------------------
int plotReducedTree::detChan(double m1eta, double m2eta) {
  // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
  if (TMath::Abs(m1eta) < fCuts[0]->etaMax && TMath::Abs(m2eta) < fCuts[0]->etaMax) return 0; 
  if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1; 
  return -1; 
}


// ----------------------------------------------------------------------
void plotReducedTree::accEffFromEffTree(string fname, numbers &a, cuts &b, int proc) {

  TFile *f = fF[fname];
  if (0 == f) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no file " << fname << " found " << endl;
    return;
  }
  TTree *t  = (TTree*)(f->Get("effTree"));
  double effFilter(1.); 
  if (!t) {
    cout << "anaBmm::accEffFromEffTree(" << a.name << "): no tree `effTree' found " << endl;
    f->pwd(); 
    return;
  } else {
    effFilter = fFilterEff[fname]; 
    cout << "anaBmm::accEffFromEffTree(" << a.name << ")" << endl
	 << " get acceptance from file " << f->GetName() 
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
	      chan = detChan(bm1eta, bm2eta); 
	      if (chan == a.index) {
		++nreco;
		if (bm1pt > b.m1pt && bm2pt > b.m2pt) {
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
    }
  }
  a.effGenFilter  = effFilter; 
  a.genFileYield  = ngen;
  a.genYield      = a.genFileYield/effFilter;
  a.genChanYield  = nchangen; 
  a.recoYield     = nreco; // reco'ed in chan, basic global reconstruction cuts 
  a.chanYield     = nchan; // reco'ed in chan, with channel-dependent (pT) cuts
  a.muidYield     = nmuid;
  a.trigYield     = nhlt;
  a.candYield     = ncand;
  
  if (a.genYield > 0) {
    a.cFrac  = a.genChanYield/a.genYield;
    a.cFracE = dEff(static_cast<int>(a.genChanYield), static_cast<int>(a.genYield));
  }
  if (a.genYield > 0) {
    a.acc = a.recoYield/a.genYield;
    a.accE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genYield));
  }
  
  if (a.genChanYield > 0) {
    a.accChan = a.recoYield/a.genChanYield;
    a.accChanE = dEff(static_cast<int>(a.recoYield), static_cast<int>(a.genChanYield));
  }
  
  if (a.recoYield > 0) {
    a.effChan = a.chanYield/a.recoYield;
    a.effChanE = dEff(static_cast<int>(a.chanYield), static_cast<int>(a.recoYield));
  }


  if (a.trigYield > 0) {
    a.effCand = a.candYield/a.trigYield;
    a.effCandE = dEff(static_cast<int>(a.candYield), static_cast<int>(a.trigYield));
  } 

  //   if (0) {
  //     if (a.chanYield > 0) {
  //       a.accMuidMC = a.muidYield/a.chanYield;
  //       a.accMuidMCE = dEff(static_cast<int>(a.muidYield), static_cast<int>(a.chanYield));
  //     } 
  //     if (a.muidYield > 0) {
  //       a.accTrigMC = a.trigYield/a.muidYield;
  //       a.accTrigMCE = dEff(static_cast<int>(a.trigYield), static_cast<int>(a.muidYield));
  //     } 
  //   }
}


// ----------------------------------------------------------------------
void plotReducedTree::bgBlind(TH1 *h, int mode, double lo, double hi) {
  
  if (0 == h) { 
    cout << "plotReducedTree::bgBlind(...): No histogram passed! mode = " << mode << endl;
    return;
  }
  
  TF1 *lF1(0);

  double histCount = h->Integral(h->FindBin(fBgLo), h->FindBin(fBgHi)-1); 
  cout << "bgBlind: histCount = " << histCount << " starting at " << h->FindBin(fBgLo) << " to " << h->FindBin(fBgHi)-1 << endl;
  fBgHist  = histCount; 
  fBgHistE  = TMath::Sqrt(histCount); 
  fBgHistExp  = histCount*(fSgHi-fSgLo)/(fBgHi-fBgLo-0.25); // FIXME fixed limits
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
void plotReducedTree::normYield(TH1 *h, int mode, double lo, double hi) {

  TF1 *lF1(0);
  
  //  lF1 = fpFunc->pol1Gauss(h); 
  if (0 == mode) { 
    fpFunc->fLo = 5.0;
    fpFunc->fHi = 5.5;
    lF1 = fpFunc->pol1ErrGauss(h, 5.27, 0.034, 5.1); 
    //    lBg = new TF1("fbg", fa_pol1_err, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 6);
  } else {
    fpFunc->fLo = 5.0;
    fpFunc->fHi = 5.5;
    lF1 = fpFunc->expoErrGauss(h, 5.27, 0.034, 5.1); 
    //    lBg = new TF1("fbg", fa_expo_err, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()), 6);
  }
  h->Fit(lF1, "rm", "", lo, hi); 

  double c  = lF1->GetParameter(0); 
  double cE = lF1->GetParError(0); 
  double ierr = lF1->IntegralError(5.15, 5.4)/h->GetBinWidth(1); 

  fNoSig = c/h->GetBinWidth(1);
  if (ierr > TMath::Sqrt(fNoSig)) {
    fNoSigE = ierr;
  } else {
    fNoSigE = cE/c*fNoSig;
  }

  cout << "N(Sig) = " << fNoSig << " +/- " << fNoSigE << endl;


  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muK} [GeV]", Form("Candidates / %3.3f GeV", h->GetBinWidth(1)), 0.05, 1.2, 1.6); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

//   for (int i = 0; i < lBg->GetNpar(); ++i) {
//     lBg->SetParameter(i, lF1->GetParameter(3+i));
//     cout << "par " << i << ": " << lBg->GetParName(i) << " = " << lBg->GetParameter(i) << endl;
//   }

//   lBg->SetLineStyle(kDashed);
//   lBg->SetLineColor(kRed);
//   lBg->SetLineWidth(3);
//   lBg->Draw("same");

  tl->SetTextSize(0.07); 
  if (0 == mode) {
    tl->DrawLatex(0.6, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.6, 0.8, "Endcap");   
  } 

  stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint) c0->SaveAs(Form("%s/norm-data-chan%d.pdf", fDirectory.c_str(), mode));
  

  //   double c  = h->GetFunction("f1")->GetParameter(0); 
  //   double cE = h->GetFunction("f1")->GetParError(0); 
  
  delete lF1; 

}


// ----------------------------------------------------------------------
void plotReducedTree::csYield(TH1 *h, int mode, double lo, double hi) {

  TF1 *lF1(0);
  
  if (0 == mode) { 
    fpFunc->fLo = 5.1;
    fpFunc->fHi = 5.6;
    lF1 = fpFunc->expoGauss(h, 5.365, 0.033); 
  } else {
    fpFunc->fLo = 5.1;
    fpFunc->fHi = 5.6;
    lF1 = fpFunc->expoGauss(h, 5.365, 0.033); 
  }
  h->Fit(lF1, "rm", "", lo, hi); 


  //   double c  = h->GetFunction("f1")->GetParameter(0); 
  //   double cE = h->GetFunction("f1")->GetParError(0); 

  double c  = lF1->GetParameter(0); 
  double cE = lF1->GetParError(0); 

  double ierr = lF1->IntegralError(5.15, 5.4)/h->GetBinWidth(1); 

  setHist(h, kBlack, 20, 1.); 
  setTitles(h, "m_{#mu#muKK} [GeV]", "Candidates/Bin", 0.05, 1.2, 1.6); 
  h->SetMinimum(0.01); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  h->Draw("e");

  fCsSig = c/h->GetBinWidth(1);
  if (ierr > TMath::Sqrt(fCsSig)) {
    fCsSigE = ierr;
  } else {
    fCsSigE = cE/c*fCsSig;
  }

  tl->SetTextSize(0.07); 
  if (0 == mode) {
    tl->DrawLatex(0.22, 0.8, "Barrel");   
  } 

  if (1 == mode) {
    tl->DrawLatex(0.22, 0.8, "Endcap");   
  } 

  stamp(0.18, "CMS, 1.14 fb^{-1}", 0.67, "#sqrt{s} = 7 TeV"); 
  if (fDoPrint) c0->SaveAs(Form("%s/cs-data-chan%d.pdf", fDirectory.c_str(), mode));

  cout << "N(Sig) = " << fCsSig << " +/- " << fCsSigE << endl;
  
  delete lF1; 

}


// ----------------------------------------------------------------------
void plotReducedTree::printNumbers(numbers &a, ostream &OUT) {
  OUT << "======================================================================" << endl;
  OUT << "numbers for \"" << a.name.c_str() << "\"" << endl;
  OUT << "fitYield     = " << a.fitYield << "+/-" << a.fitYieldE << endl;
  OUT << "genFileYield = " << a.genFileYield << endl;
  OUT << "genYield     = " << a.genYield << endl;
  OUT << "genChanYield = " << a.genChanYield << endl;
  OUT << "recoYield    = " << a.recoYield << endl;
  OUT << "chanYield    = " << a.chanYield << endl;
  OUT << "muidYield    = " << a.muidYield << endl;
  OUT << "trigYield    = " << a.trigYield << endl;
  OUT << "candYield    = " << a.candYield << endl;
  OUT << "ana0Yield    = " << a.ana0Yield << endl;
  OUT << "anaYield     = " << a.anaYield << endl;
  OUT << "anaMuYield   = " << a.anaMuonYield << endl;
  OUT << "anaTrigYield = " << a.anaTriggerYield << endl;
  OUT << "anaWmcYield  = " << a.anaWmcYield << endl;
  OUT << "mBsLo        = " << a.mBsLo << endl;
  OUT << "mBsHi        = " << a.mBsHi << endl;
  OUT << "mBdLo        = " << a.mBdLo << endl;
  OUT << "mBdHi        = " << a.mBdHi << endl;
  OUT << "PSS          = " << a.pss << endl;
  OUT << "PDS          = " << a.pds << endl;
  OUT << "PSD          = " << a.psd << endl;
  OUT << "PDD          = " << a.pdd << endl;
  OUT << "bsRare       = " << a.bsRare << endl;
  OUT << "bdRare       = " << a.bdRare << endl;
  OUT << "gen filter   = " << a.effGenFilter << endl;
  OUT << "acceptance   = " << a.acc << "+/-" << a.accE << endl;
  OUT << "accChan      = " << a.accChan << "+/-" << a.accChanE << endl; 
  OUT << "cFrac        = " << a.cFrac << "+/-" << a.cFracE << endl; 
  OUT << "effChan      = " << a.effChan << "+/-" << a.effChanE << endl; 
  OUT << "effCand      = " << a.effCand << "+/-" << a.effCandE << endl;
  OUT << "effAna       = " << a.effAna << "+/-" << a.effAnaE << endl; 
  OUT << "effMuidMC    = " << a.effMuidMC << "+/-" << a.effMuidMCE << endl;
  OUT << "effMuidPid   = " << a.effMuidPid << "+/-" << a.effMuidPidE << endl;
  OUT << "effMuidPidMC = " << a.effMuidPidMC << "+/-" << a.effMuidPidMCE << endl;
  OUT << "accMuidMC    = " << a.accMuidMC << "+/-" << a.accMuidMCE << endl;
  OUT << "effTrigMC    = " << a.effTrigMC << "+/-" << a.effTrigMCE << endl;
  OUT << "effTrigPid   = " << a.effTrigPid << "+/-" << a.effTrigPidE << endl;
  OUT << "effTrigPidMC = " << a.effTrigPidMC << "+/-" << a.effTrigPidMCE << endl;
  OUT << "accTrigMC    = " << a.accTrigMC << "+/-" << a.accTrigMCE << endl;
  OUT << "effProd(MC)  = " << a.effProdMC << endl;
  OUT << "effProd(Pid) = " << a.effProdPid << endl;
  OUT << "effProd(MC)A = " << a.effProdMC*a.acc << endl;
  OUT << "effTot       = " << a.effTot << "+/-" << a.effTotE << endl; 
  OUT << "effTotChan   = " << a.effTotChan << "+/-" << a.effTotChanE << endl; 
  OUT << "combGenYield     = " << a.combGenYield << endl; 
  OUT << "prodGenYield     = " << a.prodGenYield << endl; 
  OUT << "prodChanGenYield = " << a.chanGenYield << endl; 
  OUT.flush();

}
