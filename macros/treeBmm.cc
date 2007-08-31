#include "treeBmm.hh"

#include "TRandom.h"

#include "treeBmm.icc"


// Run with: ./bmm -c chains/bg-test -D root ; ./bmm -c chains/sg-test -D root ;

// ----------------------------------------------------------------------
void treeBmm::startAnalysis() {
  cout << "startAnalysis: ..." << endl;
}



// ----------------------------------------------------------------------
void treeBmm::readCuts(TString filename, int dump, double ptMin, double etaMax) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());
  fn.ReplaceAll("bjk", "");
  fn.ReplaceAll("bmm", "");
  fn.ReplaceAll("cuts", "");
  fn.ReplaceAll(".", "");
  fn.ReplaceAll("/", "");

  fChannel = fn;

  if (dump) {
    cout << "Setting thresholds to pT = " << ptMin <<  " and eta = " << etaMax << endl;
    cout << "====================================" << endl;
    cout << "Cut file  " << fCutFile.Data() << endl;
    cout << "Cut label " << fn.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  fMinPt = ptMin;
  fMaxEta = etaMax;

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin; 

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }

    if (!strcmp(CutName, "PTLO")) {
      PTLO = CutValue; ok = 1;
      if (dump) cout << "PTLO:           " << PTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, PTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(l) [GeV]");
    }

    if (!strcmp(CutName, "PTHI")) {
      PTHI = CutValue; ok = 1;
      if (dump) cout << "PTHI:           " << PTHI << " GeV" << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, PTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(l) [GeV]");
    }

    if (!strcmp(CutName, "ETALO")) {
      ETALO = CutValue; ok = 1;
      if (dump) cout << "ETALO:           " << ETALO << endl;
      ibin = 13;
      hcuts->SetBinContent(ibin, ETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta_{T}^{min}(l)");
    }

    if (!strcmp(CutName, "ETAHI")) {
      ETAHI = CutValue; ok = 1;
      if (dump) cout << "ETAHI:           " << ETAHI << endl;
      ibin = 14;
      hcuts->SetBinContent(ibin, ETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#eta_{T}^{max}(l)");
    }

    if (!strcmp(CutName, "RMMLO")) {
      RMMLO = CutValue; ok = 1;
      if (dump) cout << "RMMLO:           " << RMMLO << endl;
      ibin = 15;
      hcuts->SetBinContent(ibin, RMMLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "R_{#mu#mu}^{min}");
    }

    if (!strcmp(CutName, "RMMHI")) {
      RMMHI = CutValue; ok = 1;
      if (dump) cout << "RMMHI:           " << RMMHI << endl;
      ibin = 16;
      hcuts->SetBinContent(ibin, RMMHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "R_{#mu#mu}^{max}");
    }

    if (!strcmp(CutName, "TIPHI")) {
      TIPHI = CutValue; ok = 1;
      if (dump) cout << "TIPHI:           " << TIPHI << endl;
      ibin = 17;
      hcuts->SetBinContent(ibin, TIPHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "TIP(l) [cm]");
    }

    if (!strcmp(CutName, "PTBS")) {
      PTBS = CutValue; ok = 1;
      if (dump) cout << "PTBS:           " << PTBS << " GeV" << endl;
      ibin = 100; 
      hcuts->SetBinContent(ibin, PTBS);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}(B_{s}) [GeV]");
    }

    if (!strcmp(CutName, "VTXCHI")) {
      VTXCHI = CutValue; ok = 1;
      if (dump) cout << "VTXCHI:           " << VTXCHI << endl;
      ibin = 102;
      hcuts->SetBinContent(ibin, VTXCHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "#chi^2");
    }

    if (!strcmp(CutName, "L3DLO")) {
      L3DLO = CutValue; ok = 1;
      if (dump) cout << "L3DLO:           " << L3DLO << endl;
      ibin = 103;
      hcuts->SetBinContent(ibin, L3DLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{3d} [cm]");
    }

    if (!strcmp(CutName, "LXYSXYLO")) {
      LXYSXYLO = CutValue; ok = 1;
      if (dump) cout << "LXYSXYLO:           " << LXYSXYLO << endl;
      ibin = 104;
      hcuts->SetBinContent(ibin, LXYSXYLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{xy}/#sigma_{xy}");
    }

    if (!strcmp(CutName, "LXYLO")) {
      LXYLO = CutValue; ok = 1;
      if (dump) cout << "LXYLO:           " << LXYLO << endl;
      ibin = 105;
      hcuts->SetBinContent(ibin, LXYLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "l_{xy} [cm]");
    }

    if (!strcmp(CutName, "COSALPHA")) {
      COSALPHA = CutValue; ok = 1;
      if (dump) cout << "COSALPHA:        " << COSALPHA << endl;
      ibin = 106;
      hcuts->SetBinContent(ibin, COSALPHA);
      hcuts->GetXaxis()->SetBinLabel(ibin, "cos(#alpha)");
    }

    if (!strcmp(CutName, "ISOVETO")) {
      ISOVETO = int(CutValue); ok = 1;
      if (dump) cout << "ISOVETO:           " << ISOVETO << endl;
      ibin = 107;
      hcuts->SetBinContent(ibin, ISOVETO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "I_{veto}");
    }

    if (!strcmp(CutName, "ISOLATION")) {
      ISOLATION = CutValue; ok = 1;
      if (dump) cout << "ISOLATION:           " << ISOLATION << endl;
      ibin = 108;
      hcuts->SetBinContent(ibin, ISOLATION);
      hcuts->GetXaxis()->SetBinLabel(ibin, "I");
    }

    if (!strcmp(CutName, "ISOCONE")) {
      ISOCONE = CutValue; ok = 1;
      if (dump) cout << "ISOCONE:           " << Form("%3.2f", ISOCONE) << endl;
      ibin = 109;
      hcuts->SetBinContent(ibin, ISOCONE);
      hcuts->GetXaxis()->SetBinLabel(ibin, "R_{I}");
    }

//     if (!strcmp(CutName, "BMMSEL")) {
//       BMMSEL = CutValue; ok = 1;
//       if (dump) cout << "BMMSEL:           " << BMMSEL << endl;
//       ibin = 110;
//       hcuts->SetBinContent(ibin, BMMSEL);
//       hcuts->GetXaxis()->SetBinLabel(ibin, "BMMSEL");
//     }


    if (!ok) cout << "==> ERROR: Don't know about variable " << CutName << endl;
  }

}

// ----------------------------------------------------------------------
void treeBmm::decayChannel(TString ch, int dump) {

  fD0 = -1;  fD1 = -1;  fD2 = -1;

  if ( !strcmp("bd2pi", ch.Data()) ) {

    fD0 = 211;  fD1 = 211;
    fChannel = TString("bd2pi");
    cout << "Selected decay mode: bd2pi." << endl;
  }
  else if ( !strcmp("bdpik", ch.Data()) ) {

    fD0 = 321;  fD1 = 211;
    fChannel = TString("bdpik");
    cout << "Selected decay mode: bdpik." << endl;
  }
  else if ( !strcmp("bdpimunu", ch.Data()) ) {

    fD0 = 13;  fD1 = 211;
    fChannel = TString("bdpimunu");
    cout << "Selected decay mode: bdpimunu." << endl;
  }
  else if ( !strcmp("bs2pi", ch.Data()) ) {

    fD0 = 211;  fD1 = 211;
    fChannel = TString("bs2pi");
    cout << "Selected decay mode: bs2pi." << endl;
  }
  else if ( !strcmp("bskk", ch.Data()) ) {

    fD0 = 321;  fD1 = 321;
    fChannel = TString("bskk");
    cout << "Selected decay mode: bskk." << endl;
  }
  else if ( !strcmp("bskpi", ch.Data()) ) {

    fD0 = 321;  fD1 = 211;
    fChannel = TString("bskpi");
    cout << "Selected decay mode: bskpi." << endl;
  }
  else if ( !strcmp("bskmunu", ch.Data()) ) {

    fD0 = 13;  fD1 = 321;
    fChannel = TString("bskmunu");
    cout << "Selected decay mode: bskmunu." << endl;
  }
  else if ( !strcmp("lbpk", ch.Data()) ) {

    fD0 = 2212;  fD1 = 321;
    fChannel = TString("lbpk");
    cout << "Selected decay mode: lbpk." << endl;
  }
  else if ( !strcmp("lbppi", ch.Data()) ) {

    fD0 = 2212;  fD1 = 211;
    fChannel = TString("lbppi");
    cout << "Selected decay mode: lbppi." << endl;
  }
  else if ( !strcmp("qcd", ch.Data()) ) {

    fD0 = 321;  fD1 = 211;
    fChannel = TString("qcd");
    cout << "Selected decay mode: qcd." << endl;
  }
  else if ( !strcmp("bsmumug", ch.Data()) || 
	    !strcmp("bsmumup0", ch.Data()) ||
	    !strcmp("bu3munu", ch.Data()) ||
	    !strcmp("bc3munu", ch.Data()) ||
	    !strcmp("bcjpmunu", ch.Data()) ) {
    
    fD0 = 13;  fD1 = 13;
    fChannel = TString("2mu_rare");
    cout << "Selected decay mode: 2mu (from rare decays)." << endl;
  }
  else if ( !strcmp("bjk", ch.Data()) ) {

    fD0 = 13;  fD1 = 13;  fD2 = 321;
    fChannel = TString("bdjpk");
    cout << "Selected decay mode: bdjpk." << endl;
  }
  else {
    
    fD0 = 13;  fD1 = 13;
    fChannel = TString("2mu");
    cout << "Selected decay mode: 2mu." << endl;
  }
}

// ----------------------------------------------------------------------
void treeBmm::eventProcessing() {



//   TGenCand *pG;
//   cout << "----------------------------------------------------------------------" << endl;
//   for (int i = 0; i < fpEvt->nGenCands(); ++i) {
//     pG = fpEvt->getGenCand(i);
    
//     int aid = TMath::Abs(pG->fID); 
//     if ( aid == 1 || aid == 2 ||
//  	 aid == 3 || aid == 4 || 
//  	 aid == 5 || aid == 6 || 
//  	 aid == 21) {

//     }

//     cout << i << " ID: " << pG->fID << " status: " << pG->fStatus << endl;
//   }

//   return;
 
  fpHistFile->cd();
  ((TH1D*)fpHistFile->Get("runs"))->Fill(fpEvt->fRunNumber);

  fER1->Fill(0.1); 

  initVariables();
  if (fDebug & 1) cout << "==> Event: " << fEvent << endl;

  processDiscrimination();

  //   // -- Generator-level preselection
  //   kinematicSelection(0);
  
  //   if (0 == fGoodKinematics) {
  //     if (fDebug & 1) cout << " --> no good kinematics ... " << endl;
  //   }

  // -- Generator-level process discrimation
  kinematicSelection(0);
  if (0 == fGoodKinematics) {
    if (fDebug & 1) cout << " --> no good kinematics ... " << endl;
  }
  
  muonMisIdRate();

  // -- L1 trigger?
  L1Selection(100);
  if (0 == fGoodL1Trigger) {
    if (fDebug & 1) cout << " --> no good L1 trigger ... " << endl;
  }

  // -- Remove events with no candidates (no PV or not enough high-pT muons)
  cutJunk(200);
  if (0 == fGoodEvent) {
    if (fDebug & 1) cout << " --> no cand, returning ... " << endl;
    return;
  }

  // -- Candidate track quality cuts?
  trackSelection(400);
  if (fGoodTT < 1) {
    if (fDebug & 1)  cout << " --> no good tracks ... " << fGoodTT << endl;
    fGoodTT = 0; 
  }

  // -- Candidate track PID?
  pidSelection(500);
  if (0 == fGoodPid) {
    if (fDebug & 1)  cout << " --> no good PID  ... " << endl;
  }

  // -- BS Candidate?
  if ( fD2 < 0 ) {
    candidateSelection(600);
  }
  else {

    normSample(700);
    kaonPlots();
  }
  if (0 == fGoodCand) {
    if (fDebug & 1)  cout << " --> no good CAND ... " << endl;
  }
  
  // -- HLT?
  HLTSelection(300);
  if (0 == fGoodHLTTrigger) {
    if (fDebug & 1) cout << " --> no good HLT trigger ... " << endl;
  }

  fillHist(); 
  vertexResolution();

  fTree->Fill();
  //  dumpSmallTree();

}



// ----------------------------------------------------------------------
void treeBmm::initVariables() {

  fpHistFile->cd();

  fAllEvent  = 1;
  fGoodEvent = 1;
  fGoodCand  = 1;

  fSignalBox = fLoSide = fHiSide = 0;

  TAnaCand *pCand;
  fKTI.clear();
  fKCI.clear();

  //  cout << "----------------------------------------------------------------------" << endl;
  //  cout << "-->initVariables> Event: " << fEvent << " with " << fpEvt->nCands() << " candidates"
  //       << endl;

  fPtL0  = fPtL1  = fPtK  = 0.;
  fEtaL0 = fEtaL1 = fEtaK =-99.;
  fTipL0 = fTipL1 = fTipK = 100.;
  fQL0   = fQL1   = fQK   = 0.;

  // -- Choose the B candidate type and its signal tracks
  fSTI[0] = fSTI[1] = fSTI[2] = -1;

  for (int i = 0; i < fpEvt->nCands(); ++i) {

    if (fpEvt->getCand(i)->fType == 1) {      // TYPE 1 = Bs -> mu mu
      
      fpB = fpEvt->getCand(i);
      fSTI[0] = fpB->fSig1;
      fSTI[1] = fpB->fSig2;
      fpL1 = fpEvt->getRecTrack(fSTI[0]);
      fpL2 = fpEvt->getRecTrack(fSTI[1]);
      
      fPtL0 = fpL1->fPlab.Pt();
      fPtL1 = fpL2->fPlab.Pt();
      
      fDptL0 = fpL1->fPlab.Pt();
      fDptL1 = fpL2->fPlab.Pt();
      
      fEtaL0 = fpL1->fPlab.Eta();
      fEtaL1 = fpL2->fPlab.Eta();
      
      fTipL0 = fpL1->fTip;
      fTipL1 = fpL2->fTip;
      
      fQL0   = fpL1->fQ;
      fQL1   = fpL2->fQ;
      
      if (fDebug & 2) {
	
	cout << "Choosing cand " << i << " which is of type " << fpB->fType
	     << " with pt,phi,eta = "
	     << fpB->fPlab.Pt()  << ", "
	     << fpB->fPlab.Phi() << ", "
	     << fpB->fPlab.Eta()
	     << " and has signal tracks at " << endl
	     << fSTI[0]
	     << " with pt,phi,eta = "
	     << fpEvt->getRecTrack(fSTI[0])->fPlab.Pt()  << ", "
	     << fpEvt->getRecTrack(fSTI[0])->fPlab.Phi() << ", "
	     << fpEvt->getRecTrack(fSTI[0])->fPlab.Eta()
	     << " and " << endl
	     << fSTI[1]
	     << " with pt,phi,eta = "
	     << fpEvt->getRecTrack(fSTI[1])->fPlab.Pt()  << ", "
	     << fpEvt->getRecTrack(fSTI[1])->fPlab.Phi() << ", "
	     << fpEvt->getRecTrack(fSTI[1])->fPlab.Eta()
	     << endl;
      }
    } else if (fpEvt->getCand(i)->fType == 2) {     // TYPE 2 = J/Psi -> mu mu
                                                
                                                     
      fpJpsi = fpEvt->getCand(i);
      fSTI[0] = fpJpsi->fSig1;
      fSTI[1] = fpJpsi->fSig2;
      fpL1 = fpEvt->getRecTrack(fSTI[0]);
      fpL2 = fpEvt->getRecTrack(fSTI[1]);
      
      fPtL0 = fpL1->fPlab.Pt();
      fPtL1 = fpL2->fPlab.Pt();
      
      fDptL0 = fpL1->fPlab.Pt();
      fDptL1 = fpL2->fPlab.Pt();
      
      fEtaL0 = fpL1->fPlab.Eta();
      fEtaL1 = fpL2->fPlab.Eta();
      
      fTipL0 = fpL1->fTip;
      fTipL1 = fpL2->fTip;
      
      fQL0   = fpL1->fQ;
      fQL1   = fpL2->fQ;
      
      if (fDebug & 2) {
	
	cout << "Choosing cand " << i << " which is of type " << fpJpsi->fType
	     << " with pt,phi,eta = "
	     << fpJpsi->fPlab.Pt()  << ", "
	     << fpJpsi->fPlab.Phi() << ", "
	     << fpJpsi->fPlab.Eta()
	     << " and has signal tracks at " << endl
	     << fSTI[0]
	     << " with pt,phi,eta = "
	     << fpEvt->getRecTrack(fSTI[0])->fPlab.Pt()  << ", "
	     << fpEvt->getRecTrack(fSTI[0])->fPlab.Phi() << ", "
	     << fpEvt->getRecTrack(fSTI[0])->fPlab.Eta()
	     << " and " << endl
	     << fSTI[1]
	     << " with pt,phi,eta = "
	     << fpEvt->getRecTrack(fSTI[1])->fPlab.Pt()  << ", "
	     << fpEvt->getRecTrack(fSTI[1])->fPlab.Phi() << ", "
	     << fpEvt->getRecTrack(fSTI[1])->fPlab.Eta()
	     << endl;
      }
    } else if (fpEvt->getCand(i)->fType == 3) {   // TYPE 3: B -> J/Psi K

      fpB = fpEvt->getCand(i);
      fSTI[2] = fpB->fSig2;
      fpK = fpEvt->getRecTrack(fSTI[2]);
      fPtK = fpK->fPlab.Pt();
      fDptK = fpK->fPlab.Pt();
      fEtaK = fpK->fPlab.Eta();
      fTipK = fpK->fTip;
      fQK   = fpK->fQ;

      if (fDebug & 2) {

	cout << "Choosing second cand " << i << " which is of type " << fpB->fType
	     << " with pt,phi,eta = "
	     << fpB->fPlab.Pt()  << ", "
	     << fpB->fPlab.Phi() << ", "
	     << fpB->fPlab.Eta()
	     << " and has signal track at " << endl
	     << fSTI[2]
	     << " with pt,phi,eta = "
	     << fpEvt->getRecTrack(fSTI[2])->fPlab.Pt()  << ", "
	     << fpEvt->getRecTrack(fSTI[2])->fPlab.Phi() << ", "
	     << fpEvt->getRecTrack(fSTI[2])->fPlab.Eta() 
	     << endl;
      }
    } else if (fpEvt->getCand(i)->fType == 0 ) {  

      pCand = fpEvt->getCand(i);
      fKCI.push_back(i);
      fKTI.push_back(pCand->fSig2);
  
    } else {
      pCand = fpEvt->getCand(i);
      cout << "Candidate " << i << " has wrong type: " << pCand->fType << endl;
    }
  }


  // ======================================================================
  // -- Systematics: Change some variables 
  // ======================================================================
  // -- muon ID
  if (0 && gRandom->Rndm() < 0.01) {
    fpEvt->fDiMuonTriggerDecision = 1; 
  }

  // -- L1
  if (0 && fpEvt->fDiMuonTriggerDecision == 0 && gRandom->Rndm() < 0.05) {
    fpEvt->fDiMuonTriggerDecision = 1; 
  }

  // -- Tracking efficiency 2*1%
  if (0) {
    if (gRandom->Rndm() < 0.005) {
      fSTI[0] = -1; 
    }

    if (gRandom->Rndm() < 0.005) {
      fSTI[1] = -1; 
    }
  }

  // -- Tracking resolution
  if (0) {
    double k0 = 1./fPtL0;
    double k1 = 1./fPtL1;

    double smear, sigma; 
    sigma = 0.0005;  // momentum scale
    //    sigma = 0.0004;  // misalignment
    //    sigma = 0.0003;  // B field

    smear = gRandom->Gaus(0., sigma);
    k0 += smear;
    smear = gRandom->Gaus(0., sigma);
    k1 += smear;

    fPtL0 = 1./k0;
    fPtL1 = 1./k1;
  }    


}



// ----------------------------------------------------------------------
void treeBmm::bookHist() {
  TH1 *h;

  cout << "-->bookHist> " << endl;

  // -- Reduced Tree
  fTree = new TTree("events", "events");

  fTree->Branch("run",    &fRun ,"run/I");

  fTree->Branch("process",&fProcessType ,"process/I");

  fTree->Branch("goodCand",       &fGoodCand ,      "goodCand/I");
  fTree->Branch("goodPV",         &fGoodPV ,        "goodPV/I");
  fTree->Branch("goodTT",         &fGoodTT ,        "goodTT/I");
  fTree->Branch("goodPid",        &fGoodPid ,       "goodPid/I");
  fTree->Branch("goodHLT",        &fGoodHLTTrigger ,"goodHLT/I");
  fTree->Branch("goodEvent",      &fGoodEvent ,     "goodEvent/I");
  fTree->Branch("goodL1",         &fGoodL1Trigger , "goodL1/I");
  fTree->Branch("goodKinematics", &fGoodKinematics ,"goodKinematics/I");

  fTree->Branch("dVtxPerp",     &fDvtxPerp,"dVtxPerp/D");
  fTree->Branch("dVtxPar",      &fDvtxPar, "dVtxPar/D");
  fTree->Branch("dPVtxX",       &fDpVtxX,  "dPVtxX/D");
  fTree->Branch("dPVtxY",       &fDpVtxY,  "dPVtxY/D");
  fTree->Branch("dPVtxZ",       &fDpVtxZ,  "dPVtxZ/D");

  fTree->Branch("rflt",         &fFltR,  "rflt/D");
  fTree->Branch("gflt",         &fFltG,  "gflt/D");
  fTree->Branch("dFltRes",      &fFltRes,  "dFltRes/D");


  fTree->Branch("mass",     &fMass,"mass/D");
  fTree->Branch("pt",       &fPt,  "pt/D");
  fTree->Branch("p",        &fP,   "p/D");
  fTree->Branch("eta",      &fEta, "eta/D");

  fTree->Branch("ptl0",     &fPtL0,  "ptl0/D");
  fTree->Branch("etal0",    &fEtaL0, "etal0/D");
  fTree->Branch("ql0",      &fQL0,   "ql0/D");
  fTree->Branch("tipl0",    &fTipL0, "tipl0/D");

  fTree->Branch("ptl1",     &fPtL1,  "ptl1/D");
  fTree->Branch("etal1",    &fEtaL1, "etal1/D");
  fTree->Branch("ql1",      &fQL1,   "ql1/D");
  fTree->Branch("tipl1",    &fTipL1, "tipl1/D");

  fTree->Branch("muIDl0",     &fMuIDL0,     "muIDl0/D");
  fTree->Branch("muIDl1",     &fMuIDL1,     "muIDl1/D");

  fTree->Branch("genMuIDl0",  &fGenMuIDL0,  "genMuIDl0/I");
  fTree->Branch("genMuIDl1",  &fGenMuIDL1,  "genMuIDl1/I");

  fTree->Branch("chi2",     &fChi2,      "chi2/D");
  fTree->Branch("cosa3",    &fCosAngle3, "cosa3/D");
  fTree->Branch("cosa",     &fCosAngle,  "cosa/D");
  fTree->Branch("l3d",      &fL3d,       "l3d/D");
  fTree->Branch("s3d",      &fS3d,       "s3d/D");
  fTree->Branch("lxy",      &fLxy,       "lxy/D");
  fTree->Branch("sxy",      &fSxy,       "sxy/D");
  fTree->Branch("tau",      &fTau,       "tau/D");
  fTree->Branch("txy",      &fTxy,       "txy/D");
  fTree->Branch("rmm",      &fRMM,       "rmm/D");
  fTree->Branch("i05",      &fIso05,     "i05/D");
  fTree->Branch("i06",      &fIso06,     "i06/D");
  fTree->Branch("i08",      &fIso08,     "i08/D");
  fTree->Branch("i10",      &fIso10,      "i10/D");
  fTree->Branch("i12",      &fIso12,     "i12/D");
  fTree->Branch("i14",      &fIso14,     "i14/D");
  fTree->Branch("iso",      &fIso,       "iso/D");
  fTree->Branch("isoveto",  &fIsoVeto,   "isoveto/I");
  fTree->Branch("iv05",     &fIsoVeto05, "iv05/I");
  fTree->Branch("iv06",     &fIsoVeto06, "iv06/I");
  fTree->Branch("iv08",     &fIsoVeto08, "iv08/I");
  fTree->Branch("iv10",     &fIsoVeto10, "iv10/I");
  fTree->Branch("iv12",     &fIsoVeto12, "iv12/I");
  fTree->Branch("iv14",     &fIsoVeto14, "iv14/I");

  // -- J/Psi K sample
  fTree->Branch("ptk",      &fPtK,  "ptk/D");
  fTree->Branch("etak",     &fEtaK, "etak/D");
  fTree->Branch("qk",       &fQK,   "qk/D");
  fTree->Branch("tipk",     &fTipK, "tipk/D");
  fTree->Branch("muIDk",      &fMuIDK,     "muIDk/D");
  fTree->Branch("genMuIDk",   &fGenMuIDK,  "genMuIDk/I");

  fTree->Branch("rmj1",      &fRMJ1,       "rmj1/D");
  fTree->Branch("rmj2",      &fRMJ2,       "rmj2/D");
  fTree->Branch("rjk",       &fRKJ,        "rkj/D");



  h = new TH1D("runs", "runs ", 100000, 0., 100000);
  h = new TH1D("mass", "mass ", 70, 5.0, 5.7);


  fER1  = new TH1D("ER1", "Event Reduction ", 1000, 0., 1000.);    fER1->Sumw2();
  fAR1  = new TH1D("AR1", "Analysis Reduction ", 1000, 0., 1000.); fAR1->Sumw2();

  // -- Tracking plots
  h = new TH1D("teff", "Track Cut efficiencies", 20, 0., 20.); h->Sumw2();

  h = new TH1D("t100", "pt leptons",  50, 0., 25.); h->Sumw2();  setTitles(h, "p_{T}^{#mu}", "events/bin");
  h = new TH1D("t101", "pt l1",       50, 0., 25.); h->Sumw2();  setTitles(h, "p_{T}^{#mu 1}", "events/bin");
  h = new TH1D("t102", "pt l2",       50, 0., 25.); h->Sumw2();  setTitles(h, "p_{T}^{#mu 2}", "events/bin");

  h = new TH1D("t110", "eta leptons", 50, -5., 5.); h->Sumw2();  setTitles(h, "#eta^{#mu}", "events/bin");
  h = new TH1D("t111", "eta l1",      50, -5., 5.); h->Sumw2();  setTitles(h, "#eta^{#mu 1}", "events/bin");
  h = new TH1D("t112", "eta l2",      50, -5., 5.); h->Sumw2();  setTitles(h, "#eta^{#mu 2}", "events/bin");

  h = new TH1D("t120", "tip leptons", 50, 0., 0.5); h->Sumw2();  setTitles(h, "tip^{#mu}_{xy}", "events/bin");
  h = new TH1D("t121", "tip l1",      50, 0., 0.5); h->Sumw2();  setTitles(h, "tip^{#mu 1}_{xy}", "events/bin");
  h = new TH1D("t122", "tip l2",      50, 0., 0.5); h->Sumw2();  setTitles(h, "tip^{#mu 2}_{xy}", "events/bin");

  h = new TH1D("mid", "Muon ID", 200, -20., 20.); h->Sumw2();

  // -- secondary vertex plots
  h = new TH1D("v100", "chi2", 100, 0., 20.); h->Sumw2();             setTitles(h, "#chi^{2}", "events/bin");
  h = new TH1D("v101", "prob(chi2, ndof)", 40, 0., 1.); h->Sumw2(); setTitles(h, "P(#chi^{2},ndof)", "events/bin");

  h = new TH1D("v110", "lxy", 50, 0., 2.0); h->Sumw2();     setTitles(h, "l_{xy} [cm]", "events/bin");
  h = new TH1D("v111", "sxy", 50, 0., 0.05); h->Sumw2();     setTitles(h, "#sigma_{xy} [cm]", "events/bin");
  h = new TH1D("v112", "lxy/sxy", 50, 0., 50.); h->Sumw2(); setTitles(h, "l_{xy}/#sigma_{xy}", "events/bin");

  h = new TH1D("v120", "l3d", 50, 0., 2.0); h->Sumw2();     setTitles(h, "l_{3d} [cm]", "events/bin");
  h = new TH1D("v121", "s3d", 50, 0., 0.05); h->Sumw2();     setTitles(h, "#sigma_{3d} [cm]", "events/bin");
  h = new TH1D("v122", "l3d/s3d", 50, 0., 50.); h->Sumw2(); setTitles(h, "l_{3d}/#sigma_{3d}", "events/bin");

  // -- sec. vertex resolution plots
  h = new TH1D("v200", "Decay length resolution (2D)", 500, -0.05, 0.05); h->Sumw2(); 
  setTitles(h, "l_{xy}^{rec} - l_{xy}^{sim} [cm]", "events/bin");
  h = new TH1D("v201", "Decay length resolution (3D)", 500, -0.1, 0.1); h->Sumw2(); 
  setTitles(h, "l_{3D}^{rec} - l_{3D}^{sim} [cm]", "events/bin");
  h = new TH1D("v202", "Proper decay time resolution (2D)", 500, -0.05, 0.05); 
  h->Sumw2(); setTitles(h, "#tau_{xy}^{rec} - #tau_{xy}^{sim} [cm/c]", "events/bin");
  h = new TH1D("v203", "Proper decay time resolution (3D)", 500, -0.05, 0.05); 
  h->Sumw2(); setTitles(h, "#tau_{3D}^{rec} - #tau_{3D}^{sim} [cm/c]", "events/bin");
  h = new TH1D("v204", "Proper decay time resolution (2D)", 500, -1000., 1000.); 
  h->Sumw2(); setTitles(h, "#tau_{xy}^{rec} - #tau_{xy}^{sim} [fs]", "events/bin");
  h = new TH1D("v205", "Proper decay time resolution (3D)", 500, -1000., 1000.); 
  h->Sumw2(); setTitles(h, "#tau_{3D}^{rec} - #tau_{3D}^{sim} [fs]", "events/bin");

  // -- Jpsi candidate plots
  h = new TH1D("j100", "mass", 500, 0., 50.); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("j101", "mass", 60, 2.8, 3.4); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("j102", "pT", 60, 0., 30.); h->Sumw2();   setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("j103", "eta", 50, -5., 5.); h->Sumw2();  setTitles(h, "#eta", "events/bin");
  h = new TH1D("j104", "proper decay time", 100, 0., 0.5); h->Sumw2(); setTitles(h, "#tau_{B} [...]", "events/bin");

  h = new TH1D("j150", "deltaRmm", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(#mu#mu)", "events/bin");
  h = new TH1D("j151", "deltaRkj", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi K)", "events/bin");
  h = new TH1D("j152", "deltaRmj1", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{1})", "events/bin");
  h = new TH1D("j153", "deltaRmj2", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{2})", "events/bin");

  // -- B candidate plots
  h = new TH1D("b100", "mass", 500, 0., 50.); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b101", "mass", 100, 5., 6.); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b102", "pT", 60, 0., 30.); h->Sumw2();   setTitles(h, "p_{T} [GeV]", "events/bin");
  h = new TH1D("b103", "eta", 50, -5., 5.); h->Sumw2();  setTitles(h, "#eta", "events/bin");
  h = new TH1D("b104", "proper decay time", 100, 0., 0.5); h->Sumw2(); setTitles(h, "#tau_{B} [...]", "events/bin");

  h = new TH1D("b110", "cos(angle)",              100, 0.99, 1.0); h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b120", "cos(angle), pT(B) > 5",   100, 0.99, 1.0); h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b130", "cos(angle), l3d/s3d > 4", 100, 0.99, 1.0); h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b131", "cos(angle), l3d/s3d > 10",100, 0.99, 1.0); h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b140", "cos(angle), lxy/sxy > 4", 100, 0.99, 1.0); h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b141", "cos(angle), lxy/sxy > 10",100, 0.99, 1.0); h->Sumw2();   setTitles(h, "#cos(angle)(p,v)", "events/bin");
  h = new TH1D("b150", "deltaRmm", 50, 0., 5.0); h->Sumw2();   setTitles(h, "#DeltaR(#mu#mu)", "events/bin");

  h = new TH1D("b200", "isolation", 40, 0., 1.); h->Sumw2();   setTitles(h, "Isolation", "events/bin");
  h = new TH1D("b201", "isolation veto", 10, 0., 10.); h->Sumw2();   setTitles(h, "Isolation veto", "events/bin");

  // -- Generator Level
  h = new TH1D("g100", "generator pT ", 40, 0., 80.0); h->Sumw2();   setTitles(h, "pT(#mu)", "events/bin");
  h = new TH1D("g102", "generator eta ", 60, -3.0, 3.0); h->Sumw2();   setTitles(h, "eta(#mu)", "events/bin");
  h = new TH1D("g103", "generator mass ", 100, 0., 10.0); h->Sumw2();   setTitles(h, "m_{#mu#mu}", "events/bin");
  h = new TH1D("g104", "generator status", 11, -1., 10.); h->Sumw2();

  h = new TH1D("g110", "generator pT ", 80, 0., 40.); h->Sumw2();   setTitles(h, "pT(K)", "events/bin");
  h = new TH1D("g112", "generator eta ", 60, -3.0, 3.); h->Sumw2();   setTitles(h, "eta(K)", "events/bin");
  h = new TH1D("g113", "generator mass ", 100, 0., 10.); h->Sumw2();   setTitles(h, "m_{J/#Psi K}", "events/bin");
  h = new TH1D("g115", "generator pT ", 80, 0., 40.); h->Sumw2();   setTitles(h, "pT(#mu_{1})", "events/bin");
  h = new TH1D("g116", "generator pT ", 80, 0., 40.); h->Sumw2();   setTitles(h, "pT(#mu_{2})", "events/bin");

  h = new TH1D("g150", "deltaRmm", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(#mu#mu)", "events/bin");
  h = new TH1D("g151", "deltaRkj", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi K)", "events/bin");
  h = new TH1D("g152", "deltaRmj1", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{1})", "events/bin");
  h = new TH1D("g153", "deltaRmj2", 50, 0., 2.5); h->Sumw2();   setTitles(h, "#DeltaR(J/#Psi #mu_{2})", "events/bin");



  // -- Event selection
  h = new TH1D("eeff", "Event Cut efficiencies", 50, 0., 50.); h->Sumw2();

  h = new TH1D("aeff", "Event Cut efficiencies", 500, 0., 500.); h->Sumw2();


  // -- Results
  h = new TH1D("b1000", "mass", 70, 5.0, 5.7); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b1001", "mass", 70, 5.0, 5.7); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b1002", "mass", 70, 5.0, 5.7); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");

  h = new TH1D("b1010", "mass", 70, 5.0, 5.7); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b1011", "mass", 70, 5.0, 5.7); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");
  h = new TH1D("b1012", "mass", 70, 5.0, 5.7); h->Sumw2(); setTitles(h, "m_{#mu^{+}#mu^{-}} [GeV]", "events/bin");


  book(0);   // histogram(0) reco level, before any cuts
  book(1);   // histogram(1) reco level, after Kinematic cuts on generator level have been applied
  book(2);   // histogram(2) reco level, after L1 cuts
  book(3);   // histogram(3) reco level, after HLT cuts
  book(4);   // histogram(4) reco level, after all offline cuts
  book(5);   // histogram(5) reco level, after factorized offline cuts

  book2(0);   // histogram(0) reco level, before any cuts
  book2(1);   // histogram(1) reco level, after Kinematic cuts on generator level have been applied
  book2(2);   // histogram(2) reco level, after L1 cuts
  book2(3);   // histogram(3) reco level, after HLT cuts
  book2(4);   // histogram(4) reco level, after all offline cuts
  book2(5);   // histogram(5) reco level, after factorized offline cuts

  h = new TH1D("d100", "p_{T, #mu} [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("d101", "TIP_{#mu}", 50, 0., 0.2); h->Sumw2();
  h = new TH1D("d102", "#chi^{2}", 50, 0., 25.); h->Sumw2();
  h = new TH1D("d103", "l_{3d}", 50, 0., 0.25); h->Sumw2();
  h = new TH1D("d104", "m_{#mu #mu} [GeV]", 70, 5.0, 5.7); h->Sumw2();
  h = new TH1D("d105", "N_{Trk}", 150, 0., 150); h->Sumw2();
  h = new TH1D("d106", "isolation veto", 11, -1., 10.); h->Sumw2();
  h = new TH1D("d107", "isolation veto, fIsoVeto == 0", 50, -0.1, 10.); h->Sumw2();
  h = new TH1D("d108", "isolation veto, m > 4", 50, -0.1, 10.); h->Sumw2();
  h = new TH1D("d109", "isolation veto, m ~ J/psi", 50, -0.1, 10.); h->Sumw2();
  h = new TH1D("d110", "#Delta #phi", 100, -5., 5.); h->Sumw2();
  h = new TH1D("d111", "#eta_{#mu}", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("d112", "p_{T, Trk} [GeV]", 100, 0., 10.);  h->Sumw2();
  h = new TH1D("d113", "p_{T, Trk in cone} [GeV]", 100, 0., 10.);  h->Sumw2();
  h = new TH1D("d114", "Mother of 1. muon", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("d115", "Mother of 2. muon", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("d116", "Mother of kaon", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("d117", "Grand-Mother of 1. muon", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("d118", "Grand-Mother of 2. muon", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("d119", "Grand-Mother of kaon", 5000, 0., 5000.);  h->Sumw2();

  // -- Kaon candidate plots for B+ -> J/Psi K+

  h = new TH1D("keff", "number of truth-matched kaons", 100, 0., 100.);

  h = new TH1D("k100", "mass", 200, 0., 10.);
  h = new TH1D("k101", "#chi^{2}", 500, 0., 50.);
  h = new TH1D("k102", "p_{T}^{B}", 60, 0., 30.);
  h = new TH1D("k103", "#eta^{B}",  50, -5., 5.);
  h = new TH1D("k104", "p_{T}^{K}", 50, 0., 25.);
  h = new TH1D("k105", "#eta^{K}",  50, -5., 5.);

  // -- Profiles
  TProfile *p1 = new TProfile("p120","#sigma_{p_{T},#mu} vs p_{T}", 25, 0., 25., "s");
  p1 = new TProfile("p121","#sigma_{p_{T},#mu}/p_{T} vs p_{T}", 25, 0., 25., "s");
  p1 = new TProfile("p122","#sigma_{p_{T}} vs p_{T}", 25, 0., 25., "s");
  p1 = new TProfile("p123","#sigma_{p_{T}}/p_{T} vs p_{T}", 25, 0., 25., "s");

  p1 = new TProfile("p130","#sigma_{p_{T}, #mu} vs #eta", 25, 0., 2.5, "s");
  p1 = new TProfile("p131","#sigma_{p_{T}, #mu}/p_{T}  vs #eta", 25, 0., 2.5, "s");
  p1 = new TProfile("p132","#sigma_{p_{T}} vs #eta", 25, 0., 2.5, "s");
  p1 = new TProfile("p133","#sigma_{p_{T}}/p_{T}  vs #eta", 25, 0., 2.5, "s");

  h  = new TH1D("p120C", "Counter for profile p120, p121", 25, 0., 25.);  h->Sumw2();
  h  = new TH1D("p122C", "Counter for profile p122, p123", 25, 0., 25.);  h->Sumw2();
  h  = new TH1D("p130C", "Counter for profile p130, p131", 25, 0., 2.5);  h->Sumw2();
  h  = new TH1D("p132C", "Counter for profile p132, p133", 25, 0., 2.5);  h->Sumw2();

  // -- 2D-histograms
  TH2D *h2 = new TH2D("D100", "fRMM vs. fMass", 60, 0., 6., 50, 0., 10.); h2->Sumw2();
  h2 = new TH2D("D101", "fRMM vs. fMass", 60, 0., 6., 50, 0., 10.); h2->Sumw2();
  h2 = new TH2D("D102", "drmin vs. cone size", 50, 0., 10., 50, 0., 10.); h2->Sumw2();

  h2 = new TH2D("D120", "#sigma_{p_{T},#mu} vs p_{T}", 25, 0., 25., 100, -1., 1.);
  h2 = new TH2D("D121", "#sigma_{p_{T},#mu}/p_{T} vs p_{T}", 25, 0., 25., 100, -1., 1.);
  h2 = new TH2D("D122", "#sigma_{p_{T}, tracks} vs p_{T}", 25, 0., 25., 100, -1., 10.);
  h2 = new TH2D("D123", "#sigma_{p_{T}, tracks}/p_{T} vs p_{T}", 25, 0., 25., 100, -1., 10.);

  h2 = new TH2D("D130", "#sigma_{p_{T}, #mu} vs #eta", 25, 0., 2.5, 100, -1., 1.);
  h2 = new TH2D("D131", "#sigma_{p_{T}, #mu}/p_{T}  vs #eta", 25, 0., 2.5, 100, -1., 1.);
  h2 = new TH2D("D132", "#sigma_{p_{T}, tracks} vs #eta", 25, 0., 2.5, 100, -1., 10.);
  h2 = new TH2D("D133", "#sigma_{p_{T}, tracks}/p_{T}  vs #eta", 25, 0., 2.5, 100, -1., 10.);

  h2 = new TH2D("K100", "#Delta m vs. gen. PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K101", "#Delta #chi^{2} vs. gen. PDG", 5000, 0., 5000., 500, -25., 25.);
  h2 = new TH2D("K102", "#Delta m vs. gen. mother PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K103", "#Delta #chi^{2} vs. gen. mother PDG", 5000, 0., 5000., 500, -25., 25.);

  h2 = new TH2D("K200", "#Delta m vs. gen. PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K201", "#Delta #chi^{2} vs. gen. PDG", 5000, 0., 5000., 500, -25., 25.);
  h2 = new TH2D("K202", "#Delta m vs. gen. mother PDG", 5000, 0., 5000., 200, -5., 5.);
  h2 = new TH2D("K203", "#Delta #chi^{2} vs. gen. mother PDG", 5000, 0., 5000., 500, -25., 25.);

  // Mis-Id. plots
  fMisID = new TH1D("MisID", "Muon mis-ID Rate", 50, 0., 50.);      fMisID->Sumw2();

  // -- 1D-histograms
  // -- pions
  h = new TH1D("m000", "PDG id all tracks", 5000, 0., 5000.);  h->Sumw2();
  h = new TH1D("m100", "p_{T} (all #pi) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m101", "p_{T} (#pi = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m102", "p_{T} (#pi != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m110", "#eta (all #pi)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m111", "#eta (#pi = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m112", "#eta (#pi != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- kaons
  h = new TH1D("m200", "p_{T} (all K) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m201", "p_{T} (K = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m202", "p_{T} (K != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m210", "#eta (all K)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m211", "#eta (K = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m212", "#eta (K != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- protons
  h = new TH1D("m300", "p_{T} (all p) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m301", "p_{T} (p = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m302", "p_{T} (p != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m310", "#eta (all p)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m311", "#eta (p = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m312", "#eta (p != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- not identified muons
  h = new TH1D("m400", "p_{T} (all #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m401", "p_{T} (#mu = #mu) [GeV]", 50, 0., 25.);  h->Sumw2();
  h = new TH1D("m402", "p_{T} (#mu != #mu) [GeV]", 50, 0., 25.);  h->Sumw2();

  h = new TH1D("m410", "#eta (all #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m411", "#eta (#mu = #mu)", 50, -5., 5.);  h->Sumw2();
  h = new TH1D("m412", "#eta (#mu != #mu)", 50, -5., 5.);  h->Sumw2();

  // -- 2D-histograms
  // -- pions
  h2 = new TH2D("M100", "p_{T} vs. #eta (all #pi) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M101", "p_{T}  vs. #eta (#pi = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M102", "p_{T} vs. #eta  (#pi != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- kaons
  h2 = new TH2D("M200", "p_{T} vs. #eta (all K) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M201", "p_{T}  vs. #eta (K = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M202", "p_{T} vs. #eta  (K != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- protons
  h2 = new TH2D("M300", "p_{T} vs. #eta (all protons) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M301", "p_{T}  vs. #eta (p = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M302", "p_{T} vs. #eta  (p != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();

  // -- not identified muons
  h2 = new TH2D("M400", "p_{T} vs. #eta (all #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M401", "p_{T}  vs. #eta (#mu = #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();
  h2 = new TH2D("M402", "p_{T} vs. #eta  (#mu != #mu) [GeV]", 50, 0., 25., 50, -5., 5.);  h2->Sumw2();


}


// ----------------------------------------------------------------------
void treeBmm::book(int offset) {

  TH1D *h;
  h = new TH1D(Form("c%d00", offset), "pT(B)",          60, 0., 30.); h->Sumw2();     setTitles(h,  "p_{T, B} [GeV]", "events/bin");
  h = new TH1D(Form("c%d01", offset), "eta(B)",         50, -5., 5.); h->Sumw2();     setTitles(h,  "#eta_{B}", "events/bin");

  h = new TH1D(Form("c%d10", offset), "leading pT leptons", 50, 0., 25.); h->Sumw2(); setTitles(h,  "p_{T, #mu}^{max} [GeV]", "events/bin");
  h = new TH1D(Form("c%d11", offset), "eta leptons",    50, -5., 5.); h->Sumw2();     setTitles(h, "#eta_{#mu}", "events/bin");
  h = new TH1D(Form("c%d12", offset), "pT leptons",     50, 0., 25.); h->Sumw2();     setTitles(h,  "p_{T, #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d13", offset), "TIP leptons",    50, 0., 0.2); h->Sumw2();     setTitles(h,  "TIP_{#mu} [cm]", "events/bin");

  h = new TH1D(Form("c%d20", offset), "deltaRmm",       50, 0., 5.0); h->Sumw2();     setTitles(h, "#Delta R(#mu #mu)", "events/bin");
  h = new TH1D(Form("c%d21", offset), "cos(angle)",    200,0.95,1.0); h->Sumw2();     setTitles(h, "cos #alpha", "events/bin");
  h = new TH1D(Form("c%d22", offset), "lxy/sxy",        50, 0., 50.); h->Sumw2();     setTitles(h, "l_{xy}/#sigma_{xy}", "events/bin");
  h = new TH1D(Form("c%d23", offset), "lxy",           100, 0., 1.0); h->Sumw2();     setTitles(h, "l_{xy} [cm]", "events/bin");
  h = new TH1D(Form("c%d24", offset), "isolation veto", 20, 0., 20.); h->Sumw2();     setTitles(h, "I_{V}", "events/bin");
  h = new TH1D(Form("c%d25", offset), "isolation",      55, 0., 1.1); h->Sumw2();     setTitles(h, "I", "events/bin");
  h = new TH1D(Form("c%d26", offset), "l3d",            50, 0., 0.5); h->Sumw2();     setTitles(h, "l_{3D} [cm]", "events/bin");
  h = new TH1D(Form("c%d27", offset), "chi2",           50, 0., 5.0); h->Sumw2();     setTitles(h, "#chi^{2}", "events/bin");
  h = new TH1D(Form("c%d29", offset), "process",       52, 0., 52.);    h->Sumw2();     setTitles(h, "Process (GGF,FEX,GSP for t,b,c) ", "events/bin");

  if ( fD2 < 0 ) {
    h = new TH1D(Form("c%d30", offset), "mass",           100, 5.0, 6.0); h->Sumw2();     setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  }
  else {
    h = new TH1D(Form("c%d30", offset), "mass",           100, 4.9, 5.9); h->Sumw2();     setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  }

  h = new TH1D(Form("c%d31", offset), "mass",          500, 0., 50.); h->Sumw2();     setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d32", offset), "mass",           34, 4., 5.7); h->Sumw2();     setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");
  h = new TH1D(Form("c%d33", offset), "mass",          160, 2., 10.); h->Sumw2();     setTitles(h, "m_{#mu #mu} [GeV]", "events/bin");

  h = new TH1D(Form("c%d35", offset), "NTRK",          100, 0., 200.); h->Sumw2();     setTitles(h, "N_{Track}", "events/bin");

  h = new TH1D(Form("c%d40", offset), "cos(angle)",    200,-1.0,1.0); h->Sumw2();     setTitles(h, "cos #alpha", "events/bin");
  h = new TH1D(Form("c%d41", offset), "cos(angle)",    200,0.95,1.0); h->Sumw2();     setTitles(h, "cos #alpha", "events/bin");
  h = new TH1D(Form("c%d42", offset), "cos(angle)",     50,0.97,1.0); h->Sumw2();     setTitles(h, "cos #alpha", "events/bin");
  h = new TH1D(Form("c%d43", offset), "cos(angle)",    200,0.99,1.0); h->Sumw2();     setTitles(h, "cos #alpha", "events/bin");

  h = new TH1D(Form("c%d50", offset), "flight length", 100,  0.0, 0.1); h->Sumw2();   setTitles(h, "t_{rec}", "events/bin");
  h = new TH1D(Form("c%d51", offset), "flight length", 100,  0.0, 0.1); h->Sumw2();   setTitles(h, "t_{gen}", "events/bin");
  h = new TH1D(Form("c%d52", offset), "flight length", 100, -0.02, 0.02); h->Sumw2(); setTitles(h, "t_{rec} - t_{gen}", "events/bin");

}

void treeBmm::book2(int offset) {

  TH2D *h;

  h = new TH2D(Form("C%d00", offset), "dpT vs. sxy", 50, 0., 0.05, 50, 0., 0.05); h->Sumw2();  
  setTitles2(h, "#sigma_{pT}", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d01", offset), "pT vs. eta",   50, 0., 25., 50, -5., 5); h->Sumw2(); 
  setTitles2(h, "p_{T, #mu} [GeV]", "#eta_{#mu} [cm]");

  h = new TH2D(Form("C%d02", offset), "pT vs. sxy", 50, 0., 25., 50, 0., 0.05); h->Sumw2(); 
  setTitles2(h, "p_{T,#mu}", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d03", offset), "pT vs. lxy/sxy", 50, 0., 25., 50, 0., 50.); h->Sumw2();  
  setTitles2(h, "p_{T,#mu}", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d04", offset), "pT vs. TIP", 50, 0., 25., 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "TIP [mm]");

  h = new TH2D(Form("C%d05", offset), "pT vs. mass", 50, 0., 25., 500, 0., 50.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "m_{#mu #mu} [GeV]");

  h = new TH2D(Form("C%d06", offset), "pT vs. deltaR", 50, 0., 25., 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d07", offset), "pT vs. cos(angle)", 50, 0., 25., 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "cos #alpha(p,v)");

  h = new TH2D(Form("C%d08", offset), "pT vs. isolation veto", 50, 0., 25., 10, 0., 10.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "Isolation veto");

  h = new TH2D(Form("C%d09", offset), "pT vs. chi2", 50, 0., 25., 100, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "#chi^{2}");

  h = new TH2D(Form("C%d10", offset), "pT vs. pT(B)", 50, 0., 25., 60, 0., 30.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "p_{T, B} [GeV]");

  h = new TH2D(Form("C%d11", offset), "pT vs. eta(B)", 50, 0., 25., 50, -5., 5.); h->Sumw2();
  setTitles2(h, "p_{T,#mu}", "#eta_{B}");

  h = new TH2D(Form("C%d12", offset), "eta vs. sxy", 50, -5., 5, 50, 0., 0.05); h->Sumw2(); 
  setTitles2(h, "#eta_{#mu}", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d13", offset), "eta vs. lxy/sxy", 50, -5., 5,  50, 0., 50.); h->Sumw2();  
  setTitles2(h, "#eta_{#mu}", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d14", offset), "eta vs. TIP", 50, -5., 5, 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "TIP [mm]");

  h = new TH2D(Form("C%d15", offset), "eta vs. mass", 50, -5., 5, 500, 0., 50.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "m_{#mu #mu} [GeV]");

  h = new TH2D(Form("C%d16", offset), "eta vs. deltaR", 50, -5., 5, 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d17", offset), "eta vs.cos(angle)", 50, -5., 5, 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "cos #alpha(p,v)");

  h = new TH2D(Form("C%d18", offset), "eta vs. isolation veto", 50, -5., 5, 10, 0., 10.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "Isolation veto");

  h = new TH2D(Form("C%d19", offset), "eta vs. chi2", 50, -5., 5, 100, 0., 5.0); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "#chi^{2}");

  h = new TH2D(Form("C%d20", offset), "eta vs. pT(B)", 50, -5., 5, 60, 0., 30.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "p_{T, B} [GeV]");

  h = new TH2D(Form("C%d21", offset), "eta vs. eta(B)", 50, -5., 5, 50, -5., 5.); h->Sumw2();
  setTitles2(h, "#eta_{#mu}", "#eta_{B}");
  //--------------------new-----------------------------------------------------------------

  h = new TH2D(Form("C%d22", offset), "mass vs. sxy", 500, 0., 50, 50, 0., 0.05); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#sigma_{xy} [cm]");

  h = new TH2D(Form("C%d23", offset), "mass vs. lxy/sxy", 500, 0., 50, 50, 0., 50.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d24", offset), "mass vs. TIP", 500, 0., 50, 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "TIP [mm]");

  h = new TH2D(Form("C%d25", offset), "mass vs. deltaR", 500, 0., 50., 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d26", offset), "mass vs. cos(angle)", 500, 0., 50, 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "cos #alpha(p,v)");

  h = new TH2D(Form("C%d27", offset), "mass vs. isolation veto", 500, 0., 50, 10, 0., 10.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "Isolation veto");

  h = new TH2D(Form("C%d28", offset), "mass vs. chi2", 500, 0., 50,  100, 0., 5.0); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#chi^{2}");

  h = new TH2D(Form("C%d29", offset), "mass vs. pT(B)", 500, 0., 50.,  60, 0., 30.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "p_{T, B} [GeV]");

  h = new TH2D(Form("C%d30", offset), "mass vs. eta(B)", 500, 0., 50., 50, -5., 5.); h->Sumw2();
  setTitles2(h, "m_{#mu #mu} [GeV]", "#eta_{B}");
  //--------------------new-----------------------------------------------------------------

  h = new TH2D(Form("C%d31", offset), "pT(B) vs. sxy",            60, 0., 30., 50, 0., 0.05); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#sigma_{xy} [cm]");
  
  h = new TH2D(Form("C%d32", offset), "pT(B) vs. lxy/sxy",        60, 0., 30., 50, 0., 50.); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "l_{xy}/#sigma_{xy}");

  h = new TH2D(Form("C%d33", offset), "pT(B) vs. TIP",            60, 0., 30., 50, 0., 0.2); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "TIP [mm]");

  h = new TH2D(Form("C%d34", offset), "pT(B) vs. deltaR",         60, 0., 30., 50, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#Delta R(#mu #mu)");

  h = new TH2D(Form("C%d35", offset), "pT(B) vs. cos(angle)",     60, 0., 30., 100, 0.95, 1.0); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "cos #alpha(p,v) ");

  h = new TH2D(Form("C%d36", offset), "pT(B) vs. isolation veto", 60, 0., 30., 10, 0., 10.); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "Isolation veto");

  h = new TH2D(Form("C%d37", offset), "pT(B) vs. chi2",           60, 0., 30., 100, 0., 5.0); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#chi^{2}");

  h = new TH2D(Form("C%d38", offset), "pT(B) vs. eta(B)",         60, 0., 30., 50, -5., 5.); h->Sumw2();
  setTitles2(h, "p_{T, B} [GeV]", "#eta_{B}");
  
}


// ----------------------------------------------------------------------
void treeBmm::processDiscrimination() {

  TGenCand *pG;

  // documentation line partons (entries { d, u, s, c, b, t } )
  double docPartCnt[6];
  double docAntiCnt[6];
  // partons
  double parPartCnt[6];
  double parAntiCnt[6];    
    
  for (int i = 0; i < 6; i++) {
    docPartCnt[i] = 0; 
    docAntiCnt[i] = 0; 
    parPartCnt[i] = 0; 
    parAntiCnt[i] = 0; 
  }

  int aid(0);
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pG = fpEvt->getGenCand(i);

    ((TH1D*)gDirectory->Get("g104"))->Fill(pG->fStatus);
    
    aid = TMath::Abs(pG->fID); 
    if ( aid == 1 || aid == 2 ||
 	 aid == 3 || aid == 4 || 
 	 aid == 5 || aid == 6 || 
 	 aid == 21) {
      if ( pG->fStatus == 3 ) {
	// 	cout << "quark/gluon from documentation #" << i << "(ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 2 &&  TMath::Abs(pG->fID) != 21) {
	// 	cout << "decayed quark/gluon #" << i << " (ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 1 ) {
	// 	cout << "undecayed (?) quark/gluon #" << i << " (ID: " << pG->fID  << ")" << endl;
      }
    }

    for (int j = 0; j < 6; j++) {

      if ( pG->fStatus == 3 ) {
	if ( pG->fID == j+1 ) {  
	  docPartCnt[j]++;
	}
	if ( pG->fID == -(j+1) ) {  
	  docAntiCnt[j]++;
	}
      }

      if ( pG->fStatus == 2 ) {
	if ( pG->fID == j+1 ) {  
	  parPartCnt[j]++;
	}
	if ( pG->fID == -(j+1) ) {  
	  parAntiCnt[j]++;
	}
      }
    }
  }


//   cout << "documentation partons: " << endl;
//   for (int j = 0; j < 6; j++) {
//     cout << docPartCnt[j] << " ";
//   } 
//   cout << endl;

//   cout << "documentation antipartons: " << endl;
//   for (int j = 0; j < 6; j++) {
//     cout << docAntiCnt[j] << " ";
//   } 
//   cout << endl;


  fProcessType = -99;
  // GGF heavy flavour
  if (docPartCnt[5] == 1 && docAntiCnt[5] == 1) {
    fProcessType = 50; 
    //    printf("====> t: GGF (%i)\n", fProcessType);
    return;
  }
  
  if (docPartCnt[4] == 1 && docAntiCnt[4] == 1) {
    fProcessType = 40; 
    //    printf("====> b: GGF (%i)\n", fProcessType);
    return;
  } 
  
  if (docPartCnt[3] == 1 && docAntiCnt[3] == 1) {
    fProcessType = 30; 
    //    printf("====> c: GGF (%i)\n", fProcessType);
    return;
  }
  
  // FEX heavy flavour
  if ((docPartCnt[5] >= 1 && docAntiCnt[5] == 0) || (docPartCnt[5] == 0 && docAntiCnt[5] >= 1) ) {
    fProcessType = 51; 
    //    printf("====> t: FEX (%i)\n", fProcessType);
    return;
  }
  
  if ((docPartCnt[4] >= 1 && docAntiCnt[4] == 0) || (docPartCnt[4] == 0 && docAntiCnt[4] >= 1) ) {
    fProcessType = 41; 
    //    printf("====> b: FEX (%i)\n", fProcessType);
    return;
  }
  
  if ((docPartCnt[3] >= 1 && docAntiCnt[3] == 0) || (docPartCnt[3] == 0 && docAntiCnt[3] >= 1) ) {
    fProcessType = 31; 
    //    printf("====> c: FEX (%i)\n", fProcessType);
    return;
  }
  
  // GSP heavy flavour
  if (docPartCnt[5] == 0 && docAntiCnt[5] == 0 && (parPartCnt[5] >= 1 || parAntiCnt[5] >= 1)) {
    fProcessType = 52;
    //    printf("====> t: GSP (%i)\n", fProcessType); 
    return;
  }
  
  if (docPartCnt[4] == 0 && docAntiCnt[4] == 0 && (parPartCnt[4] >= 1 || parAntiCnt[4] >= 1)) {
    fProcessType = 42; 
    //    printf("====> b: GSP (%i)\n", fProcessType);
    return;
  }
  
  if (docPartCnt[3] == 0 && docAntiCnt[3] == 0 && (parPartCnt[3] >= 1 || parAntiCnt[3] >= 1)) {
    fProcessType = 32; 
    //    printf("====> c: GSP (%i)\n", fProcessType);
    return;
  }
  
  // light flavors
  if ((docPartCnt[5] == 0 && docAntiCnt[5] == 0) && (parPartCnt[5] == 0 && parAntiCnt[5] == 0)
      && (docPartCnt[4] == 0 && docAntiCnt[4] == 0) && (parPartCnt[4] == 0 && parAntiCnt[4] == 0)
      && (docPartCnt[3] == 0 && docAntiCnt[3] == 0) && (parPartCnt[3] == 0 && parAntiCnt[3] == 0)
      ) {
    fProcessType = 1; 
    //    printf("====> UDS: light flavors (%i)\n", fProcessType);
    return;
  }

  // if no process type was determined
  //  printf("====> Could not determine process type !!!\n");


}


// ----------------------------------------------------------------------
// -- This is on the GENERATOR level
void treeBmm::kinematicSelection(int o) {
  // offset: 0

  TGenCand *pG;
  TGenCand *pM;
  TGenCand *pGM;

  double pT(0.);
  double eta(0.);

  int cnt(0);  // muon cnt
  int kcnt(0); // kaon cnt

  int m0(-1), m1(-1), m2(-1); 
  int pid0(-1), pid1(-1), pidTmp(-1), pid2(-1);
  double pT0(0.), pT1(0.), pT2(0.); 
 
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {

    pG = fpEvt->getGenCand(i);
    
    if (              pG->fP.Pt() < fMinPt   ) { continue; }
    if ( TMath::Abs(pG->fP.Eta()) > fMaxEta   ) { continue; }
    
    if ( fD2 < 0 ) {
      if ( TMath::Abs(pG->fID) == fD0 || TMath::Abs(pG->fID) == fD1 || !strcmp("qcd", fChannel) ) {

	pT = pG->fP.Pt();
	eta = pG->fP.Eta();
	if ( (pT > fMinPt) && (TMath::Abs(eta) < fMaxEta) ) {
	  ++cnt;
	}
      
	// -- Get the highest pT muons sorted into m0, m1
	if (pT > pT0) {
	  if (m0 > -1) {
	    m1 = m0;
	    pid1 = pid0;
	    pT1 = pT0;
	  }
	
	  pT0 = pT;
	  m0 = i;
	  pid0 = TMath::Abs(pG->fID);
	}  else {
	  if (pT > pT1) {
	    pT1 = pT;
	    m1 = i;
	    pid1 = TMath::Abs(pG->fID);
	  }
	}
      }

      if ( fD0 != fD1  && pid0 == pid1 && strcmp("qcd", fChannel) ) {  
	
	//  cout << " **** SAME PARTICLE ****" << endl;
	
	if ( pid1 == fD0 ) { 
	  
	  pidTmp = fD1; 
	}
	if ( pid1 == fD1 ) {
	  
	  pidTmp = fD0; 
	}
	
	m1 = -1;
	pid1 = -1;
	pT1 = 0.;
	cnt = 1;
	
	for (int i = 0; i < fpEvt->nGenCands(); ++i) {
	  pG = fpEvt->getGenCand(i);
	  if (TMath::Abs(pG->fID) == pidTmp ) {
	    
	    pT = pG->fP.Pt();
	    eta = pG->fP.Eta();
	    if (pT > 3. && eta > -2.4 && eta < 2.4) {
	      ++cnt;
	    }	
	    // -- Get the second highest pT candidate for m1
	    if (pT > pT1) {
	      pT1 = pT;
	      m1 = i;
	      pid1 = TMath::Abs(pG->fID);
	      
	    }
	    
	  }
	  
	}
      }
    } else {    
      
      if ( TMath::Abs(pG->fID) == 13 ) {
	pM   = fpEvt->getGenCand(pG->fMom1);
	
	if ( TMath::Abs(pM->fID) == 443 ) {
	  pGM  = fpEvt->getGenCand(pM->fMom1);
	
	  if (TMath::Abs(pGM->fID) == 521 ) {
	
	    pT = pG->fP.Pt();
	    eta = pG->fP.Eta();
	    
	    if ( (pT > fMinPt) && (TMath::Abs(eta) < fMaxEta) ) {
	      ++cnt;
	    }
	    
	    // -- Get the highest pT muons sorted into m0, m1
	    if (pT > pT0) {
	      if (m0 > -1) {
		m1 = m0;
		pid1 = pid0;
		pT1 = pT0;
	      }
	      
	      pT0 = pT;
	      m0 = i;
	      pid0 = TMath::Abs(pG->fID);
	    }  else {
	      if (pT > pT1) {
		pT1 = pT;
		m1 = i;
		pid1 = TMath::Abs(pG->fID);
	      }   
	    }
	  }
	}
      }
      if ( TMath::Abs(pG->fID) == 321) {	
	pM   = fpEvt->getGenCand(pG->fMom1);
	
	if (TMath::Abs(pM->fID) == 521 ) {

	  pT = pG->fP.Pt();
	  eta = pG->fP.Eta();
	  
	  if ( (pT > fMinPt) && (TMath::Abs(eta) < fMaxEta) ) {
	    ++kcnt;
	  }
	  
	  // -- Get the highest pT muons sorted into m0, m1
	  if (pT > pT2) {
	    
	    pT2 = pT;
	    m2 = i;
	    pid2 = TMath::Abs(pG->fID);
	  }
	}
      }
    }
  }

  if ( m0 > 0 && m1 > 0 ) {

    TGenCand *gm0 = fpEvt->getGenCand(m0);
    TGenCand *gm1 = fpEvt->getGenCand(m1);
    
    double dphi = gm0->fP.DeltaPhi(gm1->fP);
    double deta = gm0->fP.Eta() - gm1->fP.Eta();
    double gRMM = TMath::Sqrt(dphi*dphi + deta*deta);

    ((TH1D*)gDirectory->Get("g150"))->Fill(gRMM);
    
    double gpT  = gm0->fP.Pt();
    ((TH1D*)gDirectory->Get("g100"))->Fill(gpT);
    if ( m2 > 0 ) { ((TH1D*)gDirectory->Get("g115"))->Fill(gpT); }
    gpT  = gm1->fP.Pt();
    ((TH1D*)gDirectory->Get("g100"))->Fill(gpT);
    if ( m2 > 0 ) { ((TH1D*)gDirectory->Get("g116"))->Fill(gpT); }
    
    double geta  = gm0->fP.Eta();
    ((TH1D*)gDirectory->Get("g102"))->Fill(geta);
    geta  = gm1->fP.Eta();
    ((TH1D*)gDirectory->Get("g102"))->Fill(geta);
    
    TLorentzVector gmm  = gm0->fP + gm1->fP;
    double gmass = gmm.M();
    ((TH1D*)gDirectory->Get("g103"))->Fill(gmass);  

    if ( m2 > 0 ) {

      TGenCand *gm0 = fpEvt->getGenCand(m0);
      TGenCand *gm1 = fpEvt->getGenCand(m1);
      TGenCand *gm2 = fpEvt->getGenCand(m2);
      
      TLorentzVector gjpsi  = gm0->fP + gm1->fP;
      
      double dphi = gjpsi.DeltaPhi(gm2->fP);
      double deta = gjpsi.Eta() - gm2->fP.Eta();
      double gRKJ = TMath::Sqrt(dphi*dphi + deta*deta);
      
      dphi = gjpsi.DeltaPhi(gm0->fP);
      deta = gjpsi.Eta() - gm0->fP.Eta();
      double gRMJ1 = TMath::Sqrt(dphi*dphi + deta*deta);
      
      dphi = gjpsi.DeltaPhi(gm1->fP);
      deta = gjpsi.Eta() - gm1->fP.Eta();
      double gRMJ2 = TMath::Sqrt(dphi*dphi + deta*deta);

      ((TH1D*)gDirectory->Get("g151"))->Fill(gRKJ);
      ((TH1D*)gDirectory->Get("g152"))->Fill(gRMJ1);
      ((TH1D*)gDirectory->Get("g153"))->Fill(gRMJ2);
      
      double gpT  = gm2->fP.Pt();
      ((TH1D*)gDirectory->Get("g110"))->Fill(gpT);
      
      double geta  = gm2->fP.Eta();
      ((TH1D*)gDirectory->Get("g112"))->Fill(geta);
      
      TLorentzVector gkj  = gm0->fP + gm1->fP + gm2->fP;
      double gmass = gkj.M();
      ((TH1D*)gDirectory->Get("g113"))->Fill(gmass);
    }
  }


     
  if (cnt >= 2 && fD2 < 0 ) {
    fGoodKinematics = 1;
    fER1->Fill(o+1.1); 

  } else if (cnt >= 2 && kcnt >=1 && fD2 > 0 ) {
    fGoodKinematics = 1;
    fER1->Fill(o+1.1); 

  } else {
    fGoodKinematics = 0;
    fER1->Fill(o+2.1);
  }
  
}

// ----------------------------------------------------------------------
void treeBmm::muonMisIdRate() { 

  double recMu(-1);
  int mcID(-1);
  double pT(0.), eta(0.);

  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {

    recMu = fpEvt->getRecTrack(it)->fMuID;
    mcID  = fpEvt->getRecTrack(it)->fMCID;

    pT  = fpEvt->getRecTrack(it)->fPlab.Pt();
    eta = fpEvt->getRecTrack(it)->fPlab.Eta();
    
    ((TH1D*)gDirectory->Get("m000"))->Fill(mcID);

    // -- Pions
    if ( TMath::Abs(mcID) == 211 ) { 

      if ( pT > 3. ) { 
	
	fMisID->Fill(0.1);
	
	((TH1D*)gDirectory->Get("m100"))->Fill(pT);   // all Pions
	((TH1D*)gDirectory->Get("m110"))->Fill(eta);
	((TH2D*)gDirectory->Get("M100"))->Fill(pT, eta);
	
	if ( recMu > 0.5 ) {
	  
	  fMisID->Fill(1.1);
	  
	  ((TH1D*)gDirectory->Get("m101"))->Fill(pT);   // mis-id. Pions
	  ((TH1D*)gDirectory->Get("m111"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M101"))->Fill(pT, eta);
	  
	}
	
	if ( recMu < 0.5 ) {
	  
	  fMisID->Fill(2.1);
	  
	  ((TH1D*)gDirectory->Get("m102"))->Fill(pT);   // correct-id. Pions
	  ((TH1D*)gDirectory->Get("m112"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M102"))->Fill(pT, eta);
	  
	}
      }
    }    

    // -- Kaons
    if ( TMath::Abs(mcID) == 321 ) { 

      if ( pT > 3. ) { 
	
	fMisID->Fill(10.1);

      ((TH1D*)gDirectory->Get("m200"))->Fill(pT);      // all Kaons
      ((TH1D*)gDirectory->Get("m210"))->Fill(eta);
      ((TH2D*)gDirectory->Get("M200"))->Fill(pT, eta);

	if ( recMu > 0.5 ) {

	  fMisID->Fill(11.1);
 
	  ((TH1D*)gDirectory->Get("m201"))->Fill(pT);   // mis-id. Kaons
	  ((TH1D*)gDirectory->Get("m211"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M201"))->Fill(pT, eta);

	}

	if ( recMu < 0.5 ) {

	  fMisID->Fill(12.1);

	  ((TH1D*)gDirectory->Get("m202"))->Fill(pT);   // correct-id. Kaons
	  ((TH1D*)gDirectory->Get("m212"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M202"))->Fill(pT, eta);

	}
      }
    }    

    // -- Protons
    if ( TMath::Abs(mcID) == 2212 ) { 

      if ( pT > 3. ) { 

	fMisID->Fill(20.1);

	((TH1D*)gDirectory->Get("m300"))->Fill(pT);      // all Protons
	((TH1D*)gDirectory->Get("m310"))->Fill(eta);
	((TH2D*)gDirectory->Get("M300"))->Fill(pT, eta);
	
	if ( recMu > 0.5 ) {

	  fMisID->Fill(21.1);
 
	  ((TH1D*)gDirectory->Get("m301"))->Fill(pT);   // mis-id. Protons
	  ((TH1D*)gDirectory->Get("m311"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M301"))->Fill(pT, eta);

	}

	if ( recMu < 0.5 ) {

	  fMisID->Fill(22.1);

	  ((TH1D*)gDirectory->Get("m302"))->Fill(pT);   // correct-id. Protons
	  ((TH1D*)gDirectory->Get("m312"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M302"))->Fill(pT, eta);

	}
      }
    }    

    // -- Muons
    if ( TMath::Abs(mcID) == 13 ) { 

      if ( pT > 3. ) { 

	fMisID->Fill(30.1);
	
	((TH1D*)gDirectory->Get("m400"))->Fill(pT);      // all Muons
	((TH1D*)gDirectory->Get("m410"))->Fill(eta);
	((TH2D*)gDirectory->Get("M400"))->Fill(pT, eta);
	
	if ( recMu > 0.5 ) {
	  
	  fMisID->Fill(31.1);
	  
	  ((TH1D*)gDirectory->Get("m401"))->Fill(pT);   // correct-id. Muons
	  ((TH1D*)gDirectory->Get("m411"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M401"))->Fill(pT, eta);
	  
	}

	if ( recMu < 0.5 ) {

	  fMisID->Fill(32.1);

	  ((TH1D*)gDirectory->Get("m402"))->Fill(pT);   // mis-id. Muons
	  ((TH1D*)gDirectory->Get("m412"))->Fill(eta);
	  ((TH2D*)gDirectory->Get("M402"))->Fill(pT, eta);

	}
      }
    }
  }
}

// ---------------------------------------------------------------------
void treeBmm::L1Selection(int o) {
  // offset: 100
  fGoodL1Trigger = 1;
  
  fER1->Fill(o+0.1);
  int L1 = fpEvt->fDiMuonTriggerDecision;

  if ( fD0 != 13 || fD1 != 13 ) { L1 = 1; }


  if (!L1) {
    fGoodL1Trigger = 0; 
    fER1->Fill(o+12.1);
  } else {
    fER1->Fill(o+11.1);
  }

  if (1 == fGoodKinematics) {
    if (0 == fGoodL1Trigger) {
      fER1->Fill(o+2.1);
    } else {
      fER1->Fill(o+1.1);
    }
  }
}


// ----------------------------------------------------------------------
void treeBmm::cutJunk(int o) {
  // offset: 200
  fGoodEvent = 1;
 
  fER1->Fill(o+0.1);
  if (fSTI[0] < 0) {
    if (fDebug & 2) cout << "treeBmm::cutJunk> no index for signal track 0 ";
    fGoodEvent = 0;
    fER1->Fill(o+24.1);
  } else {
    fER1->Fill(o+23.1);

  }

  if (fSTI[1] < 0) {
    if (fDebug & 2)    cout << "treeBmm::cutJunk> no index for signal track 1 " ;
    fGoodEvent = 0;
    fER1->Fill(o+26.1);
  } else {
    fER1->Fill(o+25.1);
  }

  if (fD2 > 0 && fSTI[2] < 0) {
    if (fDebug & 2)    cout << "treeBmm::cutJunk> no index for signal track 2 " ;
    fGoodEvent = 0;
    fER1->Fill(o+28.1);
  } else {
    fER1->Fill(o+27.1);
  }

  // -- Study primary vertex
  if (-99 == fpEvt->fPrimaryVertex2.fType) {
    if (fDebug & 2)    cout << "treeBmm::cutJunk> no primary vertex " ;
    fGoodEvent = 0;
    fER1->Fill(o+30.1);
  } else {
    fER1->Fill(o+29.1);
  }

  if (0 == fGoodEvent) { 
    if (fDebug & 2)    cout << endl;
    fER1->Fill(o+12.1);
  } else {
    fER1->Fill(o+11.1);
  }

  if (1 == fGoodKinematics) {
    if (0 == fGoodEvent) { 
      fER1->Fill(o+2.1);
    } else {
      fER1->Fill(o+1.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger)) {
    if (0 == fGoodEvent) { 
      fER1->Fill(o+4.1);
    } else {
      fER1->Fill(o+3.1);
    }
  }

}

// ----------------------------------------------------------------------
void treeBmm::trackSelection(int o) {
  // offset: 400

  fpHistFile->cd();
  TH1D *h;

  fGoodTT = 1;

  for (int it = 0; it < 2; ++it) {
    TAnaTrack *pT = fpEvt->getRecTrack(fSTI[it]);
    double x(0.);
    
    fER1->Fill(o + 20+ it*10. + 0.1);

    x = pT->fPlab.Pt();

    h = (TH1D*)gDirectory->Get("t100"); h->Fill(x);
    h = (TH1D*)gDirectory->Get(Form("t10%d", it+1)); h->Fill(x);
    if (x < PTLO) fGoodTT = -1;
    if (x > PTHI) fGoodTT = -2;
    if (fGoodTT > 0) fER1->Fill(o + 20 + it*10. + 1.1);
    if ((fGoodTT < 1) && (fDebug & 0x8)) {
      break;
    }

    x = pT->fPlab.Eta();
    h = (TH1D*)gDirectory->Get("t110"); h->Fill(x);
    h = (TH1D*)gDirectory->Get(Form("t11%d", it+1)); h->Fill(x);
    if (x < ETALO) fGoodTT = -3;
    if (x > ETAHI) fGoodTT = -4;
    if (fGoodTT > 0) fER1->Fill(o + 20 + it*10. + 2.1);
    if ((fGoodTT < 1) && (fDebug & 0x8)) {
      break;
    }

    x = pT->fTip;
    h = (TH1D*)gDirectory->Get("t120"); h->Fill(x);
    h = (TH1D*)gDirectory->Get(Form("t12%d", it+1)); h->Fill(x);
    if (x > TIPHI) fGoodTT = -5;
    if (fGoodTT > 0) fER1->Fill(o + 20 + it*10. + 3.1);
    if ((fGoodTT < 1) && (fDebug & 0x8)) {
      break;
    }

  }

  fER1->Fill(o + 10.1);
  if (fGoodTT > 0) {
    fER1->Fill(o + 11.1);
  } else {
    fER1->Fill(o + 12.1);
  }

  if (1 == fGoodKinematics) {
    if (fGoodTT > 0) {
      fER1->Fill(o + 1.1);
    } else {
      fER1->Fill(o + 2.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger)) {
    if (fGoodTT > 0) { 
      fER1->Fill(o + 3.1);
    } else {
      fER1->Fill(o + 4.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger) &&  (1 == fGoodEvent)) {
    if (fGoodTT > 0) { 
      fER1->Fill(o + 5.1);
    } else {
      fER1->Fill(o + 6.1);
    }
  }

  // --- HLT Trigger not checked yet ?????!!!!!!
  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger) &&  (1 == fGoodEvent) &&  (1 == fGoodHLTTrigger)) {
    if (fGoodTT > 0) { 
      fER1->Fill(o+7.1);
    } else {
      fER1->Fill(o+8.1);
    }
  }

}

// ----------------------------------------------------------------------
void treeBmm::pidSelection(int o) {
  // offset: 500

  fGoodPid = 0;

  fGenMuIDL0 = 0;
  fGenMuIDL1 = 0;

  // -- fMuID = 1 if recontructed as muons, 0 else
  fMuIDL0 = fpEvt->getRecTrack(fSTI[0])->fMuID;
  fMuIDL1 = fpEvt->getRecTrack(fSTI[1])->fMuID;

  int genID0 = fpEvt->getRecTrack(fSTI[0])->fMCID;
  int genID1 = fpEvt->getRecTrack(fSTI[1])->fMCID;

  fGenMuIDL0 = TMath::Abs(genID0);
  fGenMuIDL1 = TMath::Abs(genID1);
  
  ((TH1D*)gDirectory->Get("mid"))->Fill(fMuIDL0);
  ((TH1D*)gDirectory->Get("mid"))->Fill(fMuIDL1);

  if ( fD2 > 0 ) { 
   
    fMuIDK = fpEvt->getRecTrack(fSTI[2])->fMuID;
    int genIDK = fpEvt->getRecTrack(fSTI[2])->fMCID;
    fGenMuIDK = TMath::Abs(genIDK);
  }

  fER1->Fill(o + 20.1);

  if (fMuIDL0 > 0.) {
    fER1->Fill(o + 21.1);
  } else {
    fER1->Fill(o + 22.1);
  }

  if (fMuIDL1 > 0.) {
    fER1->Fill(o + 23.1);
  } else {
    fER1->Fill(o + 24.1);
  }  

  if (fMuIDL0 > 0. && fMuIDL1 > 0.) {
    fGoodPid = 1;
    fER1->Fill(o + 25.1);
  } else {
    fER1->Fill(o + 26.1);
  }  

  if (1 == fGoodKinematics) {
    if (fGoodPid > 0) {
      fER1->Fill(o + 1.1);
    } else {
      fER1->Fill(o + 2.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger)) {
    if (fGoodPid > 0) { 
      fER1->Fill(o + 3.1);
    } else {
      fER1->Fill(o + 4.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger) &&  (1 == fGoodEvent)) {
    if (fGoodPid > 0) { 
      fER1->Fill(o + 5.1);
    } else {
      fER1->Fill(o + 6.1);
    }
  }

  // --- HLT Trigger not checked yet ?????!!!!!!
  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger) &&  (1 == fGoodEvent) &&  (1 == fGoodHLTTrigger)) {
    if (fGoodPid > 0) { 
      fER1->Fill(o+7.1);
    } else {
      fER1->Fill(o+8.1);
    }
  }

  // --- HLT Trigger not checked yet ?????!!!!!!
  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger) &&  (1 == fGoodEvent) &&  (1 == fGoodHLTTrigger) && (1 == fGoodTT)) {
    if (fGoodPid > 0) { 
      fER1->Fill(o+9.1);
    } else {
      fER1->Fill(o+10.1);
    }
  }

  if (1 == fGoodTT) {
    fER1->Fill(o+11.1);
    if (fGoodPid > 0) { 
      fER1->Fill(o+12.1);
    } else {
      fER1->Fill(o+13.1);
    }
  }


}

// ----------------------------------------------------------------------
void treeBmm::candidateSelection(int offset) {
  // offset: 600
  fGoodCand  = 1;

  TH1D *h;

  fChi2 = fpB->fVtx.fChi2;
  fNdof = fpB->fVtx.fNdof;
  fProb = fpB->fVtx.fProb;

  fL3d  = fpB->fVtx.fD3d;
  fS3d  = fpB->fVtx.fD3dE;
  fLxy  = fpB->fVtx.fDxy;
  fSxy  = fpB->fVtx.fDxyE;

  fMass = fpB->fMass;
  if (fChainFileName.Contains("sg-001")) {
    fMass *= 5.374/5.29;
  }



  fPt   = fpB->fPlab.Pt();
  fP    = fpB->fPlab.Mag();
  fEta  = fpB->fPlab.Eta();
  fTau  = fL3d*fMass/fP;
  fTxy  = fLxy*fMass/fPt;

  // -- Opening cone of signal muons
  double dphi = fpL1->fPlab.DeltaPhi(fpL2->fPlab);
  double deta = fpL1->fPlab.Eta() - fpL2->fPlab.Eta();
  
  fRMM  = TMath::Sqrt(dphi*dphi + deta*deta);

  // -- 3D version: l is the vector from the PV to the SV
  TVector3 l = fpB->fVtx.fPoint - fpEvt->fPrimaryVertex2.fPoint;
  fCosAngle3  = TMath::Cos(fpB->fPlab.Angle(l));

  // -- 2D version in transverse plane
  TVector2 tPV = fpEvt->fPrimaryVertex2.fPoint.XYvector();
  TVector2 tSV = fpB->fVtx.fPoint.XYvector();
  TVector2 tB  = fpB->fPlab.XYvector();
  TVector2 tl  = tSV - tPV;
  fCosAngle = TMath::Cos(tB.DeltaPhi(tl));
  
  
  // -- Compute isolation with tracks
  TGenCand *pG, *pM;
  TVector3 ptv;
  int gIndex(-1);
  double muonID(-1);
  double genPt(0.), genEta(0.),trkPt(0.), res0(0.), res1(0.);
  double df2, de2;
  double dr(0.), drmin(99.);
  double sum(0.), sum05(0.), sum06(0.), sum08(0.), sum10(0.), sum12(0.), sum14(0.);
  double ptmin(0.9);

  double cone = 0.5*fRMM + 0.4;
  
  fIsoVeto = fIsoVeto05 = fIsoVeto06 = fIsoVeto08 = fIsoVeto10 = fIsoVeto12 = fIsoVeto14 = 0;
  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {

    ptv   = fpEvt->getRecTrack(it)->fPlab;
    trkPt = ptv.Pt();
    
    gIndex = fpEvt->getRecTrack(it)->fGenIndex;
    muonID = fpEvt->getRecTrack(it)->fMuID;

    if ( gIndex > 0 ) {

      pG      = fpEvt->getGenCand(gIndex);
      genPt   = pG->fP.Pt();
      genEta  = TMath::Abs(pG->fP.Eta());

      res0    = trkPt - genPt;
      res1    = (trkPt - genPt)/genPt;


      //---------------------------Resolution Plots------------------------------
 
      // pT Resolution of GlobalMuons
      if ( muonID > 0.5 ) {

	((TProfile*)gDirectory->Get("p120"))->Fill(genPt, res0, 1);
	((TProfile*)gDirectory->Get("p121"))->Fill(genPt, res1, 1);
	((TH1D*)gDirectory->Get("p120C"))->Fill(genPt);

	((TH2D*)gDirectory->Get("D120"))->Fill(genPt, res0);
	((TH2D*)gDirectory->Get("D121"))->Fill(genPt, res1);
	
	if ( genPt > 8. && genPt < 12. ) {
 
	  ((TProfile*)gDirectory->Get("p130"))->Fill(genEta, res0, 1);
	  ((TProfile*)gDirectory->Get("p131"))->Fill(genEta, res1, 1);
	  ((TH1D*)gDirectory->Get("p130C"))->Fill(genEta);
 
	  ((TH2D*)gDirectory->Get("D130"))->Fill(genEta, res0);
	  ((TH2D*)gDirectory->Get("D131"))->Fill(genEta, res1);
	}
      }
      

      // pT Resolution of all tracks
      ((TProfile*)gDirectory->Get("p122"))->Fill(genPt, res0, 1);
      ((TProfile*)gDirectory->Get("p123"))->Fill(genPt, res1, 1);
      ((TH1D*)gDirectory->Get("p122C"))->Fill(genPt);

	((TH2D*)gDirectory->Get("D122"))->Fill(genPt, res0);
	((TH2D*)gDirectory->Get("D123"))->Fill(genPt, res1);
      
      if ( genPt > 8. && genPt < 12. ) {

	((TProfile*)gDirectory->Get("p132"))->Fill(genEta, res0, 1); 
	((TProfile*)gDirectory->Get("p133"))->Fill(genEta, res1, 1); 
	((TH1D*)gDirectory->Get("p132C"))->Fill(genEta);

	((TH2D*)gDirectory->Get("D132"))->Fill(genEta, res0); 
	((TH2D*)gDirectory->Get("D133"))->Fill(genEta, res1); 
      }

      //---------------------------Resolution Plots------------------------------

      pM     = fpEvt->getGenCand(pG->fMom1);

      if (it == fSTI[0]) {
	
	((TH1D*)gDirectory->Get("d114"))->Fill(pM->fID);
      }
      if (it == fSTI[1]) {
	
	((TH1D*)gDirectory->Get("d115"))->Fill(pM->fID);
      }
    }

    if (it == fSTI[0]) { continue; }
    if (it == fSTI[1]) { continue; }

    df2 = fpB->fPlab.DeltaPhi(ptv);
    ((TH1D*)gDirectory->Get("d110"))->Fill(df2);
    de2 = (ptv.Eta() - fpB->fPlab.Eta());

    dr  = TMath::Sqrt(df2*df2 + de2*de2);

    ((TH1D*)gDirectory->Get("d112"))->Fill(ptv.Pt());

    if ((dr < cone) && (ptv.Pt() > ptmin)) {
      sum += ptv.Pt();
      ((TH1D*)gDirectory->Get("d113"))->Fill(ptv.Pt());
    }

    // -- Systematics
    double lostTrack = (gRandom->Rndm() < 0.01? 0. : 1.);
    lostTrack = 1.0; // this is equivalent to NO smearing
    
    if ((dr < 0.5) && (ptv.Pt() > ptmin)) { sum05 += lostTrack*ptv.Pt(); ++fIsoVeto05; }
    if ((dr < 0.6) && (ptv.Pt() > ptmin)) { sum06 += lostTrack*ptv.Pt(); ++fIsoVeto06; }
    if ((dr < 0.8) && (ptv.Pt() > ptmin)) { sum08 += lostTrack*ptv.Pt(); ++fIsoVeto08; }
    if ((dr < 1.0) && (ptv.Pt() > ptmin)) { sum10 += lostTrack*ptv.Pt(); ++fIsoVeto10; }
    if ((dr < 1.2) && (ptv.Pt() > ptmin)) { sum12 += lostTrack*ptv.Pt(); ++fIsoVeto12; }
    if ((dr < 1.4) && (ptv.Pt() > ptmin)) { sum14 += lostTrack*ptv.Pt(); ++fIsoVeto14; }

    if ((ptv.Pt() > ptmin) && (dr < cone)) {
      ++fIsoVeto;
      //      cout << "  ... isolation veto triggered" << endl;
    }

    if ((ptv.Pt() > ptmin) && (dr < drmin)) {
      drmin = dr;
    }

  }

  fIso05 = fPt/(fPt + sum05);
  fIso06 = fPt/(fPt + sum06);
  fIso08 = fPt/(fPt + sum08);
  fIso10 = fPt/(fPt + sum10);
  fIso12 = fPt/(fPt + sum12);
  fIso14 = fPt/(fPt + sum14);

  fIso = fPt/(fPt + sum);

  // -- Choose cone size for isolation 
  if (TMath::Abs(ISOCONE - 0.5) < 0.001) {
    fIso = fIso05;
  } else if (TMath::Abs(ISOCONE - 0.6) < 0.001) {
    fIso = fIso06;
  } else if (TMath::Abs(ISOCONE - 0.8) < 0.001) {
    fIso = fIso08;
  } else if (TMath::Abs(ISOCONE - 1.0) < 0.001) {
    fIso = fIso10;
  } else if (TMath::Abs(ISOCONE - 1.2) < 0.001) {
    fIso = fIso12;
  } else if (TMath::Abs(ISOCONE - 1.4) < 0.001) {
    fIso = fIso14;
  } else {
    cout << "don't know about cone size " << ISOCONE << endl;
    exit(1);
  }


  // -- Systematics: Vertexing 
  if (0) {
    double smear, sigma;
    sigma = 0.3*fLxy;
    
    smear = gRandom->Gaus(0., sigma);
    fLxy += smear;
  }



  ((TH2D*)gDirectory->Get("D100"))->Fill(fMass, fRMM);
  ((TH1D*)gDirectory->Get("d106"))->Fill(drmin);

  if (fIsoVeto == 0) {
    ((TH1D*)gDirectory->Get("d107"))->Fill(drmin);
    if (fMass > 4.) {
      ((TH1D*)gDirectory->Get("d108"))->Fill(drmin);
    }
    if (fMass > 2.8 && fMass < 3.2) {
      ((TH1D*)gDirectory->Get("d109"))->Fill(drmin);
    }

    ((TH2D*)gDirectory->Get("D101"))->Fill(fMass, fRMM);

    ((TH2D*)gDirectory->Get("D102"))->Fill(cone, drmin);

  }



  // -- Fill histograms
  h = (TH1D*)gDirectory->Get("v100"); h->Fill(fChi2);
  h = (TH1D*)gDirectory->Get("v101"); h->Fill(fProb);

  h = (TH1D*)gDirectory->Get("v110"); h->Fill(fLxy);
  h = (TH1D*)gDirectory->Get("v111"); h->Fill(fSxy);
  h = (TH1D*)gDirectory->Get("v112"); h->Fill(fLxy/fSxy);

  h = (TH1D*)gDirectory->Get("v120"); h->Fill(fL3d);
  h = (TH1D*)gDirectory->Get("v121"); h->Fill(fS3d);
  h = (TH1D*)gDirectory->Get("v122"); h->Fill(fL3d/fS3d);

  h = (TH1D*)gDirectory->Get("b100"); h->Fill(fMass);
  h = (TH1D*)gDirectory->Get("b101"); h->Fill(fMass);
  h = (TH1D*)gDirectory->Get("b102"); h->Fill(fPt);
  h = (TH1D*)gDirectory->Get("b103"); h->Fill(fEta);
  h = (TH1D*)gDirectory->Get("b104"); h->Fill(fTau);

  // -- Angle between B momentum and sv direction
  h = (TH1D*)gDirectory->Get("b110"); h->Fill(fCosAngle);
  if (fPt > 5) {
    h = (TH1D*)gDirectory->Get("b120"); h->Fill(fCosAngle);
  }

  if (fL3d/fS3d > 4.) {
    h = (TH1D*)gDirectory->Get("b130"); h->Fill(fCosAngle);
  }
  if (fL3d/fS3d > 10.) {
    h = (TH1D*)gDirectory->Get("b131"); h->Fill(fCosAngle);
  }


  if (fLxy/fSxy > 4.) {
    h = (TH1D*)gDirectory->Get("b140"); h->Fill(fCosAngle);
  }
  if (fLxy/fSxy > 10.) {
    h = (TH1D*)gDirectory->Get("b141"); h->Fill(fCosAngle);
  }

  h = (TH1D*)gDirectory->Get("b150"); h->Fill(fRMM);
  h = (TH1D*)gDirectory->Get("b200"); h->Fill(fIso);
  h = (TH1D*)gDirectory->Get("b201"); h->Fill(fIsoVeto);


  //   cout << "SV: " << fpB->fVtx.fPoint.X()
  //        << " mass = " << fMass
  //        << " pT =  "  << fPt
  //        << " eta =  " << fEta
  //        << " fl3d = " <<  fL3d
  //        << " fchi2 = " <<  fChi2
  //        << endl;

  // -- candidate selection for real
  ((TH1D*)gDirectory->Get("b1000"))->Fill(fMass);
  TH1D *he = (TH1D*)gDirectory->Get("eeff");
  if (fChi2 > VTXCHI) {
    fGoodCand = 0;
    he->Fill(21.1);
    return;
  } else {
    he->Fill(20.1);
    ((TH1D*)gDirectory->Get("b1001"))->Fill(fMass);
  }

  if (fLxy < LXYLO) {
    fGoodCand = 0;
    he->Fill(23.1);
    return;
  } else {
    he->Fill(22.1);
    ((TH1D*)gDirectory->Get("b1002"))->Fill(fMass);
  }


}

// ----------------------------------------------------------------------
void treeBmm::normSample(int offset) {
  // offset: 700

  fGoodCand  = 1;

  TH1D *h;

  // -- B+ Candidate
  fChi2 = fpB->fVtx.fChi2;
  fNdof = fpB->fVtx.fNdof;
  fProb = fpB->fVtx.fProb;

  fL3d  = fpB->fVtx.fD3d;
  fS3d  = fpB->fVtx.fD3dE;
  fLxy  = fpB->fVtx.fDxy;
  fSxy  = fpB->fVtx.fDxyE;

  fPt   = fpB->fPlab.Pt();
  fP    = fpB->fPlab.Mag();
  fEta  = fpB->fPlab.Eta();
  fTau  = fL3d*fMass/fP;
  fTxy  = fLxy*fMass/fPt;

  fMass = fpB->fMass;

  // -- Jpsi Candidate
  fPtJ   = fpJpsi->fPlab.Pt();
  fPJ    = fpJpsi->fPlab.Mag();
  fEtaJ  = fpJpsi->fPlab.Eta();

  //  fTauJ  = fL3d*fMass/fP;
  //  fTxyJ  = fLxy*fMass/fPt;
  fMassJ = fpJpsi->fMass;

  // -- Opening cone of signal muons    
  double dphi = fpL1->fPlab.DeltaPhi(fpL2->fPlab);
  double deta = fpL1->fPlab.Eta() - fpL2->fPlab.Eta();

  fRMM  = TMath::Sqrt(dphi*dphi + deta*deta);

  // -- 3D version: l is the vector from the PV to the SV
  TVector3 l = fpB->fVtx.fPoint - fpEvt->fPrimaryVertex2.fPoint;
  fCosAngle3  = TMath::Cos(fpB->fPlab.Angle(l));

  // -- 2D version in transverse plane
  TVector2 tPV = fpEvt->fPrimaryVertex2.fPoint.XYvector();
  TVector2 tSV = fpB->fVtx.fPoint.XYvector();
  TVector2 tB  = fpB->fPlab.XYvector();
  TVector2 tl  = tSV - tPV;
  fCosAngle = TMath::Cos(tB.DeltaPhi(tl));
  
  
  // -- Compute isolation with tracks
  TGenCand *pG, *pM, *pGM;
  TVector3 ptv;
  int gIndex(-1);
  double muonID(-1);
  double genPt(0.), genEta(0.),trkPt(0.), res0(0.), res1(0.);
  double df2, de2;
  double dr(0.), drmin(99.);
  double sum(0.), sum05(0.), sum06(0.), sum08(0.), sum10(0.), sum12(0.), sum14(0.);
  double ptmin(1.0);

  //  double cone = 0.5*fRMM + 0.4;
  double cone = 1.1;
  
  fIsoVeto = fIsoVeto05 = fIsoVeto06 = fIsoVeto08 = fIsoVeto10 = fIsoVeto12 = fIsoVeto14 = 0;
  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {

    ptv   = fpEvt->getRecTrack(it)->fPlab;
    trkPt = ptv.Pt();
    
    gIndex = fpEvt->getRecTrack(it)->fGenIndex;
    muonID = fpEvt->getRecTrack(it)->fMuID;
    
    if ( gIndex > 0 ) {

      pG     = fpEvt->getGenCand(gIndex);
      genPt  = pG->fP.Pt();
      genEta = TMath::Abs(pG->fP.Eta());

      res0   = trkPt - genPt;
      res1   = (trkPt - genPt)/genPt;


      //---------------------------Resolution Plots------------------------------

      // pT Resolution of GlobalMuons
      if ( muonID > 0.5 ) {

	((TProfile*)gDirectory->Get("p120"))->Fill(genPt, res0, 1);
	((TProfile*)gDirectory->Get("p121"))->Fill(genPt, res1, 1);
	((TH1D*)gDirectory->Get("p120C"))->Fill(genPt);

	((TH2D*)gDirectory->Get("D120"))->Fill(genPt, res0);
	((TH2D*)gDirectory->Get("D121"))->Fill(genPt, res1);
	
	if ( genPt > 8. && genPt < 12. ) {
 
	  ((TProfile*)gDirectory->Get("p130"))->Fill(genEta, res0, 1);
	  ((TProfile*)gDirectory->Get("p131"))->Fill(genEta, res1, 1);
	  ((TH1D*)gDirectory->Get("p130C"))->Fill(genEta);
 
	  ((TH2D*)gDirectory->Get("D130"))->Fill(genEta, res0);
	  ((TH2D*)gDirectory->Get("D131"))->Fill(genEta, res1);
	}
      }
      

      // pT Resolution of all tracks
      ((TProfile*)gDirectory->Get("p122"))->Fill(genPt, res0, 1);
      ((TProfile*)gDirectory->Get("p123"))->Fill(genPt, res1, 1);
      ((TH1D*)gDirectory->Get("p122C"))->Fill(genPt);

      ((TH2D*)gDirectory->Get("D122"))->Fill(genPt, res0);
      ((TH2D*)gDirectory->Get("D123"))->Fill(genPt, res1);
      
      if ( genPt > 8. && genPt < 12. ) {

	((TProfile*)gDirectory->Get("p132"))->Fill(genEta, res0, 1); 
	((TProfile*)gDirectory->Get("p133"))->Fill(genEta, res1, 1); 
	((TH1D*)gDirectory->Get("p132C"))->Fill(genEta);

	((TH2D*)gDirectory->Get("D132"))->Fill(genEta, res0); 
	((TH2D*)gDirectory->Get("D133"))->Fill(genEta, res1); 
      }

      //---------------------------Resolution Plots------------------------------


      pM     = fpEvt->getGenCand(pG->fMom1);
      pGM    = fpEvt->getGenCand(pM->fMom1);
      
      if (it == fSTI[0]) {
	
	((TH1D*)gDirectory->Get("d114"))->Fill(pM->fID);
	((TH1D*)gDirectory->Get("d117"))->Fill(pGM->fID);
      }
      if (it == fSTI[1]) {
	
	((TH1D*)gDirectory->Get("d115"))->Fill(pM->fID);
	((TH1D*)gDirectory->Get("d118"))->Fill(pGM->fID);
      }
      if (it == fSTI[2]) {
	
	((TH1D*)gDirectory->Get("d116"))->Fill(pM->fID);
	((TH1D*)gDirectory->Get("d119"))->Fill(pGM->fID);
      }
      
    }

    df2 = fpJpsi->fPlab.DeltaPhi(ptv);
    de2 = (ptv.Eta() - fpJpsi->fPlab.Eta());
    dr  = TMath::Sqrt(df2*df2 + de2*de2);
    
    if (it == fSTI[0]) { fRMJ1 = dr; continue; }
    if (it == fSTI[1]) { fRMJ2 = dr; continue; }
    if (it == fSTI[2]) { fRKJ  = dr; continue; }

    ((TH1D*)gDirectory->Get("d110"))->Fill(df2);
    ((TH1D*)gDirectory->Get("d112"))->Fill(ptv.Pt());

    if ((ptv.Pt() > ptmin) && (dr < cone)) {
      sum += ptv.Pt();
      ((TH1D*)gDirectory->Get("d113"))->Fill(ptv.Pt());
    }

    // -- Systematics
    double lostTrack = (gRandom->Rndm() < 0.01? 0. : 1.);
    lostTrack = ptmin; // this is equivalent to NO smearing
    
    if ((dr < 0.5) && (ptv.Pt() > ptmin)) { sum05 += lostTrack*ptv.Pt(); ++fIsoVeto05; }
    if ((dr < 0.6) && (ptv.Pt() > ptmin)) { sum06 += lostTrack*ptv.Pt(); ++fIsoVeto06; }
    if ((dr < 0.8) && (ptv.Pt() > ptmin)) { sum08 += lostTrack*ptv.Pt(); ++fIsoVeto08; }
    if ((dr < 1.0) && (ptv.Pt() > ptmin)) { sum10 += lostTrack*ptv.Pt(); ++fIsoVeto10; }
    if ((dr < 1.2) && (ptv.Pt() > ptmin)) { sum12 += lostTrack*ptv.Pt(); ++fIsoVeto12; }
    if ((dr < 1.4) && (ptv.Pt() > ptmin)) { sum14 += lostTrack*ptv.Pt(); ++fIsoVeto14; }
    

    if ((ptv.Pt() > ptmin) && (dr < cone)) {
      ++fIsoVeto;
      //      cout << "  ... isolation veto triggered" << endl;
    }

    if ((ptv.Pt() > ptmin) && (dr < drmin)) {
      drmin = dr;
    }

  }

  fIso05 = fPt/(fPt + sum05);
  fIso06 = fPt/(fPt + sum06);
  fIso08 = fPt/(fPt + sum08);
  fIso10 = fPt/(fPt + sum10);
  fIso12 = fPt/(fPt + sum12);
  fIso14 = fPt/(fPt + sum14);

  fIso = fPt/(fPt + sum);

  // -- Choose cone size for isolation 
  if (TMath::Abs(ISOCONE - 0.5) < 0.001) {
    fIso = fIso05;
  } else if (TMath::Abs(ISOCONE - 0.6) < 0.001) {
    fIso = fIso06;
  } else if (TMath::Abs(ISOCONE - 0.8) < 0.001) {
    fIso = fIso08;
  } else if (TMath::Abs(ISOCONE - 1.0) < 0.001) {
    fIso = fIso10;
  } else if (TMath::Abs(ISOCONE - 1.2) < 0.001) {
    fIso = fIso12;
  } else if (TMath::Abs(ISOCONE - 1.4) < 0.001) {
    fIso = fIso14;
  } else {
    cout << "don't know about cone size " << ISOCONE << endl;
    exit(1);
  }


  // -- Systematics: Vertexing 
  if (0) {
    double smear, sigma;
    sigma = 0.3*fLxy;
    
    smear = gRandom->Gaus(0., sigma);
    fLxy += smear;
  }



  ((TH2D*)gDirectory->Get("D100"))->Fill(fMass, fRMM);
  ((TH1D*)gDirectory->Get("d106"))->Fill(drmin);

  if (fIsoVeto == 0) {
    ((TH1D*)gDirectory->Get("d107"))->Fill(drmin);
    if (fMass > 4.) {
      ((TH1D*)gDirectory->Get("d108"))->Fill(drmin);
    }
    if (fMass > 2.8 && fMass < 3.2) {
      ((TH1D*)gDirectory->Get("d109"))->Fill(drmin);
    }

    ((TH2D*)gDirectory->Get("D101"))->Fill(fMass, fRMM);

    ((TH2D*)gDirectory->Get("D102"))->Fill(cone, drmin);

  }



  // -- Fill histograms
  h = (TH1D*)gDirectory->Get("v100"); h->Fill(fChi2);
  h = (TH1D*)gDirectory->Get("v101"); h->Fill(fProb);

  h = (TH1D*)gDirectory->Get("v110"); h->Fill(fLxy);
  h = (TH1D*)gDirectory->Get("v111"); h->Fill(fSxy);
  h = (TH1D*)gDirectory->Get("v112"); h->Fill(fLxy/fSxy);

  h = (TH1D*)gDirectory->Get("v120"); h->Fill(fL3d);
  h = (TH1D*)gDirectory->Get("v121"); h->Fill(fS3d);
  h = (TH1D*)gDirectory->Get("v122"); h->Fill(fL3d/fS3d);

  h = (TH1D*)gDirectory->Get("b100"); h->Fill(fMass);
  h = (TH1D*)gDirectory->Get("b101"); h->Fill(fMass);
  h = (TH1D*)gDirectory->Get("b102"); h->Fill(fPt);
  h = (TH1D*)gDirectory->Get("b103"); h->Fill(fEta);
  h = (TH1D*)gDirectory->Get("b104"); h->Fill(fTau);

  h = (TH1D*)gDirectory->Get("j100"); h->Fill(fMassJ);
  h = (TH1D*)gDirectory->Get("j101"); h->Fill(fMassJ);
  h = (TH1D*)gDirectory->Get("j102"); h->Fill(fPtJ);
  h = (TH1D*)gDirectory->Get("j103"); h->Fill(fEtaJ);
  h = (TH1D*)gDirectory->Get("j104"); h->Fill(fTauJ);

  // -- Angle between B momentum and sv direction
  h = (TH1D*)gDirectory->Get("b110"); h->Fill(fCosAngle);
  if (fPt > 5) {
    h = (TH1D*)gDirectory->Get("b120"); h->Fill(fCosAngle);
  }

  if (fL3d/fS3d > 4.) {
    h = (TH1D*)gDirectory->Get("b130"); h->Fill(fCosAngle);
  }
  if (fL3d/fS3d > 10.) {
    h = (TH1D*)gDirectory->Get("b131"); h->Fill(fCosAngle);
  }


  if (fLxy/fSxy > 4.) {
    h = (TH1D*)gDirectory->Get("b140"); h->Fill(fCosAngle);
  }
  if (fLxy/fSxy > 10.) {
    h = (TH1D*)gDirectory->Get("b141"); h->Fill(fCosAngle);
  }

  h = (TH1D*)gDirectory->Get("b150"); h->Fill(fRMM);
  h = (TH1D*)gDirectory->Get("b200"); h->Fill(fIso);
  h = (TH1D*)gDirectory->Get("b201"); h->Fill(fIsoVeto);

  h = (TH1D*)gDirectory->Get("j150"); h->Fill(fRMM);
  h = (TH1D*)gDirectory->Get("j151"); h->Fill(fRKJ);
  h = (TH1D*)gDirectory->Get("j152"); h->Fill(fRMJ1);
  h = (TH1D*)gDirectory->Get("j153"); h->Fill(fRMJ2);


  //   cout << "SV: " << fpB->fVtx.fPoint.X()
  //        << " mass = " << fMass
  //        << " pT =  "  << fPt
  //        << " eta =  " << fEta
  //        << " fl3d = " <<  fL3d
  //        << " fchi2 = " <<  fChi2
  //        << endl;

  // -- candidate selection for real
  ((TH1D*)gDirectory->Get("b1000"))->Fill(fMass);
  TH1D *he = (TH1D*)gDirectory->Get("eeff");
  if (fChi2 > VTXCHI) {
    fGoodCand = 0;
    he->Fill(21.1);
    return;
  } else {
    he->Fill(20.1);
    ((TH1D*)gDirectory->Get("b1001"))->Fill(fMass);
  }

  if (fLxy < LXYLO) {
    fGoodCand = 0;
    he->Fill(23.1);
    return;
  } else {
    he->Fill(22.1);
    ((TH1D*)gDirectory->Get("b1002"))->Fill(fMass);
  }


}

// ----------------------------------------------------------------------
void treeBmm::kaonPlots() {

  TAnaCand   *pCand;  
  TAnaTrack  *pTrack;

  TGenCand *pG;
  TGenCand *pM;

  int gIndex(-1), gPDG(-1), mPDG(-1);
  int tmKaons(0), tmkIndex(-1);

  for (unsigned int i = 0; i < fKCI.size(); i++) {

    pCand = fpEvt->getCand(fKCI[i]);
    pTrack = fpEvt->getRecTrack(fKTI[i]);
    gIndex = pTrack->fGenIndex;

    if ( gIndex > 0 ) { 

      pG  = fpEvt->getGenCand(gIndex);
      gPDG = pG->fID;

      pM  = fpEvt->getGenCand(pG->fMom1);
      mPDG = pM->fID;

      if ( gPDG == 321 && mPDG == 521 ) {
	tmkIndex = i;
	tmKaons++;
      } 
    }
  }

  ((TH1D*)gDirectory->Get("keff"))->Fill(0.1 + tmKaons);

  if ( fSTI[2] >= 0 ) {  

    ((TH1D*)gDirectory->Get("keff"))->Fill(10.1);

    if ( tmkIndex >= 0) {

      if ( tmkIndex == fSTI[2] ) {
	((TH1D*)gDirectory->Get("keff"))->Fill(11.1);
      } else {
	((TH1D*)gDirectory->Get("keff"))->Fill(12.1);
      }

      pCand  = fpEvt->getCand(fKCI[tmkIndex]);
      pTrack = fpEvt->getRecTrack(fKTI[tmkIndex]);

      ((TH1D*)gDirectory->Get("k100"))->Fill(pCand->fMass);
      
      ((TH1D*)gDirectory->Get("k101"))->Fill(pCand->fVtx.fChi2);
      ((TH1D*)gDirectory->Get("k102"))->Fill(pCand->fPlab.Pt());
      ((TH1D*)gDirectory->Get("k103"))->Fill(pCand->fPlab.Eta());
      ((TH1D*)gDirectory->Get("k104"))->Fill(pTrack->fPlab.Pt());
      ((TH1D*)gDirectory->Get("k105"))->Fill(pTrack->fPlab.Eta());
      
      gIndex = fpK->fGenIndex;
      
      if ( gIndex > 0 ) { 
	
	pG  = fpEvt->getGenCand(gIndex); 
	gPDG = pG->fID; 
	
	pM  = fpEvt->getGenCand(pG->fMom1);
	mPDG = pM->fID;
	
	if ( tmkIndex == fSTI[2] ) {

	  ((TH2D*)gDirectory->Get("K100"))->Fill(gPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K101"))->Fill(gPDG, pCand->fMass  - fpB->fMass);
	  ((TH2D*)gDirectory->Get("K102"))->Fill(mPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K103"))->Fill(mPDG, pCand->fMass  - fpB->fMass);

	} else {

	  ((TH2D*)gDirectory->Get("K200"))->Fill(gPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K201"))->Fill(gPDG, pCand->fMass  - fpB->fMass);
	  ((TH2D*)gDirectory->Get("K202"))->Fill(mPDG, pCand->fVtx.fChi2  - fpB->fVtx.fChi2);
	  ((TH2D*)gDirectory->Get("K203"))->Fill(mPDG, pCand->fMass  - fpB->fMass);
	}

      }
    }
  }
}

// ----------------------------------------------------------------------
void treeBmm::HLTSelection(int o) {
  // offset: 300

  double HLTlPT(4.0);
  double HLTlETA(2.4);
  double HLTlTIP(10000.1);
  double HLTCHI2(20.);
  double HLTL3D(0.015);
  double HLTMASSLO(5.0);
  double HLTMASSHI(6.0);

  if ( strcmp("2mu", fChannel) ) {

    HLTMASSLO = 5.0;
    HLTMASSHI = 6.0;
  }


  fER1->Fill(o+0.1);
  if (1 == fGoodKinematics 
      && 1 == fGoodL1Trigger
      && (fPtL0 > HLTlPT) && (fPtL1 > HLTlPT)
      && (fEtaL0 > -HLTlETA) && (fEtaL0 < HLTlETA)
      && (fEtaL1 > -HLTlETA) && (fEtaL1 < HLTlETA)
      && (fQL0*fQL1 < -0.5)
      && (fTipL0 < HLTlTIP) && (fTipL1 < HLTlTIP)
      && (fChi2 < HLTCHI2)
      && (fL3d > HLTL3D)
      && (fMass > HLTMASSLO) && (fMass < HLTMASSHI)
      ) {
    fGoodHLTTrigger = 1;
    fER1->Fill(o+11.1);
  } else {
    fGoodHLTTrigger = 0;
    fER1->Fill(o+12.1);
  }

  if (1 == fGoodKinematics) {
    if (0 == fGoodHLTTrigger) { 
      fER1->Fill(o+2.1);
    } else {
      fER1->Fill(o+1.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger)) {
    if (0 == fGoodHLTTrigger) { 
      fER1->Fill(o+4.1);
    } else {
      fER1->Fill(o+3.1);
    }
  }

  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger) &&  (1 == fGoodEvent)) {
    if (0 == fGoodHLTTrigger) { 
      fER1->Fill(o+6.1);
    } else {
      fER1->Fill(o+5.1);
    }
  }


  if ((1 == fGoodKinematics) && (1 == fGoodL1Trigger)) {
    fER1->Fill(o+30.1);
    if ((fPtL0 > HLTlPT) && (fPtL1 > HLTlPT)) {
      fER1->Fill(o+31.1);
      if ((fEtaL0 > -HLTlETA) && (fEtaL0 < HLTlETA) && (fEtaL1 > -HLTlETA) && (fEtaL1 < HLTlETA)) {
	fER1->Fill(o+32.1);
	if (fQL0*fQL1 < -0.5) {
	  fER1->Fill(o+33.1);
	  if ((fTipL0 < HLTlTIP) && (fTipL1 < HLTlTIP)) {
	    fER1->Fill(o+34.1);
	    if (fChi2 < HLTCHI2) {
	      fER1->Fill(o+35.1);
	      if (fL3d > HLTL3D) {
		fER1->Fill(o+36.1);
		if ((fMass > HLTMASSLO) && (fMass < HLTMASSHI)) {
		  fER1->Fill(o+37.1);
		}      
	      }
	    }
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
void treeBmm::fillHist() {

  TH1D *ha = (TH1D*)gDirectory->Get("aeff");

  if (fpL1->fPlab.Pt() > PTLO) {
    ha->Fill(101.1);
  } else {
    ha->Fill(102.1);
  }

  if (fpL2->fPlab.Pt() > PTLO) {
    ha->Fill(103.1);
  } else {
    ha->Fill(104.1);
  }

  if (fPt > PTBS) {
    ha->Fill(105.1);
  } else {
    ha->Fill(106.1);
  }

  if ((fRMM > 0.4) && (fRMM < 1.2)) {
    ha->Fill(107.1);
  } else {
    ha->Fill(108.1);
  }

  if ( fGoodKinematics
       && (fPtL0 > 4.3)
       && (fPtL1 > 4.3)
       && (fPt > 12.)
       && (fRMM > 0.4) && (fRMM < 1.2)
       ) {

    ha->Fill(109.1);
  } else {
    ha->Fill(110.1);
  }

  // -- efficiency
  fillAnalysisEff();


}

// ----------------------------------------------------------------------
void treeBmm::vertexResolution() {


  TH1D *h;

  TVector3 resPV = fpEvt->fPrimaryVertex2.fSimPoint - fpEvt->fPrimaryVertex2.fPoint;

  fDpVtxX      = resPV.X(); 
  fDpVtxY      = resPV.Y(); 
  fDpVtxZ      = resPV.Z(); 
 
  // -------------------------------------------------------------------
  TGenCand *pM;
  TGenCand *pG; // B-meson
  TGenCand *pD1;
  TGenCand *pD2;

  TVector3 simSecVtx;
  TVector3 simPrimVtx;

  int simVtx(0);
  double genPt(0.), genPlab(0.), mass(0.);

  // -- Check for two muons from B-Meson
  int truth_bmm(0), truth_bjk(0), truth_jmm(0);

  if ( fpL1->fGenIndex > 0 && fpL2->fGenIndex > 0 ) {

    TGenCand *g1 = fpEvt->getGenCand(fpL1->fGenIndex);
    TGenCand *g2 = fpEvt->getGenCand(fpL2->fGenIndex);
    
    int l1 = g1->fID;
    int l2 = g2->fID;
  
    if ( TMath::Abs(l1) == 13 && TMath::Abs(l2) == 13 &&
	 g1->fMom1 >0 && g2->fMom1 > 0 ) {
    
      TGenCand *m1    = fpEvt->getGenCand(g1->fMom1);
      TGenCand *m2    = fpEvt->getGenCand(g2->fMom2);
    
      int b1 = m1->fID;
      int b2 = m2->fID;
    
      // -- Bs -> mu mu
      if ( fD2 < 0 && TMath::Abs(b1) == 531 && TMath::Abs(b2) == 531 ) {
      
	truth_bmm = 1;
      }

      // -- B+ -> J/Psi K+ 
      if ( fD2 > 0 && TMath::Abs(b1) == 443 && TMath::Abs(b2) == 443 &&
	   m1->fMom1 >0 && m2->fMom1 > 0 ) {
      
	TGenCand *gm1    = fpEvt->getGenCand(m1->fMom1);
	TGenCand *gm2    = fpEvt->getGenCand(m2->fMom2);
    
	int bj1 = gm1->fID;
	int bj2 = gm2->fID;
     
	if ( TMath::Abs(bj1) == 521 && TMath::Abs(bj2) == 521 ) {

	  truth_jmm = 1;

	  if ( fpK->fGenIndex > 0 ) {

	    TGenCand *g3 = fpEvt->getGenCand(fpK->fGenIndex);
	    
	    int kaon = g3->fID;
    
	    if ( TMath::Abs(kaon) == 321 && g3->fMom1 >0 ) {
	      
	      TGenCand *m3    = fpEvt->getGenCand(g3->fMom1);
	      
	      int b3 = m3->fID;
	      
	      if ( TMath::Abs(b3) == 521 ) {
		cout << "truth" << endl;
		truth_bjk = 1;
	      }
	    }
	  }
	}
      }
    }
  }


  // -- Bs -> mu mu
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    
    pG   = fpEvt->getGenCand(i);
    
    if ( truth_bmm && TMath::Abs(pG->fID) == 531 &&
	 pG->fMom1 >0 && pG->fDau1 > 0 && pG->fDau2) {
      
      pM    = fpEvt->getGenCand(pG->fMom1);
      pD1   = fpEvt->getGenCand(pG->fDau1);
      pD2   = fpEvt->getGenCand(pG->fDau2);
      
      if ( TMath::Abs(pG->fDau2 - pG->fDau1) == 1 &&
	   TMath::Abs(pD1->fID) == 13 &&  TMath::Abs(pD2->fID) == 13 ) {
	
	if ( fpL1->fPlab.Pt() > 3. && fpL2->fPlab.Pt() > 3. && fLxy/fSxy > 2. ) { 
	  
	  simPrimVtx.SetXYZ( 0.1*pM->fV.x(),
			     0.1*pM->fV.y(),
			     0.1*pM->fV.z() );
	  
	  simSecVtx.SetXYZ( 0.1*pG->fV.x(),
			    0.1*pG->fV.y(),
			    0.1*pG->fV.z() );
	  
	  genPt    = pG->fP.Perp();
	  genPlab  = pG->fP.P();
	  mass     = pG->fMass;
	  
	  simVtx++;
	}
      }
    }
    

    // -- B+ -> J/Psi K+ 
    if ( truth_jmm && TMath::Abs(pG->fID) == 521 &&
	 pG->fMom1 >0 && pG->fDau1 > 0 && pG->fDau2 ) {
      
      pM    = fpEvt->getGenCand(pG->fMom1);
      pD1   = fpEvt->getGenCand(pG->fDau1);
      pD2   = fpEvt->getGenCand(pG->fDau2);
      
      if ( TMath::Abs(pG->fDau2 - pG->fDau1) == 1 &&
	   (TMath::Abs(pD1->fID) == 321 &&  TMath::Abs(pD2->fID) == 443) || 
	   (TMath::Abs(pD1->fID) == 443 &&  TMath::Abs(pD2->fID) == 321) ) {
	
	if ( fpL1->fPlab.Pt() > 3. && fpL2->fPlab.Pt() > 3. && fpK->fPlab.Pt() > 3. && fLxy/fSxy > 2. ) { 
	  
	  simPrimVtx.SetXYZ( 0.1*pM->fV.x(),
			     0.1*pM->fV.y(),
			     0.1*pM->fV.z() );
	  
	  simSecVtx.SetXYZ( 0.1*pG->fV.x(),
			    0.1*pG->fV.y(),
			    0.1*pG->fV.z() );
	  
	  genPt    = pG->fP.Perp();
	  genPlab  = pG->fP.P();
	  mass     = pG->fMass;
	  
	  simVtx++;
	}
      }
    }
  }
  // -------------------------------------------------------------------

  // calcuate the resultion of the decay vertex (
  //**   TVector3 Res = fpB->fVtx.fSimPoint -  fpB->fVtx.fPoint; 
  TVector3 Res = simSecVtx -  fpB->fVtx.fPoint; 
  double theta = fpB->fPlab.Angle(Res);
 
  if (fDebug & 3) {

    cout << "SecVtx (gen) " << " " << simSecVtx.x() << " " << simSecVtx.y() << " " << simSecVtx.z() << endl;
    cout << "SecVtx (rec) " << " " << fpB->fVtx.fPoint.x() 
                            << " " << fpB->fVtx.fPoint.y() 
                            << " " << fpB->fVtx.fPoint.z() << endl;
  }

  fDvtxPerp    = Res.Mag()*TMath::Sin(theta);  
  fDvtxPar     = Res.Mag()*TMath::Cos(theta);

  // -- proper time resolution
  //**   TVector3 fltGen = fpB->fVtx.fSimPoint - fpEvt->fPrimaryVertex2.fSimPoint;
  //**   TVector3 fltGen = simSecVtx            - fpEvt->fPrimaryVertex2.fSimPoint;
  TVector3 fltGen = simSecVtx - simPrimVtx;
  TVector3 fltRec = fpB->fVtx.fPoint    - fpEvt->fPrimaryVertex2.fPoint;

  if (fDebug & 3) {

    cout << "PrimVtx (gen) " << " " << simPrimVtx.x() << " " << simPrimVtx.y() << " " << simPrimVtx.z() << endl;
    cout << "PrimVtx (rec) " << " " << fpEvt->fPrimaryVertex2.fPoint.x() 
	                     << " " << fpEvt->fPrimaryVertex2.fPoint.y() 
	                     << " " << fpEvt->fPrimaryVertex2.fPoint.z() << endl;
  }

  fFltR   = fltRec.Perp()*fMass/fPt;
  fFltG   = fltGen.Perp()*fMass/fPt;
  fFltRes = fFltR - fFltG;

  if (fDebug & 3) {

    cout <<" FltRec: " << fFltR << endl;
    cout <<" FltGen: " << fFltG << endl;
    cout <<" FltRes: " << fFltRes << endl;
  }

  if ( fltGen.Mag() == 0 || fltRec.Mag() == 0 ) { return; }

  double DLxy = fltRec.Perp() - fltGen.Perp();
  double DL3d = fltRec.Mag() - fltGen.Mag();
  double DTxy = fltRec.Perp()*fMass/fPt - fltGen.Perp()*mass/genPt;
  double DT3d = fltRec.Mag()*fMass/fP - fltGen.Mag()*mass/genPlab;
  double DTxy_s = (0.01*1.E15/299792458.)*(fltRec.Perp()*fMass/fPt - fltGen.Perp()*mass/genPt);
  double DT3d_s = (0.01*1.E15/299792458.)*(fltRec.Mag()*fMass/fP - fltGen.Mag()*mass/genPlab);

  if (fDebug & 3) {

    cout << " DLxy = " << DLxy
	 << " DL3d = " << DL3d << endl;
    cout << " DTxy = " << DTxy
	 << " DT3d = " << DT3d << endl;
  }

  if ( simVtx ) {
     
      h = (TH1D*)gDirectory->Get("v200"); h->Fill(DLxy);
      h = (TH1D*)gDirectory->Get("v201"); h->Fill(DL3d);
      h = (TH1D*)gDirectory->Get("v202"); h->Fill(DTxy);
      h = (TH1D*)gDirectory->Get("v203"); h->Fill(DT3d);
      h = (TH1D*)gDirectory->Get("v204"); h->Fill(DTxy_s);
      h = (TH1D*)gDirectory->Get("v205"); h->Fill(DT3d_s);
  } 
}



// ----------------------------------------------------------------------
void treeBmm::boxes() {
  //  cout << "  boxes " << endl;
}


// ----------------------------------------------------------------------
void treeBmm::histogram(int offset) {

  ((TH1D*)gDirectory->Get(Form("c%d00", offset)))->Fill(fPt);
  ((TH1D*)gDirectory->Get(Form("c%d01", offset)))->Fill(fEta);
  
  ((TH1D*)gDirectory->Get(Form("c%d10", offset)))->Fill(fPtL0);
  ((TH1D*)gDirectory->Get(Form("c%d11", offset)))->Fill(fEtaL0);
  ((TH1D*)gDirectory->Get(Form("c%d11", offset)))->Fill(fEtaL1);
  ((TH1D*)gDirectory->Get(Form("c%d12", offset)))->Fill(fPtL0);
  ((TH1D*)gDirectory->Get(Form("c%d12", offset)))->Fill(fPtL1);
  ((TH1D*)gDirectory->Get(Form("c%d13", offset)))->Fill(fTipL0);
  ((TH1D*)gDirectory->Get(Form("c%d13", offset)))->Fill(fTipL1);

  ((TH1D*)gDirectory->Get(Form("c%d20", offset)))->Fill(fRMM);
  ((TH1D*)gDirectory->Get(Form("c%d21", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d22", offset)))->Fill(fLxy/fSxy);
  ((TH1D*)gDirectory->Get(Form("c%d23", offset)))->Fill(fLxy);
  ((TH1D*)gDirectory->Get(Form("c%d24", offset)))->Fill(fIsoVeto);
  ((TH1D*)gDirectory->Get(Form("c%d25", offset)))->Fill(fIso);
  ((TH1D*)gDirectory->Get(Form("c%d26", offset)))->Fill(fL3d);
  ((TH1D*)gDirectory->Get(Form("c%d27", offset)))->Fill(fChi2);
  ((TH1D*)gDirectory->Get(Form("c%d29", offset)))->Fill(fProcessType);
  ((TH1D*)gDirectory->Get(Form("c%d30", offset)))->Fill(fMass);
  ((TH1D*)gDirectory->Get(Form("c%d31", offset)))->Fill(fMass);
  ((TH1D*)gDirectory->Get(Form("c%d32", offset)))->Fill(fMass);
  ((TH1D*)gDirectory->Get(Form("c%d33", offset)))->Fill(fMass);

  ((TH1D*)gDirectory->Get(Form("c%d35", offset)))->Fill(fpEvt->nRecTracks());

  ((TH1D*)gDirectory->Get(Form("c%d40", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d41", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d42", offset)))->Fill(fCosAngle);
  ((TH1D*)gDirectory->Get(Form("c%d43", offset)))->Fill(fCosAngle);

  ((TH1D*)gDirectory->Get(Form("c%d50", offset)))->Fill(fFltR);
  ((TH1D*)gDirectory->Get(Form("c%d51", offset)))->Fill(fFltG);
  ((TH1D*)gDirectory->Get(Form("c%d52", offset)))->Fill(fFltRes);


  ((TH2D*)gDirectory->Get(Form("C%d00", offset)))->Fill(fDptL0, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d01", offset)))->Fill(fPtL0, fEtaL0);
  ((TH2D*)gDirectory->Get(Form("C%d01", offset)))->Fill(fPtL1, fEtaL1);
  ((TH2D*)gDirectory->Get(Form("C%d02", offset)))->Fill(fPtL0, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d03", offset)))->Fill(fPtL0, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d04", offset)))->Fill(fPtL0, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d04", offset)))->Fill(fPtL1, fTipL1);
  ((TH2D*)gDirectory->Get(Form("C%d05", offset)))->Fill(fPtL0, fMass);
  ((TH2D*)gDirectory->Get(Form("C%d06", offset)))->Fill(fPtL0, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d07", offset)))->Fill(fPtL0, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d08", offset)))->Fill(fPtL0, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d09", offset)))->Fill(fPtL0, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d10", offset)))->Fill(fPtL0, fPt);
  ((TH2D*)gDirectory->Get(Form("C%d11", offset)))->Fill(fPtL0, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d12", offset)))->Fill(fEtaL0, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d13", offset)))->Fill(fEtaL0, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d14", offset)))->Fill(fEtaL0, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d15", offset)))->Fill(fEtaL0, fMass);
  ((TH2D*)gDirectory->Get(Form("C%d16", offset)))->Fill(fEtaL0, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d17", offset)))->Fill(fEtaL0, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d18", offset)))->Fill(fEtaL0, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d19", offset)))->Fill(fEtaL0, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d20", offset)))->Fill(fEtaL0, fPt);
  ((TH2D*)gDirectory->Get(Form("C%d21", offset)))->Fill(fEtaL0, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d22", offset)))->Fill(fMass, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d23", offset)))->Fill(fMass, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d24", offset)))->Fill(fMass, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d25", offset)))->Fill(fMass, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d26", offset)))->Fill(fMass, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d27", offset)))->Fill(fMass, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d28", offset)))->Fill(fMass, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d29", offset)))->Fill(fMass, fPt);
  ((TH2D*)gDirectory->Get(Form("C%d30", offset)))->Fill(fMass, fEta);

  ((TH2D*)gDirectory->Get(Form("C%d31", offset)))->Fill(fPt, fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d32", offset)))->Fill(fPt, fLxy/fSxy);
  ((TH2D*)gDirectory->Get(Form("C%d33", offset)))->Fill(fPt, fTipL0);
  ((TH2D*)gDirectory->Get(Form("C%d34", offset)))->Fill(fPt, fRMM);
  ((TH2D*)gDirectory->Get(Form("C%d35", offset)))->Fill(fPt, fCosAngle);
  ((TH2D*)gDirectory->Get(Form("C%d36", offset)))->Fill(fPt, fIsoVeto);
  ((TH2D*)gDirectory->Get(Form("C%d37", offset)))->Fill(fPt, fChi2);
  ((TH2D*)gDirectory->Get(Form("C%d38", offset)))->Fill(fPt, fEta);
}


// ----------------------------------------------------------------------
void treeBmm::fillAnalysisEff() {

  if (1) {
    histogram(0);
  }

  if (1 == fGoodKinematics) {  
    fAR1->Fill(1.1);  histogram(1); 
  } else {
    fAR1->Fill(2.1);
  }

  if (1 == fGoodKinematics && 1 == fGoodL1Trigger) {
    fAR1->Fill(11.1);  histogram(2);
  } else {
    fAR1->Fill(12.1);
  }

  if (1 == fGoodKinematics && 1 == fGoodL1Trigger  && 1 == fGoodHLTTrigger) {
    fAR1->Fill(21.1); histogram(3);
  } else {
    fAR1->Fill(22.1);
  }

  // -- cumulative cut effiencies
  if (1 == fGoodKinematics && 1 == fGoodL1Trigger  && 1 == fGoodHLTTrigger) {
  //  if (1 == fGoodKinematics) {

    if ((fPtL0 > PTLO) && (fPtL1 > PTLO)) {
      fAR1->Fill(100.1);

      if ((RMMLO < fRMM) && (fRMM < RMMHI)) {
	fAR1->Fill(101.1);
	
	if (fPt > PTBS) {
	  fAR1->Fill(120.1);
	  
	  if (ETALO < fEta && fEta < ETAHI) {
	    fAR1->Fill(121.1);
	    
	    if (fCosAngle > COSALPHA) {
	      fAR1->Fill(122.1);
	      
	      if (fLxy/fSxy > LXYSXYLO) {
		fAR1->Fill(123.1);
		histogram(5);
	  
		
		if (fIso > ISOLATION) {
		  fAR1->Fill(124.1);
		  
		  if (fChi2 < VTXCHI) {
		    fAR1->Fill(125.1);  histogram(4);
		  }
		}
	      }
	    }
	  }
	}

	
	// -- somewhat tightened preselection for the 'factorizing' cuts
	if (
	    fLxy/fSxy > 7
	    // && fCosAngle > 0.9
	    ) {
	  
	  if (fIso > ISOLATION) {
	    fAR1->Fill(224.1);
	  } else {
	    fAR1->Fill(324.1);
	  }
	  
	  if (fChi2 < VTXCHI) {
	    fAR1->Fill(226.1);
	  } else {	  
	    fAR1->Fill(326.1);
	  }
	}
      }
    }    
  }
}

// ----------------------------------------------------------------------
void treeBmm::setTitles(TH1 *h, const char *sx, const char *sy, float size,
			float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy);
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

void treeBmm::setTitles2(TH2 *h, const char *sx, const char *sy, float size,
			float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy);
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}
 

// void treeBmm::dumpSmallTree() {

//   if  (fSTI[0] < 0 || fSTI[1] < 0) return;



//   TDirectory *old = gDirectory;

//   static int first(1);

//   double px1 = fpL1->fPlab.X();
//   double py1 = fpL1->fPlab.Y();
//   double pz1 = fpL1->fPlab.Z();

//   TGenCand *g1 = fpEvt->getGenCand(fpL1->fGenIndex);

//   if (fpL1->fGenIndex != fGm0 || fpL2->fGenIndex != fGm1) {
//     if (fpL1->fGenIndex != fGm0) {
//       cout << "mismatch in generator pointer 1" << endl;
//     }

//     if (fpL2->fGenIndex != fGm1) {
//       cout << "mismatch in generator pointer 2" << endl;
//     }
    
//     cout << " ---------------------------------------------------------------------- " << endl;
//   }

//   double gx1 = g1->fP.X();
//   double gy1 = g1->fP.Y();
//   double gz1 = g1->fP.Z();

//   double px2 = fpL2->fPlab.X();
//   double py2 = fpL2->fPlab.Y();
//   double pz2 = fpL2->fPlab.Z();

//   TGenCand *g2 = fpEvt->getGenCand(fpL2->fGenIndex);
//   double gx2 = g2->fP.X();
//   double gy2 = g2->fP.Y();
//   double gz2 = g2->fP.Z();

//   TLorentzVector gL1, gL2;
//   gL1.SetXYZM(gx1, gy1, gz1, 0.105);
//   gL2.SetXYZM(gx2, gy2, gz2, 0.105);

//   TLorentzVector gB = gL1 + gL2; 
//   double gMass = gB.M();

//   if (first) {
//     first = 0; 
    
//     fFSarah = new TFile("sarah.root", "RECREATE");

//     fTSarah = new TTree("events", "events");

// //     t->Branch("run",    &run,"run/I");
//     fTSarah->Branch("mass",   &fMass,"mass/D");
//     fTSarah->Branch("gmass",  &gMass,"gmass/D");

//     fTSarah->Branch("px1",     &px1,"px1/D");
//     fTSarah->Branch("py1",     &py1,"py1/D");
//     fTSarah->Branch("pz1",     &pz1,"pz1/D");

//     fTSarah->Branch("gx1",     &gx1,"gx1/D");
//     fTSarah->Branch("gy1",     &gy1,"gy1/D");
//     fTSarah->Branch("gz1",     &gz1,"gz1/D");

//     fTSarah->Branch("px2",     &px2,"px2/D");
//     fTSarah->Branch("py2",     &py2,"py2/D");
//     fTSarah->Branch("pz2",     &pz2,"pz2/D");
    
//     fTSarah->Branch("gx2",     &gx2,"gx2/D");
//     fTSarah->Branch("gy2",     &gy2,"gy2/D");
//     fTSarah->Branch("gz2",     &gz2,"gz2/D");
//   }
  
//   fTSarah->Fill();

//   old->cd();


// }
