// -- The big global
TAna00Event *pEvent;
TAna00Event *pEvent2;
TAna00Event *pTempEvt;


// ----------------------------------------------------------------------
void v1(const char *filename = "bmm-test.root") {
  
  TH1D *h0 = new TH1D("h0", " ", 100, 5.0, 5.6);
  TH1D *h1 = new TH1D("h1", " ", 100, 5.0, 5.6);
  TH1D *h2 = new TH1D("h2", " ", 100, 5.0, 5.6);
  TH2D *h3 = new TH2D("h3", " ", 20, 0.0, 20., 20, 0.0, 20.);
  TH2D *h4 = new TH2D("h4", " ", 20, 0.0, 20., 20, 0.0, 20.);
  setTitles(h1, "p^{#mu}_{T} [GeV]", "Entries/Bin", 0.06);
  setTitles(h2, "#eta^{#mu} ", "Entries/Bin", 0.06);
  setFilledHist(h1);
  setFilledHist(h2);

  TH1D *h11 = new TH1D("h11", " ", 30, 0.0, 30.);
  TH1D *h12 = new TH1D("h12", " ", 50, -10.0, 10.0);
  setTitles(h11, "p^{B}_{T} [GeV]", "Entries/Bin", 0.06);
  setTitles(h12, "#eta^{B} ", "Entries/Bin", 0.06);
  setFilledHist(h11);
  setFilledHist(h12);

  TChain chain("T1");
  chain.Add(filename);

  TLorentzVector p4B, p4m0, p4m1; 
  TLorentzVector *p4m; 

  // -- Set up for reading
  Int_t nentries(0), nb(0);
  Int_t iEvent(0), ic(0);
  
  pEvent = new TAna00Event(0);
  TAnaTrack *pTrack;
  TAnaCand  *pCand;
  TGenCand  *pGen;
  
  chain.SetBranchAddress("TAna00Event", &pEvent);
  nentries = chain.GetEntries();
  
  cout << "Found " << nentries << " entries in the chain" << endl;
  
  for (iEvent = 0; iEvent < nentries; iEvent++) {
    pEvent->Clear();
    nb += chain.GetEvent(iEvent); 
    
    for (ic = 0; ic < pEvent->nCands(); ++ic) {
      pCand = pEvent->getCand(ic);
      pCand->dump();
      if (pCand->fVtx.fType == 1) h1->Fill(pCand->fMass);
      if (pCand->fVtx.fType == 2) h2->Fill(pCand->fMass);
      if (pCand->fVtx.fType == 3) h3->Fill(pCand->fMass);
    }
    
  }
  
}



// ----------------------------------------------------------------------
void plotMass(const char *filename = "bmm-test.root") {
  
  TH1D *h1 = new TH1D("mll1", "Generator PID", 40, 5.0, 5.8);
  TH1D *h2 = new TH1D("mll2", "Reco PID", 40, 5.0, 5.8);

  TH1D *h10 = new TH1D("h10", "track pT", 100, 0., 20.);
  TH1D *h11 = new TH1D("h11", "muon  pT", 100, 0., 20.);
  
  setTitles(h1, "m_{#mu#mu} [GeV]", "Entries/Bin", 0.06);
  setTitles(h2, "m_{#mu#mu} [GeV]", "Entries/Bin", 0.06);
  setFilledHist(h1);
  setFilledHist(h2);

  TChain chain("T1");
  chain.Add(filename);

  TLorentzVector p4B, p4m0, p4m1; 
  TVector3 m0p3, m1p3;  

  // -- Set up for reading
  Int_t nentries(0), nb(0);
  Int_t iEvent(0), it(0);
  
  pEvent = new TAna00Event(0);
  TAnaTrack *pTrack;
  TAnaCand  *pCand;
  TGenCand  *pGen;
  
  chain.SetBranchAddress("TAna00Event", &pEvent);
  nentries = chain.GetEntries();
  
  cout << "Found " << nentries << " entries in the chain" << endl;
  
  for (iEvent = 0; iEvent < nentries; iEvent++) {
    pEvent->Clear();
    nb += chain.GetEvent(iEvent); 

    int signalTracks(0), bt0(-1), bt1(-1);
    for (it = 0; it < pEvent->nRecTracks(); ++it) {
      pTrack = pEvent->getRecTrack(it);
      
      h10->Fill(pTrack->fPlab.Pt());

      if (trkIsDescendant(pTrack, 531) > 0) {
	++signalTracks;
 	cout << pTrack->fGenIndex << "  " 
 	     << pEvent->getGenCand(pTrack->fGenIndex)->fID << " " 
 	     << pTrack->fMCID << endl;

	h11->Fill(pTrack->fPlab.Pt());

	if (bt0 < 0) {
	  bt0 = it;
	} else {
	  if (bt1 < 0) {
	    bt1 = it;
	  }
	  else {
	    cout << "XXXXXXXXXX too many Bs tracks" << endl;
	  }
	}
	
      }
      
    }
    
    
    if (signalTracks == 2) {
      if (bt0 > -1 && bt1 > -1) {
    	// cout << "Using tracks " << bt0 << " and " << bt1 << " as Bs tracks" << endl;
	m0p3 = pEvent->getRecTrack(bt0)->fPlab;
	m1p3 = pEvent->getRecTrack(bt1)->fPlab;
  	p4m0.SetVectM(m0p3, 0.105);
  	p4m1.SetVectM(m1p3, 0.105);
 	p4B = p4m0 + p4m1;
      }
      
      h1->Fill(p4B.M());
      //       cout << " --> Found Bs with mass = " << p4B.M()
      // 	   << endl; 
    }
    
    if (pEvent->nCands() > 0) {
      pCand = pEvent->getCand(0);
      h2->Fill(pCand->fMass);
    }

  }
  
  // -- Display

  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);

  gStyle->SetOptStat(0);  
  c0.Clear();
  c0.Divide(1,2);
  c0.cd(1); shrinkPad(0.15);
  h1->Draw();
  tl->DrawLatex(0.7, 0.8, Form("Entries: %4d", int(h1->GetEntries())));
  tl->DrawLatex(0.7, 0.7, Form("RMS: %4.3f", double(h1->GetRMS())));

  c0.cd(2); shrinkPad(0.15);
  h2->Draw();
  tl->DrawLatex(0.7, 0.8, Form("Entries: %4d", int(h2->GetEntries())));
  tl->DrawLatex(0.7, 0.7, Form("RMS: %4.3f", double(h2->GetRMS())));

  c0.SaveAs("recplots-0.eps");
  c0.SaveAs("recplots-0.ps");

}


// ----------------------------------------------------------------------
int trkIsDescendant(TAnaTrack *pTrk, int ID, int matchCharge = 0) {
  
  int result(-2);

  int genIndex = pTrk->fGenIndex;
  if ((genIndex < 0) || (genIndex > pEvent->nGenCands())) {
    return result;
  }

  result = 0; 

  TGenCand *pGen = pEvent->getGenCand(genIndex);

  int mother = pGen->fMom1;
  while (mother > 0) {
    pGen = pEvent->getGenCand(mother);
    if (pGen->fID == ID) {
      return 1;
    }

    if (matchCharge) {
      if (pGen->fID == ID) {
	return 1;
      } 
    } else {
      if (TMath::Abs(pGen->fID) == ID) {
	return 1;
      } 
    }      
    
    mother = pGen->fMom1;
  }

  return result; 
}


// ----------------------------------------------------------------------
int genIsDescendant(TGenCand *pGen, int ID, int matchCharge = 0) {
  
  int result(-2);

  int mother = pGen->fMom1;
  while (mother > 0) {
    pGen = pEvent->getGenCand(mother);
    if (pGen->fID == ID) {
      return 1;
    }

    if (matchCharge) {
      if (pGen->fID == ID) {
	return 1;
      } 
    } else {
      if (TMath::Abs(pGen->fID) == ID) {
	return 1;
      } 
    }      
    
    mother = pGen->fMom1;
  }

  return result; 
}


// ----------------------------------------------------------------------
void genPlots(const char *filename = "bmm-test.root") {

  TH1D *h0 = new TH1D("h0", " ", 600, 0.0, 600.);
  TH1D *h1 = new TH1D("h1", " ", 30, 0.0, 30.);
  TH1D *h2 = new TH1D("h2", " ", 50, -10.0, 10.0);
  TH2D *h3 = new TH2D("h3", " ", 20, 0.0, 20., 20, 0.0, 20.);
  TH2D *h4 = new TH2D("h4", " ", 20, 0.0, 20., 20, 0.0, 20.);
  setTitles(h1, "p^{#mu}_{T} [GeV]", "Entries/Bin", 0.06);
  setTitles(h2, "#eta^{#mu} ", "Entries/Bin", 0.06);
  setFilledHist(h1);
  setFilledHist(h2);

  TH1D *h11 = new TH1D("h11", " ", 30, 0.0, 30.);
  TH1D *h12 = new TH1D("h12", " ", 50, -10.0, 10.0);
  setTitles(h11, "p^{B}_{T} [GeV]", "Entries/Bin", 0.06);
  setTitles(h12, "#eta^{B} ", "Entries/Bin", 0.06);
  setFilledHist(h11);
  setFilledHist(h12);

  TChain chain("T1");
  chain.Add(filename);

  TLorentzVector p4B, p4m0, p4m1; 
  TLorentzVector *p4m; 
  TVector3 m0p3, m1p3;  

  // -- Set up for reading
  Int_t nentries(0), nb(0);
  Int_t iEvent(0), it(0);
  
  pEvent = new TAna00Event(0);
  TAnaTrack *pTrack;
  TAnaCand  *pCand;
  TGenCand  *pGen;
  
  chain.SetBranchAddress("TAna00Event", &pEvent);
  nentries = chain.GetEntries();
  
  cout << "Found " << nentries << " entries in the chain" << endl;
  
  //  for (iEvent = 0; iEvent < nentries; iEvent++) {
  for (iEvent = 0; iEvent < 100; iEvent++) {
    pEvent->Clear();
    nb += chain.GetEvent(iEvent); 

    int mID(0); 
    int amID(0);

    //    cout << "Found " << pEvent->nGenCands() << " gen tracks in event" << endl;
    int signalTracks(0), bt0(-1), bt1(-1);
    double pt0(-1.), pt1(-1.); 
    for (it = 0; it < pEvent->nGenCands(); ++it) {
      pGen = pEvent->getGenCand(it);

      mID = pGen->fID;
      amID = TMath::Abs(mID);

      h0->Fill(amID);
      
      // -- Signal B 
      if (amID == 531) {
	p4m = pGen->fP;
	h11->Fill(p4m->Pt());
	h12->Fill(p4m->Eta());
      }

      // -- Muons that can be traced back to a Bs (also cascade muons)
      if ((amID == 13) && (genIsDescendant(pGen, 531) > 0)) {
	p4m = pGen->fP;
	h1->Fill(p4m->Pt());
	h2->Fill(p4m->Eta());

	if (pt0 < 0) {
	  pt0 = p4m->Pt();
	} else {
	  pt1 = p4m->Pt();
	}
      }
    }

    if (pt0 > pt1) {
      h3->Fill(pt0, pt1);
    } else{
      h3->Fill(pt1, pt0);
    }

    if (pt0 > 7.0 && pt1 > 7.0) {
      h4->Fill(pt0, pt1);
    }
  }

  // -- Display
  gStyle->SetOptStat(0);  
  c0.Clear();
  c0.Divide(2,2);
  c0.cd(1); shrinkPad(0.15);
  h11->Draw();
  c0.cd(2); shrinkPad(0.15);
  gPad->SetLogy(1);
  h0->Draw();
  h0->SetMinimum(0.5);

  c0.cd(3); shrinkPad(0.15);
  h1->Draw();
  c0.cd(4); shrinkPad(0.15);
  h3->Draw("box");

  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->DrawLatex(0.5, 0.8, Form("p_{i} > 7 GeV: %4d", int(h4->GetEntries())));

  
  c0.SaveAs("genplots-0.eps");
  c0.SaveAs("genplots-0.ps");
}

// ----------------------------------------------------------------------
void yield(int ch = 1, int add = 1) {


  const char *filename1[200];
  const char *filename2[200];

  switch(ch)
  {

  case 1:
    sprintf(filename1, "test/test-no-hlt-bad.root");
    sprintf(filename2, "test/test-tk-no-hlt-bad.root");
    break;

  case 2:
    sprintf(filename1, "test/test-no-hlt.root");
    sprintf(filename2, "test/test-tk-no-hlt.root");
    break;

  case 3:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt.root");
    break;
      
  case 4:
    sprintf(filename1, "test/test-with-hlt-bad_1.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    break;
      
  case 5:
    sprintf(filename1, "test/test-with-hlt-bad_6.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_6.root");
    break;
      
  case 6:
    sprintf(filename1, "test/test-no-hlt.root");
    sprintf(filename2, "test/test-with-hlt.root");
    break;
      
  case 7:
    sprintf(filename1, "test/test-tk-no-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt.root");
    break;
      
  case 8:
    sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_cscs/output_6.root");
    sprintf(filename2, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_fnal/output_6.root");
    break;

  case 9:
    sprintf(filename1, "test/test-tk-no-hlt-byChi2.root");
    sprintf(filename2, "test/test-tk-no-hlt-byHits.root");
    break;

  case 10:
    sprintf(filename1, "test/test-tk-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    break;

  case 11:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-with-hlt-bad_1.root");
    break;

  case 12:
    sprintf(filename1, "test/test-tk-no-hlt.root");
    sprintf(filename2, "test/test-tk-no-hlt-byHits.root");
    break;
      
  default: 
    cout << " Enter a number < 13 !!!" << endl;
    return;
    
  };


  TFile *f1 = new TFile(Form("%s", filename1));
  TH2D *h1 =  (TH2D*)gDirectory->Get("glb");

  TFile *f2 = new TFile(Form("%s", filename2));
  TH2D *h2 =  (TH2D*)gDirectory->Get("glb");

  if ( add ) {

    h2->Add(h1, -1);
    h2->GetZaxis()->SetRangeUser(-4,4);
    
    double tot(0), nbins(1);
    for (int i = 0; i < h2->GetXaxis()->GetNbins()+1; i++ ) {
      for (int j = 0; j < h2->GetYaxis()->GetNbins()+1; j++ ) {
	
	if (h2->GetBinContent(i,j) == 0) {
	  h2->SetBinContent(i,j,-10);
	} else {
	  tot += h2->GetBinContent(i,j);
	  nbins++;
	}
      }
    }
    
    h2->SetTitle(Form ("Global muon yield: N_{w/o reco} - N_{w/ reco} #rightarrow mean = %4.4f"
		     , tot/nbins));
  } else {
    
    h2->Divide(h1);
    h2->GetZaxis()->SetRangeUser(0,4);

    double tot(1), nbins(1);
    for (int i = 0; i < h2->GetXaxis()->GetNbins()+1; i++ ) {
      for (int j = 0; j < h2->GetYaxis()->GetNbins()+1; j++ ) {
	
	if (h2->GetBinContent(i,j) != 0) {
	  tot += h2->GetBinContent(i,j);
	  nbins++;
	}
      }
    }

    h2->SetTitle(Form ("Global muon yield: #frac{N_{w/o reco}}{N_{w/ reco}} #rightarrow mean = %4.4f"
		     , tot/nbins));
  }



  gStyle->SetOptStat(0);  
  c0.Clear();

  c0.cd(); 

  gPad->SetLogy(0);
  h2->Draw("colz");

  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.03);
  tl->DrawLatex(0.56, 0.86, Form("%s", filename1));
  tl->DrawLatex(0.56, 0.80, Form("%s", filename2));

  c0.SaveAs(Form("yield%i_%i.pdf",ch, add));

}
// ----------------------------------------------------------------------
void muonId() { 

  const char *filename1[200];

  for (int ch = 1; ch < 13; ch++) {

    switch(ch)
    {
	
      case 1:
	sprintf(filename1, "test/test-no-hlt.root");
	break;

      case 2:
	sprintf(filename1, "test/test-with-hlt.root");
	break;
      
      case 3:
	sprintf(filename1, "test/test-with-hlt-bad_1.root");
	break;
      
      case 4:
	sprintf(filename1, "test/test-with-hlt-bad_6.root");
	break;
	
      case 5:
	sprintf(filename1, "test/test-tk-no-hlt.root");
	break;
      
      case 6:
	sprintf(filename1, "test/test-tk-with-hlt.root");
	break;
      
      case 7:
	sprintf(filename1, "test/test-tk-with-hlt-bad_1.root");
	break;
      
      case 8:
	sprintf(filename1, "test/test-tk-with-hlt-bad_6.root");
	break;

      case 9:
	sprintf(filename1, "test/test-tk-no-hlt-byChi2.root");
	break;
      
      case 10:
	sprintf(filename1, "test/test-tk-no-hlt-byHits.root");
	break;
      
      case 11:
	sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_cscs/output_6.root");
	break;
      
      case 12:
	sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_fnal/output_6.root");
	break;
      
      case 13:
	sprintf(filename1, "grid/analysed_Spring07_noHLT_new/Bs2MuMu/output_2.root");
	break;

      default: 
	cout << " Enter a number < 14 !!!" << endl;
	return;
    
    };
  
 
    TChain chain("T1");
    chain.Add(Form("%s", filename1));
      
    // -- Set up for reading
    Int_t nentries(0), nb(0);
    Int_t iEvent(0), it(0), idx(0);
    
    delete pEvent; 
    
    pEvent = new TAna00Event(0);
    
    TAnaTrack *pTrack;
    TAnaCand  *pCand;
    TGenCand  *pGen;
    
    chain.SetBranchAddress("TAna00Event", &pEvent);
    nentries = chain.GetEntries();
    
    double pT(0.);
    int totMu(0),   totMu8(0),  recMu(0), lostMu(0), misMu(0), recPar(0); 
    int totGlb(0), totGlb8(0),  totMuG(0);  
    for (iEvent = 0; iEvent < nentries; iEvent++) {

      pEvent->Clear();    
      nb += chain.GetEvent(iEvent); 

      for (it = 0; it < pEvent->nSigTracks(); ++it) {
	
	pTrack = pEvent->getSigTrack(it);
	pT     = pTrack->fPlab.Pt();
	idx    = pTrack->fMuType;

	if ( idx !=  3 ) { 
	  continue; 
	}

	totGlb++;

	if ( pT > 8. ) {
	  totGlb8++;
	}
      }
	
      for (it = 0; it < pEvent->nRecTracks(); ++it) {
	
	pTrack = pEvent->getRecTrack(it);	
	pT = pTrack->fPlab.Pt();
	
	if (pTrack->fMuID >= 0.) {
	  totMuG++;
	}

	if (TMath::Abs(pTrack->fMCID) == 13) {
	  totMu++;
	}

	if ( pT < 8.) {
	  continue;
	}

	if (TMath::Abs(pTrack->fMCID) == 13) {
	  totMu8++;
	  if (pTrack->fMuID < 0.0) {
	    lostMu++;
	  } else {
	    recMu++;
	  }
	} else if (TMath::Abs(pTrack->fMCID) == 211 ||
		   TMath::Abs(pTrack->fMCID) == 321 ||
		   TMath::Abs(pTrack->fMCID) == 2212 ) {
	  if (pTrack->fMuID < 0.0) {
	    recPar++;
	  } else {
	    misMu++;
	  }
	}
      }
    }

    cout  << Form("%s",filename1) << ": "  << nentries << " events "  << endl
	  << "  " << totMu << " (resp. " << totMu8 << ") muons, " 
	  << totGlb << " (resp. " << totGlb8 << ") global muons  ----  " 
	  << totMuG << " assoc. muons " << endl;
    cout << "  muon-eff  = " << Form("%4.1f", 100*recMu/(1.*totMu8)) << " %"
	 << "  ( muon-lost = " << Form("%4.1f", 100*lostMu/(1.*totMu8)) << " %)" << endl;
    cout << "  muon-mis  = " << Form("%4.1f", 100*misMu/(1.*totMu8)) << " %" << endl << endl;
  }
}

// ----------------------------------------------------------------------
void tracks(int ch = 3) { 

  const char *filename1[200];
  const char *filename2[200];

  int start(0), stop(4500), offset(0);

  switch(ch)
  {

  case 1:
    sprintf(filename1, "test/test-no-hlt-bad.root");
    sprintf(filename2, "test/test-tk-no-hlt-bad.root");
    break;

  case 2:
    sprintf(filename1, "test/test-no-hlt.root");
    sprintf(filename2, "test/test-tk-no-hlt.root");
    break;

  case 3:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt.root");
    break;
      
  case 4:
    sprintf(filename1, "test/test-with-hlt-bad_1.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    start = 0; stop = 10000;
    break;
      
  case 5:
    sprintf(filename1, "test/test-with-hlt-bad_6.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_6.root");
    start = 0; stop = 10000;
    break;
      
  case 6:
    sprintf(filename1, "test/test-no-hlt.root");
    sprintf(filename2, "test/test-with-hlt.root");
    break;
      
  case 7:
    sprintf(filename1, "test/test-tk-no-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt.root");
    break;
      
  case 8:
    sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_cscs/output_6.root");
    sprintf(filename2, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_fnal/output_6.root");
    start = 0; stop = 10000;
    break;

  case 9:
    sprintf(filename1, "test/test-tk-no-hlt-byChi2.root");
    sprintf(filename2, "test/test-tk-no-hlt-byHits.root");
    break;

  case 10:
    sprintf(filename1, "test/test-tk-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    start = 1500; stop = 3000;
    break;

  case 11:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-with-hlt-bad_1.root");
    start = 1500; stop = 3000;
    break;

  case 12:
    sprintf(filename1, "test/test-tk-no-hlt.root");
    sprintf(filename2, "test/test-tk-no-hlt-byHits.root");
    break;

  case 13:
    sprintf(filename1, "test/test-tk-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    start = 3000; stop = 4500; offset = 3000;
    break;

  case 14:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-with-hlt-bad_1.root");
    start = 3000; stop = 4500; offset = 3000;
    break;
      
  case 15:
    sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu/output_6.root");
    sprintf(filename2, "grid/analysed_Spring07_noHLT_new/Bs2MuMu/output_2.root");
    start = 2541; stop = 4041; offset = -500;
    break;

  default: 
    cout << " Enter a number < 16 !!!" << endl;
    return;
    
  };

  // -- Plots for sample 1 
  TH1D *h0 = new TH1D("h0", "mass ", 500, 0.0, 50.);
  TH1D *h1 = new TH1D("h1", "pT ", 3000, 0.0, 30.);
  TH1D *h2 = new TH1D("h2", "eta ", 5000, -10.0, 10.0);

  setTitles(h0, "m_{#mu #mu}  [GeV]", "Entries/Bin", 0.06);
  setTitles(h1, "p^{#mu}_{T, rec} [GeV]", "Entries/Bin", 0.06);
  setTitles(h2, "#eta^{#mu}_{rec} ", "Entries/Bin", 0.06);

  setFilledHist(h0);
  setFilledHist(h1);
  setFilledHist(h2);

  TH1D *h11 = new TH1D("h11", " ", 30, 0.0, 30.);
  TH1D *h12 = new TH1D("h12", " ", 50, -10.0, 10.0);

  setTitles(h11, "p^{#mu}_{T, gen} [GeV]", "Entries/Bin", 0.06);
  setTitles(h12, "#eta^{#mu}_{gen} ", "Entries/Bin", 0.06);
  
  setFilledHist(h11);
  setFilledHist(h12);

  // -- Correlation plots (pT, mass)
  TH2D *h4 = new TH2D("h4", "Global muon matching, p_{T} > 3.", 400, -20., 20., 400, -20., 20.);
  TH2D *h5 = new TH2D("h5", "Truth-matching, p_{T} > 3. (Associator)", 400, -20., 20., 400, -20., 20.);
  TH2D *h6 = new TH2D("h6", "Generator particles, p_{T} > 3.", 400, -20., 20., 400, -20., 20.);
  TH2D *h7 = new TH2D("h7", "Inv. mass of candidate with best #chi^{2}", 220, -1., 10., 220, -1., 10.);

  h4->GetXaxis()->SetTitle("p_{T, rec}^{1} [GeV]"); h4->GetYaxis()->SetTitle("p_{T, rec}^{2} [GeV]");
  h5->GetXaxis()->SetTitle("p_{T, rec}^{1} [GeV]"); h5->GetYaxis()->SetTitle("p_{T, rec}^{2} [GeV]");
  h6->GetXaxis()->SetTitle("p_{T, gen}^{1} [GeV]"); h6->GetYaxis()->SetTitle("p_{T, gen}^{2} [GeV]");
  h7->GetXaxis()->SetTitle("m_{#mu #mu}^{1} [GeV]"); h7->GetYaxis()->SetTitle("m_{#mu #mu}^{2} [GeV]");

  // -- Correlation plots (numbers)
  TH2D *h03 = new TH2D("h03", " ", 200, 0.0, 200., 200, 0.0, 200.);
  TH2D *h13 = new TH2D("h13", " ", 1000, 0.0, 1000., 1000, 0.0, 1000.);
  TH2D *h23 = new TH2D("h23", " ", 30, 0.0, 30., 30, 0.0, 30.);  
  
  h03->GetXaxis()->SetTitle("N_{1}^{rec}"); h03->GetYaxis()->SetTitle("N_{2}^{rec}");
  h13->GetXaxis()->SetTitle("N_{1}^{gen}"); h13->GetYaxis()->SetTitle("N_{2}^{gen}");
  h23->GetXaxis()->SetTitle("N_{1}^{cand}"); h23->GetYaxis()->SetTitle("N_{2}^{cand}"); 

  TChain chain("T1");
  chain.Add(Form("%s", filename1));
  
  TChain chain2("T1");
  chain2.Add(Form("%s", filename2));
  
  TLorentzVector p4B, p4m0, p4m1; 
  TLorentzVector *p4m, *p4m2; 
  TVector3 m0p3, m1p3;  
  
  // -- Set up for reading
  Int_t nentries(0), nentries2(0), nb(0), nb2(0);
  Int_t iEvent(0), jEvent(0), it(0);
  
  delete pEvent; 
  delete pEvent2;
  
  pEvent = new TAna00Event(0);
  pEvent2 = new TAna00Event(0);
  
  TAnaTrack *pTrack, *pTrack2;
  TAnaCand  *pCand, *pCand2;
  TGenCand  *pGen, *pGen2;
  
  chain.SetBranchAddress("TAna00Event", &pEvent);
  nentries = chain.GetEntries();
  
  chain2.SetBranchAddress("TAna00Event", &pEvent2);
  nentries2 = chain2.GetEntries();
  
  cout << "Found " << nentries << " entries in the 1. chain" << endl;
  cout << "Found " << nentries2 << " entries in the 2. chain" << endl;
  
  int ella(0), tmp(0), match(0), nmatch(0);
  // for (iEvent = 0; iEvent < nentries; iEvent++) {
  for (iEvent = start; iEvent < stop; iEvent++) {
    
    if ( int((iEvent-start)/100) > ella ) {
      cout << "... events " << iEvent << " of " << nentries << " (matches " << nmatch << ")" << endl;
      ella++;
    }

    pEvent->Clear();    
    nb += chain.GetEvent(iEvent); 

    match = 0;

    pEvent2->Clear();
    tmp = chain2.GetEvent(iEvent+offset);

    if ( (pEvent->fRunNumber == pEvent2->fRunNumber)
	 && (pEvent->fEventNumber == pEvent2->fEventNumber) ) {
      
      pEvent2->Clear();
      nb2 += chain2.GetEvent(iEvent+offset);
      match = 1; nmatch++;

    } else {

      // for (jEvent = 0; jEvent < nentries2; jEvent++) {
      for (jEvent = start + offset; jEvent < stop + offset; jEvent++) {
	
	pEvent2->Clear();
	tmp = chain2.GetEvent(jEvent);
	
	if ( (pEvent->fRunNumber == pEvent2->fRunNumber)
	     && (pEvent->fEventNumber == pEvent2->fEventNumber) ) {
	  
	  pEvent2->Clear();
	  nb2 += chain2.GetEvent(jEvent);
	  match = 1; nmatch++;
	  break;
	}
      }
    }
   
    if ( match == 0 ) { 
      
      continue; 
    }
    
 
    //    cout << "Found " << pEvent->nGenCands() << " gen tracks in event" << endl;
    int signalTracks(0), bt0(-1), bt1(-1);
    double pt0(-1.), pt1(-1.);
    
    double pT1(0.), pT2(0.);
    double gen_pT1(0.), gen_pT2(0.);
    double mass1(0.), mass2(0.);
    
    int id(0), id2(0);
    int bmmsel(0), bmmsel2(0);
    double chi_min(99999), chi_min2(99999.);
    
    h03->Fill(pEvent->nRecTracks(), pEvent2->nRecTracks());
    h13->Fill(pEvent->nGenCands(), pEvent2->nGenCands());
    h23->Fill(pEvent->nCands(), pEvent2->nCands());
    
    for (it = 0; it < TMath::Min(pEvent->nRecTracks(), pEvent2->nRecTracks()); ++it) {

      pTrack = pEvent->getRecTrack(it);
      pTrack2 = pEvent2->getRecTrack(it);
      
      pT1 = pTrack->fPlab.Pt();
      pT2 = pTrack2->fPlab.Pt();

      if ( pT1 < 3. || pT2 < 3. ) {
	continue;
      }
      
      if ( pTrack->fMuID >= 0. && pTrack2->fMuID >= 0.) {
	h4->Fill(pTrack->fPlab.Pt(),pTrack2->fPlab.Pt());
      } else if ( pTrack->fMuID >= 0. && !(pTrack2->fMuID >= 0.) ) {
	h4->Fill(pTrack->fPlab.Pt(),-pTrack2->fPlab.Pt());
      } else if ( !(pTrack->fMuID >= 0.) && pTrack2->fMuID >= 0.) {
	h4->Fill(-pTrack->fPlab.Pt(),pTrack2->fPlab.Pt());
      } else if ( !(pTrack->fMuID >= 0.) && !(pTrack2->fMuID >= 0.) ) {
	h4->Fill(-pTrack->fPlab.Pt(),-pTrack2->fPlab.Pt());
      } 

      if ( TMath::Abs(pTrack->fMCID) == 13 && TMath::Abs(pTrack2->fMCID) == 13) {
	h5->Fill(pTrack->fPlab.Pt(),pTrack2->fPlab.Pt());
      } else if ( TMath::Abs(pTrack->fMCID) == 13 && !(TMath::Abs(pTrack2->fMCID) == 13) ) {
	h5->Fill(pTrack->fPlab.Pt(),-pTrack2->fPlab.Pt());
      } else if ( !(TMath::Abs(pTrack->fMCID) == 13) && TMath::Abs(pTrack2->fMCID) == 13) {
	h5->Fill(-pTrack->fPlab.Pt(),pTrack2->fPlab.Pt());
      } else if ( !(TMath::Abs(pTrack->fMCID) == 13) && !(TMath::Abs(pTrack2->fMCID) == 13) ) {
	h5->Fill(-pTrack->fPlab.Pt(),-pTrack2->fPlab.Pt());
      } 

      h1->Fill(pTrack->fPlab.Pt());
      h2->Fill(pTrack->fPlab.Eta());
    }

    for (it = 2; it < TMath::Min(pEvent->nGenCands(),pEvent2->nGenCands()); ++it) {

      pGen = pEvent->getGenCand(it);
      pGen2 = pEvent2->getGenCand(it);

      p4m = pGen->fP;
      p4m2 = pGen2->fP;

      gen_pT1 = p4m->Pt();
      gen_pT2 = p4m->Pt();

      if ( gen_pT1 < 3. || gen_pT2 < 3. ) {
	continue;
      }

      if ( TMath::Abs(pGen->fID) == 13 && TMath::Abs(pGen2->fID) == 13) {
	h6->Fill(pTrack->fPlab.Pt(),pTrack2->fPlab.Pt());
      } else if ( TMath::Abs(pGen->fID) == 13 && !(TMath::Abs(pGen2->fID) == 13) ) {
	h6->Fill(pTrack->fPlab.Pt(),-pTrack2->fPlab.Pt());
      } else if ( !(TMath::Abs(pGen->fID) == 13) && TMath::Abs(pGen2->fID) == 13) {
	h6->Fill(-pTrack->fPlab.Pt(),pTrack2->fPlab.Pt());
      } else if ( !(TMath::Abs(pGen->fID) == 13) && !(TMath::Abs(pGen2->fID) == 13) ) {
	h6->Fill(-pTrack->fPlab.Pt(),-pTrack2->fPlab.Pt());
      } 

      h11->Fill(p4m->Pt());
      h12->Fill(p4m->Eta());
    }
    
    
    for (it = 0; it < pEvent->nCands(); ++it) {
      
      pCand = pEvent->getCand(it);
      
      id     = int(pCand->fType/10);
      bmmsel = pCand->fType - 10*id;
      
      if ( id == 531 && bmmsel == 2 ) {
	
	if ( pCand->fVtx.fChi2 < chi_min ) {
	  mass1 = pCand->fMass;
	  chi_min = pCand->fVtx.fChi2;
	}
      }
    }

    h0->Fill(mass1);
    
    for (it = 0; it < pEvent2->nCands(); ++it) {
      
      pCand2 = pEvent2->getCand(it);

      id2     = int(pCand2->fType/10);
      bmmsel2 = pCand2->fType - 10*id2;
      
      if ( id2 == 531 && bmmsel2 == 2 ) {
	
	if ( pCand2->fVtx.fChi2 < chi_min2 ) {
	  mass2 = pCand2->fMass;
	  chi_min2 = pCand2->fVtx.fChi2;
	}
      }
    }

    h7->Fill(mass1, mass2);
  }
 

  // -- Display
  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.04);

  gStyle->SetOptStat(0);  

  c0.Clear(); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h03->Draw("colz");
  tl->DrawLatex(0.26, 0.86, Form("%s (2)",filename2));
  tl->DrawLatex(0.36, 0.80, Form("vs. %s (1)", filename1));
  c0.SaveAs(Form("sample%i_ntracks.pdf",ch));

  c0.Clear();
  c0.Divide(2,2);

  c0.cd(1); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h4->Draw("colz");
  c0.cd(2); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h5->Draw("colz");

  c0.cd(3); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h6->Draw("colz");
  c0.cd(4); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h7->Draw("colz");

  c0.cd(1);
  tl->DrawLatex(0.26, 0.86, Form("%s (2)",filename2));
  tl->DrawLatex(0.36, 0.80, Form("vs. %s (1)", filename1));
  
  c0.SaveAs(Form("sample%i_tracks.pdf",ch));

  delete h0;
  delete h1;
  delete h2;
  delete h4; 
  delete h5; 
  delete h6; 
  delete h7;
  delete h11; 
  delete h12; 
  delete h03;
  delete h13;
  delete h23;
  delete tl;
}

// ----------------------------------------------------------------------
void muons(int ch = 3) { 

  const char *filename1[200];
  const char *filename2[200];

  int start(0), stop(4500), offset(0);

  switch(ch)
  {

  case 1:
    sprintf(filename1, "test/test-no-hlt-bad.root");
    sprintf(filename2, "test/test-tk-no-hlt-bad.root");
    break;

  case 2:
    sprintf(filename1, "test/test-no-hlt.root");
    sprintf(filename2, "test/test-tk-no-hlt.root");
    break;

  case 3:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt.root");
    break;
      
  case 4:
    sprintf(filename1, "test/test-with-hlt-bad_1.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    start = 0; stop = 10000;
    break;
      
  case 5:
    sprintf(filename1, "test/test-with-hlt-bad_6.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_6.root");
    start = 0; stop = 10000;
    break;
      
  case 6:
    sprintf(filename1, "test/test-no-hlt.root");
    sprintf(filename2, "test/test-with-hlt.root");
    break;
      
  case 7:
    sprintf(filename1, "test/test-tk-no-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt.root");
    break;
      
  case 8:
    sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_cscs/output_6.root");
    sprintf(filename2, "grid/analysed_Spring07_withHLT_new/Bs2MuMu_fnal/output_6.root");
    start = 0; stop = 10000;
    break;

  case 9:
    sprintf(filename1, "test/test-tk-no-hlt-byChi2.root");
    sprintf(filename2, "test/test-tk-no-hlt-byHits.root");
    break;

  case 10:
    sprintf(filename1, "test/test-tk-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    start = 1500; stop = 3000;
    break;

  case 11:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-with-hlt-bad_1.root");
    start = 1500; stop = 3000;
    break;

  case 12:
    sprintf(filename1, "test/test-tk-no-hlt.root");
    sprintf(filename2, "test/test-tk-no-hlt-byHits.root");
    break;

  case 13:
    sprintf(filename1, "test/test-tk-with-hlt.root");
    sprintf(filename2, "test/test-tk-with-hlt-bad_1.root");
    start = 3000; stop = 4500; offset = 3000;
    break;

  case 14:
    sprintf(filename1, "test/test-with-hlt.root");
    sprintf(filename2, "test/test-with-hlt-bad_1.root");
    start = 3000; stop = 4500; offset = 3000;
    break;
      
  case 15:
    sprintf(filename1, "grid/analysed_Spring07_withHLT_new/Bs2MuMu/output_6.root");
    sprintf(filename2, "grid/analysed_Spring07_noHLT_new/Bs2MuMu/output_2.root");
    start = 2541; stop = 4041; offset = -500;
    break;

  default: 
    cout << " Enter a number < 16 !!!" << endl;
    return;
    
  };


  // -- Correlation plots (pT, mass)
  TH2D *h4 = new TH2D("h4", "Stand-alone muons (two leading #mu)", 210, -1., 20., 210, -1., 20.);
  TH2D *h5 = new TH2D("h5", "Tracker muons (two leading #mu)", 210, -1., 20., 210, -1., 20.);
  TH2D *h6 = new TH2D("h6", "Combined/global muons (two leading #mu)", 210, -1., 20., 210, -1., 20.);
  TH2D *h7 = new TH2D("h7", "L1 muons (two leading #mu)", 210, -1., 20., 210, -1., 20.);

  h4->GetXaxis()->SetTitle("p_{T, std}^{1} [GeV]"); h4->GetYaxis()->SetTitle("p_{T, std}^{2} [GeV]");
  h5->GetXaxis()->SetTitle("p_{T, trk}^{1} [GeV]"); h5->GetYaxis()->SetTitle("p_{T, trk}^{2} [GeV]");
  h6->GetXaxis()->SetTitle("p_{T, cmb}^{1} [GeV]"); h6->GetYaxis()->SetTitle("p_{T, cmb}^{2} [GeV]");
  h7->GetXaxis()->SetTitle("p_{T, l1}^{1} [GeV]"); h7->GetYaxis()->SetTitle("p_{T, l1}^{2} [GeV]");

  TH2D *h03 = new TH2D("h03", " ", 20, 0.0, 20., 20, 0.0, 20.);
  h03->GetXaxis()->SetTitle("N_{1}^{glb #mu}"); h03->GetYaxis()->SetTitle("N_{2}^{glb #mu}");

  TChain chain("T1");
  chain.Add(Form("%s", filename1));
  
  TChain chain2("T1");
  chain2.Add(Form("%s", filename2));
  
  // -- Set up for reading
  Int_t nentries(0), nentries2(0), nb(0), nb2(0);
  Int_t iEvent(0), jEvent(0), it(0);
  
  delete pEvent; 
  delete pEvent2;

  pEvent = new TAna00Event(0);
  pEvent2 = new TAna00Event(0);
  
  TAnaTrack *pTrack, *pTrack2;
  
  chain.SetBranchAddress("TAna00Event", &pEvent);
  nentries = chain.GetEntries();
  
  chain2.SetBranchAddress("TAna00Event", &pEvent2);
  nentries2 = chain2.GetEntries();
  
  cout << "Found " << nentries << " entries in the 1. chain" << endl;
  cout << "Found " << nentries2 << " entries in the 2. chain" << endl;
  
  int ella(0), tmp(0), match(0), nmatch(0);
  // for (iEvent = 0; iEvent < nentries; iEvent++) {
  for (iEvent = start; iEvent < stop; iEvent++) {
    
    if ( int((iEvent-start)/100) > ella ) {
      cout << "... events " << iEvent << " of " << nentries << " (matches " << nmatch << ")" << endl;
      ella++;
    }

    pEvent->Clear();    
    nb += chain.GetEvent(iEvent); 

    match = 0;

    pEvent2->Clear();
    tmp = chain2.GetEvent(iEvent+offset);

    if ( (pEvent->fRunNumber == pEvent2->fRunNumber)
	 && (pEvent->fEventNumber == pEvent2->fEventNumber) ) {
      
      pEvent2->Clear();
      nb2 += chain2.GetEvent(iEvent+offset);
      match = 1; nmatch++;

    } else {

      // for (jEvent = 0; jEvent < nentries2; jEvent++) {
      for (jEvent = start + offset; jEvent < stop + offset; jEvent++) {
	
	pEvent2->Clear();
	tmp = chain2.GetEvent(jEvent);
	
	if ( (pEvent->fRunNumber == pEvent2->fRunNumber)
	     && (pEvent->fEventNumber == pEvent2->fEventNumber) ) {
	  
	  pEvent2->Clear();
	  nb2 += chain2.GetEvent(jEvent);
	  match = 1; nmatch++;
	  break;
	}
      }
    }
   
    if ( match == 0 ) { 
      
      continue; 
    }
 
    //    cout << "Found " << pEvent->nGenCands() << " gen tracks in event" << endl;
    
    double pT1(0.), pT2(0.);
    double lead_pT1[4], lead_pT2[4];
    double scnd_pT1[4], scnd_pT2[4];

    int n1[4], n2[4];
    int lead_i1[4], lead_i2[4];  
    int scnd_i1[4], scnd_i2[4];  
  
    for (int i = 0; i < 4; i++) {

      n1[i] = n2[i] = 0;
      lead_i1[i] = lead_i2[i] = scnd_i1[i] = scnd_i2[i] = -1;
      lead_pT1[i] = lead_pT2[i] = scnd_pT1[i] = scnd_pT2[i] = 0.;
    }

    //h03->Fill(pEvent->nSigTracks(), pEvent2->nSigTracks());

    int idx(0);
    for (it = 0; it < pEvent->nSigTracks(); ++it) {
      
      pTrack = pEvent->getSigTrack(it);
      pT1    = pTrack->fPlab.Pt();
      idx    = pTrack->fMuType - 1;
      
      if ( idx > -1 ) {
	
	n1[idx]++;

	if (pT1 > lead_pT1[idx]) {

	  if (lead_i1[idx]  > -1) {

	    scnd_i1[idx]  = lead_i1[idx];
	    scnd_pT1[idx] = lead_pT1[idx];
	  }
	  
	  lead_i1[idx] = it;
	  lead_pT1[idx] = pT1;

	} else {

	  if (pT1 > scnd_pT1[idx]) {

	    scnd_i1[idx]  = it;
	    scnd_pT1[idx] = pT1;
	  }
	}
      }     
    }
      

    for (it = 0; it < pEvent2->nSigTracks(); ++it) {
      
      pTrack2 = pEvent2->getSigTrack(it);
      pT2     = pTrack2->fPlab.Pt();
      idx     = pTrack2->fMuType - 1;
      
      if ( idx > -1 ) {
	
	n2[idx]++;

	if (pT2 > lead_pT2[idx]) {

	  if (lead_i2[idx]  > -1) {

	    scnd_i2[idx]  = lead_i2[idx];
	    scnd_pT2[idx] = lead_pT2[idx];
	  }
	  
	  lead_i2[idx] = it;
	  lead_pT2[idx] = pT2;

	} else {

	  if (pT2 > scnd_pT2[idx]) {

	    scnd_i2[idx]  = it;
	    scnd_pT2[idx] = pT2;
	  }
	}
      }     
    }
    
    h03->Fill(n1[2], n2[2]);

    if ( lead_pT1[0] || lead_pT2[0] ) {  h4->Fill(lead_pT1[0], lead_pT2[0]); }
    if ( scnd_pT1[0] || scnd_pT2[0] ) {  h4->Fill(scnd_pT1[0], scnd_pT2[0]); }

    if ( lead_pT1[1] || lead_pT2[1] ) {  h5->Fill(lead_pT1[1], lead_pT2[1]); }
    if ( scnd_pT1[1] || scnd_pT2[1] ) {  h5->Fill(scnd_pT1[1], scnd_pT2[1]); }
    
    if ( lead_pT1[2] || lead_pT2[2] ) {  h6->Fill(lead_pT1[2], lead_pT2[2]); }
    if ( scnd_pT1[2] || scnd_pT2[2] ) {  h6->Fill(scnd_pT1[2], scnd_pT2[2]); }
    
    if ( lead_pT1[3] || lead_pT2[3] ) {  h7->Fill(lead_pT1[3], lead_pT2[3]); }
    if ( scnd_pT1[3] || scnd_pT2[3] ) {  h7->Fill(scnd_pT1[3], scnd_pT2[3]); }
    
  }
  

  // -- Display
  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.04);

  gStyle->SetOptStat(0);  

  c0.Clear(); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h03->Draw("colz");
  tl->DrawLatex(0.26, 0.86, Form("%s (2)",filename2));
  tl->DrawLatex(0.36, 0.80, Form("vs. %s (1)", filename1));
  c0.SaveAs(Form("sample%i_nmuons.pdf",ch));

  c0.Clear();
  c0.Divide(2,2);

  c0.cd(1); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h4->Draw("colz");
  c0.cd(2); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h5->Draw("colz");

  c0.cd(3); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h6->Draw("colz");
  c0.cd(4); shrinkPad(0.15, 0.1, 0.15, 0.1); gPad->SetLogz(1);
  h7->Draw("colz");

  c0.cd(1);
  tl->DrawLatex(0.26, 0.86, Form("%s (2)",filename2));
  tl->DrawLatex(0.36, 0.80, Form("vs. %s (1)", filename1));

  c0.SaveAs(Form("sample%i_muons.pdf",ch));

  delete h4; 
  delete h5; 
  delete h6; 
  delete h7;
  delete h03;
  delete tl;

}


// ======================================================================
#if 1


TCanvas c0("c0", "", 615,0,656,700);

// ----------------------------------------------------------------------
void setTitles(TH1 *h, const char *sx, const char *sy, float size = 0.05, 
               float xoff = 1.2, float yoff = 1.3, float lsize = 0.05, int font = 132) {
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

// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}

// ----------------------------------------------------------------------
void setFilledHist(TH1 *h, int lcol = kBlack, int fcol = kYellow, int fstyle = 1000, int width = 2) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(lcol);   h->SetLineWidth(width);   
  h->SetFillStyle(fstyle); h->SetFillColor(fcol);
}



// ----------------------------------------------------------------------
void shrinkPad(double b = 0.15, double l = 0.15, double r = 0.1, double t = 0.1) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}
