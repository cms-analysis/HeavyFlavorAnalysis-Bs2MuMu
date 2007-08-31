// -- The big global
TAna00Event *pEvent;


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
  
  for (iEvent = 0; iEvent < nentries; iEvent++) {
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
