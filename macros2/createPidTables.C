// ----------------------------------------------------------------------
void makeAll(const char *dir = "../macros/pidtables/130104") {

  // -- single file cases
  read1File(Form("%s/TnP_L3.root", dir), "",          Form("%s/L3_data_all.dat", dir)); 
  read1File(Form("%s/TnP_L3.root", dir), "_seagulls", Form("%s/L3_data_seagulls.dat", dir)); 
  read1File(Form("%s/TnP_L3.root", dir), "_cowboys",  Form("%s/L3_data_cowboys.dat", dir)); 
  
  read1File(Form("%s/TnP_L3.root", dir), "_mc",          Form("%s/L3_mc_all.dat", dir)); 
  read1File(Form("%s/TnP_L3.root", dir), "_seagulls_mc", Form("%s/L3_mc_seagulls.dat", dir)); 
  read1File(Form("%s/TnP_L3.root", dir), "_cowboys_mc",  Form("%s/L3_mc_cowboys.dat", dir)); 

  // -- two file cases
  read2Files(Form("%s/TnP_L1L2_Track2.root", dir), Form("%s/TnP_L1L2_Track7.root", dir), 
	    "",          Form("%s/L1L2_data_all.dat", dir)); 
  read2Files(Form("%s/TnP_L1L2_Track2.root", dir), Form("%s/TnP_L1L2_Track7.root", dir), 
	    "_seagulls", Form("%s/L1L2_data_seagulls.dat", dir)); 
  read2Files(Form("%s/TnP_L1L2_Track2.root", dir), Form("%s/TnP_L1L2_Track7.root", dir), 
	    "_cowboys",  Form("%s/L1L2_data_cowboys.dat", dir)); 

  read2Files(Form("%s/TnP_MuonID_Track2.root", dir), Form("%s/TnP_MuonID_Track7.root", dir), 
	    "",          Form("%s/MuonID_data_all.dat", dir)); 
  read2Files(Form("%s/TnP_MuonID_Track2.root", dir), Form("%s/TnP_MuonID_Track7.root", dir), 
	    "_seagulls", Form("%s/MuonID_data_seagulls.dat", dir)); 
  read2Files(Form("%s/TnP_MuonID_Track2.root", dir), Form("%s/TnP_MuonID_Track7.root", dir), 
	    "_cowboys",  Form("%s/MuonID_data_cowboys.dat", dir)); 

  read2Files(Form("%s/TnP_L1L2_Track2.root", dir), Form("%s/TnP_L1L2_Track7.root", dir), 
	    "_mc",          Form("%s/L1L2_mc_all.dat", dir)); 
  read2Files(Form("%s/TnP_L1L2_Track2.root", dir), Form("%s/TnP_L1L2_Track7.root", dir), 
	    "_seagulls_mc", Form("%s/L1L2_mc_seagulls.dat", dir)); 
  read2Files(Form("%s/TnP_L1L2_Track2.root", dir), Form("%s/TnP_L1L2_Track7.root", dir)
	    , "_cowboys_mc",  Form("%s/L1L2_mc_cowboys.dat", dir)); 

  read2Files(Form("%s/TnP_MuonID_Track2.root", dir), Form("%s/TnP_MuonID_Track7.root", dir), 
	    "_mc",          Form("%s/MuonID_mc_all.dat", dir)); 
  read2Files(Form("%s/TnP_MuonID_Track2.root", dir), Form("%s/TnP_MuonID_Track7.root", dir), 
	    "_seagulls_mc", Form("%s/MuonID_mc_seagulls.dat", dir)); 
  read2Files(Form("%s/TnP_MuonID_Track2.root", dir), Form("%s/TnP_MuonID_Track7.root", dir), 
	    "_cowboys_mc",  Form("%s/MuonID_mc_cowboys.dat", dir)); 

}



// ----------------------------------------------------------------------
void read1File(const char *fname = "TnP_L3.root", const char *pfix = "_seagulls", const char *oname = "L3_seagulls.dat") {

  TFile *f = TFile::Open(fname); 

  TH2D *hc = ((TH2D*)gFile->Get(Form("eff2D%s", pfix)))->Clone("hc"); 
  TH2D *hu = ((TH2D*)gFile->Get(Form("errhigh%s", pfix)))->Clone("hu"); 
  TH2D *hl = ((TH2D*)gFile->Get(Form("errlow%s", pfix)))->Clone("hl"); 

  TH2D *h1 = ((TH2D*)gFile->Get(Form("eff2D%s", pfix)))->Clone("h1"); 
  h1->Clear();
  h1->SetTitle("h1");

  cout << hc->GetEntries() << endl;

  zone(2,2);
  hc->DrawCopy("colz");

  double error(0.); 
  for (int ix = 1; ix <= hc->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hc->GetNbinsY(); ++iy) {
      error = 0.5*(hu->GetBinContent(ix, iy) - hl->GetBinContent(ix, iy));
      hc->SetBinError(ix, iy, error); 
      cout << hc->GetBinContent(ix, iy) << "+/-" << hc->GetBinError(ix, iy) << endl;
    }
  }

  PidTable a;
  // -- read in table from histogram
  a.readFromEffHist(gDirectory, "hc"); 
    
  // -- cast table into histogram
  a.eff2d(h1);
  c0->cd(2);
  h1->DrawCopy("colz");

  // -- subtract the original from the new histogam
  h1->Add(hc, -1.); 

  // -- this should be zero
  c0->cd(3);
  h1->DrawCopy("colz");
  
  a.shiftPmax(24., 500.);   
  a.shiftTmax(2.39, 2.5); 
  a.printAll();

  TString pname(oname); 
  a.dumpToFile(pname.Data()); 
}


// ----------------------------------------------------------------------
void read2Files(const char *f1name = "TnP_L1L2_Track2.root", const char *f2name = "TnP_L1L2_Track7.root", 
		const char *pfix = "_seagulls", const char *oname = "L1L2_seagulls.dat") {

  // -- file 1
  TFile *f1 = TFile::Open(f1name); 

  TH2D *hc1 = ((TH2D*)gFile->Get(Form("eff2D%s", pfix)))->Clone("hc1"); 
  TH2D *hu1 = ((TH2D*)gFile->Get(Form("errhigh%s", pfix)))->Clone("hu1"); 
  TH2D *hl1 = ((TH2D*)gFile->Get(Form("errlow%s", pfix)))->Clone("hl1"); 

  TH2D *h1 = ((TH2D*)gFile->Get(Form("eff2D%s", pfix)))->Clone("h1"); 
  h1->Clear();
  h1->SetTitle("h1");

  TArrayD *xbins1 = h1->GetXaxis()->GetXbins();
  int nx1 = xbins1->GetSize();
  TArrayD *ybins1 = h1->GetYaxis()->GetXbins();
  int ny1 = ybins1->GetSize();
  cout << "hist 1: " << nx1 << " " << ny1 << " -> " << xbins1->GetAt(1) << endl;

  cout << hc1->GetEntries() << endl;

  double error(0.); 
  for (int ix = 1; ix <= hc1->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hc1->GetNbinsY(); ++iy) {
      error = 0.5*(hu1->GetBinContent(ix, iy) - hl1->GetBinContent(ix, iy));
      hc1->SetBinError(ix, iy, error); 
      cout << hc1->GetBinContent(ix, iy) << "+/-" << hc1->GetBinError(ix, iy) << endl;
    }
  }

  zone(2,2);
  hc1->DrawCopy("colz");

  PidTable a;
  // -- read in table from histogram
  a.readFromEffHist(gDirectory, "hc1"); 



  // -- file 2
  TFile *f2 = TFile::Open(f2name); 

  TH2D *hc2 = ((TH2D*)gFile->Get(Form("eff2D%s", pfix)))->Clone("hc2"); 
  TH2D *hu2 = ((TH2D*)gFile->Get(Form("errhigh%s", pfix)))->Clone("hu2"); 
  TH2D *hl2 = ((TH2D*)gFile->Get(Form("errlow%s", pfix)))->Clone("hl2"); 

  TH2D *h2 = ((TH2D*)gFile->Get(Form("eff2D%s", pfix)))->Clone("h2"); 
  h2->Clear();
  h2->SetTitle("h2");

  cout << hc2->GetEntries() << endl;

  for (int ix = 1; ix <= hc2->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hc2->GetNbinsY(); ++iy) {
      error = 0.5*(hu2->GetBinContent(ix, iy) - hl2->GetBinContent(ix, iy));
      hc2->SetBinError(ix, iy, error); 
      cout << ix << " " << iy << " " << hc2->GetBinContent(ix, iy) << "+/-" << hc2->GetBinError(ix, iy) << endl;
    }
  }

  c0->cd(2);
  hc2->DrawCopy("colz");

  PidTable b;
  // -- read in table from histogram
  b.readFromEffHist(gDirectory, "hc2"); 

  TArrayD *xbins2 = h2->GetXaxis()->GetXbins();
  int nx2 = xbins2->GetSize();
  TArrayD *ybins2 = h2->GetYaxis()->GetXbins();
  int ny2 = ybins2->GetSize();
  cout << "hist 2: " << nx2 << " " << ny2 << " -> " << xbins2->GetAt(1) << endl;

  int nx(0), ny(0); 
  double xbins[100]; 
  double ybins[100]; 
  
  for (int i = 0; i < nx1; ++i) {
    xbins[i] = xbins1->GetAt(i); 
    ++nx;
  }
  xbins[nx] = 3.0;
  nx += 1;

  for (int i = 0; i < ny1; ++i) {
    ybins[i] = ybins1->GetAt(i); 
    ++ny;
  }

  for (int i = 1; i < ny2; ++i) {
    ybins[ny] = ybins2->GetAt(i); 
    ++ny;
  }
  ybins[ny] = 100.;
  ny += 1;

  for (int i = 0; i < nx; ++i) cout << "x " << i << " -> " << xbins[i] << endl;
  for (int i = 0; i < ny; ++i) cout << "y " << i << " -> " << ybins[i] << endl;

  TH2D *h3 = new TH2D("h3", "h3", nx-1, xbins, ny-1, ybins); 

  // -- combine b with a
  PidTable cc;
  cc.fillTable(a); 
  cc.fillTable(b); 
  cc.shiftPmax(24., 500.);   
  cc.shiftTmax(2.39, 2.5); 
  cc.printAll();
  
  cc.eff2d(h3);
  c0->cd(3);
  h3->DrawCopy("colz");
  
  
  TString pname(oname); 
  cc.dumpToFile(pname.Data()); 
}


