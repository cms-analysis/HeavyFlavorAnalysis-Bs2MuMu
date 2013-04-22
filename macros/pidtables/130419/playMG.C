
// ----------------------------------------------------------------------
void mkTH2D(const char *ofile = "mg-0.root") {

  string dir = "all_and_pass_histograms_tmva";

  std::vector<string> samples; 
  samples.push_back("pion"); 
  samples.push_back("kaon"); 

  std::vector<string> charge; 
  charge.push_back("Plus"); 
  charge.push_back("Minus"); 

  TFile *fout = TFile::Open(ofile, "RECREATE"); 

  TFile *f0all(0), *f1all(0); 
  TFile *f0pass(0), *f1pass(0); 
  TH1D *h0all(0), *h1all(0), *h0pass(0), *h1pass(0);

  string name1, name2; 
  for (int i = 0; i < samples.size(); ++i) {
    for (int iq = 0; iq < charge.size(); ++iq) {
    
      f0all = TFile::Open(Form("%s/hAllGenParticle%sGenPt_barrel_%s.root", dir.c_str(), charge[iq].c_str(), samples[i].c_str())); 
      f0pass = TFile::Open(Form("%s/hAssocGenParticle%sGenPt_barrel_%s.root", dir.c_str(), charge[iq].c_str(), samples[i].c_str())); 
      
      f1all = TFile::Open(Form("%s/hAllGenParticle%sGenPt_endcap_%s.root", dir.c_str(), charge[iq].c_str(), samples[i].c_str())); 
      f1pass = TFile::Open(Form("%s/hAssocGenParticle%sGenPt_endcap_%s.root", dir.c_str(), charge[iq].c_str(), samples[i].c_str())); 
      
      name1 = Form("hAllGenParticle%sGenPt_Eta_0_1.4", charge[iq].c_str()); 
      name2 = Form("hAssocGenParticle%sGenPt_Eta_0_1.4", charge[iq].c_str());
      
      h0all   = (TH1D*)((TH1D*)f0all->Get(name1.c_str()))->Clone(Form("%s_%s", name1.c_str(), samples[i].c_str())); 
      h0pass  = (TH1D*)((TH1D*)f0pass->Get(name2.c_str()))->Clone(Form("%s_%s", name2.c_str(), samples[i].c_str())); 
      
      name1 = Form("hAllGenParticle%sGenPt_Eta_1.4_2.4", charge[iq].c_str()); 
      name2 = Form("hAssocGenParticle%sGenPt_Eta_1.4_2.4", charge[iq].c_str());
      
      h1all   = (TH1D*)((TH1D*)f1all->Get(name1.c_str()))->Clone(Form("%s_%s", name1.c_str(), samples[i].c_str())); 
      h1pass  = (TH1D*)((TH1D*)f1pass->Get(name2.c_str()))->Clone(Form("%s_%s", name2.c_str(), samples[i].c_str())); 
      
      h0all->SetDirectory(fout);
      h0pass->SetDirectory(fout);
      h1all->SetDirectory(fout);
      h1pass->SetDirectory(fout);
      
      fout->cd();
      TH2D *h2all = new TH2D(Form("all_%s_%s", samples[i].c_str(), charge[iq].c_str() ), Form("all_%s_%s", samples[i].c_str(), charge[iq].c_str()), 
			     2, 0., 2.8, h0all->GetNbinsX(), h0all->GetBinLowEdge(1), h0all->GetBinLowEdge(h0all->GetNbinsX()+1)); h2all->Sumw2(); 
      TH2D *h2pass = new TH2D(Form("pass_%s_%s", samples[i].c_str(), charge[iq].c_str()), Form("pass_%s_%s", samples[i].c_str(), charge[iq].c_str()),
			      2, 0., 2.8, h0all->GetNbinsX(), h0all->GetBinLowEdge(1), h0all->GetBinLowEdge(h0all->GetNbinsX()+1)); h2pass->Sumw2(); 
      
      for (int ibin = 1; ibin <= h0all->GetNbinsX(); ++ibin) {
	h2all->SetBinContent(1, ibin, h0all->GetBinContent(ibin)); h2all->SetBinError(1, ibin, h0all->GetBinError(ibin)); 
	h2all->SetBinContent(2, ibin, h1all->GetBinContent(ibin)); h2all->SetBinError(2, ibin, h1all->GetBinError(ibin)); 
	
	h2pass->SetBinContent(1, ibin, h0pass->GetBinContent(ibin)); h2pass->SetBinError(1, ibin, h0pass->GetBinError(ibin)); 
	h2pass->SetBinContent(2, ibin, h1pass->GetBinContent(ibin)); h2pass->SetBinError(2, ibin, h1pass->GetBinError(ibin)); 
      }      
      
    }
  }
  
  fout->Write(); 
  fout->Close(); 
}

 
// ----------------------------------------------------------------------
void mkPidTables(const char *filename) {

  TFile *f1 = TFile::Open(filename);

  PidTable a("fakeTemplate.dat");
  PidTable b;

  double mean(0.), sf(0.); 

  std::vector<string> charge; 
  charge.push_back("Plus"); 
  charge.push_back("Minus"); 
  
  for (int iq = 0; iq < charge.size(); ++iq) {

    a.flush(); 
    b.clear(); 
    b.readFromHist(f1, Form("pass_pion_%s", charge[iq].c_str()), Form("all_pion_%s", charge[iq].c_str())); 
    a.fillEff(b); 
    a.dumpToFile(Form("pion%sFakeRate-mvaMuon.dat", charge[iq].c_str())); 
  
    a.flush(); 
    b.clear(); 
    b.readFromHist(f1, Form("pass_kaon_%s", charge[iq].c_str()), Form("all_kaon_%s", charge[iq].c_str())); 
    cout << " ==> b.printAll(): after read in" << endl;
    b.printAll(); 
    cout << " ==> a.printAll(): first time" << endl;
    a.fillEff(b); 
    cout << " ==> a.printAll(): after filling" << endl;
    a.printAll(); 
    a.dumpToFile(Form("kaon%sFakeRate-mvaMuon.dat", charge[iq].c_str())); 

    a.flush(); 
    b.clear(); 
    b.readFromHist(f1, Form("pass_pion_%s", charge[iq].c_str()), Form("all_pion_%s", charge[iq].c_str())); 
    a.fillEff(b); 
    mean = a.getMean(); 
    sf = 0.00015/mean; 
    a.scale(sf); 
    a.dumpToFile(Form("proton%sFakeRate-mvaMuon.dat", charge[iq].c_str())); 

  }
}


// ----------------------------------------------------------------------
void drawPidTables(const char *qual = "mvaMuon") {

  PidTable a;
  std::vector<string> charge; 
  charge.push_back("Plus"); 
  charge.push_back("Minus"); 

  std::vector<string> sample; 
  sample.push_back("pion"); 
  sample.push_back("kaon"); 
  sample.push_back("proton"); 

  std::vector<string> ssample; 
  ssample.push_back("#pi"); 
  ssample.push_back("K"); 
  ssample.push_back("p"); 

  string name; 

  double xbins[] = {0., 1.4, 2.4};
  double ybins[] = {0., 4., 5., 7., 10., 15., 20., 50.}; 
  TH2D *h2 = new TH2D("h2", "", 2, xbins, 7, ybins); 
  h2->SetMinimum(0.0); 
  h2->SetMaximum(0.001); 
  
  gStyle->SetOptStat(0); 

  for (int is = 0; is < sample.size(); ++is) {
    for (int iq = 0; iq < charge.size(); ++iq) {

      shrinkPad(0.15, 0.15, 0.2); 
      a.clear(); 
      name = Form("%s%sFakeRate-%s.dat", sample[is].c_str(), charge[iq].c_str(), qual); 
      h2->Reset();
      setTitles(h2, Form("|#eta(%s^{%s} )|", ssample[is].c_str(), (iq==0?"+":"-")), Form("p_{T}(%s^{%s} )", ssample[is].c_str(), (iq==0?"+":"-"))); 

      a.readFromFile(name.c_str()); 
      a.eff2d(h2); 
      h2->Draw("colztext");
      c0->SaveAs(Form("%s%sFakeRate-%s.pdf", sample[is].c_str(), charge[iq].c_str(), qual)); 
    }
  }


}
