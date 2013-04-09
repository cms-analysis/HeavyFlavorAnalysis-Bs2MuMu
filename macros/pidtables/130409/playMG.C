
// ----------------------------------------------------------------------
void mkTH2D(const char *ofile = "mg-0.root") {

  string dir = "/shome/ursl/mg/all_and_pass_histograms";

  std::vector<string> samples; 
  samples.push_back("pion"); 
  samples.push_back("kaon"); 

  std::vector<string> region; 
  region.push_back("barrel"); 
  region.push_back("endcap"); 

  std::vector<string> etarange; 
  etarange.push_back("0_1.4"); 
  etarange.push_back("1.4_2.4"); 
  
  TFile *fout = TFile::Open(ofile, "RECREATE"); 

  TFile *f0all(0), *f1all(0); 
  TFile *f0pass(0), *f1pass(0); 
  TH1D *h0all(0), *h1all(0), *h0pass(0), *h1pass(0);

  string name1, name2; 
  for (int i = 0; i < samples.size(); ++i) {
 
    f0all = TFile::Open(Form("%s/hAllGenParticleGenPt_barrel_%s.root", dir.c_str(), samples[i].c_str())); 
    f0pass = TFile::Open(Form("%s/hAssocGenParticleGenPt_barrel_%s.root", dir.c_str(), samples[i].c_str())); 
    
    f1all = TFile::Open(Form("%s/hAllGenParticleGenPt_endcap_%s.root", dir.c_str(), samples[i].c_str())); 
    f1pass = TFile::Open(Form("%s/hAssocGenParticleGenPt_endcap_%s.root", dir.c_str(), samples[i].c_str())); 
    
    name1 = Form("hAllGenParticleGenPt_Eta_0_1.4"); 
    name2 = Form("hAssocGenParticleGenPt_Eta_0_1.4");
    
    h0all   = (TH1D*)((TH1D*)f0all->Get(name1.c_str()))->Clone(Form("%s_%s", name1.c_str(), samples[i].c_str())); 
    h0pass  = (TH1D*)((TH1D*)f0pass->Get(name2.c_str()))->Clone(Form("%s_%s", name2.c_str(), samples[i].c_str())); 

    name1 = Form("hAllGenParticleGenPt_Eta_1.4_2.4"); 
    name2 = Form("hAssocGenParticleGenPt_Eta_1.4_2.4");
    
    h1all   = (TH1D*)((TH1D*)f1all->Get(name1.c_str()))->Clone(Form("%s_%s", name1.c_str(), samples[i].c_str())); 
    h1pass  = (TH1D*)((TH1D*)f1pass->Get(name2.c_str()))->Clone(Form("%s_%s", name2.c_str(), samples[i].c_str())); 
    
    h0all->SetDirectory(fout);
    h0pass->SetDirectory(fout);
    h1all->SetDirectory(fout);
    h1pass->SetDirectory(fout);

    fout->cd();
    TH2D *h2all = new TH2D(Form("all_%s", samples[i].c_str()), Form("all_%s", samples[i].c_str()), 
			   2, 0., 2.8, h0all->GetNbinsX(), h0all->GetBinLowEdge(1), h0all->GetBinLowEdge(h0all->GetNbinsX()+1)); h2all->Sumw2(); 
    TH2D *h2pass = new TH2D(Form("pass_%s", samples[i].c_str()), Form("pass_%s", samples[i].c_str()),
			    2, 0., 2.8, h0all->GetNbinsX(), h0all->GetBinLowEdge(1), h0all->GetBinLowEdge(h0all->GetNbinsX()+1)); h2pass->Sumw2(); 
    
    for (int ibin = 1; ibin < h0all->GetNbinsX(); ++ibin) {
      h2all->SetBinContent(1, ibin, h0all->GetBinContent(ibin)); h2all->SetBinError(1, ibin, h0all->GetBinError(ibin)); 
      h2all->SetBinContent(2, ibin, h1all->GetBinContent(ibin)); h2all->SetBinError(2, ibin, h1all->GetBinError(ibin)); 

      h2pass->SetBinContent(1, ibin, h0pass->GetBinContent(ibin)); h2pass->SetBinError(1, ibin, h0pass->GetBinError(ibin)); 
      h2pass->SetBinContent(2, ibin, h1pass->GetBinContent(ibin)); h2pass->SetBinError(2, ibin, h1pass->GetBinError(ibin)); 
    }      

    zone(3,2); 
    h0all->Draw();
    c0->cd(2);
    h1all->Draw();
    c0->cd(3);
    h2all->Draw("colz");

    c0->cd(4);
    h0pass->Draw();
    c0->cd(5);
    h1pass->Draw();
    c0->cd(6);
    h2pass->Draw("colz");

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

  a.flush(); 
  b.flush(); 
  b.readFromHist(f1, "pass_pion", "all_pion"); 
  a.fillEff(b); 
  a.dumpToFile("pionFakeRate-tightMuon.dat"); 
  mean = a.getMean(); 
  sf = 0.311; 
  a.scale(sf); 
  a.dumpToFile("pionFakeRate-mvaTmMuon.dat"); 
  

  a.flush(); 
  b.flush(); 
  b.readFromHist(f1, "pass_kaon", "all_kaon"); 
  a.fillEff(b); 
  a.dumpToFile("kaonFakeRate-tightMuon.dat"); 
  sf = 0.311; 
  a.scale(sf); 
  a.dumpToFile("kaonFakeRate-mvaTmMuon.dat"); 

  a.flush(); 
  b.flush(); 
  b.readFromHist(f1, "pass_pion", "all_pion"); 
  a.fillEff(b); 
  mean = a.getMean(); 
  sf = 0.00015/mean; 
  a.scale(sf); 
  a.dumpToFile("protonFakeRate-tightMuon.dat"); 

  sf = 0.311;
  a.scale(sf); 
  a.dumpToFile("protonFakeRate-mvaTmMuon.dat"); 
  
}
