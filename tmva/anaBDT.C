// ----------------------------------------------------------------------
void makeAll(const char *location="tmva2", int imax = 10) {
  fillHistograms(location, imax);
  plotProfiles(location); 
}


// ----------------------------------------------------------------------
void plotProfiles(const char *location="tmva1") {

  TFile *f = TFile::Open(Form("output-%s.root", location)); 

  double miny(0.), maxy(2.5); 
  
  hpl(p_ntree_SSB1, Form("(0.:1000.,%f:%f)", miny, maxy));
  c0.SaveAs(Form("%s_ntree_prof.pdf", location));

  hpl(hNTrees, "fillyellow"); 
  c0.SaveAs(Form("%s_ntree_hist.pdf", location));

  hpl(p_ncuts_SSB1, Form("(0.:500.,%f:%f)", miny, maxy));
  c0.SaveAs(Form("%s_ncuts_prof.pdf", location));
  hpl(hnCuts, "fillyellow"); 
  c0.SaveAs(Form("%s_ncuts_hist.pdf", location));

  hpl(p_nevts_SSB1, Form("(0.:1000.,%f:%f)", miny, maxy));
  c0.SaveAs(Form("%s_nevts_prof.pdf", location));
  hpl(hNevts, "fillyellow"); 
  c0.SaveAs(Form("%s_nevts_hist.pdf", location));

  hpl(p_depth_SSB1, Form("(0.:11.,%f:%f)", miny, maxy));
  c0.SaveAs(Form("%s_depth_prof.pdf", location));
  hpl(hMaxDepth, "fillyellow");
  c0.SaveAs(Form("%s_depth_hist.pdf", location));

  hpl(p_nodes_SSB1, Form("(1.e4:1.e6,%f:%f)", miny, maxy));
  c0.SaveAs(Form("%s_nodes_prof.pdf", location));
  hpl(hnodesmax, "fillyellow");
  c0.SaveAs(Form("%s_nodes_hist.pdf", location));

  hpl(p_beta_SSB1, Form("(0.:1.,%f:%f)", miny, maxy));
  c0.SaveAs(Form("%s_beta_prof.pdf", location));
  hpl(hbeta, "fillyellow");
  c0.SaveAs(Form("%s_beta_hist.pdf", location));
}


// ----------------------------------------------------------------------
TH1D* findParHistogram(const char *name, int n, TFile *f) {
  TH1D *h1 = (TH1D*)f->Get(Form("%s_%d", name, n)); 
  int nbins(60), lo(0.), hi(3.); 
  if (0 == h1) {
    f->cd();
    h1 = new TH1D(Form("%s_%d", name, n), Form("%s_%d", name, n), nbins, lo, hi);
    h1->SetDirectory(f);
  }
  return h1;
}

// ----------------------------------------------------------------------
void fillHistograms(const char *location = "tmva1", int imax = 10) {
  zone(1);
  string based = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/ursl/tmva"; 

  TFile *fout = TFile::Open(Form("output-%s.root", location), "RECREATE");
  TFile *f; 
  TH1D *h, *hssb, *hProb; 
  TH1D *hNTrees   = new TH1D("hNTrees", "NTrees", 100, 0., 1000.); 
  TH1D *hNevts    = new TH1D("hNevts", "nEvtsMin", 100, 0., 1000.); 
  TH1D *hMaxDepth = new TH1D("hMaxDepth", "MaxDepth", 15, 0., 15.); 
  TH1D *hnCuts    = new TH1D("hnCuts", "nCuts", 101, 0., 101.); 
  TH1D *hbeta     = new TH1D("hbeta", "beta", 10, 0., 1.); 
  TH1D *hnodesmax = new TH1D("hnodesmax", "nodesmax", 100, 1.e4, 1.e6); 

  TH1D *hSSB      = new TH1D("hSSB", "S/sqrt(S+B)", 4000, 0., 4000.); 
  TH1D *hSSB1     = new TH1D("hSSB1", "S/sqrt(S+B), KS", 4000, 0., 4000.); hSSB1->SetLineColor(kBlue); 

  TProfile *p_ntree_SSB0 = new TProfile("p_ntree_SSB0", "SSB0  vs ntree", 50, 0., 1000., "S"); 
  TProfile *p_ntree_SSB1 = new TProfile("p_ntree_SSB1", "SSB1  vs ntree", 50, 0., 1000., "S"); 
  TProfile *p_ntree_KSSG = new TProfile("p_ntree_KSSG", "KS Sg vs ntree", 50, 0., 1000., "S"); 
  TProfile *p_ntree_KSBG = new TProfile("p_ntree_KSBG", "KS Bg vs ntree", 50, 0., 1000., "S"); 

  TProfile *p_nevts_SSB0 = new TProfile("p_nevts_SSB0", "SSB0  vs nevts", 50, 0., 1000., "S"); 
  TProfile *p_nevts_SSB1 = new TProfile("p_nevts_SSB1", "SSB1  vs nevts", 50, 0., 1000., "S"); 
  TProfile *p_nevts_KSSG = new TProfile("p_nevts_KSSG", "KS Sg vs nevts", 50, 0., 1000., "S"); 
  TProfile *p_nevts_KSBG = new TProfile("p_nevts_KSBG", "KS Bg vs nevts", 50, 0., 1000., "S"); 

  TProfile *p_depth_SSB0 = new TProfile("p_depth_SSB0", "SSB0  vs depth", 15, 0., 15., "S"); 
  TProfile *p_depth_SSB1 = new TProfile("p_depth_SSB1", "SSB1  vs depth", 15, 0., 15., "S"); 
  TProfile *p_depth_KSSG = new TProfile("p_depth_KSSG", "KS Sg vs depth", 15, 0., 15., "S"); 
  TProfile *p_depth_KSBG = new TProfile("p_depth_KSBG", "KS Bg vs depth", 15, 0., 15., "S"); 

  TProfile *p_ncuts_SSB0 = new TProfile("p_ncuts_SSB0", "SSB0  vs ncuts", 50, 0., 100., "S"); 
  TProfile *p_ncuts_SSB1 = new TProfile("p_ncuts_SSB1", "SSB1  vs ncuts", 50, 0., 100., "S"); 
  TProfile *p_ncuts_KSSG = new TProfile("p_ncuts_KSSG", "KS Sg vs ncuts", 50, 0., 100., "S"); 
  TProfile *p_ncuts_KSBG = new TProfile("p_ncuts_KSBG", "KS Bg vs ncuts", 50, 0., 100., "S"); 

  TProfile *p_beta_SSB0 = new TProfile("p_beta_SSB0", "SSB0  vs beta", 10, 0., 1., "S"); 
  TProfile *p_beta_SSB1 = new TProfile("p_beta_SSB1", "SSB1  vs beta", 10, 0., 1., "S"); 
  TProfile *p_beta_KSSG = new TProfile("p_beta_KSSG", "KS Sg vs beta", 10, 0., 1., "S"); 
  TProfile *p_beta_KSBG = new TProfile("p_beta_KSBG", "KS Bg vs beta", 10, 0., 1., "S"); 

  TProfile *p_nodes_SSB0 = new TProfile("p_nodes_SSB0", "SSB0  vs nodes", 50, 1.e4, 1.e6, "S"); 
  TProfile *p_nodes_SSB1 = new TProfile("p_nodes_SSB1", "SSB1  vs nodes", 50, 1.e4, 1.e6, "S"); 
  TProfile *p_nodes_KSSG = new TProfile("p_nodes_KSSG", "KS Sg vs nodes", 50, 1.e4, 1.e6, "S"); 
  TProfile *p_nodes_KSBG = new TProfile("p_nodes_KSBG", "KS Bg vs nodes", 50, 1.e4, 1.e6, "S"); 

  int ntree, nevts, depth, ncuts, nnodesmax; 
  double beta; 
  double ksSg, ksBg; 
  int intrees, inevts, imaxdepth, incuts, ibeta, inodesmax; 
  int first(1), failedFiles(0); 
  for (int i = 0; i < imax; ++i) {
    cout << Form("%s/%s/root/TMVA-even-%d.root", based.c_str(), location, i) << endl;
    f = TFile::Open(Form("%s/%s/root/TMVA-even-%d.root", based.c_str(), location, i)); 
    if (!f) {
      cout << "error opening file! Skipping ... " << endl;
      ++failedFiles;
      continue;
    }
    // -- get parameter distribution
    h = (TH1D*)f->Get("hSetup");
    if (first) {
      first = 0; 
      for (int ibin = 1; ibin < h->GetNbinsX(); ++ibin) {
	if (!strcmp("NTrees", h->GetXaxis()->GetBinLabel(ibin))) intrees = ibin; 
	if (!strcmp("nEventsMin", h->GetXaxis()->GetBinLabel(ibin)))  inevts = ibin; 
	if (!strcmp("MaxDepth", h->GetXaxis()->GetBinLabel(ibin))) imaxdepth = ibin; 
	if (!strcmp("nCuts", h->GetXaxis()->GetBinLabel(ibin)))  incuts = ibin; 
	if (!strcmp("AdaBoostBeta", h->GetXaxis()->GetBinLabel(ibin))) ibeta = ibin; 
	if (!strcmp("NNodesMax", h->GetXaxis()->GetBinLabel(ibin))) inodesmax = ibin; 
      }
    }

    hProb = (TH1D*)f->Get("res_ks"); 
    ksSg = hProb->GetBinContent(1);
    ksBg = hProb->GetBinContent(2);

    // -- extract maximum S/SQRT(S+B)
    hssb = (TH1D*)f->Get("res_ssb");
    double ssbMax = hssb->GetBinContent(1); 

    hSSB->SetBinContent(i+1, ssbMax); 
    if (ksSg > 0.10 && ksBg > 0.10) {
      hSSB1->SetBinContent(i+1, ssbMax); 

      hNTrees->Fill(h->GetBinContent(intrees));   
      hNevts->Fill(h->GetBinContent(inevts));   
      hMaxDepth->Fill(h->GetBinContent(imaxdepth));   
      hnCuts->Fill(h->GetBinContent(incuts));   
      hbeta->Fill(h->GetBinContent(ibeta));   
      hnodesmax->Fill(h->GetBinContent(inodesmax));
    }
    if (ssbMax > 1.5) { 
      cout << " i = " << i << " ssbMax = " << ssbMax << endl;
      cout << "   ntrees   = " << h->GetBinContent(intrees) << endl;
      cout << "   nevts    = " << h->GetBinContent(inevts) << endl;
      cout << "   maxdepth = " << h->GetBinContent(imaxdepth) << endl;
      cout << "   ncuts    = " << h->GetBinContent(incuts) << endl;
      cout << "   beta     = " << h->GetBinContent(ibeta) << endl;
      cout << "   KS prob SG = " << ksSg << " BG = " << ksBg << endl;
    }

    // -- fill histograms per setting
    ntree = h->GetBinContent(intrees);
    p_ntree_KSSG->Fill(ntree, ksSg);
    p_ntree_KSBG->Fill(ntree, ksBg);
    p_ntree_SSB0->Fill(ntree, ssbMax);
    findParHistogram("ksSg_ntree", ntree, fout)->Fill(ksSg);
    findParHistogram("ksBg_ntree", ntree, fout)->Fill(ksBg);
    findParHistogram("ssb0_ntree", ntree, fout)->Fill(ssbMax);
    if (ksSg > 0.10 && ksBg > 0.10) {
      findParHistogram("ssb1_ntree", ntree, fout)->Fill(ssbMax);
      p_ntree_SSB1->Fill(ntree, ssbMax);
    }

    nevts = h->GetBinContent(inevts);
    p_nevts_KSSG->Fill(nevts, ksSg);
    p_nevts_KSBG->Fill(nevts, ksBg);
    p_nevts_SSB0->Fill(nevts, ssbMax);
    findParHistogram("ksSg_nevts", nevts, fout)->Fill(ksSg);
    findParHistogram("ksBg_nevts", nevts, fout)->Fill(ksBg);
    findParHistogram("ssb0_nevts", nevts, fout)->Fill(ssbMax);
    if (ksSg > 0.10 && ksBg > 0.10) {
      findParHistogram("ssb1_nevts", nevts, fout)->Fill(ssbMax);
      p_nevts_SSB1->Fill(nevts, ssbMax);
    }

    depth = h->GetBinContent(imaxdepth);
    p_depth_KSSG->Fill(depth, ksSg);
    p_depth_KSBG->Fill(depth, ksBg);
    p_depth_SSB0->Fill(depth, ssbMax);
    findParHistogram("ksSg_depth", depth, fout)->Fill(ksSg);
    findParHistogram("ksBg_depth", depth, fout)->Fill(ksBg);
    findParHistogram("ssb0_depth", depth, fout)->Fill(ssbMax);
    if (ksSg > 0.10 && ksBg > 0.10) {
      findParHistogram("ssb1_depth", depth, fout)->Fill(ssbMax);
      p_depth_SSB1->Fill(depth, ssbMax);
    }

    ncuts = h->GetBinContent(incuts);
    p_ncuts_KSSG->Fill(ncuts, ksSg);
    p_ncuts_KSBG->Fill(ncuts, ksBg);
    p_ncuts_SSB0->Fill(ncuts, ssbMax);
    findParHistogram("ksSg_ncuts", ncuts, fout)->Fill(ksSg);
    findParHistogram("ksBg_ncuts", ncuts, fout)->Fill(ksBg);
    findParHistogram("ssb0_ncuts", ncuts, fout)->Fill(ssbMax);
    if (ksSg > 0.10 && ksBg > 0.10) {
      findParHistogram("ssb1_ncuts", ncuts, fout)->Fill(ssbMax);
      p_ncuts_SSB1->Fill(ncuts, ssbMax);
    }

    beta = h->GetBinContent(ibeta);
    p_beta_KSSG->Fill(beta, ksSg);
    p_beta_KSBG->Fill(beta, ksBg);
    p_beta_SSB0->Fill(beta, ssbMax);
    findParHistogram("ksSg_beta", beta, fout)->Fill(ksSg);
    findParHistogram("ksBg_beta", beta, fout)->Fill(ksBg);
    findParHistogram("ssb0_beta", beta, fout)->Fill(ssbMax);
    if (ksSg > 0.10 && ksBg > 0.10) {
      findParHistogram("ssb1_beta", beta, fout)->Fill(ssbMax);
      p_beta_SSB1->Fill(beta, ssbMax);
    }

    nnodesmax = h->GetBinContent(inodesmax);
    p_nodes_KSSG->Fill(nnodesmax, ksSg);
    p_nodes_KSBG->Fill(nnodesmax, ksBg);
    p_nodes_SSB0->Fill(nnodesmax, ssbMax);
    findParHistogram("ksSg_nodes", nnodesmax, fout)->Fill(ksSg);
    findParHistogram("ksBg_nodes", nnodesmax, fout)->Fill(ksBg);
    findParHistogram("ssb0_nodes", nnodesmax, fout)->Fill(ssbMax);
    if (ksSg > 0.10 && ksBg > 0.10) {
      findParHistogram("ssb1_nodes", beta, fout)->Fill(ssbMax);
      p_nodes_SSB1->Fill(nnodesmax, ssbMax);
    }


    f->Close();
  }

  fout->Write();
  fout->Close();
  
  cout << "Failed to open " << failedFiles << " files " << endl;
}


// ----------------------------------------------------------------------
void histogramSSB(const char *location) {

  static TH1D *h1(0), *h2(0); 
  if (0 == h1) {h1 = new TH1D("h1", "h1", 100, 0., 2.5); setFilledHist(h1, kBlue, kBlue, 3356); setTitles(h1, "S/#sqrt{S+B}", "# BDTs");}
  if (0 == h2) {h2 = new TH1D("h2", "h2", 100, 0., 2.5); setFilledHist(h2, kRed,  kRed,  3356); }
  

  char directory[1000]; 
  system("/bin/rm -f d.plotDAC");
  //  system(Form("/bin/ls /shome/ursl/tmva/%s/log/TMVA-??.log > d.files", location));
  system(Form("/bin/ls /shome/ursl/tmva/%s/log/*.log > d.files", location));
  
  ifstream is("d.files");
  vector<string> lines;
  char  buffer[200];
  string fname; 
  double ssb; 
  int type; 
  while (is.getline(buffer, 200, '\n')) {
    fname = buffer;
    parseLogFile(fname, ssb, type); 
    if (type > 0) cout << " type = " << type << " ssb = " << ssb << " " << fname << endl;
    if (1 == type) h1->Fill(ssb);
    if (2 == type) h2->Fill(ssb);
  }
  is.close();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  h1->Draw();
  h2->Draw("samehist");

  tl->DrawLatex(0.25, 0.92, Form("sel: %4.3f", h2->GetMean())); 
  tl->DrawLatex(0.65, 0.92, Form("all: %4.3f", h1->GetMean())); 

  c0.SaveAs(Form("histogramSSB-%s.pdf", location)
 );
  
}

// ----------------------------------------------------------------------
double parseLogFile(string file, double &ssb, int &type) {
  ifstream is(file.c_str());
  char  buffer[1000];
  string sbuffer; 
  ssb = -1.;
  type = -1; 
  while (is.getline(buffer, 1000, '\n')) {
    sbuffer = buffer; 
    // ===> REMOVED TMVA-9.root, kssg = 0.780601 ksbg = 0.718995 SSB = 1.264315
    // ===> KEEP TMVA-17.root, kssg = 0.999343 ksbg = 0.681100 SSB = 1.338510
    
    //    cout << sbuffer << endl;
    if (string::npos != sbuffer.find("===> REMOVED TMVA")) {
      string::size_type m1 = sbuffer.find("SSB ="); 
      string stype = sbuffer.substr(m1+5, sbuffer.size()-m1); 
      type = 1; 
      ssb = atof(stype.c_str());
      //      cout << file << " removed, ssb = " << ssb << endl;
      break;
    }

    if (string::npos != sbuffer.find("===> KEEP TMVA")) {
      string::size_type m1 = sbuffer.find("SSB ="); 
      string stype = sbuffer.substr(m1+5, sbuffer.size()-m1); 
      type = 2; 
      ssb = atof(stype.c_str());
      //      cout << file << " keep, ssb = " << ssb << endl;
      break;
    }
  }  

  return ssb;
}
