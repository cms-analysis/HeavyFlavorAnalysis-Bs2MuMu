// ----------------------------------------------------------------------
void parameters(int ch=0) {

  vector<string> dolist; vector<int> color; vector<int> symbol; vector<int> symbol; vector<string> legend; 
  dolist.push_back(Form("TMVA-%d-combined.root", 20+ch)); color.push_back(kBlack);  symbol.push_back(20); legend.push_back("TMVA default");
  dolist.push_back(Form("TMVA-%d-combined.root", 104+ch)); color.push_back(kCyan);  symbol.push_back(23); legend.push_back("nNodesMax=5");
  dolist.push_back(Form("TMVA-%d-combined.root", 106+ch)); color.push_back(kBlue);  symbol.push_back(24); legend.push_back("+nTrees=800+maxDepth=2");
  dolist.push_back(Form("TMVA-%d-combined.root", 108+ch)); color.push_back(kRed);   symbol.push_back(25); legend.push_back("this analysis");

  shrinkPad(0.1, 0.15); 
  TH2F* frame(0); 
  if (0 == ch) {
    frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.65, 100, 0.999, 1.0001);
  } else {
    frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.5, 100, 0.999, 1.0001);
  }

  frame->GetXaxis()->SetTitle(" #epsilon_{S}");
  frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");
  frame->Draw();  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  double x(0.15), y0(0.4), y(0);
  float rocInt(0.); 
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  tl->DrawLatex(x, y0+0.05, "Setup");
  tl->DrawLatex(x+0.45, y0+0.05, "integral");
  tl->SetTextSize(0.04);
  for (unsigned i = 0; i < dolist.size(); ++i) {
    TFile *f = TFile::Open(dolist[i].c_str()); 
    TGraph *g = (TGraph*)f->Get("groc"); 
    sscanf(g->GetTitle(), "integral = %f", &rocInt);
    g->SetMarkerStyle(20+i);
    g->SetMarkerSize(1);
    g->SetMarkerColor(color[i]); 
    g->SetLineColor(color[i]); 
    g->SetLineWidth(2); 
    g->Draw("l"); 
    y = y0 - i*0.05;
    tl->SetTextColor(color[i]); 
    tl->DrawLatex(x, y, legend[i].c_str());
    tl->DrawLatex(x+0.5, y, Form("(%4.3f)", rocInt));
  }

  c0->SaveAs(Form("overlayROC-parameters-%d.pdf", ch));
}

// ----------------------------------------------------------------------
void variables(int ch=0) {

  vector<string> dolist; vector<int> color; vector<int> symbol; vector<int> symbol; vector<string> legend; 
  dolist.push_back(Form("TMVA-%d-combined.root", 24+ch)); color.push_back(kCyan);  symbol.push_back(24); legend.push_back("pvips docatrk");
  dolist.push_back(Form("TMVA-%d-combined.root", 22+ch)); color.push_back(kBlue);  symbol.push_back(23); legend.push_back("fls3d alpha");
  dolist.push_back(Form("TMVA-%d-combined.root", 26+ch)); color.push_back(20);   symbol.push_back(25); legend.push_back(" + pvips docatrk");
  dolist.push_back(Form("TMVA-%d-combined.root", 28+ch)); color.push_back(25);     symbol.push_back(21); legend.push_back(" + pvip iso closetrk");
  dolist.push_back(Form("TMVA-%d-combined.root", 30+ch)); color.push_back(30);     symbol.push_back(21); legend.push_back(" + maxdoca");
  dolist.push_back(Form("TMVA-%d-combined.root", 32+ch)); color.push_back(35);     symbol.push_back(21); legend.push_back(" +  pt eta");
  dolist.push_back(Form("TMVA-%d-combined.root", 20+ch)); color.push_back(kBlack); symbol.push_back(20); legend.push_back("all ");
  dolist.push_back(Form("TMVA-%d-combined.root", 108+ch)); color.push_back(kRed); symbol.push_back(20); legend.push_back("this analysis ");

  TH2F* frame(0); 
  if (0 == ch) {
    frame = new TH2F("frame", "BDT output distributions", 100, 0.6, 0.9, 100, 0.8, 1.01);
  } else {
    frame = new TH2F("frame", "BDT output distributions", 100, 0.5, 0.75, 100, 0.8, 1.01);
  }
  frame->GetXaxis()->SetTitle(" #epsilon_{S}");
  frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");
  frame->Draw();  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  double x(0.15), y0(0.5), y(0);
  float rocInt(0.); 
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  tl->SetTextSize(0.05);
  tl->DrawLatex(x, y0+0.05, "Variables");
  tl->DrawLatex(x+0.25, y0+0.05, "integral");
  tl->SetTextSize(0.04);
  for (unsigned i = 0; i < dolist.size(); ++i) {
    TFile *f = TFile::Open(dolist[i].c_str()); 
    TGraph *g = (TGraph*)f->Get("groc"); 
    sscanf(g->GetTitle(), "integral = %f", &rocInt);
    g->SetMarkerStyle(20+i);
    g->SetMarkerSize(1);
    g->SetMarkerColor(color[i]); 
    g->SetLineColor(color[i]); 
    g->SetLineWidth(2); 
    g->Draw("l"); 
    y = y0 - i*0.05;
    tl->SetTextColor(color[i]); 
    tl->DrawLatex(x, y, legend[i].c_str());
    tl->DrawLatex(x+0.3, y, Form("(%4.3f)", rocInt));
  }


  c0->SaveAs(Form("overlayROC-variables-%d.pdf", ch));
}


// ----------------------------------------------------------------------
void preselections(int ch=0) {

  vector<string> dolist; vector<int> color; vector<int> symbol; vector<int> symbol; vector<string> legend; 
  dolist.push_back(Form("TMVA-%d-combined.root", 10+ch)); color.push_back(kBlack);  symbol.push_back(24); legend.push_back("tight");
  dolist.push_back(Form("TMVA-%d-combined.root", 20+ch)); color.push_back(kBlue);  symbol.push_back(23); legend.push_back("tighter");
  dolist.push_back(Form("TMVA-%d-combined.root", 108+ch)); color.push_back(kRed); symbol.push_back(20); legend.push_back("this analysis ");

  TH2F* frame(0); 
  if (0 == ch) {
    frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.6, 100, 0.9995, 1.0001);
  } else {
    frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.5, 100, 0.9995, 1.0001);
  }
  frame->GetXaxis()->SetTitle(" #epsilon_{S}");
  frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");
  frame->Draw();  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  double x(0.20), y0(0.5), y(0);
  float rocInt(0.); 
  tl->SetTextColor(kBlack);
  tl->SetTextSize(0.05);
  tl->SetTextSize(0.05);
  tl->DrawLatex(x, y0+0.05, "Preselection");
  tl->DrawLatex(x+0.35, y0+0.05, "integral");
  tl->SetTextSize(0.04);
  for (unsigned i = 0; i < dolist.size(); ++i) {
    TFile *f = TFile::Open(dolist[i].c_str()); 
    TGraph *g = (TGraph*)f->Get("groc"); 
    sscanf(g->GetTitle(), "integral = %f", &rocInt);
    g->SetMarkerStyle(20+i);
    g->SetMarkerSize(1);
    g->SetMarkerColor(color[i]); 
    g->SetLineColor(color[i]); 
    g->SetLineWidth(2); 
    g->Draw("l"); 
    y = y0 - i*0.05;
    tl->SetTextColor(color[i]); 
    tl->DrawLatex(x, y, legend[i].c_str());
    tl->DrawLatex(x+0.4, y, Form("(%4.3f)", rocInt));
  }


  c0->SaveAs(Form("overlayROC-preselection-%d.pdf", ch));
}
