int kCstCol = 103;
void plot(const char *filName1="treebmm/csg-004.default.root"
	  , const char *filName2="treebmm/csg-005.default.root")
{

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetTitle(0);

	gStyle->SetStatFont(132);
	gStyle->SetTextFont(132);
	gStyle->SetLabelFont(132, "X");
	gStyle->SetLabelFont(132, "Y");
	gStyle->SetTitleFont(132);

	gROOT->ForceStyle();

	TLatex *tl = new TLatex;
	tl->SetNDC(kTRUE);
	tl->SetTextSize(0.08);	

	leggL = new TLegend(0.2,0.7,0.6,0.89);
	leggL->SetFillStyle(0); leggL->SetBorderSize(0); leggL->SetTextSize(0.05);  leggL->SetFillColor(0); 

	TFile *f  = new TFile(Form("%s", filName1));
	TFile *f2 = new TFile(Form("%s", filName2));

	TH1D *p1;
	TH1D *p2;

	TH2D *h1;
	TH2D *h2;


	TCanvas *c1 = new TCanvas("c1", "", 1200, 1000);
	c1->Clear();
	c1->Divide(2,2);

	TCanvas *c2 = new TCanvas("c2", "", 700, 900);
	c2->Clear();
	c2->Divide(1,2);

//---------------------------------------------------------------------------------- 
	c1->cd(1);

 	h1 = (TH2D*)f->Get("D121");
	h1->FitSlicesY();
	TH1D *h1_2 = (TH1D*)gDirectory->Get("D121_2"); 	
	h1_2->SetName("D121_2_1");
	TH1D *h1_0 = (TH1D*)gDirectory->Get("D121_0"); 	
	h1_0->SetName("D121_0_1");
	TH1D *h1_1 = (TH1D*)gDirectory->Get("D121_1"); 	
	h1_1->SetName("D121_1_1");
	TH1D *h1_chi2 = (TH1D*)gDirectory->Get("D121_chi2"); 	
	h1_chi2->SetName("D121_chi2_1");

 	h2 = (TH2D*)f2->Get("D121");
	h2->FitSlicesY();
	TH1D *h2_2 = (TH1D*)gDirectory->Get("D121_2"); 	
	h2_2->SetName("D121_2_2");
	TH1D *h2_0 = (TH1D*)gDirectory->Get("D121_0"); 	
	h2_0->SetName("D121_0_2");
	TH1D *h2_1 = (TH1D*)gDirectory->Get("D121_1"); 	
	h2_1->SetName("D121_1_2");
	TH1D *h2_chi2 = (TH1D*)gDirectory->Get("D121_chi2"); 	
	h2_chi2->SetName("D121_chi2_2");

  	setHist(h1_2, kBlue,  8, 0.9);
  	setHist(h2_2, kBlack, 24, 0.9);	

 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");


  	leggL->Clear();
	//	leggLe = leggL->AddEntry(h1_2, Form("perfect al. (B_{s} #rightarrow #mu #mu)"), "p"); leggLe->SetTextColor(kBlue);
    	leggLe = leggL->AddEntry(h1_2, Form("perfect al. (B^{+} #rightarrow J/#Psi K^{+})"), "p"); leggLe->SetTextColor(kBlue);
	leggLe = leggL->AddEntry(h2_2, Form("short-term al. (B^{+} #rightarrow J/#Psi K^{+})"), "p"); leggLe->SetTextColor(kBlack);

 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

	c2->cd(1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");
 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

//----------------------------------------------------------------------------------

	c1->cd(2);
	c1_2->Divide(3,2);


	c1_2->cd(4);
	h1->ProjectionX("h121_1_x");
	h2->ProjectionX("h121_2_x");

	setHist(h121_1_x, kBlue,  8, 0.4);
  	setHist(h121_2_x, kBlack, 24, 0.4);	

	h121_1_x->SetTitle("p_{T} distribution");
	h121_1_x->DrawCopy("E");  gPad->SetLogy(1); 
 	h121_2_x->DrawCopy("ESAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

	c1_2->cd(1);

	setHist(h1_0, kBlue,  8, 0.4);
  	setHist(h2_0, kBlack, 24, 0.4);	

	h1_0->GetYaxis()->SetRangeUser(0.1*h1_0->GetMinimum(1E-20), 10*h1_0->GetMaximum());
 	h1_0->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_0->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(2);

	setHist(h1_1, kBlue,  8, 0.4);
  	setHist(h2_1, kBlack, 24, 0.4);	

	h1_1->GetYaxis()->SetRangeUser(0.1*h1_1->GetMinimum(1E-20), 10*h1_1->GetMaximum());
 	h1_1->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_1->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(3);

	setHist(h1_chi2, kBlue,  8, 0.4);
  	setHist(h2_chi2, kBlack, 24, 0.4);	

	h1_chi2->GetYaxis()->SetRangeUser(0.1*h1_chi2->GetMinimum(1E-20), 10*h1_chi2->GetMaximum());
 	h1_chi2->DrawCopy("P"); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_chi2->DrawCopy("PSAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(5);
  	setHist(h1_2, kBlue,  8, 0.4);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->SetTitle(Form("%s vs. RMS",h1_2->GetTitle()));
	h1_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(6);
  	setHist(h2_2, kBlack, 24, 0.4);	
 	h2_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h2_2->SetTitle(Form("%s vs. RMS",h2_2->GetTitle()));
	h2_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");	

	// --------------- RMS ------------------------
 	p1 = (TH1D*)f->Get("p121");
  	p2 = (TH1D*)f2->Get("p121");
 	
 	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);
 	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);

  	setHist(hnew1, kBlue,  20, 0.4);
  	setHist(hnew2, kBlack, 20, 0.4);

	for (int i = 0; i < p1->GetNbinsX(); i++ ) {
	  
	  hnew1->SetBinContent(i, p1->GetBinError(i));
  	  hnew2->SetBinContent(i, p2->GetBinError(i));
	}

	c1_2->cd(5);
	hnew1->DrawCopy("HISTSAME");
	c1_2->cd(6);
  	hnew2->DrawCopy("HISTSAME");

	delete hnew1;
	delete hnew2;
	
//----------------------------------------------------------------------------------
 
 	h1 = (TH2D*)f->Get("D131");
	h1->FitSlicesY(0, -1, -1, 0);
	TH1D *h1_2 = (TH1D*)gDirectory->Get("D131_2"); 	
	h1_2->SetName("D131_2_1");
	TH1D *h1_0 = (TH1D*)gDirectory->Get("D131_0"); 	
	h1_0->SetName("D131_0_1");
	TH1D *h1_chi2 = (TH1D*)gDirectory->Get("D131_chi2"); 	
	h1_chi2->SetName("D131_chi2_1");

 	h2 = (TH2D*)f2->Get("D131");
	h2->FitSlicesY();
	TH1D *h2_2 = (TH1D*)gDirectory->Get("D131_2"); 	
	h1_2->SetName("D131_2_2");
	TH1D *h2_0 = (TH1D*)gDirectory->Get("D131_0"); 	
	h2_0->SetName("D131_0_2");
	TH1D *h2_chi2 = (TH1D*)gDirectory->Get("D131_chi2"); 	
	h2_chi2->SetName("D131_chi2_2");

  	setHist(h1_2, kBlue,  8, 0.9);
  	setHist(h2_2, kBlack, 24, 0.9);	

	c1->cd(3);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");

 	leggL->Draw();

	//	tl->DrawLatex(0.05, 0.83, "_{[%]}");
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");

	c2->cd(2);
 	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");

 	leggL->Draw();
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");
//----------------------------------------------------------------------------------

	c1->cd(4);

	c1_4->Divide(3,2);

	c1_4->cd(4);
	h1->ProjectionX("h131_1_x");
	h2->ProjectionX("h131_2_x");

	setHist(h131_1_x, kBlue,  8, 0.4);
  	setHist(h131_2_x, kBlack, 24, 0.4);	

	h131_1_x->SetTitle("#eta distribution");
	h131_1_x->DrawCopy("E");  gPad->SetLogy(1); 
 	h131_2_x->DrawCopy("ESAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");

	c1_4->cd(1);

	setHist(h1_0, kBlue,  8, 0.4);
  	setHist(h2_0, kBlack, 24, 0.4);	

	h1_0->GetYaxis()->SetRangeUser(0.1*h1_0->GetMinimum(1E-20), 10*h1_0->GetMaximum());
 	h1_0->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_0->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(2);

	setHist(h1_1, kBlue,  8, 0.4);
  	setHist(h2_1, kBlack, 24, 0.4);	

	h1_1->GetYaxis()->SetRangeUser(0.1*h1_1->GetMinimum(1E-20), 10*h1_1->GetMaximum());
 	h1_1->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_1->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(3);

	setHist(h1_chi2, kBlue,  8, 0.4);
  	setHist(h2_chi2, kBlack, 24, 0.4);	

	h1_chi2->GetYaxis()->SetRangeUser(0.1*h1_chi2->GetMinimum(1E-20), 10*h1_chi2->GetMaximum());
 	h1_chi2->DrawCopy("P"); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_chi2->DrawCopy("PSAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(5);
  	setHist(h1_2, kBlue,  8, 0.4);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->SetTitle(Form("%s vs. RMS",h1_2->GetTitle()));
	h1_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(6);
  	setHist(h2_2, kBlack, 24, 0.4);	
 	h2_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h2_2->SetTitle(Form("%s vs. RMS",h2_2->GetTitle()));
	h2_2->DrawCopy(); gPad->SetLogy(1);

	tl->DrawLatex(0.85, 0.05, "_{#eta}"); 	

	// --------------- RMS ------------------------
 	p1 = (TH1D*)f->Get("p131");
  	p2 = (TH1D*)f2->Get("p131");
 	
 	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 2.5);
 	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 2.5);

  	setHist(hnew1, kBlue,  20, 0.4);
  	setHist(hnew2, kBlack, 20, 0.4);

	for (int i = 0; i < p1->GetNbinsX(); i++ ) {
	  
	  hnew1->SetBinContent(i, p1->GetBinError(i));
  	  hnew2->SetBinContent(i, p2->GetBinError(i));
	}

	c1_4->cd(5);
	hnew1->DrawCopy("HISTSAME");
	c1_4->cd(6);
  	hnew2->DrawCopy("HISTSAME");

	delete hnew1;
	delete hnew2;

  	setHist(h1_2, kBlue,  8, 0.9);
  	setHist(h2_2, kBlack, 24, 0.9);	

 	c1->SaveAs("pT_res_fits.pdf");
 	c2->SaveAs("pT_res.pdf");

	delete c2;


}

//=================================================================================================
void plot2(const char *filName1="treebmm/csg-004.default.root"
	  , const char *filName2="treebmm/csg-003.default.root")
{

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetTitle(0);

	gStyle->SetStatFont(132);
	gStyle->SetTextFont(132);
	gStyle->SetLabelFont(132, "X");
	gStyle->SetLabelFont(132, "Y");
	gStyle->SetTitleFont(132);

	gROOT->ForceStyle();

	TLatex *tl = new TLatex;
	tl->SetNDC(kTRUE);
	tl->SetTextSize(0.08);	

	leggL = new TLegend(0.2,0.7,0.6,0.89);
	leggL->SetFillStyle(0); leggL->SetBorderSize(0); leggL->SetTextSize(0.05);  leggL->SetFillColor(0); 

	TFile *f  = new TFile(Form("%s", filName1));
	TFile *f2 = new TFile(Form("%s", filName2));

	TH1D *p1;
	TH1D *p2;

	TH2D *h1;
	TH2D *h2;


	TCanvas *c1 = new TCanvas("c1", "", 1200, 1000);
	c1->Clear();
	c1->Divide(2,2);

	TCanvas *c2 = new TCanvas("c2", "", 700, 900);
	c2->Clear();
	c2->Divide(1,2);

//---------------------------------------------------------------------------------- 
	c1->cd(1);

 	h1 = (TH2D*)f->Get("D121");
	h1->FitSlicesY();
	TH1D *h1_2 = (TH1D*)gDirectory->Get("D121_2"); 	
	h1_2->SetName("D121_2_1");
	TH1D *h1_0 = (TH1D*)gDirectory->Get("D121_0"); 	
	h1_0->SetName("D121_0_1");
	TH1D *h1_1 = (TH1D*)gDirectory->Get("D121_1"); 	
	h1_1->SetName("D121_1_1");
	TH1D *h1_chi2 = (TH1D*)gDirectory->Get("D121_chi2"); 	
	h1_chi2->SetName("D121_chi2_1");

 	h2 = (TH2D*)f2->Get("D121");
	h2->FitSlicesY();
	TH1D *h2_2 = (TH1D*)gDirectory->Get("D121_2"); 	
	h2_2->SetName("D121_2_2");
	TH1D *h2_0 = (TH1D*)gDirectory->Get("D121_0"); 	
	h2_0->SetName("D121_0_2");
	TH1D *h2_1 = (TH1D*)gDirectory->Get("D121_1"); 	
	h2_1->SetName("D121_1_2");
	TH1D *h2_chi2 = (TH1D*)gDirectory->Get("D121_chi2"); 	
	h2_chi2->SetName("D121_chi2_2");

  	setHist(h1_2, kBlue,  8, 0.9);
  	setHist(h2_2, kCstCol, 24, 0.9);	

 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");


  	leggL->Clear();
	//	leggLe = leggL->AddEntry(h1_2, Form("perfect al. (B_{s} #rightarrow #mu #mu)"), "p"); leggLe->SetTextColor(kBlue);
    	leggLe = leggL->AddEntry(h1_2, Form("perfect al. (B^{+} #rightarrow J/#Psi K^{+})"), "p"); leggLe->SetTextColor(kBlue);
	leggLe = leggL->AddEntry(h2_2, Form("short-term al. (B^{+} #rightarrow J/#Psi K^{+})"), "p"); leggLe->SetTextColor(kCstCol);

 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

	c2->cd(1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");
 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

//----------------------------------------------------------------------------------

	c1->cd(2);
	c1_2->Divide(3,2);


	c1_2->cd(4);
	h1->ProjectionX("h121_1_x");
	h2->ProjectionX("h121_2_x");

	setHist(h121_1_x, kBlue,  8, 0.4);
  	setHist(h121_2_x, kCstCol, 24, 0.4);	

	h121_1_x->SetTitle("p_{T} distribution");
	h121_1_x->DrawCopy("E");  gPad->SetLogy(1); 
 	h121_2_x->DrawCopy("ESAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

	c1_2->cd(1);

	setHist(h1_0, kBlue,  8, 0.4);
  	setHist(h2_0, kCstCol, 24, 0.4);	

	h1_0->GetYaxis()->SetRangeUser(0.1*h1_0->GetMinimum(1E-20), 10*h1_0->GetMaximum());
 	h1_0->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_0->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(2);

	setHist(h1_1, kBlue,  8, 0.4);
  	setHist(h2_1, kCstCol, 24, 0.4);	

	h1_1->GetYaxis()->SetRangeUser(0.1*h1_1->GetMinimum(1E-20), 10*h1_1->GetMaximum());
 	h1_1->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_1->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(3);

	setHist(h1_chi2, kBlue,  8, 0.4);
  	setHist(h2_chi2, kCstCol, 24, 0.4);	

	h1_chi2->GetYaxis()->SetRangeUser(0.1*h1_chi2->GetMinimum(1E-20), 10*h1_chi2->GetMaximum());
 	h1_chi2->DrawCopy("P"); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_chi2->DrawCopy("PSAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(5);
  	setHist(h1_2, kBlue,  8, 0.4);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->SetTitle(Form("%s vs. RMS",h1_2->GetTitle()));
	h1_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(6);
  	setHist(h2_2, kCstCol, 24, 0.4);	
 	h2_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h2_2->SetTitle(Form("%s vs. RMS",h2_2->GetTitle()));
	h2_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");	

	// --------------- RMS ------------------------
 	p1 = (TH1D*)f->Get("p121");
  	p2 = (TH1D*)f2->Get("p121");
 	
 	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);
 	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);

  	setHist(hnew1, kBlue,  20, 0.4);
  	setHist(hnew2, kCstCol, 20, 0.4);

	for (int i = 0; i < p1->GetNbinsX(); i++ ) {
	  
	  hnew1->SetBinContent(i, p1->GetBinError(i));
  	  hnew2->SetBinContent(i, p2->GetBinError(i));
	}

	c1_2->cd(5);
	hnew1->DrawCopy("HISTSAME");
	c1_2->cd(6);
  	hnew2->DrawCopy("HISTSAME");

	delete hnew1;
	delete hnew2;
	
//----------------------------------------------------------------------------------
 
 	h1 = (TH2D*)f->Get("D131");
	h1->FitSlicesY(0, -1, -1, 0);
	TH1D *h1_2 = (TH1D*)gDirectory->Get("D131_2"); 	
	h1_2->SetName("D131_2_1");
	TH1D *h1_0 = (TH1D*)gDirectory->Get("D131_0"); 	
	h1_0->SetName("D131_0_1");
	TH1D *h1_chi2 = (TH1D*)gDirectory->Get("D131_chi2"); 	
	h1_chi2->SetName("D131_chi2_1");

 	h2 = (TH2D*)f2->Get("D131");
	h2->FitSlicesY();
	TH1D *h2_2 = (TH1D*)gDirectory->Get("D131_2"); 	
	h1_2->SetName("D131_2_2");
	TH1D *h2_0 = (TH1D*)gDirectory->Get("D131_0"); 	
	h2_0->SetName("D131_0_2");
	TH1D *h2_chi2 = (TH1D*)gDirectory->Get("D131_chi2"); 	
	h2_chi2->SetName("D131_chi2_2");

  	setHist(h1_2, kBlue,  8, 0.9);
  	setHist(h2_2, kCstCol, 24, 0.9);	

	c1->cd(3);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");

 	leggL->Draw();

	//	tl->DrawLatex(0.05, 0.83, "_{[%]}");
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");

	c2->cd(2);
 	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");

 	leggL->Draw();
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");
//----------------------------------------------------------------------------------

	c1->cd(4);

	c1_4->Divide(3,2);

	c1_4->cd(4);
	h1->ProjectionX("h131_1_x");
	h2->ProjectionX("h131_2_x");

	setHist(h131_1_x, kBlue,  8, 0.4);
  	setHist(h131_2_x, kCstCol, 24, 0.4);	

	h131_1_x->SetTitle("#eta distribution");
	h131_1_x->DrawCopy("E");  gPad->SetLogy(1); 
 	h131_2_x->DrawCopy("ESAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");

	c1_4->cd(1);

	setHist(h1_0, kBlue,  8, 0.4);
  	setHist(h2_0, kCstCol, 24, 0.4);	

	h1_0->GetYaxis()->SetRangeUser(0.1*h1_0->GetMinimum(1E-20), 10*h1_0->GetMaximum());
 	h1_0->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_0->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(2);

	setHist(h1_1, kBlue,  8, 0.4);
  	setHist(h2_1, kCstCol, 24, 0.4);	

	h1_1->GetYaxis()->SetRangeUser(0.1*h1_1->GetMinimum(1E-20), 10*h1_1->GetMaximum());
 	h1_1->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_1->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(3);

	setHist(h1_chi2, kBlue,  8, 0.4);
  	setHist(h2_chi2, kCstCol, 24, 0.4);	

	h1_chi2->GetYaxis()->SetRangeUser(0.1*h1_chi2->GetMinimum(1E-20), 10*h1_chi2->GetMaximum());
 	h1_chi2->DrawCopy("P"); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_chi2->DrawCopy("PSAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(5);
  	setHist(h1_2, kBlue,  8, 0.4);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->SetTitle(Form("%s vs. RMS",h1_2->GetTitle()));
	h1_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(6);
  	setHist(h2_2, kCstCol, 24, 0.4);	
 	h2_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h2_2->SetTitle(Form("%s vs. RMS",h2_2->GetTitle()));
	h2_2->DrawCopy(); gPad->SetLogy(1);

	tl->DrawLatex(0.85, 0.05, "_{#eta}"); 	

	// --------------- RMS ------------------------
 	p1 = (TH1D*)f->Get("p131");
  	p2 = (TH1D*)f2->Get("p131");
 	
 	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 2.5);
 	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 2.5);

  	setHist(hnew1, kBlue,  20, 0.4);
  	setHist(hnew2, kCstCol, 20, 0.4);

	for (int i = 0; i < p1->GetNbinsX(); i++ ) {
	  
	  hnew1->SetBinContent(i, p1->GetBinError(i));
  	  hnew2->SetBinContent(i, p2->GetBinError(i));
	}

	c1_4->cd(5);
	hnew1->DrawCopy("HISTSAME");
	c1_4->cd(6);
  	hnew2->DrawCopy("HISTSAME");

	delete hnew1;
	delete hnew2;

  	setHist(h1_2, kBlue,  8, 0.9);
  	setHist(h2_2, kCstCol, 24, 0.9);	

 	c1->SaveAs("pT_res1_fits.pdf");
 	c2->SaveAs("pT_res1.pdf");

	delete c2;


}

// ======================================================================================

void plot3(const char *filName1="treebmm/csg-003.default.root"
	  , const char *filName2="treebmm/csg-005.default.root")
{

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetTitle(0);

	gStyle->SetStatFont(132);
	gStyle->SetTextFont(132);
	gStyle->SetLabelFont(132, "X");
	gStyle->SetLabelFont(132, "Y");
	gStyle->SetTitleFont(132);

	gROOT->ForceStyle();

	TLatex *tl = new TLatex;
	tl->SetNDC(kTRUE);
	tl->SetTextSize(0.08);	

	leggL = new TLegend(0.2,0.7,0.6,0.89);
	leggL->SetFillStyle(0); leggL->SetBorderSize(0); leggL->SetTextSize(0.05);  leggL->SetFillColor(0); 

	TFile *f  = new TFile(Form("%s", filName1));
	TFile *f2 = new TFile(Form("%s", filName2));

	TH1D *p1;
	TH1D *p2;

	TH2D *h1;
	TH2D *h2;

	TCanvas *c1 = new TCanvas("c1", "", 1200, 1000);
	c1->Clear();
	c1->Divide(2,2);

	TCanvas *c2 = new TCanvas("c2", "", 700, 900);
	c2->Clear();
	c2->Divide(1,2);

//---------------------------------------------------------------------------------- 
	c1->cd(1);

 	h1 = (TH2D*)f->Get("D121");
	h1->FitSlicesY();
	TH1D *h1_2 = (TH1D*)gDirectory->Get("D121_2"); 	
	h1_2->SetName("D121_2_1");
	TH1D *h1_0 = (TH1D*)gDirectory->Get("D121_0"); 	
	h1_0->SetName("D121_0_1");
	TH1D *h1_1 = (TH1D*)gDirectory->Get("D121_1"); 	
	h1_1->SetName("D121_1_1");
	TH1D *h1_chi2 = (TH1D*)gDirectory->Get("D121_chi2"); 	
	h1_chi2->SetName("D121_chi2_1");

 	h2 = (TH2D*)f2->Get("D121");
	h2->FitSlicesY();
	TH1D *h2_2 = (TH1D*)gDirectory->Get("D121_2"); 	
	h2_2->SetName("D121_2_2");
	TH1D *h2_0 = (TH1D*)gDirectory->Get("D121_0"); 	
	h2_0->SetName("D121_0_2");
	TH1D *h2_1 = (TH1D*)gDirectory->Get("D121_1"); 	
	h2_1->SetName("D121_1_2");
	TH1D *h2_chi2 = (TH1D*)gDirectory->Get("D121_chi2"); 	
	h2_chi2->SetName("D121_chi2_2");

  	setHist(h1_2, kCstCol,  8, 0.9);
  	setHist(h2_2, kBlack, 24, 0.9);	

 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");


  	leggL->Clear();
	leggLe = leggL->AddEntry(h1_2, Form("perfect al. (B_{s} #rightarrow #mu #mu)"), "p"); leggLe->SetTextColor(kCstCol);
    
	leggLe = leggL->AddEntry(h2_2, Form("short-term al. (B^{+} #rightarrow J/#Psi K^{+})"), "p"); leggLe->SetTextColor(kBlack);

 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

	c2->cd(1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");
 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

//----------------------------------------------------------------------------------

	c1->cd(2);
	c1_2->Divide(3,2);


	c1_2->cd(4);
	h1->ProjectionX("h121_1_x");
	h2->ProjectionX("h121_2_x");

	setHist(h121_1_x, kCstCol,  8, 0.4);
  	setHist(h121_2_x, kBlack, 24, 0.4);	

	h121_1_x->SetTitle("p_{T} distribution");
	h121_1_x->DrawCopy("E");  gPad->SetLogy(1); 
 	h121_2_x->DrawCopy("ESAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");

	c1_2->cd(1);

	setHist(h1_0, kCstCol,  8, 0.4);
  	setHist(h2_0, kBlack, 24, 0.4);	

	h1_0->GetYaxis()->SetRangeUser(0.1*h1_0->GetMinimum(1E-20), 10*h1_0->GetMaximum());
 	h1_0->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_0->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(2);

	setHist(h1_1, kCstCol,  8, 0.4);
  	setHist(h2_1, kBlack, 24, 0.4);	

	h1_1->GetYaxis()->SetRangeUser(0.1*h1_1->GetMinimum(1E-20), 10*h1_1->GetMaximum());
 	h1_1->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_1->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(3);

	setHist(h1_chi2, kCstCol,  8, 0.4);
  	setHist(h2_chi2, kBlack, 24, 0.4);	

	h1_chi2->GetYaxis()->SetRangeUser(0.1*h1_chi2->GetMinimum(1E-20), 10*h1_chi2->GetMaximum());
 	h1_chi2->DrawCopy("P"); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_chi2->DrawCopy("PSAME");

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(5);
  	setHist(h1_2, kCstCol,  8, 0.4);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->SetTitle(Form("%s vs. RMS",h1_2->GetTitle()));
	h1_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");

	c1_2->cd(6);
  	setHist(h2_2, kBlack, 24, 0.4);	
 	h2_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h2_2->SetTitle(Form("%s vs. RMS",h2_2->GetTitle()));
	h2_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");	

	// --------------- RMS ------------------------
 	p1 = (TH1D*)f->Get("p121");
  	p2 = (TH1D*)f2->Get("p121");
 	
 	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);
 	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);

  	setHist(hnew1, kCstCol,  20, 0.4);
  	setHist(hnew2, kBlack, 20, 0.4);

	for (int i = 0; i < p1->GetNbinsX(); i++ ) {
	  
	  hnew1->SetBinContent(i, p1->GetBinError(i));
  	  hnew2->SetBinContent(i, p2->GetBinError(i));
	}

	c1_2->cd(5);
	hnew1->DrawCopy("HISTSAME");
	c1_2->cd(6);
  	hnew2->DrawCopy("HISTSAME");

	delete hnew1;
	delete hnew2;
	
//----------------------------------------------------------------------------------
 
 	h1 = (TH2D*)f->Get("D131");
	h1->FitSlicesY(0, -1, -1, 0);
	TH1D *h1_2 = (TH1D*)gDirectory->Get("D131_2"); 	
	h1_2->SetName("D131_2_1");
	TH1D *h1_0 = (TH1D*)gDirectory->Get("D131_0"); 	
	h1_0->SetName("D131_0_1");
	TH1D *h1_chi2 = (TH1D*)gDirectory->Get("D131_chi2"); 	
	h1_chi2->SetName("D131_chi2_1");

 	h2 = (TH2D*)f2->Get("D131");
	h2->FitSlicesY();
	TH1D *h2_2 = (TH1D*)gDirectory->Get("D131_2"); 	
	h1_2->SetName("D131_2_2");
	TH1D *h2_0 = (TH1D*)gDirectory->Get("D131_0"); 	
	h2_0->SetName("D131_0_2");
	TH1D *h2_chi2 = (TH1D*)gDirectory->Get("D131_chi2"); 	
	h2_chi2->SetName("D131_chi2_2");

  	setHist(h1_2, kCstCol,  8, 0.9);
  	setHist(h2_2, kBlack, 24, 0.9);	

	c1->cd(3);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");

 	leggL->Draw();

	//	tl->DrawLatex(0.05, 0.83, "_{[%]}");
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");

	c2->cd(2);
 	h1_2->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_2->DrawCopy("SAME");

 	leggL->Draw();
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");
//----------------------------------------------------------------------------------

	c1->cd(4);
	c1_4->Divide(3,2);


	c1_4->cd(4);
	h1->ProjectionX("h131_1_x");
	h2->ProjectionX("h131_2_x");

	setHist(h131_1_x, kCstCol,  8, 0.4);
  	setHist(h131_2_x, kBlack, 24, 0.4);	

	h131_1_x->SetTitle("#eta distribution");
	h131_1_x->DrawCopy("E");  gPad->SetLogy(1); 
 	h131_2_x->DrawCopy("ESAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");

	c1_4->cd(1);

	setHist(h1_0, kCstCol,  8, 0.4);
  	setHist(h2_0, kBlack, 24, 0.4);	

	h1_0->GetYaxis()->SetRangeUser(0.1*h1_0->GetMinimum(1E-20), 10*h1_0->GetMaximum());
 	h1_0->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_0->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(2);

	setHist(h1_1, kCstCol,  8, 0.4);
  	setHist(h2_1, kBlack, 24, 0.4);	

	h1_1->GetYaxis()->SetRangeUser(0.1*h1_1->GetMinimum(1E-20), 10*h1_1->GetMaximum());
 	h1_1->DrawCopy(); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_1->DrawCopy("SAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(3);

	setHist(h1_chi2, kCstCol,  8, 0.4);
  	setHist(h2_chi2, kBlack, 24, 0.4);	

	h1_chi2->GetYaxis()->SetRangeUser(0.1*h1_chi2->GetMinimum(1E-20), 10*h1_chi2->GetMaximum());
 	h1_chi2->DrawCopy("P"); gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	h2_chi2->DrawCopy("PSAME");

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(5);
  	setHist(h1_2, kCstCol,  8, 0.4);
 	h1_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h1_2->SetTitle(Form("%s vs. RMS",h1_2->GetTitle()));
	h1_2->DrawCopy(); gPad->SetLogy(1); 

	tl->DrawLatex(0.85, 0.05, "_{#eta}");

	c1_4->cd(6);
  	setHist(h2_2, kBlack, 24, 0.4);	
 	h2_2->GetYaxis()->SetRangeUser(5E-3, 1E-1);
	h2_2->SetTitle(Form("%s vs. RMS",h2_2->GetTitle()));
	h2_2->DrawCopy(); gPad->SetLogy(1);

	tl->DrawLatex(0.85, 0.05, "_{#eta}"); 	

	// --------------- RMS ------------------------
 	p1 = (TH1D*)f->Get("p131");
  	p2 = (TH1D*)f2->Get("p131");
 	
 	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 2.5);
 	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 2.5);

  	setHist(hnew1, kCstCol,  20, 0.4);
  	setHist(hnew2, kBlack, 20, 0.4);

	for (int i = 0; i < p1->GetNbinsX(); i++ ) {
	  
	  hnew1->SetBinContent(i, p1->GetBinError(i));
  	  hnew2->SetBinContent(i, p2->GetBinError(i));
	}

	c1_4->cd(5);
	hnew1->DrawCopy("HISTSAME");
	c1_4->cd(6);
  	hnew2->DrawCopy("HISTSAME");

	delete hnew1;
	delete hnew2;

  	setHist(h1_2, kCstCol,  8, 0.9);
  	setHist(h2_2, kBlack, 24, 0.9);	

 	c1->SaveAs("pT_res2_fits.pdf");
 	c2->SaveAs("pT_res2.pdf");

	delete c2;

}

// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color = 1, Int_t symbol = 20, Double_t size = 1.6, Double_t width = 1.) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}
