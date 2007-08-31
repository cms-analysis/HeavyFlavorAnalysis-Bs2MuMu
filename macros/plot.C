void plot(const char *filName1, const char *filName2="", const char *opt1 = "", const char *opt2="")
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

	TFile *f = new TFile(Form("%s", filName1));
	if(opt1=="SAME") TFile *f2 = new TFile(Form("%s", filName2));

	TH1D *h1;
	TH1D *h2;

	TH1D *cnt1;
	TH1D *cnt2;

	TCanvas *c1 = new TCanvas("c1", "", 1000, 500);
	c1->Clear();
	if (0 != strcmp(opt2, "single")) {
	  c1->Divide(4,3);
	}

	TCanvas *c2 = new TCanvas("c2", "", 300, 300);
	c2->Clear();

	TCanvas *c3 = new TCanvas("c3", "", 300, 300);
	c3->Clear();

	TCanvas *c4a = new TCanvas("c4a", "", 600, 300);
	c4a->Clear();
	c4a->Divide(2,1);

	TCanvas *c4b = new TCanvas("c4b", "", 600, 300);
	c4b->Clear();
	c4b->Divide(2,1);

	char histname[20];
	sprintf(histname, "c1"); // before HLT cuts
	sprintf(histname, "c2"); // after HLT cuts
	sprintf(histname, "c0"); // before any cuts

//----------------------------------------------------------------------------------
	if(0 != strcmp(opt2, "single")) c1->cd(1);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "00"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "00"));
	//	scaleHisto(h1, h2, h1, h1, 1.02, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("m");
	
 	if (0 == strcmp(opt1, "SAME")) {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "pT(B_{s}) [GeV/c]");
	setTitles(h1, h2, "p_{T} [GeV]", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1);
	if (0 == strcmp(opt2, "single")) c1->SaveAs("pT_Bs.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(2);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "10"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "10"));
	//	scaleHisto(h1, h2, h1, h1, 1.02, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e");

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "pT^{#mu}_{min} [GeV/c]");
	setTitles(h1, h2, "p_{T}^{#mu} [GeV]", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	if(opt2=="single") c1->SaveAs("pT_mu.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(3);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "11"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "11"));
	//	scaleHisto(h1, h2, h1, h1, 1.02, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e");

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "#eta leptons");
	setTitles(h1, h2, "#eta_{#mu}", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	if(opt2=="single") c1->SaveAs("eta_mu.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(4);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "13"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "13"));
	//	scaleHisto(h1, h2, h1, h1, 1.2, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "TIP(#mu) [cm]");
	setTitles(h1, h2, "#TIP(#mu)", "events/bin", 0.7, 0.85, 0.48, 0.85, opt1, 0);
	if(opt2=="single") c1->SaveAs("TIP.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(5);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "25"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "25"));
	scaleHisto(h1, h2, h1 , h1, 1.2, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "Isolation");
	setTitles(h1, h2, "Isolation", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
        if(opt2=="single") c1->SaveAs("mass.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(6);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "26"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "26"));
	scaleHisto(h1, h2, h1, h1, 1.2, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "|V_{0} - V_{B_{S}}| [cm]");
	setTitles(h1, h2, "|V_{0} - V_{B_{S}}| [cm]", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	if(opt2=="single") c1->SaveAs("Vertex.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(7);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "27"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "27"));
	scaleHisto(h1, h2, h1, h1, 1.2, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "#chi^{2}");
	setTitles(h1, h2, "#chi^{2}", "events/bin", 0.7, 0.85, 0.48, 0.85, opt1, 0);
	if(opt2=="single") c1->SaveAs("Chi2.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(8);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "22"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "22"));
	scaleHisto(h1, h2, h1, h1, 1.2,0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "L_{xy}/#sigma_{xy}");
	setTitles(h1, h2, "L_{xy}/#sigma_{xy}", "events/bin", 0.7, 0.85, 0.48, 0.85, opt1, 0);
	if(opt2=="single") c1->SaveAs("Lxy_sigma.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(9);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "23"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "23"));
	scaleHisto(h1, h2, h1, h1, 1.2, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "#sigma_{xy} [cm]");
	setTitles(h1, h2, "#sigma_{xy} [cm]", "events/bin", 0.7, 0.85, 0.48, 0.85, opt1, 0);
	if(opt2=="single") c1->SaveAs("sigmaxy.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(10);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "21"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "21"));
	scaleHisto(h1, h2, h1, h1, 1.2 , 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e"); gPad->SetLogy(1);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "cos(#alpha)");
	setTitles(h1, h2, "cos(angle) (p,v)", "events/bin", 0.7, 0.85, 0.48, 0.85, opt1, 0);
	if(opt2=="single") c1->SaveAs("cosAngle.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(11);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "20"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "20"));
	//	scaleHisto(h1, h2, h1 , h1, 1.002, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e");gPad->SetLogy(0);

	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "#DeltaR(#mu#mu)");
	setTitles(h1, h2, "#DeltaR(#mu#mu)", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	if(opt2=="single") c1->SaveAs("Rmm.pdf");
//----------------------------------------------------------------------------------
	if(opt2!="single") c1->cd(12);
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "24"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "24"));
	//	scaleHisto(h1, h2, h1, h1, 1.002, 0);
	setHist(h1, kBlue, 20);
	h1->Draw("e");
	if (opt1=="SAME") {
	  setHist(h2, kBlack);
	  h2->Draw("SAMEHIST");
	}

	tl->DrawLatex(0.15, 0.93, "Isolation Veto");
	setTitles(h1, h2, "IsoVeto", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	if(opt2=="single") c1->SaveAs("IsoVeto.pdf");
//----------------------------------------------------------------------------------

	if(opt2!="single") { c1->SaveAs("summaryPlot.pdf");  }

	legg = new TLegend(0.5,0.7,0.8,0.89);
	legg->SetFillStyle(0); legg->SetBorderSize(0); legg->SetTextSize(0.05);  legg->SetFillColor(0);

	leggL = new TLegend(0.2,0.7,0.6,0.89);
	leggL->SetFillStyle(0); leggL->SetBorderSize(0); leggL->SetTextSize(0.05);  leggL->SetFillColor(0); 
//----------------------------------------------------------------------------------
	c2->cd();
	h1 = (TH1D*)f->Get(Form("%s%s", histname, "30"));
	h2 = (TH1D*)f2->Get(Form("%s%s", histname, "30"));
	//	scaleHisto(h1, h2, h1, h1, 1., 0);
	setHist(h1, kBlue);
	setHist(h2, kBlack, 20, 0.1); 
	h1->GetYaxis()->SetRangeUser(0., 1.2*h1->GetMaximum());
	h1->Draw("HIST");gPad->SetLogy(0);
	h2->Draw("ESAME");

	legg->Clear();
	legge = legg->AddEntry(h1, Form("perfect al."), "l"); legge->SetTextColor(kBlue);
	legge = legg->AddEntry(h2, Form("short-term al."), "p"); legge->SetTextColor(kBlack);
	legg->Draw();

	tl->DrawLatex(0.15, 0.93, "Mass B [GeV/c^{2}]");
	setTitles(h1, h2, "Mass m_{B} [GeV]", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	c2->SaveAs("mass.pdf");

//----------------------------------------------------------------------------------
	c3->cd();
	h1 = (TH1D*)f->Get("j101");
	h2 = (TH1D*)f2->Get("j101");
	//	scaleHisto(h1, h2, h1, h1, 0.9, 0);
	setHist(h1, kBlue);
	setHist(h2, kBlack, 20, 0.1);
	h1->GetYaxis()->SetRangeUser(0., 1.2*h1->GetMaximum());
	h1->Draw("HIST");gPad->SetLogy(0);
	h2->Draw("ESAME");


	leggLe = leggL->AddEntry(h1, Form("perfect al."), "l"); leggLe->SetTextColor(kBlue);
	leggLe = leggL->AddEntry(h2, Form("short-term al."), "p"); leggLe->SetTextColor(kBlack);
	leggL->Draw();

	tl->DrawLatex(0.15, 0.93, "Mass [GeV/c^{2}]");
	setTitles(h1, h2, "Mass m_{#mu+}m_{#mu-} [GeV]", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
	c3->SaveAs("mass_jpsi.pdf");



//----------------------------------------------------------------------------------
 	c4a->cd(1);
 	h1 = (TH1D*)f->Get("p120");
 	h2 = (TH1D*)f2->Get("p120");

 	cnt1 = (TH1D*)f->Get("p120C");
	cnt2 = (TH1D*)f2->Get("p120C");
 	
	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T})/bin p_{T}", 25, 0., 25);
	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T})/bin p_{T}", 25, 0., 25);

 	setHist(hnew1, kBlue,  20, 0.1);
 	setHist(hnew2, kBlack, 20, 0.1);

        for (int i = 0; i < h1->GetNbinsX(); i++ ) {

 	  hnew1->SetBinContent(i, h1->GetBinError(i)/(i+0.5));
	  if (cnt1->GetBinContent(i))  hnew1->SetBinError(i, 0.01/TMath::Sqrt(cnt1->GetBinContent(i)));

	  hnew2->SetBinContent(i, h2->GetBinError(i)/(i+0.5));
	  if (cnt2->GetBinContent(i))  hnew2->SetBinError(i, 0.01/TMath::Sqrt(cnt2->GetBinContent(i)));
        }

 	hnew1->GetYaxis()->SetRangeUser(5E-3, 1E-1);
 	hnew1->Draw("HIST");gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	hnew2->Draw("ESAME");

 	leggL->Clear();
 	leggLe = leggL->AddEntry(hnew1, Form("perfect al."), "l"); leggLe->SetTextColor(kBlue);
 	leggLe = leggL->AddEntry(hnew2, Form("short-term al."), "p"); leggLe->SetTextColor(kBlack);
 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
 	setTitles(h1, h2, "#frac{#sigma_{p_{T}}}{p_{T}} vs p_{T}", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);

 	TCanvas *c5 = new TCanvas("c5", "", 300, 300);
 	c5->Clear();
 	setHist(cnt1, kBlue,  20, 0.1);
 	setHist(cnt2, kBlack, 20, 0.1);
 	cnt1->GetYaxis()->SetRangeUser(0, 1.2*cnt1->GetMaximum());
 	cnt1->Draw("HIST");gPad->SetLogy(0); gPad->SetTopMargin(0.12);
 	cnt2->Draw("HISTSAME");

//----------------------------------------------------------------------------------
 	c4a->cd(2);
 	h1 = (TH1D*)f->Get("p121");
 	h2 = (TH1D*)f2->Get("p121");
 	
	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);
	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0., 25);

 	setHist(hnew1, kBlue,  20, 0.1);
 	setHist(hnew2, kBlack, 20, 0.1);

        for (int i = 0; i < h1->GetNbinsX(); i++ ) {

 	  hnew1->SetBinContent(i, h1->GetBinError(i));
	  if (cnt1->GetBinContent(i))  hnew1->SetBinError(i, 0.01/TMath::Sqrt(cnt1->GetBinContent(i)));

	  hnew2->SetBinContent(i, h2->GetBinError(i));
	  if (cnt2->GetBinContent(i))  hnew2->SetBinError(i, 0.01/TMath::Sqrt(cnt2->GetBinContent(i)));
        }

 	hnew1->GetYaxis()->SetRangeUser(5E-3, 1E-1);
 	hnew1->Draw("HIST");gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	hnew2->Draw("ESAME");

 	leggL->Clear();
 	leggLe = leggL->AddEntry(hnew1, Form("perfect al."), "l"); leggLe->SetTextColor(kBlue);
 	leggLe = leggL->AddEntry(hnew2, Form("short-term al."), "p"); leggLe->SetTextColor(kBlack);
 	leggL->Draw();

	tl->DrawLatex(0.85, 0.05, "_{p_{T}}");
	tl->DrawLatex(0.75, 0.95, "_{global muons}");
 	setTitles(h1, h2, "#frac{#sigma_{p_{T}}}{p_{T}} vs p_{T}", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
 	c4a->SaveAs("pT_res.pdf");

//----------------------------------------------------------------------------------
 	c4b->cd(1);
 	h1 = (TH1D*)f->Get("p130");
 	h2 = (TH1D*)f2->Get("p130");

 	cnt1 = (TH1D*)f->Get("p130C");
 	cnt2 = (TH1D*)f2->Get("p130C");
 
	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T})/10 GeV", 25, 0., 2.5);
	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T})/10 GeV", 25, 0., 2.5);

 	setHist(hnew1, kBlue,  20, 0.1);
 	setHist(hnew2, kBlack, 20, 0.1);

        for (int i = 0; i < h1->GetNbinsX(); i++ ) {

 	  hnew1->SetBinContent(i, 100*h1->GetBinError(i)/10.);
	  if (cnt1->GetBinContent(i))  hnew1->SetBinError(i, 1/TMath::Sqrt(cnt1->GetBinContent(i)));

 	  hnew2->SetBinContent(i, 100*h2->GetBinError(i)/10.);
	  if (cnt2->GetBinContent(i))  hnew2->SetBinError(i, 1/TMath::Sqrt(cnt2->GetBinContent(i)));
        }

 	hnew1->GetYaxis()->SetRangeUser(1e-1, 10);
 	hnew1->Draw("HIST");gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	hnew2->Draw("ESAME");

 	leggL->Clear();
 	leggLe = leggL->AddEntry(hnew1, Form("perfect al."), "l"); leggLe->SetTextColor(kBlue);
 	leggLe = leggL->AddEntry(hnew2, Form("short-term al."), "p"); leggLe->SetTextColor(kBlack);
 	leggL->Draw();

	tl->DrawLatex(0.05, 0.83, "_{[%]}");
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
 	setTitles(h1, h2, "#frac{#sigma_{p_{T}}}{p_{T}} vs #eta", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
 
 	TCanvas *c51 = new TCanvas("c51", "", 300, 300);
 	c51->Clear();
 	setHist(cnt1, kBlue,  20, 0.1);
 	setHist(cnt2, kBlack, 20, 0.1);
 	cnt1->GetYaxis()->SetRangeUser(0, 1.2*cnt1->GetMaximum());
 	cnt1->Draw("HIST");gPad->SetLogy(0); gPad->SetTopMargin(0.12);
 	cnt2->Draw("HISTSAME");

//----------------------------------------------------------------------------------
 	c4b->cd(2);
 	h1 = (TH1D*)f->Get("p131");
 	h2 = (TH1D*)f2->Get("p131");

	TH1D *hnew1 = new TH1D("pTres1","#sigma(#delta p_{T}/p_{T})", 25, 0.,2.5);
	TH1D *hnew2 = new TH1D("pTres2","#sigma(#delta p_{T}/p_{T})", 25, 0.,2.5);

 	setHist(hnew1, kBlue,  20, 0.1);
 	setHist(hnew2, kBlack, 20, 0.1);

        for (int i = 0; i < h1->GetNbinsX(); i++ ) {

 	  hnew1->SetBinContent(i, 100*h1->GetBinError(i));
	  if (cnt1->GetBinContent(i))  hnew1->SetBinError(i, 1/TMath::Sqrt(cnt1->GetBinContent(i)));

 	  hnew2->SetBinContent(i, 100*h2->GetBinError(i));
	  if (cnt2->GetBinContent(i))  hnew2->SetBinError(i, 1/TMath::Sqrt(cnt2->GetBinContent(i)));
        }

 	hnew1->GetYaxis()->SetRangeUser(1e-1, 10);
 	hnew1->Draw("HIST");gPad->SetLogy(1); gPad->SetTopMargin(0.12);
 	hnew2->Draw("ESAME");

 	leggL->Clear();
 	leggLe = leggL->AddEntry(hnew1, Form("perfect al."), "l"); leggLe->SetTextColor(kBlue);
 	leggLe = leggL->AddEntry(hnew2, Form("short-term al."), "p"); leggLe->SetTextColor(kBlack);
 	leggL->Draw();

	tl->DrawLatex(0.05, 0.83, "_{[%]}");
	tl->DrawLatex(0.85, 0.05, "_{#eta}");
	tl->DrawLatex(0.60, 0.95, "_{p_{T}^{glb. #mu} = 8 - 12 GeV}");
 	setTitles(h1, h2, "#frac{#sigma_{p_{T}}}{p_{T}} vs #eta", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
 	c4b->SaveAs("eta_res.pdf");


// 	TCanvas *c6 = new TCanvas("c6", "", 600, 600);
// 	c6->Clear();
// 	c6->Divide(2,2);

//  //----------------------------------------------------------------------------------
//  	for (int pad = 0; pad < 4; pad++ ) {	

//  	  c6->cd(pad+1);
//  	  h1 = (TH1D*)f->Get(Form("j15%i", pad));
//  	  h2 = (TH1D*)f2->Get(Form("j15%i", pad));

//  	  //scaleHisto(h1, h2, h1, h1, 10.,0);
//  	  setHist(h1, kBlue);
//  	  setHist(h2, kBlack, 20, 0.1);
//  	  h1->GetYaxis()->SetRangeUser(0., 1.2*h1->GetMaximum());
//  	  h1->Draw("HIST");gPad->SetLogy(0); gPad->SetTopMargin(0.12);
//  	  h2->Draw("ESAME");

//  	  legg->Draw();

//  	  if ( pad==0 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(#mu #mu)");
//  	  if ( pad==1 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(J/#Psi K)");
//  	  if ( pad==2 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(J/#Psi #mu_{1})");
//  	  if ( pad==3 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(J/#Psi #mu_{2})");

//  	  setTitles(h1, h2, "#DeltaR", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
//  	}

//  	c6->SaveAs("deltaR.pdf");
//  	gPad->SetLogy(0);

// 	TCanvas *c7 = new TCanvas("c7", "", 600, 600);
// 	c7->Clear();
// 	c7->Divide(2,2);

//  //----------------------------------------------------------------------------------
//  	for (int pad = 0; pad < 4; pad++ ) {	

//  	  c7->cd(pad+1);
//  	  h1 = (TH1D*)f->Get(Form("g15%i", pad));
//  	  h2 = (TH1D*)f2->Get(Form("g15%i", pad));

//  	  //scaleHisto(h1, h2, h1, h1, 10.,0);
//  	  setHist(h1, kBlue);
//  	  setHist(h2, kBlack, 20, 0.1);
//  	  h1->GetYaxis()->SetRangeUser(0., 1.2*h1->GetMaximum());
//  	  h1->Draw("HIST");gPad->SetLogy(0); gPad->SetTopMargin(0.12);
//  	  h2->Draw("ESAME");

//  	  legg->Draw();

//  	  if ( pad==0 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(#mu #mu)");
//  	  if ( pad==1 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(J/#Psi K)");
//  	  if ( pad==2 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(J/#Psi #mu_{1})");
//  	  if ( pad==3 ) tl->DrawLatex(0.15, 0.93, "#DeltaR(J/#Psi #mu_{2})");

//  	  setTitles(h1, h2, "#DeltaR", "events/bin", 0.7, 0.85, 0.7, 0.68, opt1, 0);
//  	}

//  	c7->SaveAs("deltaR_gen.pdf");
//  	gPad->SetLogy(0);


}

//----------------------------------------------------------------------------------
void scaleHisto(TH1 *ha, TH1 *hb, TH1 *hc, TH1 *hd, double scale, const int &neg) {

  double extra(0.);

  double hmax[4] = {0,0,0,0};                         double max(0);
  double hmin[4] = {1000000,1000000,1000000,1000000}; double min(1000000);

  if (ha) hmax[0] = findMax(ha);
  if (hb) hmax[1] = findMax(hb);
  if (hc) hmax[2] = findMax(hc);
  if (hd) hmax[3] = findMax(hd);

  if (ha) hmin[0] = findMin(ha);
  if (hb) hmin[1] = findMin(hb);
  if (hc) hmin[2] = findMin(hc);
  if (hd) hmin[3] = findMin(hd);

  for ( int i = 0; i < 4; i++ ) {

    if ( hmax[i] > max ) { max = hmax[i]; }
    if ( hmin[i] < min ) { min = hmin[i]; }

  }

  extra = scale*(max - min);

  if( neg ) {

    ha->GetYaxis()->SetRangeUser(min-extra, max+extra);
    hb->GetYaxis()->SetRangeUser(min-extra, max+extra);
    hc->GetYaxis()->SetRangeUser(min-extra, max+extra);
    hd->GetYaxis()->SetRangeUser(min-extra, max+extra);
  }
  else {

    ha->GetYaxis()->SetRangeUser(0.1, max+extra);
    hb->GetYaxis()->SetRangeUser(0.1, max+extra);
    hc->GetYaxis()->SetRangeUser(0.1, max+extra);
    hd->GetYaxis()->SetRangeUser(0.1, max+extra);
  }
}

// ----------------------------------------------------------------------
double findMax(TH1 *h) {

    double max = -2000.;
    double tmp = 0.;

    for (int i = 0; i < h->GetNbinsX(); ++i) {

      tmp = h->GetBinContent(i);
      if ( tmp > max && tmp != 0 ) max = tmp;

    }

    return max;
}
//-----------------------------------------------------------------------------------------
double findMin(TH1 *h) {

    double min = 2000.;
    double tmp = 0.;

    for (int i = 0; i < h->GetNbinsX(); ++i) {

      tmp = h->GetBinContent(i);
      if ( tmp < min  && tmp != 0 ) min = tmp;

    }

    return min;
}

//-----------------------------------------------------------------------------------------
void setTitles(TH1 *ha, TH1 *hb, const char *sx, const char *sy, float x1, float y1, float x2, float y2,
                                 const char *DrawOpt1, const int *stat = 1) {

  if (ha == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    ha->SetXTitle(sx);                  ha->SetYTitle(sy);
    ha->SetTitleOffset(0.05, "x");      ha->SetTitleOffset(0.9, "y");
    ha->SetTitleSize(1.2, "x");         ha->SetTitleSize(1.2, "y");
    ha->SetLabelSize(0.06, "x");        ha->SetLabelSize(0.06, "y");
    ha->SetLabelFont(132, "x");         ha->SetLabelFont(132, "y");
    ha->GetXaxis()->SetTitleFont(132);  ha->GetYaxis()->SetTitleFont(132);
    ha->SetNdivisions(508, "X");        ha->GetYaxis()->CenterTitle();
    ha->SetMarkerColor(1);              ha->SetMarkerColor(4);
    ha->SetTitle("");
    ha->SetMarkerStyle(1);
    //ha->SetMarkerSize(1);

    if(DrawOpt1 == "SAME") {
      hb->Scale(ha->GetSumOfWeights()/hb->GetSumOfWeights());
      hb->SetTitle("");
    }

    if ( stat ) {
      TLatex *ts = new TLatex;
      ts->SetNDC(kTRUE);
      float fontsize = 0.05;
      ts->SetTextSize(fontsize);
      
      ts->SetTextColor(kBlue);
      ts->DrawLatex(x1, y1, Form("Entries %4.0f", ha->GetEntries()));
      y1 -= fontsize*1.1;
      ts->DrawLatex(x1, y1, Form("Mean %4.2f", ha->GetMean()));
      y1 -= fontsize*1.1;
      ts->DrawLatex(x1, y1, Form("RMS %4.2f", ha->GetRMS()));

      if (DrawOpt1 == "SAME") {
	ts->SetTextColor(kBlack);
	ts->DrawLatex(x2, y2, Form("Entries %4.0f", hb->GetEntries()));
	y2 -= fontsize*1.1;
	ts->DrawLatex(x2, y2, Form("Mean %4.2f", hb->GetMean()));
	y2 -= fontsize*1.1;
	ts->DrawLatex(x2, y2, Form("RMS %4.2f", hb->GetRMS()));
      }
    }
  }
}

// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color = 1, Int_t symbol = 20, Double_t size = 1.6, Double_t width = 1.) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(0); h->SetFillColor(color);
}



void massPlot() {
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  
  gStyle->SetStatFont(132);
  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132, "X");
  gStyle->SetLabelFont(132, "Y");
  gStyle->SetTitleFont(132);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  c028->SetTitleOffset(1.02, "x");      
  c028->SetTitleOffset(1.20, "y");
  c028->SetAxisRange(0, 20, "x");
  c028->SetXTitle("m_{#mu^{+}, #mu^{-}} [GeV/c^{2}]");
  c028->Draw();

  c0->SaveAs("mass.pdf");
}
