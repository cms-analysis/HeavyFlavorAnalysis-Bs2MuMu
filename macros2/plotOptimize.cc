#include "plotOptimize.hh"

#include "../macros/AnalysisDistribution.hh"
#include "../../../AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#include "../interface/HFMasses.hh"
#include "../macros/bayesianlimit.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "THStack.h"

#include <iomanip>
#include <string>
#include <list>

using namespace std; 
using std::string; 

ClassImp(plotOptimize)


// ----------------------------------------------------------------------
void plotOptimize::optimizeULs(int nruns, int seed) {
  int version(-1); 
  fDoPrint = false; 
  ofstream OUT(Form("optimizeUL-%d.txt", seed)); 

  int NCUTS(12);
  //                  0       1      2     3        4       5        6          7      8           9      10      11
  //string cuts[] = {"m2pt", "m1pt","pt", "alpha", "chi2", "fls3d", "docatrk", "iso", "closetrk", "lip", "lips", "mwindow"}; 
  double loCuts[] = {4.0,    4.0,   5.0,  0.01,    1.2,      5,       0.0,      0.70,  1,         0.0,   1.,     0.020 };
  double hiCuts[] = {6.0,    7.0,   11.,  0.10,    2.5,     25,       0.1,      0.95,  4,         1.0,   3.5,    0.100};

  // to add: 
  // LIP(S)2

  if (seed > 0) {
    cout << "Setting random number seed " << seed << endl;
    OUT  << "Setting random number seed " << seed << endl;
    gRandom->SetSeed(seed);
  }

  
  string cutline; 
  double cut; 
  for (int j = 0; j < nruns; ++j) {
    ++version; 
    for (int i = 0; i < NCUTS; ++i) {
      cut = gRandom->Rndm()*(hiCuts[i]-loCuts[i]) + loCuts[i];
      for (int ic = 0; ic < 2; ++ic) {
	if (0 == i) fCuts[ic]->m2pt = cut; 
	if (1 == i) fCuts[ic]->m1pt = cut; 
	if (2 == i) fCuts[ic]->pt = cut; 
	if (3 == i) fCuts[ic]->alpha = cut; 
	if (4 == i) fCuts[ic]->chi2dof = cut; 
	if (5 == i) fCuts[ic]->fls3d = cut; 
	if (6 == i) fCuts[ic]->docatrk = cut; 
	if (7 == i) fCuts[ic]->iso = cut;
	if (8 == i) fCuts[ic]->closetrk = static_cast<int>(cut);
	if (9 == i) fCuts[ic]->pvlip = cut;
	if (10== i) fCuts[ic]->pvlips = cut;
	if (11== i) fCuts[ic]->mBsLo = 5.370 - cut; 
	if (11== i) fCuts[ic]->mBsHi = 5.370 + cut; 
	if (11== i) fCuts[ic]->mBsLo = 5.30;
	if (11== i) fCuts[ic]->mBsHi = 5.45;
      }
    }

	

    if (gRandom->Rndm() > 0.5) { 
      fDoApplyCowboyVetoAlsoInSignal = true; 
    } else {
      fDoApplyCowboyVetoAlsoInSignal = false; 
    }

    printCuts(OUT); 

    cout << "--> loopTree: signal MC" << endl;
    loopTree(0);  // signal eff
    c0->Modified(); c0->Update();
    loopTree(1);  // Bd2MuMu eff
    c0->Modified(); c0->Update();
    cout << "--> loopTree: signal data" << endl;
    loopTree(5);  // data signal
    c0->Modified(); c0->Update();
    
    
    // -- simple blimits in two channels
    int nobs(-1); 
    for (int ichan = 0; ichan < 2; ++ichan) {
	fhMassWithAllCuts[ichan]->Draw();
	double scale   = fLumi["SgData"]/39.4;
	double nbs     = 2.0e9*(1.0-0.12)*scale;
	double siglo   = fNumbersBs[ichan]->mBsLo; 
	double sighi   = fNumbersBs[ichan]->mBsHi; 
	double mFactor = (sighi-siglo)/(fBgHi-fBgLo-(5.45-5.20)); 
	double lscale  = fLumi["SgData"]/fLumi["SgMc"];
    
	double signal  = fNumbersBs[ichan]->anaWmcYield*lscale;
	fBgHistExp     = fNumbersBs[ichan]->bgObs*mFactor;
	
	fBgExp = fBgHistExp;
	fBgExpE = 0.2*fBgExp;
	nobs = static_cast<int>(fBgExp + 0.5 + signal); 
	double nulbayes  = blimit(0.95, nobs, 1.0, 0.2, fBgExp, fBgExpE, 1);
	double effTot = fNumbersBs[ichan]->effTot*fNumbersBs[0]->pss;
	double ulbayes  = nulbayes/(effTot*nbs);
	cout << "==> effTot:     " << effTot << " chan=" << ichan << " version=" << version << endl;
	cout << "==> nObs:       " << fNumbersBs[ichan]->bgObs << " chan=" << ichan << " version=" << version  << endl;
	cout << "==> nExp:       " << fBgHistExp << " chan=" << ichan << " version=" << version << endl;
	cout << "==> ul(Bayes):  " << ulbayes << " chan=" << ichan << " version=" << version << endl;
	cout << "==> Signal:     " << signal << " chan=" << ichan << " version=" << version << endl;
	cout << "==> SSB:        " << signal/TMath::Sqrt(signal+fBgHistExp) << " chan=" << ichan << " version=" << version << endl;
	
	OUT << "==> effTot = " << effTot << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> nObs   = " << fNumbersBs[ichan]->bgObs << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> nExp   = " << fBgHistExp << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> UL     = " << ulbayes << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> Signal = " << signal << " chan=" << ichan << " version=" << version << endl;
	OUT << "==> SSB    = " << signal/TMath::Sqrt(signal+fBgHistExp) << " chan=" << ichan << " version=" << version << endl;
      }
  }

}


struct bla{
  double ul, ssb, ssb1, ssb2;
  double nobs, nexp; 
  double sig, eff; 
  double mlo, mhi; 
  double m1pt, m2pt, pt; 
  double chi2dof, iso, alpha, fls3d, docatrk; 
  double lip, lips; 
  int closetrk, cowboyVeto; 
};


// ----------------------------------------------------------------------
void plotOptimize::bestUL(const char *fname, int mode) {
  // mode = 0 UL
  //        1 SM S/sqrt(S+B) 
  //        2 SM sig
  TFile *f = TFile::Open(fname); 
  
  TTree *t = (TTree*)f->Get("t");

  int nsettings(20); 

  int chan, file, run; 
  float mlo, mhi, pt, m1pt, m2pt, iso, chi2dof, alpha, fls3d, docatrk; 
  float ul, ssb, ssb1, ssb2, nobs, nexp, sig, eff; 
  int closetrk, cowboyVeto; 
  float lip, lips; 

  t->SetBranchAddress("chan", &chan);
  t->SetBranchAddress("file", &file);
  t->SetBranchAddress("run", &run);
  t->SetBranchAddress("ul", &ul);
  t->SetBranchAddress("ssb", &ssb);
  t->SetBranchAddress("ssb1", &ssb1);
  t->SetBranchAddress("ssb2", &ssb2);

  t->SetBranchAddress("nobs", &nobs);
  t->SetBranchAddress("nexp", &nexp);
  t->SetBranchAddress("sig", &sig);
  t->SetBranchAddress("eff", &eff);

  t->SetBranchAddress("mlo", &mlo);
  t->SetBranchAddress("mhi", &mhi);
  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("m1pt", &m1pt);
  t->SetBranchAddress("m2pt", &m2pt);
  t->SetBranchAddress("iso", &iso);
  t->SetBranchAddress("chi2dof", &chi2dof);
  t->SetBranchAddress("alpha", &alpha);
  t->SetBranchAddress("fls3d", &fls3d);
  t->SetBranchAddress("docatrk", &docatrk);
  t->SetBranchAddress("closetrk", &closetrk);
  t->SetBranchAddress("pvlip", &lip);
  t->SetBranchAddress("pvlips", &lips);
  t->SetBranchAddress("cowboyVeto", &cowboyVeto);

  // -- initialize to some minimum values beyond which neither ul nor ssb are interesting
  bla ini = {2.5e-8, 1.4, 1.2, 1.2, 
	     0., 0., 
	     4.1, 0.,
	     0., 0., 0., 
	     0., 0., 0., 0., 0.,
	     0., 0., 0, 0
  };

  list<bla> bestList0(1, ini) ;

  ini.ul = 3.0e-8; 
  ini.ssb = 0.8;
  ini.sig = 2.0;
  list<bla> bestList1(1, ini) ;

  int nb(0); 
  int nentries = Int_t(t->GetEntries());
  cout << "Searching for best fom in " << nentries << " cut settings" << endl;
  for (int jentry = 0; jentry < nentries; jentry++) {
    nb = t->GetEntry(jentry);

    if (jentry %20000 == 0) cout << "setting " << jentry 
				<< " bestList0.size() = " << bestList0.size()
				<< " bestList1.size() = " << bestList1.size()
				<< endl;
    
    ini.ul = ul; 
    ini.ssb  = ssb; 
    ini.ssb1 = ssb1; 
    ini.ssb2 = ssb2; 
    ini.nobs = nobs; 
    ini.nexp = nexp;
    ini.sig = sig; 
    ini.eff = eff; 
    ini.mlo = mlo; 
    ini.mhi = mhi; 
    ini.pt = pt; 
    ini.m1pt = m1pt; 
    ini.m2pt = m2pt; 
    ini.iso = iso; 
    ini.chi2dof = chi2dof; 
    ini.alpha = alpha; 
    ini.fls3d = fls3d; 
    ini.docatrk = docatrk; 
    ini.closetrk = closetrk; 
    ini.lip = lip; 
    ini.lips = lips; 
    ini.cowboyVeto = cowboyVeto; 

    if (chan == 0) {
      for (list<bla>::iterator i = bestList0.begin(); i != bestList0.end(); ++i) {
	if (0 == mode) {
	  if (ini.ul < i->ul) { 
	    bestList0.insert(i, ini); 
	    break;
	  }
	}
	if (1 == mode) {
	  if (ini.ssb > i->ssb) { 
	    bestList0.insert(i, ini); 
	    break;
	  }
	}
	
	if (2 == mode) {
	  if (ini.sig > i->sig) { 
	    bestList0.insert(i, ini); 
	    break;
	  }
	}
	
      }
    } else {
      for (list<bla>::iterator i = bestList1.begin(); i != bestList1.end(); ++i) {
	if (0 == mode) {
	  if (ini.ul < i->ul) { 
	    bestList1.insert(i, ini); 
	    break;
	  }
	}
	if (1 == mode) {
	  if (ini.ssb > i->ssb) { 
	    bestList1.insert(i, ini); 
	    break;
	  }
	}
	
	if (2 == mode) {
	  if (ini.sig > i->sig) { 
	    bestList1.insert(i, ini); 
	    break;
	  }
	}
	
      }


    }

  }

  int icnt(0); 
  cout << "ssb:  S/sqrt(S+B) for SM BF" << endl;
  cout << "nobs: number of observed events in histogram" << endl;
  cout << "S/B: number of sg/bg events expected in Bs signal window" << endl;
  cout << "e:    efficiency" << endl;
  cout << "channel 0, mode " << mode << endl;
  for (list<bla>::iterator i = bestList0.begin(); i != bestList0.end(); ++i) {
    ++icnt;
    if (icnt > nsettings) break;
    ini = *i;
    cout << Form("ul=%2.1e ssb=%3.2f S/B=%3.2f/%3.2f e=%5.4f nobs=%2.0f ", ini.ul, ini.ssb,
		 ini.sig, ini.nexp, ini.eff, ini.nobs)
	 << Form("%4.3f<m<%4.3f pT=%3.1f ",  ini.mlo,  ini.mhi,  ini.pt)
	 << Form("pt1=%3.1f pt2=%3.1f I=%3.2f c2=%3.2f a=%4.3f f=%3.1f d=%4.3f ", 
		 ini.m1pt, ini.m2pt, ini.iso, ini.chi2dof, ini.alpha, ini.fls3d, ini.docatrk)
	 << Form("l=%4.3f/%3.2f n=%d v=%d", ini.lip, ini.lips, ini.closetrk, ini.cowboyVeto)
	 << endl;
  }

  cout << "ssb:  S/sqrt(S+B) for SM BF" << endl;
  cout << "nobs: number of observed events in histogram" << endl;
  cout << "S/B: number of sg/bg events expected in Bs signal window" << endl;
  cout << "e:    efficiency" << endl;
  cout << "channel 1, mode " << mode << endl;
  icnt = 0; 
  for (list<bla>::iterator i = bestList1.begin(); i != bestList1.end(); ++i) {
    ++icnt;
    if (icnt > nsettings) break;
    ini = *i;
    cout << Form("ul=%2.1e ssb=%3.2f S/B=%3.2f/%3.2f e=%5.4f nobs=%2.0f ", ini.ul, ini.ssb,
		 ini.sig, ini.nexp, ini.eff, ini.nobs)
	 << Form("%4.3f<m<%4.3f pT=%3.1f ",  ini.mlo,  ini.mhi,  ini.pt)
	 << Form("pt1=%3.1f pt2=%3.1f I=%3.2f c2=%3.2f a=%4.3f f=%3.1f d=%4.3f ", 
		 ini.m1pt, ini.m2pt, ini.iso, ini.chi2dof, ini.alpha, ini.fls3d, ini.docatrk)
	 << Form("l=%4.3f/%3.2f n=%d v=%d", ini.lip, ini.lips, ini.closetrk, ini.cowboyVeto)
	 << endl;
  }

}



// ----------------------------------------------------------------------
void plotOptimize::readOptimize(int nfiles) {
  TFile f("optimize.root", "RECREATE"); 

  TTree *t = new TTree("t","t");
  

  t->Branch("chan", &_chan ,"chan/I");
  t->Branch("file", &_file ,"file/I");
  t->Branch("run", &_run ,"run/I");
  t->Branch("ul", &_ul, "ul/F");
  t->Branch("nobs", &_nobs, "nobs/F");
  t->Branch("nexp", &_nexp, "nexp/F");
  t->Branch("sig", &_sig, "sig/F");
  t->Branch("ssb", &_ssb, "ssb/F");
  t->Branch("ssb1", &_ssb1, "ssb1/F");
  t->Branch("ssb2", &_ssb2, "ssb2/F");
  t->Branch("eff", &_eff, "eff/F");
  t->Branch("mlo", &_mlo, "mlo/F");
  t->Branch("mhi", &_mhi, "mhi/F");
  t->Branch("pt", &_pt, "pt/F");
  t->Branch("m1pt", &_m1pt, "m1pt/F");
  t->Branch("m2pt", &_m2pt, "m1pt/F");
  t->Branch("iso", &_iso, "iso/F");
  t->Branch("chi2dof", &_chi2dof, "chi2dof/F");
  t->Branch("alpha", &_alpha, "alpha/F");
  t->Branch("fls3d", &_fls3d, "fls3d/F");
  t->Branch("docatrk", &_docatrk, "docatrk/F");
  t->Branch("closetrk", &_closetrk, "closetrk/I");
  t->Branch("pvlip", &_pvlip, "pvlip/F");
  t->Branch("pvlips", &_pvlips, "pvlips/F");
  t->Branch("cowboyVeto", &_cowboyVeto, "cowboyveto/I");
  
  for (int i = 0; i <= nfiles; ++i) {
    _file = i; readFile(Form("optjobs/optimizeUL-%d.txt", i-1), t);
  }

  t->Write(); 
  f.Close();

}


// ----------------------------------------------------------------------
void plotOptimize::readFile(const char *fname, TTree *t) {

  char  buffer[200];
  ifstream is(fname);
  if (!is.is_open()) {
    cout << "skipping file " << fname << endl;
    return;
  }
  string line;
  double nu(0.); 

  cout << "reading file  " << fname << endl;
  while (is.getline(buffer, 200, '\n')) {
    line = buffer; 
    //    cout << line << endl;
    if (string::npos != line.find("mBsLo")) {
      sscanf(buffer, "mBsLo %f", &_mlo); 
      _ul  = +1.; // reset 
    }
    if (string::npos != line.find("mBsHi")) sscanf(buffer, "mBsHi %f", &_mhi); 
    if (string::npos != line.find("pt")) sscanf(buffer, "pt %f", &_pt); 
    if (string::npos != line.find("m1pt")) sscanf(buffer, "m1pt %f", &_m1pt); 
    if (string::npos != line.find("m2pt")) sscanf(buffer, "m2pt %f", &_m2pt); 
    if (string::npos != line.find("iso")) sscanf(buffer, "iso %f", &_iso); 
    if (string::npos != line.find("chi2dof")) sscanf(buffer, "chi2dof %f", &_chi2dof); 
    if (string::npos != line.find("alpha")) sscanf(buffer, "alpha %f", &_alpha); 
    if (string::npos != line.find("fls3d")) sscanf(buffer, "fls3d %f", &_fls3d); 
    if (string::npos != line.find("docatrk")) sscanf(buffer, "docatrk %f", &_docatrk); 
    if (string::npos != line.find("closetrk")) sscanf(buffer, "closetrk %i", &_closetrk); 
    if (string::npos != line.find("pvlip")) sscanf(buffer, "pvlip %f", &_pvlip); 
    if (string::npos != line.find("pvlips")) sscanf(buffer, "pvlips %f", &_pvlips); 
    if (string::npos != line.find("fDoApplyCowboyVetoAlsoInSignal")) sscanf(buffer, "fDoApplyCowboyVetoAlsoInSignal %i", &_cowboyVeto); 

    if (string::npos != line.find("chan=0")) {
      sscanf(buffer, "==> effTot = %f chan=0 version=%d", &_eff, &_run); 
      sscanf(buffer, "==> nObs   = %f chan=0 version=%d", &_nobs, &_run); 
      sscanf(buffer, "==> nExp   = %f chan=0 version=%d", &_nexp, &_run); 
      sscanf(buffer, "==> UL     = %f chan=0 version=%d", &_ul, &_run); 
      sscanf(buffer, "==> Signal = %f chan=0 version=%d", &_sig, &_run); 
      sscanf(buffer, "==> SSB    = %f chan=0 version=%d", &_ssb, &_run); 
      if (_ssb > 0) {
	_chan = 0; 
	nu = 1.0; 
	_ssb0 = nu*_sig/TMath::Sqrt(nu*_sig + _nexp);
	nu = 1.0/0.38; 
	_ssb1 = nu*_sig/TMath::Sqrt(nu*_sig + _nexp);
	nu = 2.0/0.38; 
	_ssb2 = nu*_sig/TMath::Sqrt(nu*_sig + _nexp);
	
// 	cout << "eff: " << _eff << " ul = " << _ul << " version = " << _run 
// 	     << " " << _ssb << " " << _ssb0 << " " << _ssb1 << " " << _ssb2 << endl;
	t->Fill();
	_ssb = -1;
      }
    }

    if (string::npos != line.find("chan=1")) {
      sscanf(buffer, "==> effTot = %f chan=1 version=%d", &_eff, &_run); 
      sscanf(buffer, "==> nObs   = %f chan=1 version=%d", &_nobs, &_run); 
      sscanf(buffer, "==> nExp   = %f chan=1 version=%d", &_nexp, &_run); 
      sscanf(buffer, "==> UL     = %f chan=1 version=%d", &_ul, &_run); 
      sscanf(buffer, "==> Signal = %f chan=1 version=%d", &_sig, &_run); 
      sscanf(buffer, "==> SSB    = %f chan=1 version=%d", &_ssb, &_run); 
      if (_ssb > 0) {
	_chan = 1; 
	nu = 1.0; 
	_ssb0 = nu*_sig/TMath::Sqrt(nu*_sig + _nexp);
	nu = 1.0/0.38; 
	_ssb1 = nu*_sig/TMath::Sqrt(nu*_sig + _nexp);
	nu = 2.0/0.38; 
	_ssb2 = nu*_sig/TMath::Sqrt(nu*_sig + _nexp);

// 	cout << "eff: " << _eff << " ul = " << _ul << " version = " << _run
// 	     << " " << _ssb << " " << _ssb0 << " " << _ssb1 << " " << _ssb2 << endl;
	t->Fill();
	_ssb = -1;
      }
    }
  }
  is.close();

}







// ----------------------------------------------------------------------
plotOptimize::plotOptimize(const char *files, const char *cuts, const char *dir, int mode) : plotClass(files, cuts, dir, mode) { 

  fDoPrint = true; 

  fNumbersFileName = fDirectory + "/anaBmm.plotOptimize." + fSuffix + ".tex";
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  fTEX.open(fNumbersFileName.c_str(), ios::app);

}

// ----------------------------------------------------------------------
plotOptimize::~plotOptimize() {
  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotOptimize::makeAll(int nfiles, int mode) {

  readOptimize(nfiles);
  bestUL("optimize.root", mode);

}


