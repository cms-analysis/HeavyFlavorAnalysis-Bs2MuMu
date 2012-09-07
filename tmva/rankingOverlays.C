#include "TH1.h"

// ----------------------------------------------------------------------
void rankingOverlays(string fname, int n = 2) {

  TH1D* vHist[15];
  vector<string> vNames;
  vNames.push_back("m1pt"); 
  vNames.push_back("m2pt"); 
  vNames.push_back("m1eta");
  vNames.push_back("m2eta");
  vNames.push_back("pt");   
  vNames.push_back("eta");  
  vNames.push_back("fls3d");
  vNames.push_back("alpha");
  vNames.push_back("maxdoca");
  vNames.push_back("pvip"); 
  vNames.push_back("pvips");
  vNames.push_back("iso");
  vNames.push_back("closetrk");
  vNames.push_back("docatrk");
  vNames.push_back("chi2dof");

  string hsName;
  TH1D* h;
  for (int i = 0; i < vNames.size(); ++i) {
    h = new TH1D(vNames[i].c_str(), vNames[i].c_str(), 100, 0., 0.5); 
    vHist[i] = h;
  }

  int seed(100), offset(0);
  string slabel;
  double w8; 
  for (int i = 0; i < n; ++i) {
    offset = seed + i; 
    h = getRanking(fname, "BDT", offset); 
    for (int ibin = 1; ibin < h->GetNbinsX(); ++ibin) {
      slabel = h->GetXaxis()->GetBinLabel(ibin);  
      if (!strcmp("", slabel.c_str())) break;
      w8 = h->GetBinContent(ibin); 
      cout << slabel << " " << w8 << endl;
      for (int iv = 0; iv < vNames.size(); ++iv) {
	if (slabel == vNames[iv]) {
	  cout << "filling " << w8 << " with name " << slabel << " into hist " << vNames[iv] << endl;
	  vHist[iv]->Fill(w8);
	}
      }
    }
  }

  c0.Clear();
  c0.Divide(4,4);
  for (int ic = 1; ic <= 15; ++ic) {
    c0->cd(ic); 
    if (vHist[ic-1]) {
      vHist[ic-1]->Draw();
    } else {
      cout << "vHist[" << ic-1 << "] not found" << endl;
    }
  }
  c0.SaveAs("rankingOverlays.pdf");
  cout << "done" << endl;

  c0.Clear();
  h->Draw();
}




// ----------------------------------------------------------------------
TH1D* getRanking(string fname, string prefix, int seed) {
  TH1D *h1 = new TH1D(Form("rank_%d", seed), Form("rank_%d", seed), 100, 0., 100.);
  // -- read in variable ranking from logfile
  vector<string> allLines; 
  char  buffer[2000];
  cout << "getRanking: open file " << Form("%s", fname.c_str()) << endl;
  ifstream is(Form("%s", fname.c_str()));
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  string::size_type m1, m2;
  string varn, vars, varp, tvars; 
  int bail(0), istart(0); 
  string after = Form("==> CREATE TOY DATA for seed = %d", seed);
  cout << "after: ->" << after << "<-" << endl;


  for (unsigned int i = 0; i < allLines.size(); ++i) {
    if (string::npos != allLines[i].find(after)) {
      istart = i; 
      break;
    }
  }

  cout << "start at line " << istart << endl;
  for (unsigned int i = istart; i < allLines.size(); ++i) {
    // -- method unspecific classification
    if ((string::npos != allLines[i].find(Form("--- %s", prefix.c_str())))
	&& (string::npos != allLines[i].find(": Rank : Variable "))
	) {
      bail = 0; 
      for (unsigned int j = i+2; j < i+100; ++j) {
	if (string::npos != allLines[j].find(": ---------------------------------")) {
	  bail = 1;
	  cout << "  -> breaking out " << endl;
	  break;
	}
	
	m1 = allLines[j].find(":"); 
	m2 = allLines[j].find(":", m1+1);
	varn = allLines[j].substr(m1+2, m2-m1-2); 
	m1 = m2; 
	m2 = allLines[j].find(":", m1+1);
	vars = allLines[j].substr(m1+2, m2-m1-2); 
	tvars = trim_right(vars, "\t");
	tvars = trim_right(tvars, " ");
	m1 = m2; 
	m2 = allLines[j].find(":", m1+1);
	varp = allLines[j].substr(m1+2, m2-m1-2); 
	cout << varn << "-> " << vars << " ->" << tvars << "<- " << " -> " << varp << endl;
	int ibin = atoi(varn.c_str()); 
	h1->GetXaxis()->SetBinLabel(ibin, tvars.c_str());
	h1->SetBinContent(ibin, atof(varp.c_str()));
      }
      if (1 == bail) break;
    }
  }

  return h1;
}    

inline std::string trim_right(const std::string &source , const std::string& t = " ") {
  std::string str = source;
  return str.erase( str.find_last_not_of(t) + 1);
}

inline std::string trim_left( const std::string& source, const std::string& t = " ") {
  std::string str = source;
  return str.erase(0 , source.find_first_not_of(t) );
}

inline std::string trim(const std::string& source, const std::string& t = " ") {
  std::string str = source;
  return trim_left( trim_right( str , t) , t );
} 
