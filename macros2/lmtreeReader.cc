#include "lmtreeReader.hh"

using std::cout;
using std::endl;
using std::string;
using std::vector;

lmtreeReader::lmtreeReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> lmtreeReader: constructor..." << endl;

  fForceJson = false;
  SAVECANDS = 0;

  Reco_B2MM_size = 0;
  Reco_B2JpsiK_size = 0;
  Reco_B2JpsiPhi_size = 0;
  Reco_Jpsi_size = 0;
  Reco_Muon_size = 0;

  Gen_Cand_size = 0;
  Gen_Muon_size = 0;

  trigger_name.push_back("HLT_DoubleMu3_Jpsi");
  trigger_name.push_back("HLT_Dimuon6p5_Jpsi");
  trigger_name.push_back("HLT_Dimuon6p5_Jpsi_Displaced");
  trigger_name.push_back("HLT_Dimuon7_Jpsi_Displaced");
  trigger_name.push_back("HLT_DoubleMu3p5_Jpsi_Displaced");
  trigger_name.push_back("HLT_DoubleMu4_Jpsi_Displaced");

  trigger_name.push_back("HLT_DoubleMu3_Bs");
  trigger_name.push_back("HLT_DoubleMu2_Bs");
  trigger_name.push_back("HLT_Dimuon4_Bs_Barrel");
  trigger_name.push_back("HLT_Dimuon6_Bs");
  trigger_name.push_back("HLT_DoubleMu4_Dimuon6_Bs");
  trigger_name.push_back("HLT_DoubleMu4_Dimuon4_Bs_Barrel");

  HLT = new std::vector <std::pair <std::string,bool> >;
  HLT->resize(trigger_name.size());
}

lmtreeReader::~lmtreeReader() {
  cout << "==> lmtreeReader: destructor..." << endl;
}

void lmtreeReader::readCuts(TString filename, int dump) {
  char buffer[1024];
  char cutName[128];
  float cut;
  int ok;
  FILE *cutFile = fopen(filename.Data(), "r");
  if (!cutFile) {
    cout << "cut file not existing: " << filename.Data() << " aborting..." << endl;
    abort();
  }
  if (dump) {
    cout << "===================================================" << endl;
    cout << "==> lmtreeReader: Cut file " << filename.Data() << endl;
    cout << "---------------------------------------------------" << endl;
  }

  while (fgets(buffer, sizeof(buffer), cutFile)) {
    if (buffer[strlen(buffer)-1] == '\n') buffer[strlen(buffer)-1] = '\0';
	// skip comment line
	if (buffer[0] == '#') continue;
	ok = sscanf(buffer, "%s %f", cutName, &cut);
	if(ok < 1) {
	  cout << "==> massReader: Cut not parse input line:" << endl;
	  cout << buffer << endl;
	  cout << "==> massReader: abort..." << endl;
	  abort();
	}
	if (!parseCut(cutName, cut)) {
	  cout << "==> lmtreeReader: Error parsing variable " << cutName << ". abort..." << endl;
	  abort();
	}
  }
  if (dump) cout << "---------------------------------------------------" << endl;
  if (cutFile) fclose(cutFile);
} // readCuts()

bool lmtreeReader::parseCut(char *cutName, float cut, int dump) {
  // parse the options...
  if (!strcmp(cutName,"SAVEALWAYS")) {
    SAVEALWAYS = (int)cut;
    if (dump) cout << "SAVEALWAYS: " << SAVEALWAYS << endl;
    return true;
  }
  if (!strcmp(cutName,"SAVECANDS")) {
    SAVECANDS = (int)cut;
    if (dump) cout << "SAVECANDS: " << SAVECANDS << endl;
    return true;
  }
  if (!strcmp(cutName,"MC")) {
    MC = (int)cut;
    if (dump) cout << "MC: " << MC << endl;
    return true;
  }
  if (!strcmp(cutName,"B2MMTYPE")) {
    B2MMTYPE = (int)cut;
    if (dump) cout << "B2MMTYPE: " << B2MMTYPE << endl;
    return true;
  }
  if (!strcmp(cutName,"B2JpsiKTYPE")) {
    B2JpsiKTYPE = (int)cut;
    if (dump) cout << "B2JpsiKTYPE: " << B2JpsiKTYPE << endl;
    return true;
  }
  if (!strcmp(cutName,"B2JpsiPhiTYPE")) {
    B2JpsiPhiTYPE = (int)cut;
    if (dump) cout << "B2JpsiPhiTYPE: " << B2JpsiPhiTYPE << endl;
    return true;
  }
  if (!strcmp(cutName,"NCANDS")) {
    NCANDS = (int)cut;
    if (dump) cout << "NCANDS: " << NCANDS << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDPT")) {
    CANDPT = cut;
    if (dump) cout << "CANDPT: " << CANDPT << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDETA")) {
    CANDETA = cut;
    if (dump) cout << "CANDETA: " << CANDETA << endl;
    return true;
  }
  if (!strcmp(cutName,"MUONPT")) {
    MUONPT = cut;
    if (dump) cout << "MUONPT: " << MUONPT << endl;
    return true;
  }
  if (!strcmp(cutName,"MUONETA")) {
    MUONETA = cut;
    if (dump) cout << "MUONETA: " << MUONETA << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDMINMASS")) {
    CANDMINMASS = cut;
    if (dump) cout << "CANDMINMASS: " << CANDMINMASS << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDMAXMASS")) {
    CANDMAXMASS = cut;
    if (dump) cout << "CANDMAXMASS: " << CANDMAXMASS << endl;
    return true;
  }
  if (!strcmp(cutName,"PVAVEW8")) {
    PVAVEW8 = cut;
    if (dump) cout << "PVAVEW8: " << PVAVEW8 << endl;
    return true;
  }
  if (!strcmp(cutName,"PVLIP")) {
    PVLIP = cut;
    if (dump) cout << "PVLIP: " << PVLIP << endl;
    return true;
  }
  if (!strcmp(cutName,"PVLIPS")) {
    PVLIPS = cut;
    if (dump) cout << "PVLIPS: " << PVLIPS << endl;
    return true;
  }
  if (!strcmp(cutName,"CANDCLOSETRK")) {
    CANDCLOSETRK = cut;
    if (dump) cout << "CANDCLOSETRK: " << CANDCLOSETRK << endl;
    return true;
  }
//  if (!strcmp(CutName, "JSON")) {
//    JSONFILE = cut;
//    if (dump) cout << "JSON FILE: " << JSONFILE << endl;
//    return true;
//  }

  return false;
} // parseCut()

void lmtreeReader::bookHist() {

  // create the tree
  T1 = new TTree("T1", "Reduced tree for Bs2MuMu studies by lm");

  T1->Branch("run", &fRun, "run/I");
  T1->Branch("LS", &fLS, "LS/I");
  T1->Branch("event", &fEvt, "event/I");
  T1->Branch("Reco_NPV", &Reco_NPV, "Reco_NPV/I");

  for (unsigned int i = 0; i < trigger_name.size(); ++i){
    T1->Branch(trigger_name[i].c_str(), &(HLT->at(i).second), (trigger_name[i] + "/O").c_str());
  }

  Reco_Muon_4mom  = new TClonesArray("TLorentzVector");
  Reco_Muon_PosM2 = new TClonesArray("TVector3");
  T1->Branch("Reco_Muon_size", &Reco_Muon_size, "Reco_Muon_size/I");
  T1->Branch("Reco_Muon_4mom", "TClonesArray", &Reco_Muon_4mom, 32000, 0);
  T1->Branch("Reco_Muon_PosM2", "TClonesArray", &Reco_Muon_PosM2, 32000, 0);
  T1->Branch("Reco_Muon_Charge", Reco_Muon_Charge, "Reco_Muon_Charge[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_Chi2AntiKink", Reco_Muon_Chi2AntiKink, "Reco_Muon_Chi2AntiKink[Reco_Muon_size]/D");
  T1->Branch("Reco_Muon_ID", Reco_Muon_ID, "Reco_Muon_ID[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_pdgId", Reco_Muon_pdgId, "Reco_Muon_pdgId[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_ValidHits", Reco_Muon_ValidHits, "Reco_Muon_ValidHits[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_PxValidHits", Reco_Muon_PxValidHits, "Reco_Muon_PxValidHits[Reco_Muon_size]/I");
  T1->Branch("Reco_Muon_Index", Reco_Muon_Index, "Reco_Muon_Index[Reco_Muon_size]/I");

  if (SAVECANDS) {
    Reco_B2MM_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_B2MM_size", &Reco_B2MM_size, "Reco_B2MM_size/I");
    T1->Branch("Reco_B2MM_MCTruth", Reco_B2MM_MCTruth, "Reco_B2MM_MCTruth[Reco_B2MM_size]/O");
    T1->Branch("Reco_B2MM_4mom", "TClonesArray", &Reco_B2MM_4mom, 32000, 0);
    T1->Branch("Reco_B2MM_Type", Reco_B2MM_Type, "Reco_B2MM_Type[Reco_B2MM_size]/I");
    T1->Branch("Reco_B2MM_Dxy", Reco_B2MM_Dxy, "Reco_B2MM_Dxy[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_DxyE", Reco_B2MM_DxyE, "Reco_B2MM_Dxyz[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Dxyz", Reco_B2MM_Dxyz, "Reco_B2MM_Dxyz[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_DxyzE", Reco_B2MM_DxyzE, "Reco_B2MM_DxyzE[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Alpha", Reco_B2MM_Alpha, "Reco_B2MM_Alpha[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Alphaxy", Reco_B2MM_Alphaxy, "Reco_B2MM_Alphaxy[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_VxRChi2", Reco_B2MM_VxRChi2, "Reco_B2MM_VxRChi2[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_DOCAtrk", Reco_B2MM_DOCAtrk, "Reco_B2MM_DOCAtrk[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Iso1_pt09", Reco_B2MM_Iso1_pt09, "Reco_B2MM_Iso1_pt09[Reco_B2MM_size]/D");
    T1->Branch("Reco_B2MM_Mu1_index", Reco_B2MM_Mu1_index, "Reco_B2MM_Mu1_index[Reco_B2MM_size]/I");
    T1->Branch("Reco_B2MM_Mu2_index", Reco_B2MM_Mu2_index, "Reco_B2MM_Mu2_index[Reco_B2MM_size]/I");

    Reco_B2JpsiK_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiK_Tk1_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiK_Jpsi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_B2JpsiK_size", &Reco_B2JpsiK_size, "Reco_B2JpsiK_size/I");
    T1->Branch("Reco_B2JpsiK_MCTruth", Reco_B2JpsiK_MCTruth, "Reco_B2JpsiK_MCTruth[Reco_B2JpsiK_size]/O");
    T1->Branch("Reco_B2JpsiK_4mom", "TClonesArray", &Reco_B2JpsiK_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiK_Type", Reco_B2JpsiK_Type, "Reco_B2JpsiK_Type[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Dxy", Reco_B2JpsiK_Dxy, "Reco_B2JpsiK_Dxy[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_DxyE", Reco_B2JpsiK_DxyE, "Reco_B2JpsiK_Dxyz[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Dxyz", Reco_B2JpsiK_Dxyz, "Reco_B2JpsiK_Dxyz[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_DxyzE", Reco_B2JpsiK_DxyzE, "Reco_B2JpsiK_DxyzE[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Alpha", Reco_B2JpsiK_Alpha, "Reco_B2JpsiK_Alpha[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Alphaxy", Reco_B2JpsiK_Alphaxy, "Reco_B2JpsiK_Alphaxy[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_VxRChi2", Reco_B2JpsiK_VxRChi2, "Reco_B2JpsiK_VxRChi2[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_DOCAtrk", Reco_B2JpsiK_DOCAtrk, "Reco_B2JpsiK_DOCAtrk[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Iso1_pt09", Reco_B2JpsiK_Iso1_pt09, "Reco_B2JpsiK_Iso1_pt09[Reco_B2JpsiK_size]/D");
    T1->Branch("Reco_B2JpsiK_Tk1_4mom", "TClonesArray", &Reco_B2JpsiK_Tk1_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiK_Tk1_Charge", Reco_B2JpsiK_Tk1_Charge, "Reco_B2JpsiK_Tk1_Charge[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Jpsi_4mom", "TClonesArray", &Reco_B2JpsiK_Jpsi_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiK_Mu1_index", Reco_B2JpsiK_Mu1_index, "Reco_B2JpsiK_Mu1_index[Reco_B2JpsiK_size]/I");
    T1->Branch("Reco_B2JpsiK_Mu2_index", Reco_B2JpsiK_Mu2_index, "Reco_B2JpsiK_Mu2_index[Reco_B2JpsiK_size]/I");

    Reco_B2JpsiPhi_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Tk1_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Tk2_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Jpsi_4mom = new TClonesArray("TLorentzVector");
    Reco_B2JpsiPhi_Phi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Reco_B2JpsiPhi_size", &Reco_B2JpsiPhi_size, "Reco_B2JpsiPhi_size/I");
    T1->Branch("Reco_B2JpsiPhi_MCTruth", Reco_B2JpsiPhi_MCTruth, "Reco_B2JpsiPhi_MCTruth[Reco_B2JpsiPhi_size]/O");
    T1->Branch("Reco_B2JpsiPhi_4mom", "TClonesArray", &Reco_B2JpsiPhi_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Type", Reco_B2JpsiPhi_Type, "Reco_B2JpsiPhi_Type[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Dxy", Reco_B2JpsiPhi_Dxy, "Reco_B2JpsiPhi_Dxy[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_DxyE", Reco_B2JpsiPhi_DxyE, "Reco_B2JpsiPhi_Dxyz[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Dxyz", Reco_B2JpsiPhi_Dxyz, "Reco_B2JpsiPhi_Dxyz[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_DxyzE", Reco_B2JpsiPhi_DxyzE, "Reco_B2JpsiPhi_DxyzE[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Alpha", Reco_B2JpsiPhi_Alpha, "Reco_B2JpsiPhi_Alpha[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Alphaxy", Reco_B2JpsiPhi_Alphaxy, "Reco_B2JpsiPhi_Alphaxy[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_VxRChi2", Reco_B2JpsiPhi_VxRChi2, "Reco_B2JpsiPhi_VxRChi2[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_DOCAtrk", Reco_B2JpsiPhi_DOCAtrk, "Reco_B2JpsiPhi_DOCAtrk[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Iso1_pt09", Reco_B2JpsiPhi_Iso1_pt09, "Reco_B2JpsiPhi_Iso1_pt09[Reco_B2JpsiPhi_size]/D");
    T1->Branch("Reco_B2JpsiPhi_Mu1_index", Reco_B2JpsiPhi_Mu1_index, "Reco_B2JpsiPhi_Mu1_index[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Mu2_index", Reco_B2JpsiPhi_Mu2_index, "Reco_B2JpsiPhi_Mu2_index[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Tk1_4mom", "TClonesArray", &Reco_B2JpsiPhi_Tk1_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Tk1_Charge", Reco_B2JpsiPhi_Tk1_Charge, "Reco_B2JpsiPhi_Tk1_Charge[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Tk2_4mom", "TClonesArray", &Reco_B2JpsiPhi_Tk2_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Tk2_Charge", Reco_B2JpsiPhi_Tk2_Charge, "Reco_B2JpsiPhi_Tk2_Charge[Reco_B2JpsiPhi_size]/I");
    T1->Branch("Reco_B2JpsiPhi_Jpsi_4mom", "TClonesArray", &Reco_B2JpsiPhi_Jpsi_4mom, 32000, 0);
    T1->Branch("Reco_B2JpsiPhi_Phi_4mom", "TClonesArray", &Reco_B2JpsiPhi_Phi_4mom, 32000, 0);

    Reco_Jpsi_4mom = new TClonesArray("TLorentzVector");
    /*
    T1->Branch("Reco_Jpsi_size", &Reco_Jpsi_size, "Reco_Jpsi_size/I");
    T1->Branch("Reco_Jpsi_MCTruth", Reco_Jpsi_MCTruth, "Reco_Jpsi_MCTruth[Reco_Jpsi_size]/O");
    T1->Branch("Reco_Jpsi_4mom", "TClonesArray", &Reco_Jpsi_4mom, 32000, 0);
 //   T1->Branch("Reco_Jpsi_Type", Reco_Jpsi_Type, "Reco_Jpsi_Type[Reco_Jpsi_size]/I");
    T1->Branch("Reco_Jpsi_Dxy", Reco_Jpsi_Dxy, "Reco_Jpsi_Dxy[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_DxyE", Reco_Jpsi_DxyE, "Reco_Jpsi_Dxyz[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Dxyz", Reco_Jpsi_Dxyz, "Reco_Jpsi_Dxyz[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_DxyzE", Reco_Jpsi_DxyzE, "Reco_Jpsi_DxyzE[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Alpha", Reco_Jpsi_Alpha, "Reco_Jpsi_Alpha[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Alphaxy", Reco_Jpsi_Alphaxy, "Reco_Jpsi_Alphaxy[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_VxRChi2", Reco_Jpsi_VxRChi2, "Reco_Jpsi_VxRChi2[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_DOCAtrk", Reco_Jpsi_DOCAtrk, "Reco_Jpsi_DOCAtrk[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Iso1_pt09", Reco_Jpsi_Iso1_pt09, "Reco_Jpsi_Iso1_pt09[Reco_Jpsi_size]/D");
    T1->Branch("Reco_Jpsi_Mu1_index", Reco_Jpsi_Mu1_index, "Reco_Jpsi_Mu1_index[Reco_Jpsi_size]/I");
    T1->Branch("Reco_Jpsi_Mu2_index", Reco_Jpsi_Mu2_index, "Reco_Jpsi_Mu2_index[Reco_Jpsi_size]/I");
*/
  }
  if (MC) {
    Gen_Cand_4mom  = new TClonesArray("TLorentzVector");
    Gen_Cand_Tk1_4mom  = new TClonesArray("TLorentzVector");
    Gen_Cand_Tk2_4mom  = new TClonesArray("TLorentzVector");
    T1->Branch("Gen_Cand_size", &Gen_Cand_size, "Gen_Cand_size/I");
    T1->Branch("Gen_Cand_4mom", "TClonesArray", &Gen_Cand_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Tk1_4mom", "TClonesArray", &Gen_Cand_Tk1_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Tk2_4mom", "TClonesArray", &Gen_Cand_Tk2_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Mu1_index", Gen_Cand_Mu1_index, "Gen_Cand_Mu1_index[Gen_Cand_size]/I");
    T1->Branch("Gen_Cand_Mu2_index", Gen_Cand_Mu2_index, "Gen_Cand_Mu2_index[Gen_Cand_size]/I");
    T1->Branch("Gen_Cand_proc", Gen_Cand_proc, "Gen_Cand_proc[Gen_Cand_size]/I");

    Gen_Muon_4mom  = new TClonesArray("TLorentzVector");
    Gen_Cand_Jpsi_4mom = new TClonesArray("TLorentzVector");
    Gen_Cand_Phi_4mom = new TClonesArray("TLorentzVector");
    T1->Branch("Gen_Muon_size", &Gen_Muon_size, "Gen_Muon_size/I");
    T1->Branch("Gen_Muon_4mom", "TClonesArray", &Gen_Muon_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Jpsi_4mom", "TClonesArray", &Gen_Cand_Jpsi_4mom, 32000, 0);
    T1->Branch("Gen_Cand_Phi_4mom", "TClonesArray", &Gen_Cand_Phi_4mom, 32000, 0);
    T1->Branch("Gen_Muon_Charge", Gen_Muon_Charge, "Gen_Muon_Charge[Gen_Muon_size]/I");
    T1->Branch("Gen_Muon_Index", Gen_Muon_Index, "Gen_Muon_Index[Gen_Muon_size]/I");
  }

} // massReader::bookHist()

void lmtreeReader::eventProcessing() {

  // Fill Generator Block candidates
  if (MC) fillGen();

  // fill trigger infos
  fillTrigger();

  // Fill Reconstructed Block candidates
  fillMuons();

  // Fill only if there are enough cands
  if (SAVECANDS) {
    fillRecoCand();
  }

  if (SAVEALWAYS || Reco_B2MM_size >= NCANDS || Reco_B2JpsiK_size  >= NCANDS || Reco_B2JpsiPhi_size >= NCANDS) T1->Fill();

  clearVariables();
} // massReader::eventProcessing()

void lmtreeReader::fillGen() {
  int nc = fpEvt->nGenCands();
  for (int j = 0; j < nc && j < array_size; j++) {
    TGenCand* Cand = fpEvt->getGenCand(j);
    if (abs(Cand->fID) == 531 || abs(Cand->fID) == 521) {
      if (abs(Cand->fID) == 531) {
        for (int k = Cand->fDau1; k <= Cand->fDau2; k++) {
          for (int l = k+1; l <= Cand->fDau2; l++) {
            TGenCand* Mu1 = fpEvt->getGenCand(k);
            TGenCand* Mu2 = fpEvt->getGenCand(l);
            if ( (Mu1->fID == 13 && Mu2->fID == -13) || (Mu1->fID == -13 && Mu2->fID == 13)) {
              new ((*Gen_Cand_4mom)[Gen_Cand_size])TLorentzVector(Cand->fP);
              Gen_Cand_pdgId[Gen_Cand_size] = Cand->fID;
              Gen_Cand_Mu1_index[Gen_Cand_size] = Mu1->fNumber;
              Gen_Cand_Mu2_index[Gen_Cand_size] = Mu2->fNumber;
              new ((*Gen_Cand_Tk1_4mom)[Gen_Cand_size])TLorentzVector();
              new ((*Gen_Cand_Tk2_4mom)[Gen_Cand_size])TLorentzVector();
              new ((*Gen_Cand_Jpsi_4mom)[Gen_Cand_size])TLorentzVector();
              new ((*Gen_Cand_Phi_4mom)[Gen_Cand_size])TLorentzVector();
              Gen_Cand_proc[Gen_Cand_size] = 5311313;
              Gen_Cand_size++;
            }
          }
        }
        for (int k = Cand->fDau1; k <= Cand->fDau2; k++) {
          for (int l = Cand->fDau1; l <= Cand->fDau2; l++) {
            TGenCand* jpsi = fpEvt->getGenCand(k);
            TGenCand* phi = fpEvt->getGenCand(l);
            if (abs(jpsi->fID) == 443 && abs(phi->fID) == 333) {
              for (int m = jpsi->fDau1; m <= jpsi->fDau2; m++) {
                for (int n = m+1; n <= jpsi->fDau2; n++) {
                  TGenCand* Mu1 = fpEvt->getGenCand(m);
                  TGenCand* Mu2 = fpEvt->getGenCand(n);
                  if ( (Mu1->fID == 13 && Mu2->fID == -13) || (Mu1->fID == -13 && Mu2->fID == 13)) {
                    for (int o = phi->fDau1; o <= phi->fDau2; o++) {
                      for (int p = o+1; p <= phi->fDau2; p++) {
                        TGenCand* ka1 = fpEvt->getGenCand(o);
                        TGenCand* ka2 = fpEvt->getGenCand(p);
                        if ( (ka1->fID == 321 && ka2->fID == -321) || (ka1->fID == -321 && ka2->fID == 321)) {
                          new ((*Gen_Cand_4mom)[Gen_Cand_size])TLorentzVector(Cand->fP);
                          Gen_Cand_pdgId[Gen_Cand_size] = Cand->fID;
                          Gen_Cand_Mu1_index[Gen_Cand_size] = Mu1->fNumber;
                          Gen_Cand_Mu2_index[Gen_Cand_size] = Mu2->fNumber;
                          new ((*Gen_Cand_Tk1_4mom)[Gen_Cand_size])TLorentzVector(ka1->fP);
                          new ((*Gen_Cand_Tk2_4mom)[Gen_Cand_size])TLorentzVector(ka2->fP);
                          new ((*Gen_Cand_Jpsi_4mom)[Gen_Cand_size])TLorentzVector(jpsi->fP);
                          new ((*Gen_Cand_Phi_4mom)[Gen_Cand_size])TLorentzVector(phi->fP);
                          Gen_Cand_proc[Gen_Cand_size] = 531443333;
                          Gen_Cand_size++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (abs(Cand->fID) == 521) {
        for (int k = Cand->fDau1; k <= Cand->fDau2; k++) {
          for (int l = Cand->fDau1; l <= Cand->fDau2; l++) {
            TGenCand* jpsi = fpEvt->getGenCand(k);
            TGenCand* kaon = fpEvt->getGenCand(l);
            if (abs(jpsi->fID) == 443 && abs(kaon->fID) == 321) {
              for (int m = jpsi->fDau1; m <= jpsi->fDau2; m++) {
                for (int n = m+1; n <= jpsi->fDau2; n++) {
                  TGenCand* Mu1 = fpEvt->getGenCand(m);
                  TGenCand* Mu2 = fpEvt->getGenCand(n);
                  if ( (Mu1->fID == 13 && Mu2->fID == -13) || (Mu1->fID == -13 && Mu2->fID == 13)) {
                    new ((*Gen_Cand_4mom)[Gen_Cand_size])TLorentzVector(Cand->fP);
                    Gen_Cand_pdgId[Gen_Cand_size] = Cand->fID;
                    Gen_Cand_Mu1_index[Gen_Cand_size] = Mu1->fNumber;
                    Gen_Cand_Mu2_index[Gen_Cand_size] = Mu2->fNumber;
                    new ((*Gen_Cand_Tk1_4mom)[Gen_Cand_size])TLorentzVector(kaon->fP);
                    new ((*Gen_Cand_Tk2_4mom)[Gen_Cand_size])TLorentzVector();
                    new ((*Gen_Cand_Jpsi_4mom)[Gen_Cand_size])TLorentzVector(jpsi->fP);
                    new ((*Gen_Cand_Phi_4mom)[Gen_Cand_size])TLorentzVector();
                    Gen_Cand_proc[Gen_Cand_size] = 521443321;
                    Gen_Cand_size++;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (abs(Cand->fID) == 13) {
      new ((*Gen_Muon_4mom)[Gen_Muon_size])TLorentzVector(Cand->fP);
      Gen_Muon_Index[Gen_Muon_size] = Cand->fNumber;
      if (Cand->fNumber != j) cout << "non e' lo stesso numero: fNumber = " << Cand->fNumber << " e j = " << j << endl;
      Gen_Muon_size++;
    }
  }
}

void lmtreeReader::fillTrigger() {

  size_t found;
  bool exists;
  int lenght = sizeof(fpEvt->fHLTNames) / sizeof (TString);
  for (unsigned int k = 0; k < trigger_name.size(); k++) {
    exists = false;
    for (int j = 0; j < lenght; j++) {
      string htl_name(fpEvt->fHLTNames[j].Data());
      found = htl_name.find(trigger_name[k]);
      if (found != string::npos) {
        std::pair < std::string, bool> fire(trigger_name[k], fpEvt->fHLTResult[j]);
        HLT->at(k) = fire;
        exists = true;
        break;
      }
    }
    if (exists == false) {
      std::pair < std::string, bool> fire(trigger_name[k], false);
      HLT->at(k) = fire;
    }
  }
  if (false) {
    cout << lenght << " triggers are" << endl;
    for (int j = 0; j < lenght; j++) {
      string htl_name(fpEvt->fHLTNames[j].Data());
      cout << htl_name << endl;
    }
  }
}

void lmtreeReader::fillRecoCand() {

  Reco_NPV = fpEvt->nPV();
  int nc = fpEvt->nCands();

  for (int j = 0; j < nc; j++) {
    TAnaCand* Cand = fpEvt->getCand(j);
    //if ( j != Cand->fIndex) cout << j << " = " << Cand->fIndex << " ?" << endl;
    if (selCand(Cand)) {
      if (Cand->fType == B2MMTYPE && Reco_B2MM_size < array_size)  fillB2MMCand(Cand);
      if (Cand->fType == B2JpsiKTYPE && Reco_B2JpsiK_size < array_size) fillB2JpsiKCand(Cand);
      if (Cand->fType == B2JpsiPhiTYPE && Reco_B2JpsiPhi_size < array_size) fillB2JpsiPhiCand(Cand);

      if (Reco_B2MM_size >= array_size) {
        cout << "event "<< fEvt <<" B2MM array too small: "<< array_size << " nCands = " << Reco_B2MM_size << endl;
        abort();
      }
    }

    //if (Cand->fType == 300443) {
    //  fillJpsiCand(Cand);
    //  jpsi++;
    //}
    //cout << "reco jpsi = " << jpsi << endl;

  }
  return;
}

void lmtreeReader::fillB2MMCand(TAnaCand* Cand) {

  TAnaTrack *mu0;
  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);
  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    mu0 = fpEvt->getSigTrack(it);
    if (TMath::Abs(mu0->fMCID) == 13) {
      if (mu1s == 0) mu1s = mu0;
      else mu2s = mu0;
    }
  }
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;
  if (mu1r->fQ + mu2r->fQ != 0) return;
  double mass = Cand->fMass;
  if (BLIND && mass > 5.200 && mass < 5.450) return;
  double pt   = Cand->fPlab.Pt();
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_B2MM_4mom)[Reco_B2MM_size])TLorentzVector(quadriMom);
  Reco_B2MM_Type[Reco_B2MM_size]   = Cand->fType;
  Reco_B2MM_Dxy[Reco_B2MM_size]    = Cand->fVtx.fDxy;
  Reco_B2MM_DxyE[Reco_B2MM_size]   = Cand->fVtx.fDxyE;
  Reco_B2MM_Dxyz[Reco_B2MM_size]   = Cand->fVtx.fD3d;
  Reco_B2MM_DxyzE[Reco_B2MM_size]  = Cand->fVtx.fD3dE;
  Reco_B2MM_VxRChi2[Reco_B2MM_size] = Cand->fVtx.fChi2/Cand->fVtx.fNdof;
  double fCandDocaTrk = 99.;
  if (Cand->fNstTracks.size() != 0) {
    fCandDocaTrk = Cand->fNstTracks[0].second.first;
  }
  Reco_B2MM_DOCAtrk[Reco_B2MM_size]   = fCandDocaTrk;
  double iso = isoClassicWithDOCA(Cand, 0.05, 0.7, 0.9); // 500um DOCA cut
  Reco_B2MM_Iso1_pt09[Reco_B2MM_size] = iso;
  int pvidx = (Cand->fPvIdx > -1? Cand->fPvIdx : 0);
  TVector3 svpv(Cand->fVtx.fPoint - fpEvt->getPV(pvidx)->fPoint);
  double alpha = svpv.Angle(Cand->fPlab);
  Reco_B2MM_Alpha[Reco_B2MM_size] = alpha;

  bool mc = false;
  int bs_mu1 = -1, bs_mu2 = -2;
  if (mu1r->fGenIndex >= 0 && mu2r->fGenIndex >= 0) {
    TGenCand* mu1_gen = fpEvt->getGenCand(mu1r->fGenIndex);
    TGenCand* mu2_gen = fpEvt->getGenCand(mu2r->fGenIndex);
    if ( (mu1_gen->fID == 13 && mu2_gen->fID == -13) || (mu1_gen->fID == -13 && mu2_gen->fID == 13)) {
      for (int i = mu1_gen->fMom1; i <= mu1_gen->fMom2 && i >= 0; i++) {
        TGenCand* Bs_gen = fpEvt->getGenCand(i);
        if (abs(Bs_gen->fID) == 531) {
          bs_mu1 = i;
          break;
        }
      }
      for (int i = mu2_gen->fMom1; i <= mu2_gen->fMom2 && i >= 0; i++) {
        TGenCand* Bs_gen = fpEvt->getGenCand(i);
        if (abs(Bs_gen->fID) == 531) {
          bs_mu2 = i;
          break;
        }
      }
    }
    if (bs_mu1 == bs_mu2) mc = true;
  }
  Reco_B2MM_MCTruth[Reco_B2MM_size] = mc;

  Reco_B2MM_Mu1_index[Reco_B2MM_size] = mu1_index;
  Reco_B2MM_Mu2_index[Reco_B2MM_size] = mu2_index;

  Reco_B2MM_size++;

  return;
}

void lmtreeReader::fillB2JpsiKCand(TAnaCand* Cand) {

  bool children[3] = {false, false, false}; // k, mu+, mu-
  int counter = 0;
  int zerocounter = 0;
  TAnaTrack *sigtrack0;
  TAnaTrack *kappas(0), *kappar(0);
  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);

  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    zerocounter++;
    sigtrack0 = fpEvt->getSigTrack(it);
    if (abs(sigtrack0->fMCID) == 321) {
      kappas = sigtrack0;
      counter++;
      children[0] = true;
    }
    if (abs(sigtrack0->fMCID) == 13) {
      if (mu1s == 0) mu1s = sigtrack0;
      else mu2s = sigtrack0;
      counter++;
      if (sigtrack0->fQ == 1) children[1] = true;
      if (sigtrack0->fQ == -1) children[2] = true;
    }
  }
  if (counter != 3) {cout << "counter " << counter << "zer0counter " << counter << endl; abort();}
  if (MC) {
    for (int i = 0; i < 3; i++) {
      if (children[i] == false) {
        cout << "skipping not true Bu2JpsiK, event " << fEvt << endl;
        return;
      }
    }
  }
  kappar = fpEvt->getRecTrack(kappas->fIndex);
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);

  if (!selTrack(kappar)) return;
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;
  if (mu1r->fQ + mu2r->fQ != 0) return;
  
  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_B2JpsiK_4mom)[Reco_B2JpsiK_size])TLorentzVector(quadriMom);
  Reco_B2JpsiK_Type[Reco_B2JpsiK_size]   = Cand->fType;
  Reco_B2JpsiK_Dxy[Reco_B2JpsiK_size]    = Cand->fVtx.fDxy;
  Reco_B2JpsiK_DxyE[Reco_B2JpsiK_size]   = Cand->fVtx.fDxyE;
  Reco_B2JpsiK_Dxyz[Reco_B2JpsiK_size]   = Cand->fVtx.fD3d;
  Reco_B2JpsiK_DxyzE[Reco_B2JpsiK_size]  = Cand->fVtx.fD3dE;
  Reco_B2JpsiK_VxRChi2[Reco_B2JpsiK_size] = Cand->fVtx.fChi2/Cand->fVtx.fNdof;
  double fCandDocaTrk = 99.;
  if (Cand->fNstTracks.size() != 0) {
    fCandDocaTrk = Cand->fNstTracks[0].second.first;
  }
  Reco_B2JpsiK_DOCAtrk[Reco_B2JpsiK_size]   = fCandDocaTrk;
  double iso  = isoClassicWithDOCA(Cand, 0.05, 0.5, 0.5); // 500um DOCA cut
  Reco_B2JpsiK_Iso1_pt09[Reco_B2JpsiK_size] = iso;
  int pvidx = (Cand->fPvIdx > -1? Cand->fPvIdx : 0);
  TVector3 svpv(Cand->fVtx.fPoint - fpEvt->getPV(pvidx)->fPoint);
  double alpha = svpv.Angle(Cand->fPlab);
  Reco_B2JpsiK_Alpha[Reco_B2JpsiK_size] = alpha;
  Reco_B2JpsiK_Mu1_index[Reco_B2JpsiK_size] = mu1_index;
  Reco_B2JpsiK_Mu2_index[Reco_B2JpsiK_size] = mu2_index;

  pt = kappar->fPlab.Pt();
  eta = kappar->fPlab.Eta();
  phi = kappar->fPlab.Phi();
  quadriMom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  new ((*Reco_B2JpsiK_Tk1_4mom)[Reco_B2JpsiK_size])TLorentzVector(quadriMom);
  Reco_B2JpsiK_Tk1_Charge[Reco_B2JpsiK_size] = kappar->fQ;

  pt   = mu1r->fPlab.Pt();
  eta  = mu1r->fPlab.Eta();
  phi  = mu1r->fPlab.Phi();
  TLorentzVector m1_4mom;
  m1_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  pt   = mu2r->fPlab.Pt();
  eta  = mu2r->fPlab.Eta();
  phi  = mu2r->fPlab.Phi();
  TLorentzVector m2_4mom;
  m2_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  TLorentzVector jpsi;
  jpsi = m1_4mom + m2_4mom;
  new ((*Reco_B2JpsiK_Jpsi_4mom)[Reco_B2JpsiK_size])TLorentzVector(jpsi);

  bool mc = false;
  int jpsi_mu1 = -1, jpsi_mu2 = -2, bu_jpsi = -3, bu_k = -4;
  if (mu1r->fGenIndex >= 0 && mu2r->fGenIndex >= 0 && kappar->fGenIndex >=0 ) {
    TGenCand* mu1_gen = fpEvt->getGenCand(mu1r->fGenIndex);
    TGenCand* mu2_gen = fpEvt->getGenCand(mu2r->fGenIndex);
    if ( (mu1_gen->fID == 13 && mu2_gen->fID == -13) || (mu1_gen->fID == -13 && mu2_gen->fID == 13)) {
      for (int i = mu1_gen->fMom1; i <= mu1_gen->fMom2 && i >= 0; i++) {
        TGenCand* Jpsi_gen = fpEvt->getGenCand(i);
        if (abs(Jpsi_gen->fID) == 443) {
          jpsi_mu1 = i;
          break;
        }
      }
      for (int i = mu2_gen->fMom1; i <= mu2_gen->fMom2 && i >= 0; i++) {
        TGenCand* Jpsi_gen = fpEvt->getGenCand(i);
        if (abs(Jpsi_gen->fID) == 443) {
          jpsi_mu2 = i;
          break;
        }
      }
      if (jpsi_mu1 == jpsi_mu2) {
        TGenCand* Jpsi_gen = fpEvt->getGenCand(jpsi_mu1);
        for (int i = Jpsi_gen->fMom1; i <= Jpsi_gen->fMom2 && i >= 0; i++) {
          TGenCand* Bu_gen = fpEvt->getGenCand(i);
          if (abs(Bu_gen->fID) == 521) {
            bu_jpsi = i;
            break;
          }
        }
        TGenCand* k_gen = fpEvt->getGenCand(kappar->fGenIndex);
        if ( fabs(k_gen->fID) == 321) {
          for (int i = k_gen->fMom1; i <= k_gen->fMom2 && i >= 0; i++) {
            TGenCand* Bu_gen = fpEvt->getGenCand(i);
            if (abs(Bu_gen->fID) == 521) {
              if (k_gen->fQ == Bu_gen->fQ) {
                bu_k = i;
                break;
              }
            }
          }
        }
      }
      if (bu_jpsi == bu_k) mc = true;
    }
  }
  Reco_B2JpsiK_MCTruth[Reco_B2JpsiK_size] = mc;

  Reco_B2JpsiK_size++;

  return;
}

void lmtreeReader::fillB2JpsiPhiCand(TAnaCand* Cand) {

  int counter = 0;
  int zerocounter = 0;
  TAnaTrack *sigtrack0;
  TAnaTrack *kappa1s(0), *kappa1r(0);
  TAnaTrack *kappa2s(0), *kappa2r(0);

  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);

  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    zerocounter++;
    sigtrack0 = fpEvt->getSigTrack(it);
    if (abs(sigtrack0->fMCID) == 321) {
      if (kappa1s == 0) kappa1s = sigtrack0;
      else kappa2s = sigtrack0;
      counter++;
    }
    if (abs(sigtrack0->fMCID) == 13) {
      if (mu1s == 0) mu1s = sigtrack0;
      else mu2s = sigtrack0;
      counter++;
    }
  }
  if (counter != 4) {cout << "Jpsi Phi counter " << counter << "zer0counter " << counter << endl; abort();}
  kappa1r = fpEvt->getRecTrack(kappa1s->fIndex);
  kappa2r = fpEvt->getRecTrack(kappa2s->fIndex);
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);

  if (!selTrack(kappa1r)) return;
  if (!selTrack(kappa2r)) return;
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;
  if (mu1r->fQ + mu2r->fQ != 0) return;
  if (kappa1r->fQ + kappa2r->fQ != 0) return;

  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_B2JpsiPhi_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(quadriMom);
  Reco_B2JpsiPhi_Type[Reco_B2JpsiPhi_size]   = Cand->fType;
  Reco_B2JpsiPhi_Dxy[Reco_B2JpsiPhi_size]    = Cand->fVtx.fDxy;
  Reco_B2JpsiPhi_DxyE[Reco_B2JpsiPhi_size]   = Cand->fVtx.fDxyE;
  Reco_B2JpsiPhi_Dxyz[Reco_B2JpsiPhi_size]   = Cand->fVtx.fD3d;
  Reco_B2JpsiPhi_DxyzE[Reco_B2JpsiPhi_size]  = Cand->fVtx.fD3dE;
  Reco_B2JpsiPhi_VxRChi2[Reco_B2JpsiPhi_size] = Cand->fVtx.fChi2/Cand->fVtx.fNdof;
  double fCandDocaTrk = 99.;
  if (Cand->fNstTracks.size() != 0) {
    fCandDocaTrk = Cand->fNstTracks[0].second.first;
  }
  Reco_B2JpsiPhi_DOCAtrk[Reco_B2JpsiPhi_size]   = fCandDocaTrk;
  double iso  = isoClassicWithDOCA(Cand, 0.05, 0.5, 0.5); // 500um DOCA cut
  Reco_B2JpsiPhi_Iso1_pt09[Reco_B2JpsiPhi_size] = iso;
  int pvidx = (Cand->fPvIdx > -1? Cand->fPvIdx : 0);
  TVector3 svpv(Cand->fVtx.fPoint - fpEvt->getPV(pvidx)->fPoint);
  double alpha = svpv.Angle(Cand->fPlab);
  Reco_B2JpsiPhi_Alpha[Reco_B2JpsiPhi_size] = alpha;
  Reco_B2JpsiPhi_Mu1_index[Reco_B2JpsiPhi_size] = mu1_index;
  Reco_B2JpsiPhi_Mu2_index[Reco_B2JpsiPhi_size] = mu2_index;

  pt = kappa1r->fPlab.Pt();
  eta = kappa1r->fPlab.Eta();
  phi = kappa1r->fPlab.Phi();
  quadriMom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  new ((*Reco_B2JpsiPhi_Tk1_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(quadriMom);
  Reco_B2JpsiPhi_Tk1_Charge[Reco_B2JpsiPhi_size] = kappa1r->fQ;
  pt = kappa2r->fPlab.Pt();
  eta = kappa2r->fPlab.Eta();
  phi = kappa2r->fPlab.Phi();
  quadriMom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  new ((*Reco_B2JpsiPhi_Tk2_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(quadriMom);
  Reco_B2JpsiPhi_Tk2_Charge[Reco_B2JpsiPhi_size] = kappa2r->fQ;

  pt   = mu1r->fPlab.Pt();
  eta  = mu1r->fPlab.Eta();
  phi  = mu1r->fPlab.Phi();
  TLorentzVector m1_4mom;
  m1_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  pt   = mu2r->fPlab.Pt();
  eta  = mu2r->fPlab.Eta();
  phi  = mu2r->fPlab.Phi();
  TLorentzVector m2_4mom;
  m2_4mom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  TLorentzVector jpsi;
  jpsi = m1_4mom + m2_4mom;
  new ((*Reco_B2JpsiPhi_Jpsi_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(jpsi);
  pt   = kappa1r->fPlab.Pt();
  eta  = kappa1r->fPlab.Eta();
  phi  = kappa1r->fPlab.Phi();
  TLorentzVector kappa1_4mom;
  kappa1_4mom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  pt   = kappa2r->fPlab.Pt();
  eta  = kappa2r->fPlab.Eta();
  phi  = kappa2r->fPlab.Phi();
  TLorentzVector kappa2_4mom;
  kappa2_4mom.SetPtEtaPhiM(pt, eta, phi, KAON_MASS);
  TLorentzVector Phi;
  Phi = kappa1_4mom + kappa2_4mom;
  new ((*Reco_B2JpsiPhi_Phi_4mom)[Reco_B2JpsiPhi_size])TLorentzVector(Phi);

  bool mc = false;
  int jpsi_mu1 = -1, jpsi_mu2 = -2, phi_k1 = -3, phi_k2 = -4, bs_jpsi = -5, bs_phi = -6;
  if (mu1r->fGenIndex >= 0 && mu2r->fGenIndex >= 0 && kappa1r->fGenIndex >=0 && kappa2r->fGenIndex >=0) {
    TGenCand* mu1_gen = fpEvt->getGenCand(mu1r->fGenIndex);
    TGenCand* mu2_gen = fpEvt->getGenCand(mu2r->fGenIndex);
    if ( (mu1_gen->fID == 13 && mu2_gen->fID == -13) || (mu1_gen->fID == -13 && mu2_gen->fID == 13)) {
      for (int i = mu1_gen->fMom1; i <= mu1_gen->fMom2 && i >= 0; i++) {
        TGenCand* Jpsi_gen = fpEvt->getGenCand(i);
        if (abs(Jpsi_gen->fID) == 443) {
          jpsi_mu1 = i;
          break;
        }
      }
      for (int i = mu2_gen->fMom1; i <= mu2_gen->fMom2 && i >= 0; i++) {
        TGenCand* Jpsi_gen = fpEvt->getGenCand(i);
        if (abs(Jpsi_gen->fID) == 443) {
          jpsi_mu2 = i;
          break;
        }
      }
      if (jpsi_mu1 == jpsi_mu2) {
        TGenCand* k1_gen = fpEvt->getGenCand(kappa1r->fGenIndex);
        TGenCand* k2_gen = fpEvt->getGenCand(kappa2r->fGenIndex);
        if ( (k1_gen->fID == 321 && k2_gen->fID == -321) || (k1_gen->fID == -321 && k2_gen->fID == 321)) {
          for (int i = k1_gen->fMom1; i <= k1_gen->fMom2 && i >= 0; i++) {
            TGenCand* phi_gen = fpEvt->getGenCand(i);
            if (abs(phi_gen->fID) == 333) {
              phi_k1 = i;
              break;
            }
          }
          for (int i = k2_gen->fMom1; i <= k2_gen->fMom2 && i >= 0; i++) {
            TGenCand* phi_gen = fpEvt->getGenCand(i);
            if (abs(phi_gen->fID) == 333) {
              phi_k2 = i;
              break;
            }
          }
          if (phi_k1 == phi_k2) {
            TGenCand* Jpsi_gen = fpEvt->getGenCand(jpsi_mu1);
            for (int i = Jpsi_gen->fMom1; i <= Jpsi_gen->fMom2 && i >= 0; i++) {
              TGenCand* Bs_gen = fpEvt->getGenCand(i);
              if (abs(Bs_gen->fID) == 531) {
                bs_jpsi = i;
                break;
              }
            }
            TGenCand* phi_gen = fpEvt->getGenCand(phi_k1);
            for (int i = phi_gen->fMom1; i <= phi_gen->fMom2 && i >= 0; i++) {
              TGenCand* Bs_gen = fpEvt->getGenCand(i);
              if (abs(Bs_gen->fID) == 531) {
                bs_phi = i;
                break;
              }
            }
          }
        }
        if (bs_jpsi == bs_phi) mc = true;
      }
    }
  }
  Reco_B2JpsiPhi_MCTruth[Reco_B2JpsiPhi_size] = mc;

  Reco_B2JpsiPhi_size++;

}

void lmtreeReader::fillJpsiCand(TAnaCand* Cand) {

  TAnaTrack *mu0;
  TAnaTrack *mu1s(0), *mu1r(0);
  TAnaTrack *mu2s(0), *mu2r(0);
  for (int it = Cand->fSig1; it <= Cand->fSig2; ++it) {
    mu0 = fpEvt->getSigTrack(it);
    if (TMath::Abs(mu0->fMCID) == 13) {
      if (mu1s == 0) mu1s = mu0;
      else mu2s = mu0;
    }
  }
  int mu1_index = mu1s->fIndex;
  int mu2_index = mu2s->fIndex;
  mu1r = fpEvt->getRecTrack(mu1_index);
  mu2r = fpEvt->getRecTrack(mu2_index);
  if (!selMuon(mu1r)) return;
  if (!selMuon(mu2r)) return;

  double mass = Cand->fMass;
  double pt   = Cand->fPlab.Pt();// cout << pt << endl;
  double eta  = Cand->fPlab.Eta();
  double phi  = Cand->fPlab.Phi();
  TLorentzVector quadriMom;
  quadriMom.SetPtEtaPhiM(pt, eta, phi, mass);
  new ((*Reco_Jpsi_4mom)[Reco_Jpsi_size])TLorentzVector(quadriMom);
  //Reco_Jpsi_Type[Reco_Jpsi_size]   = Cand->fType;
  Reco_Jpsi_Dxy[Reco_Jpsi_size]    = Cand->fVtx.fDxy;
  Reco_Jpsi_DxyE[Reco_Jpsi_size]   = Cand->fVtx.fDxyE;
  Reco_Jpsi_Dxyz[Reco_Jpsi_size]   = Cand->fVtx.fD3d;
  Reco_Jpsi_DxyzE[Reco_Jpsi_size]  = Cand->fVtx.fD3dE;
  Reco_Jpsi_VxRChi2[Reco_Jpsi_size] = Cand->fVtx.fChi2/Cand->fVtx.fNdof;
  double fCandDocaTrk = 99.;
  if (Cand->fNstTracks.size() != 0) {
    fCandDocaTrk = Cand->fNstTracks[0].second.first;
  }
  Reco_Jpsi_DOCAtrk[Reco_Jpsi_size]   = fCandDocaTrk;
  Reco_Jpsi_Mu1_index[Reco_Jpsi_size] = mu1_index;
  Reco_Jpsi_Mu2_index[Reco_Jpsi_size] = mu2_index;

  Reco_Jpsi_size++;

  return;
}

void lmtreeReader::fillMuons() {

  int nm = fpEvt->nRecTracks();
  for (int j = 0; j < nm && Reco_Muon_size < Muon_array_size; j++) {
    TAnaTrack* Track = fpEvt->getRecTrack(j);
    if (selMuon(Track)) {
      Reco_Muon_Index[Reco_Muon_size] = j;
      Reco_Muon_ValidHits[Reco_Muon_size] = Track->fValidHits;
      Reco_Muon_PxValidHits[Reco_Muon_size] = numberOfPixLayers(Track);
      double pt   = Track->fPlab.Pt();
      double eta  = Track->fPlab.Eta();
      double phi  = Track->fPlab.Phi();
      TLorentzVector quadriMom;
      quadriMom.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
      new ((*Reco_Muon_4mom)[Reco_Muon_size])TLorentzVector(quadriMom);

      Reco_Muon_ID[Reco_Muon_size] = Track->fMuID;
      Reco_Muon_Charge[Reco_Muon_size] = Track->fQ;
      Reco_Muon_pdgId[Reco_Muon_size] = Track->fMCID;

      TVector3 triVec;
      triVec.SetXYZ(0.0, 0.0, 0.0);
      double muonchi2 = 999999.9;
      if (Track->fMuIndex >= 0 && Track->fMuID >= 0) {
        TAnaMuon* Muon = fpEvt->getMuon(Track->fMuIndex);
        double x = Muon->fPositionAtM2.X();
        double y = Muon->fPositionAtM2.Y();
        double z = Muon->fPositionAtM2.Z();
        triVec.SetXYZ(x, y, z);
        muonchi2 = Muon->fMuonChi2;
      }
      new ((*Reco_Muon_PosM2)[Reco_Muon_size])TVector3(triVec);
      Reco_Muon_Chi2AntiKink[Reco_Muon_size] = muonchi2;

      Reco_Muon_size++;
      if (Reco_Muon_size >= Muon_array_size) {
        cout << "muon array too small: "<< Reco_Muon_size << " nMuons = " << nm << endl;
        abort();
      }
    }
  }
  return;
}

bool lmtreeReader::selCand(TAnaCand* Candi) {
  if (! ( (Candi->fType == B2MMTYPE) || (Candi->fType == B2JpsiKTYPE) || (Candi->fType == B2JpsiPhiTYPE))) return false;
  if (! (Candi->fPlab.Pt() > CANDPT)) return false;
  if (! (fabs(Candi->fPlab.Eta()) < CANDETA)) return false;
  if (! (Candi->fMass > CANDMINMASS && Candi->fMass < CANDMAXMASS) ) return false;
  TAnaVertex *pv = fpEvt->getPV(Candi->fPvIdx);
  double PvNtrk = pv->getNtracks();
  double PvNdof = pv->fNdof;
  double PvAveW8 = ((PvNdof+2.)/2.)/PvNtrk;
  if (! (PvAveW8 > PVAVEW8)) return false;
  if (! (fabs(Candi->fPvLip) < PVLIP)) return false;
  if (! (fabs(Candi->fPvLip /Candi->fPvLipE) < PVLIPS)) return false;
  if (! (nCloseTracks(Candi, 0.03, 0.5) < CANDCLOSETRK)) return false;
  return true;
}

bool lmtreeReader::selMuon(TAnaTrack* Muon) {
  if (! ((Muon->fTrackQuality & 0x1<<2) == 0x1<<2)) return false;  // track is HP
  //if (! (Muon->fMuIndex >= 0 && Muon->fMuID >= 0)) return false;  // track is a muon too
  //if (! ((Muon->fMuID & 6) == 6)) return false;  // muon is global and tracker
  if (! (Muon->fPlab.Pt() > MUONPT)) return false;
  if (! (fabs(Muon->fPlab.Eta()) < MUONETA)) return false;
  return true;
}

bool lmtreeReader::selTrack(TAnaTrack* Track) {

  // check the candidate
  if (! ((Track->fTrackQuality & 0x1<<2) == 0x1<<2)) return false;  // track is HP
  if (! (Track->fPlab.Pt() > 0.5 )) return false; // track pt > 0.5 GeV/c
  if (! (Track->fTip < 1.0)) return false; // transverse impact parameter
  if (! (Track->fLip < 25.0)) return false; // longitudinal impact parameter

  return true;
}

void lmtreeReader::clearVariables() {

  Reco_NPV = 0;

  if (MC) {
    if (Gen_Cand_size > 0) {
      Gen_Cand_4mom->Clear();
      Gen_Cand_Phi_4mom->Clear();
      Gen_Cand_Jpsi_4mom->Clear();
      Gen_Cand_Tk1_4mom->Clear();
      Gen_Cand_Tk2_4mom->Clear();
    }
    if (Gen_Muon_size > 0) Gen_Muon_4mom->Clear();
    Gen_Cand_size = 0;
    Gen_Muon_size = 0;
  }

  if (SAVECANDS) {
    if (Reco_B2MM_size > 0) Reco_B2MM_4mom->Clear();
    if (Reco_B2JpsiK_size > 0) {
      Reco_B2JpsiK_4mom->Clear();
      Reco_B2JpsiK_Tk1_4mom->Clear();
      Reco_B2JpsiK_Jpsi_4mom->Clear();
    }
    if (Reco_B2JpsiPhi_size > 0) {
      Reco_B2JpsiPhi_4mom->Clear();
      Reco_B2JpsiPhi_Tk1_4mom->Clear();
      Reco_B2JpsiPhi_Tk2_4mom->Clear();
      Reco_B2JpsiPhi_Jpsi_4mom->Clear();
      Reco_B2JpsiPhi_Phi_4mom->Clear();
    }
    if (Reco_Jpsi_size > 0) Reco_Jpsi_4mom->Clear();
    Reco_B2MM_size = 0;
    Reco_B2JpsiK_size = 0;
    Reco_B2JpsiPhi_size = 0;
    Reco_Jpsi_size = 0;
  }

  if (Reco_Muon_size > 0) {
    Reco_Muon_4mom->Clear();
    Reco_Muon_PosM2->Clear();
  }
  Reco_Muon_size = 0;

}

void lmtreeReader::closeHistFile() {

  fpHistFile = T1->GetCurrentFile();
  treeReader01::closeHistFile();

} // massReader::closeHistFile()

double lmtreeReader::isoClassicWithDOCA(TAnaCand *pC, double docaCut, double r, double ptmin) {
  const double ptCut(ptmin), coneSize(r);
  const bool verbose(false);

  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.);
  TAnaTrack *pT;
  vector<int> cIdx, pIdx;
  int pvIdx = pC->fPvIdx;
  //int pvIdx2= nearestPV(pvIdx, 0.1);
  //if (TMath::Abs(fCandPvLipS2) > 2) pvIdx2 = -1;

  bool sameOrCloseVertex(false);

  //fCandI0trk = 0;
  //fCandI1trk = 0;
  //fCandI2trk = 0;

  if (verbose) cout << "Looking at cand " << pC->fType << " with " << pC->fSig2 - pC->fSig1 + 1 << " sig tracks"
            << " from PV = " << pvIdx
            << endl;
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i);
    if (verbose) cout << " track idx = " << pT->fIndex << " with ID = " << pT->fMCID << endl;
    cIdx.push_back(pT->fIndex);
    candPtScalar += pT->fPlab.Perp();
    if (verbose) {
      int tIdx = fpEvt->getRecTrack(pT->fIndex)->fPvIdx;
      if (pvIdx != tIdx) {
        cout << "Signal track pointing to PV " << tIdx << " instead of " << pvIdx << endl;
      }
    }
  }

  candPt = pC->fPlab.Perp();

  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nRecTracks(); ++i) {
    pT = fpEvt->getRecTrack(i);
    if (verbose) {
      cout << "   track " << i
           << " with pT = " << pT->fPlab.Perp()
           << " eta = " << pT->fPlab.Eta()
           << " pointing at PV " << pT->fPvIdx;
    }

//     // -- check that any track associated with a definitive vertex is from the same PV
//     if (pT->fPvIdx > -1) {
//       sameOrCloseVertex = (pT->fPvIdx == pvIdx);
//       //       sameOrCloseVertex = (pT->fPvIdx == pvIdx) || (pT->fPvIdx == pvIdx2);
//       if (!sameOrCloseVertex) continue;
//     }

    if (pT->fPvIdx != pvIdx) {
      if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
      continue;
    }


    pt = pT->fPlab.Perp();
    if (pt < ptCut) {
      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  {
      if (verbose) cout << " skipped because it is a sig track " << endl;
      continue;
    }
    if (pT->fPlab.DeltaR(pC->fPlab) < coneSize) {
      pIdx.push_back(i);
      //++fCandI0trk;
      sumPt += pt;
      if (verbose) cout << endl;
    }
    else {
      if (verbose) cout << " skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
    }
  }

  // -- Now consider the DOCA tracks
  int nsize = pC->fNstTracks.size();
  if (nsize>0) {
    for(int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;
      // double docaE = pC->fNstTracks[i].second.second;

      if(doca > docaCut) continue; // check the doca cut

      pT = fpEvt->getRecTrack(trkId);

//       // -- check that any track associated with a definitive vertex is from the same or the closest other PV
//       if (pT->fPvIdx > -1) {
// 	sameOrCloseVertex = (pT->fPvIdx == pvIdx);
// 	//	sameOrCloseVertex = (pT->fPvIdx == pvIdx) || (pT->fPvIdx == pvIdx2);
// 	if (!sameOrCloseVertex) continue;
//       }

	  if ((pT->fPvIdx > -1) && (pT->fPvIdx != pvIdx)) {
	if (verbose) cout << " doca track " << trkId << " skipped because it is from a different PV " << pT->fPvIdx <<endl;
	continue;
	  }

      pt = pT->fPlab.Perp();
      if (pt < ptCut) {
    if (verbose) cout << " doca track " << trkId << " skipped because of pt = " << pt << endl;
    continue;
      }

	  if (pT->fPlab.DeltaR(pC->fPlab) > coneSize) {
	if (verbose) cout << " doca track " << trkId << " skipped because of deltaR = " << pT->fPlab.DeltaR(pC->fPlab) << endl;
	continue;
	  }

      // -- Skip tracks already included above
      if (pIdx.end() != find(pIdx.begin(), pIdx.end(), trkId))  continue;
      if (cIdx.end() != find(cIdx.begin(), cIdx.end(), trkId))  continue;

      //++fCandI1trk;
      sumPt += pt;
      if (verbose) cout << " doca track " << trkId << " included "<<doca<<" "<<pt<<endl;

    } // for loop over tracks
  } // end if

  //fCandI2trk = fCandI0trk + fCandI1trk;

  iso = candPt/(candPt + sumPt);

  //   if (verbose) cout << "--> iso = " << candPt << " .. " << sumPt << " = " << iso << endl;
  //   if (verbose) cout << "--> iso = " << pC->fPlab.Perp() << " .. " << sumPt << " = " << pC->fPlab.Perp()/(pC->fPlab.Perp() + sumPt) << endl;

  return iso;

}

int lmtreeReader::nCloseTracks(TAnaCand *pC, double dcaCut, double ptCut) {
  int cnt(0);
  int nsize = pC->fNstTracks.size();
  int pvIdx = pC->fPvIdx;

  TAnaTrack *pT;
  double pt(0.);
  if (nsize > 0) {
    for (int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;

      if (doca > dcaCut) continue; // check the doca cut

      pT = fpEvt->getRecTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((pT->fPvIdx > -1) && (pT->fPvIdx != pvIdx)) continue;

      pt = pT->fPlab.Perp();
      if (pt < ptCut) continue;

      ++cnt;
    }
  }
  return cnt;
}
