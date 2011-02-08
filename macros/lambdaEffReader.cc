#include "lambdaEffReader.hh"
#include "../interface/HFMasses.hh"

#include <cmath>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;

// ----------------------------------------------------------------------
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/lambdaEffReader.default.cuts
//           ./runTreeReaders -f test.root
// ----------------------------------------------------------------------

/*! Strips off trailing characters from a string
  /param instring String where characters should be cut off
  /param symbol Symbol after which text should be cut off, including the symbol
  /return remaining string
  */
std::string lambdaEffReader::stripOff(std::string instring, char symbol)
{
    std::string::iterator iter = std::find(instring.begin(), instring.end(), symbol);
    std::string outstring;
    std::copy(instring.begin(), iter, std::back_inserter(outstring));
    return outstring;
}

//==========================================================================

lambdaEffReader::lambdaEffReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName)
{
    cout << "==> lambdaEffReader: constructor..." << endl;
}

// ----------------------------------------------------------------------
lambdaEffReader::~lambdaEffReader()
{
    cout << "==> lambdaEffReader: destructor..." << endl;
}

// ----------------------------------------------------------------------
void lambdaEffReader::startAnalysis()
{
    cout << "==> lambdaEffReader: Starting analysis..." << endl;
    // book reduced trees
    bookReducedTree();
}

// ----------------------------------------------------------------------
void lambdaEffReader::endAnalysis()
{
}

// ----------------------------------------------------------------------
void lambdaEffReader::eventProcessing()
{
    //cout << "------ next event " << endl;

    ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());
    ((TH1D*)fpHistFile->Get("h2"))->Fill(fpEvt->nCands());

    TAnaCand *tacCur;

    const int nGenCands(fpEvt->nGenCands());
    for(int gci=0; gci!=nGenCands; gci++)
    {
	// initialise tree variables
	fMom5122 = fJpsiL0 = fJpsiL0excl = false;
	fL0PrPi = fL0PrPiexcl = false;
	fJpsiMuMu = fJpsiMuMuexcl = false;
	fL0matched = false;
	fgenL0vtxR = fgenL0vtxZ = 9999;
	fgenprpt = fgenpreta = fgenprphi = fgenpry = 9999;
	fgenpipt = fgenpieta = fgenpiphi = fgenpiy = 9999;
	fgenL0pt = fgenL0eta = fgenL0phi = fgenL0y = 9999;
	frecprpt = frecpreta = frecprphi = frecpry = 9999;
	frecpipt = frecpieta = frecpiphi = frecpiy = 9999;
	frecL0pt = frecL0eta = frecL0phi = frecL0y = 9999;
	// get current candidate
	TGenCand *gcCur = fpEvt->getGenCand(gci);
	if (3122==abs(gcCur->fID))
	{
	    // check if the Lambda0 comes from a LambdaB and detect if it is in the proper channel Jpsi+Lambda
	    const int mom5122id = getMom(gcCur, 5122);
	    fMom5122 = (mom5122id > 1);
	    if(fMom5122)
	    {
		TGenCand *gcMom5122 = fpEvt->getGenCand(mom5122id);
		if (getDau(gcMom5122, 3122)>-1 && getDau(gcMom5122, 443))
		{
		    fJpsiL0 = true;
		    if (2==getNDau(gcMom5122))
			fJpsiL0excl = true;
		}
		// now check if the Jpsi decays in two muons
		if(getDau(gcMom5122, 443) > -1)
		{
		    TGenCand *gcMom443 = fpEvt->getGenCand(getDau(gcMom5122, 443));
		    if (getDau(gcMom443,13,1)>-1 && getDau(gcMom443,13,2)>-1)
		    {
		        fJpsiMuMu = true;
			if(2==getNDau(gcMom443)) fJpsiMuMuexcl = true;
		    }
		}
	    }
	    // now check if the Lambda decays in pr+pi
	    if(getDau(gcCur, 211) > -1 && getDau(gcCur, 2212) > -1)
	    {
		fL0PrPi = true;
		if (2==getNDau(gcCur)) fL0PrPiexcl = true;
		// and now we try to find the reconstructed event
		//const TVector3 genVtx = gcCur->fV; would be wrong as this is the production vtx of the L0
		const TGenCand *gcPr = fpEvt->getGenCand(getDau(gcCur, 2212));
		const TVector3 genVtx = gcPr->fV;
		// Take the radius and z of the L0 vertex in absolut coordinates as this will be used
		// for detector acceptance cuts, which are defined by the geometry of the tracker
		fgenL0vtxR = genVtx.XYvector().Mod();
		fgenL0vtxZ = genVtx.Z();
		const TLorentzVector tlvGenPr = gcPr->fP;
		const TGenCand *gcPi = fpEvt->getGenCand(getDau(gcCur, 211));
		const TLorentzVector tlvGenPi = gcPi->fP;
		setKinematicVariables(gcPr->fP, fgenprpt, fgenpreta, fgenprphi, fgenpry);
		setKinematicVariables(gcPi->fP, fgenpipt, fgenpieta, fgenpiphi, fgenpiy);
		setKinematicVariables(gcCur->fP, fgenL0pt, fgenL0eta, fgenL0phi, fgenL0y);

		for (int rci=0; rci!=fpEvt->nCands(); rci++)
		{
		    tacCur = fpEvt->getCand(rci);
		    if(tacCur->fType == 503122)
		    {
			const int iPr = getSigId(tacCur,2212,1);
			const int iPi = getSigId(tacCur,211,1);

			if(iPr > -1 && iPi > -1)
			{
			    const TVector3 recVtx = tacCur->fVtx.fPoint;

			    const TAnaTrack *tatPr = fpEvt->getSigTrack(iPr);
			    TLorentzVector tlvSigPr;
			    tlvSigPr.SetVectM(tatPr->fPlab,MPROTON);

			    const TAnaTrack *tatPi = fpEvt->getSigTrack(iPi);
			    TLorentzVector tlvSigPi;
			    tlvSigPi.SetVectM(tatPi->fPlab,MPROTON);

			    //if(recVtx.DeltaR(genVtx) < .3 && tlvSigPr.DeltaR(tlvGenPr) < .1 && tlvSigPi.DeltaR(tlvGenPi) < .2)
			    const double rRatio = recVtx.Mag() /  genVtx.Mag();
			    if(tlvSigPr.DeltaR(tlvGenPr) < .1 && tlvSigPi.DeltaR(tlvGenPi) < .2 && rRatio > 0.77 && rRatio < 1.3)
			    {
				// set tree variables
				fL0matched = true;
				setKinematicVariables(tlvSigPr, frecprpt, frecpreta, frecprphi, frecpry);
				setKinematicVariables(tlvSigPi, frecpipt, frecpieta, frecpiphi, frecpiy);
				TLorentzVector tlvSigL0;
				tlvSigL0.SetVectM(tacCur->fPlab, MLAMBDA_0);
				setKinematicVariables(tlvSigL0, frecL0pt, frecL0eta, frecL0phi, frecL0y);

				// fill histos
				((TH1D*)fpHistFile->Get("htm_vtxdeltaR"))->Fill(recVtx.DeltaR(genVtx));
				((TH1D*)fpHistFile->Get("htm_vtxdR"))->Fill(recVtx.Mag() -  genVtx.Mag());
				((TH1D*)fpHistFile->Get("htm_vtxratioR"))->Fill(recVtx.Mag() /  genVtx.Mag());
				((TH1D*)fpHistFile->Get("htm_vtx2d"))->Fill((recVtx.XYvector() - genVtx.XYvector()).Mod());
				((TH1D*)fpHistFile->Get("htm_vtx3d"))->Fill((recVtx - genVtx).Mag());

				((TH1D*)fpHistFile->Get("htm_prdpt"))->Fill(tlvSigPr.Pt() - tlvGenPr.Pt());
				((TH1D*)fpHistFile->Get("htm_prdphi"))->Fill(tlvSigPr.Phi() - tlvGenPr.Phi());
				((TH1D*)fpHistFile->Get("htm_prdeta"))->Fill(tlvSigPr.Eta() - tlvGenPr.Eta());
				((TH1D*)fpHistFile->Get("htm_prdeltaR"))->Fill(tlvSigPr.DeltaR(tlvGenPr));

				((TH1D*)fpHistFile->Get("htm_pidpt"))->Fill(tlvSigPi.Pt() - tlvGenPi.Pt());
				((TH1D*)fpHistFile->Get("htm_pidphi"))->Fill(tlvSigPi.Phi() - tlvGenPi.Phi());
				((TH1D*)fpHistFile->Get("htm_pideta"))->Fill(tlvSigPi.Eta() - tlvGenPi.Eta());
				((TH1D*)fpHistFile->Get("htm_pideltaR"))->Fill(tlvSigPi.DeltaR(tlvGenPi));
			    }
			}
		    }
		}
		fTree->Fill();
	    }
	}
    }
}

// ----------------------------------------------------------------------
// get index of a mother particle with certain id.
// returns -1 if no such particle present
// pos=1 (the default) looks for the first such particle
// pos=2 looks for the second and so on
int lambdaEffReader::getMom(TGenCand *cand, int id, int pos)
{
    int ret(-1);
    const int nGenCands(fpEvt->nGenCands());
    if(cand->fMom1 > -1 && cand->fMom2 > -1 && cand->fMom1 < nGenCands && cand->fMom2 < nGenCands)
    {
	for(int i=cand->fMom1; i<=cand->fMom2; i++)
	{
	    TGenCand *gcCur = fpEvt->getGenCand(i);
	    if (id == abs(gcCur->fID))
	    {
		pos--;
		if (0==pos) ret = i;
	    }
	}
    }
    return ret;
}

// ----------------------------------------------------------------------
// get index of a daughter particle with certain id.
// returns -1 if no such particle present
// pos=1 (the default) looks for the first such particle
// pos=2 looks for the second and so on
int lambdaEffReader::getDau(TGenCand *cand, int id, int pos)
{
    int ret(-1);
    const int nGenCands(fpEvt->nGenCands());
    if(cand->fDau1 > -1 && cand->fDau2 > -1 && cand->fDau1 < nGenCands && cand->fDau2 < nGenCands)
    {
	for(int i=cand->fDau1; i<=cand->fDau2; i++)
	{
	    TGenCand *gcCur = fpEvt->getGenCand(i);
	    if (id == abs(gcCur->fID))
	    {
		// check if this is a direct daughter of the particle
		if (getMom(gcCur, abs(cand->fID)) > -1)
		{
		    pos--;
		    if (0==pos) ret = i;
		}
	    }
	}
    }
    return ret;
}

// ----------------------------------------------------------------------
// returns the number of direct daughters of a particle in gen block
// for particles decayed by PYTHIA, a simple difference would suffice
// particles decayed by GEANT have indirect daughters within a daughter block,
// so we need to loop, check and count
int lambdaEffReader::getNDau(TGenCand* cand)
{
    int ret(0);
    const int nGenCands(fpEvt->nGenCands());
    if(cand->fDau1 > -1 && cand->fDau2 > -1 && cand->fDau1 < nGenCands && cand->fDau2 < nGenCands)
    {
	for (int i=cand->fDau1; i<=cand->fDau2; i++)
	{
	    TGenCand *gcCur = fpEvt->getGenCand(i);
	    if (getMom(gcCur, abs(cand->fID)) > -1)
		ret++;
	}
    }
    return ret;
}

// ----------------------------------------------------------------------
int lambdaEffReader::getCheckedGenIndex(const TAnaCand *tac, const int id, int pos)
{
    return getCheckedGenIndex(getSigId(tac,id,pos));
}

int lambdaEffReader::getCheckedGenIndex(const int sigid)
{
    if(sigid < 0 || sigid >= fpEvt->nSigTracks()) return -2;
    const TAnaTrack * tat = fpEvt->getSigTrack(sigid);
    int ret = tat->fGenIndex;
    if (ret < 0) return -1;
    if (ret >= fpEvt->nGenCands()) return -3;
    return ret;
}

// ----------------------------------------------------------------------
int lambdaEffReader::getCandId(const TAnaCand *tac, const int id, int pos)
{
    if (tac->fDau1<0 || tac->fDau2<0) return -1;
    int ret(-1);
    for(int i = tac->fDau1; i <= tac->fDau2; i++)
    {
        const TAnaCand* tacDau = fpEvt->getCand(i);
        const int tacDauId = tacDau->fType - int(tacDau->fType / 10000) * 10000;
        if(TMath::Abs(tacDauId) == id) pos--;
        if(0 == pos)
        {
            ret = i;
            break;
        }
    }
    return ret;
}

// ----------------------------------------------------------------------
int lambdaEffReader::getSigId(const TAnaCand *tac, const int id, int pos)
{
    if (tac->fSig1<0 || tac->fSig2<0) return -1;
    int ret(-1);
    for(int i = tac->fSig1; i <= tac->fSig2; i++)
    {
        const TAnaTrack* tat = fpEvt->getSigTrack(i);
        if(TMath::Abs(tat->fMCID) == id) pos--;
        if(0 == pos)
        {
            ret = i;
            break;
        }
    }
    return ret;
}

// ----------------------------------------------------------------------
void lambdaEffReader::setKinematicVariables(TLorentzVector tlv, double &pt, double &eta, double &phi, double &rapidity)
{
    pt = tlv.Pt();
    eta = tlv.Eta();
    phi = tlv.Phi();
    rapidity = tlv.Rapidity();
}

// ----------------------------------------------------------------------
void lambdaEffReader::doL1stuff()
{
    fL1TDMu0 = fL1TDMu3 = false;
    fL1TMuBH = fL1TMu0 = fL1TMu3 = fL1TMu5 = fL1TMu7 = fL1TMu10 = fL1TMu14 = fL1TMu20 = false;
    for(int i=0; i!=NL1T; i++)
    {
	const std::string name = fpEvt->fL1TNames[i].Data();
	//const double prescale = fpEvt->fL1TPrescale[i];
	const bool result = fpEvt->fL1TResult[i];
	const bool error = fpEvt->fL1TError[i];
	//cout << name << " - prescale: " << prescale << " result: " << result << " error: " << error << endl;
	if (result && !error)
	{
	    if ("L1_DoubleMuOpen"== name) fL1TDMu0 = true;
	    if ("L1_DoubleMu3"== name) fL1TDMu3 = true;
	    if ("L1_SingleMuBeamHalo"== name) fL1TMuBH = true;
	    if ("L1_SingleMu0"== name) fL1TMu0 = true;
	    if ("L1_SingleMu3"== name) fL1TMu3 = true;
	    if ("L1_SingleMu5"== name) fL1TMu5 = true;
	    if ("L1_SingleMu7"== name) fL1TMu7 = true;
	    if ("L1_SingleMu10"== name) fL1TMu10 = true;
	    if ("L1_SingleMu14"== name) fL1TMu14 = true;
	    if ("L1_SingleMu20"== name) fL1TMu20 = true;
	}
    }
}

void lambdaEffReader::doHLTstuff()
{
    fHLTqrk = fHLTqrkLS = false;
    fHLTDMuOp = fHLTDMu0 = fHLTDMu3 = fHLTDMu5 = false;
    fHLTMu3t3jp = fHLTMu3t5jp = fHLTMu5t0jp = false;
    fHLTMu0jp = fHLTMu0jpT = fHLTMu3jp = fHLTMu3jpT = fHLTMu5jp = fHLTMu5jpT = false;
    fHLTL1DMu0 = fHLTL2DMu0 = fHLTL2Mu0 = false;
    for(int i=0; i!=NHLT; i++)
    {
	const std::string name = fpEvt->fHLTNames[i].Data();
	if (name == "") break; // if name is empty we reached the end of the list (according to Urs)
	const double prescale = fpEvt->fHLTPrescale[i];
	const bool wasrun = fpEvt->fHLTWasRun[i];
	const bool result = fpEvt->fHLTResult[i];
	const bool error = fpEvt->fHLTError[i];
	//cout << name << " - wasrun: " << wasrun << " prescale: " << prescale << " result: " << result << " error: " << error << endl;
	if (result && wasrun &&!error)
	{
	    // Quarkonium triggers
	    if ("HLT_DoubleMu0_Quarkonium_v1" == name)    fHLTqrk = true;
	    if ("HLT_DoubleMu0_Quarkonium_LS_v1" == name) fHLTqrkLS = true;
	    // DoubleMu triggers (usually prescaled)
	    if ("HLT_DoubleMuOpen"       == name) fHLTDMuOp = true;
	    if ("HLT_DoubleMu0"          == name) fHLTDMu0 = true;
	    if ("HLT_DoubleMu3_v2"       == name) fHLTDMu3 = true;
	    if ("HLT_DoubleMu5_v1"       == name) fHLTDMu5 = true;
	    // Jpsi triggers
	    if ("HLT_Mu3_Track3_Jpsi"    == name) fHLTMu3t3jp = true;
	    if ("HLT_Mu3_Track3_Jpsi_v2" == name) fHLTMu3t3jp = true;
	    if ("HLT_Mu3_Track3_Jpsi_v3" == name) fHLTMu3t3jp = true;
	    if ("HLT_Mu3_Track5_Jpsi_v1" == name) fHLTMu3t5jp = true;
	    if ("HLT_Mu3_Track5_Jpsi_v2" == name) fHLTMu3t5jp = true;
	    if ("HLT_Mu3_Track5_Jpsi_v3" == name) fHLTMu3t5jp = true;
	    if ("HLT_Mu5_Track0_Jpsi"    == name) fHLTMu5t0jp = true;
	    if ("HLT_Mu5_Track0_Jpsi_v2" == name) fHLTMu5t0jp = true;
	    // Jpsi with OST
	    if ("HLT_Mu0_TkMu0_OST_Jpsi" == name) fHLTMu0jp = true;
	    if ("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1" == name) fHLTMu0jpT = true;
	    if ("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2" == name) fHLTMu0jpT = true;
	    if ("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3" == name) fHLTMu0jpT = true;
	    if ("HLT_Mu3_TkMu0_OST_Jpsi" == name) fHLTMu3jp = true;
	    if ("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2" == name) fHLTMu3jpT = true;
	    if ("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3" == name) fHLTMu3jpT = true;
	    if ("HLT_Mu5_TkMu0_OST_Jpsi" == name) fHLTMu5jp = true;
	    if ("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1" == name) fHLTMu5jpT = true;
	    if ("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2" == name) fHLTMu5jpT = true;
	    // L1 passthrough triggers
	    if ("HLT_L1DoubleMuOpen"     == name) fHLTL1DMu0 = true;
	    if ("HLT_L2DoubleMu0"        == name) fHLTL2DMu0 = true;
	    if ("HLT_Mu5_L2Mu0"          == name) fHLTL2Mu0 = true;
	}
    }
}

// ----------------------------------------------------------------------
void lambdaEffReader::initVariables()
{

}

// ----------------------------------------------------------------------
void lambdaEffReader::fillHist()
{
}

// ----------------------------------------------------------------------
void lambdaEffReader::bookHist()
{
    cout << "==> lambdaEffReader: bookHist " << endl;

    TH1D *h;
    h = new TH1D("h1", "Ntrk", 500, 0., 1000.);
    h = new TH1D("h2", "NCand", 20, 0., 20.);
    h = new TH1D("h3", "cand ID", 1000100, -100., 1000000.);

    h = new TH1D("htm_vtxdeltaR", "vtx #Delta R"       ,  100, 0. , 0.32);
    h = new TH1D("htm_vtxdR",     "vtx r_{rec}-r_{gen}",  100, -3. ,3.);
    h = new TH1D("htm_vtxratioR", "vtx r_{rec}/r_{gen}",  100, 0.5, 1.5);
    h = new TH1D("htm_vtx2d",     "vtx 2d distance",      100, 0. ,10.);
    h = new TH1D("htm_vtx3d",     "vtx 3d distance",      100, 0. , 5.);

    h = new TH1D("htm_prdpt",     "p: #Delta p_{T}",      100, -1. ,1.);
    h = new TH1D("htm_prdphi",    "p: #Delta #phi",       100, -0.2,0.2);
    h = new TH1D("htm_prdeta",    "p: #Delta #eta",       100, -0.2,0.2);
    h = new TH1D("htm_prdeltaR",  "p: #Delta R"   ,       100, 0. ,0.12);

    h = new TH1D("htm_pidpt",     "#pi: #Delta p_{T}",    100, -1. ,1.);
    h = new TH1D("htm_pidphi",    "#pi: #Delta #phi",     100, -0.2,0.2);
    h = new TH1D("htm_pideta",    "#pi: #Delta #eta",     100, -0.2,0.2);
    h = new TH1D("htm_pideltaR",  "#pi: #Delta R"   ,     100, 0. ,0.22);

    return;
}


// ----------------------------------------------------------------------
void lambdaEffReader::bookReducedTree()
{
    cout << "==> lambdaEffReader: bookReducedTree" << endl;

    // create the events tree ======================================
    fTree = new TTree("events", "events");

    // run info
    fTree->Branch("run",     &fRun,     "run/I");
    fTree->Branch("event",   &fEvt,   "event/I");
    fTree->Branch("LS",      &fLS,      "LS/I");

    fTree->Branch("Mom5122", &fMom5122, "Mom5122/O");
    fTree->Branch("JpsiL0", &fJpsiL0, "JpsiL0/O");
    fTree->Branch("JpsiL0excl", &fJpsiL0excl, "JpsiL0excl/O");

    fTree->Branch("L0PrPi", &fL0PrPi, "L0PrPi/O");
    fTree->Branch("L0PrPiexcl", &fL0PrPiexcl, "L0PrPiexcl/O");
    fTree->Branch("JpsiMuMu", &fJpsiMuMu, "JpsiMuMu/O");
    fTree->Branch("JpsiMuMuexcl", &fJpsiMuMuexcl, "JpsiMuMuexcl/O");
    fTree->Branch("L0matched", &fL0matched, "L0matched/O");

    fTree->Branch("genL0vtxR", &fgenL0vtxR, "genL0vtxR/D");
    fTree->Branch("genL0vtxZ", &fgenL0vtxZ, "genL0vtxZ/D");

    fTree->Branch("genprpt", &fgenprpt, "genprpt/D");
    fTree->Branch("genpreta", &fgenpreta, "genpreta/D");
    fTree->Branch("genprphi", &fgenprphi, "genprphi/D");
    fTree->Branch("genpry", &fgenpry, "genpry/D");

    fTree->Branch("genpipt", &fgenpipt, "genpipt/D");
    fTree->Branch("genpieta", &fgenpieta, "genpieta/D");
    fTree->Branch("genpiphi", &fgenpiphi, "genpiphi/D");
    fTree->Branch("genpiy", &fgenpiy, "genpiy/D");

    fTree->Branch("genL0pt", &fgenL0pt, "genL0pt/D");
    fTree->Branch("genL0eta", &fgenL0eta, "genL0eta/D");
    fTree->Branch("genL0phi", &fgenL0phi, "genL0phi/D");
    fTree->Branch("genL0y", &fgenL0y, "genL0y/D");

    fTree->Branch("recprpt", &frecprpt, "recprpt/D");
    fTree->Branch("recpreta", &frecpreta, "recpreta/D");
    fTree->Branch("recprphi", &frecprphi, "recprphi/D");
    fTree->Branch("recpry", &frecpry, "recpry/D");

    fTree->Branch("recpipt", &frecpipt, "recpipt/D");
    fTree->Branch("recpieta", &frecpieta, "recpieta/D");
    fTree->Branch("recpiphi", &frecpiphi, "recpiphi/D");
    fTree->Branch("recpiy", &frecpiy, "recpiy/D");

    fTree->Branch("recL0pt", &frecL0pt, "recL0pt/D");
    fTree->Branch("recL0eta", &frecL0eta, "recL0eta/D");
    fTree->Branch("recL0phi", &frecL0phi, "recL0phi/D");
    fTree->Branch("recL0y", &frecL0y, "recL0y/D");
    
}

// ----------------------------------------------------------------------
void lambdaEffReader::readCuts(TString filename, int dump)
{
    // read cut file
    fCutFile = filename;
    if (dump) cout << "==> lambdaEffReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
    ifstream infile(fCutFile.Data());

    TString fn(fCutFile.Data());

    if (dump)
    {
        cout << "====================================" << endl;
        cout << "==> lambdaEffReader: Cut file  " << fCutFile.Data() << endl;
        cout << "------------------------------------" << endl;
    }

    TH1D *hcuts = new TH1D("hcuts", "", 40, 0., 40.);
    hcuts->GetXaxis()->SetBinLabel(1, fn.Data());

    // loop over lines of input file
    std::string strLine;
    while (std::getline(infile,strLine))
    {
        // strip off everything after a hash or a slash
        strLine = stripOff(strLine, '#');
        strLine = stripOff(strLine, '/');
        std::istringstream iss(strLine, std::istringstream::in);
        // intermediate storage of the entries
        std::string key, value;
        if (iss >> key >> value)
        {
            if (dump) cout << key << ": " << value << endl;
            // a switch statement doesn't work as C++ does not allow to take strings as switch argument
            // http://www.codeguru.com/cpp/cpp/cpp_mfc/article.php/c4067
            if("CUTLbCandidate" == key) // The candidate to retrieve signal values
            {
                setCut(CUTLbCandidate, value);
                //setCut(CUTLbCandidate, value, hcuts, 1, "Candidate");
                continue;
            }

            cout << "==> lambdaEffReader: ERROR in cutfile: Don't know what to do with " << key << "=" << value << endl;
        }
    }

    if (dump)  cout << "------------------------------------" << endl;

}
