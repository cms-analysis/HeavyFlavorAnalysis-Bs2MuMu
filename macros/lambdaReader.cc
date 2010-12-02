#include "lambdaReader.hh"
#include "../interface/HFMasses.hh"

#include "TRandom.h"
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
// Run with: ./runTreeReaders -c chains/bg-test -D root -C cuts/lambdaReader.default.cuts
//           ./runTreeReaders -f test.root
// ----------------------------------------------------------------------

/*! Converts simple types to strings
  \param i variable to be converted
  \return string
  */
template <typename T>
std::string toString(T i)
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
}

/*! Strips off trailing characters from a string
  /param instring String where characters should be cut off
  /param symbol Symbol after which text should be cut off, including the symbol
  /return remaining string
  */
std::string stripOff(std::string instring, char symbol)
{
    std::string::iterator iter = std::find(instring.begin(), instring.end(), symbol);
    std::string outstring;
    std::copy(instring.begin(), iter, std::back_inserter(outstring));
    return outstring;
}

//==========================================================================
decayMap::decayMap() { };

decayMap::pos_t decayMap::getPos(std::string key)
{
    pos_t ret;
    map_t::iterator pos = myMap.find(key);
    if(pos!=myMap.end())
    {
        ret=pos->second;
    }
    else
    {
        ret=myMap.size();
        myMap[key]=ret;
    }
    return ret;
}

void decayMap::printMap()
{
    createRevMap();
    cout << "map contains " << myMap.size() << " entries" << endl;
    cout << "revmap contains " << myRevMap.size() << " entries" << endl;
    for(revmap_t::const_iterator it = myRevMap.begin(); it!=myRevMap.end(); it++)
    {
        cout << it->first << " - " << it->second << endl;
    }
}

void decayMap::createRevMap()
{
    myRevMap.empty();
    for(map_t::const_iterator it = myMap.begin(); it!=myMap.end(); it++)
    {
        myRevMap[it->second]=it->first;
    }
}

unsigned int decayMap::readFile(std::string filename)
{
    std::ifstream infile(filename.c_str());
    if(!infile.good())
    {
	throw std::invalid_argument("File " + filename + " not found.");
    }
    std::string s;
    unsigned int counter(0);
    while (getline(infile,s))
    {
	getPos(s);
	counter++;
    }
    return counter;
}

//==========================================================================
class findDaughters
{
public:
    findDaughters() : nDaughters(0) {};
    findDaughters(TGenCand* curParticle, TAna01Event* fpEvt) : nDaughters(0)
    {
        strDaughters = "";
        strDauGrdau = toString(abs(curParticle->fID)) + " -> ";
        //strDauGrdau = "";
        const int nGenCands(fpEvt->nGenCands());
        mother = curParticle;
        if(mother->fDau1 > -1 && mother->fDau1 < nGenCands && mother->fDau2 > -1 && mother->fDau2 < nGenCands && mother->fDau1 < mother->fDau2)
        {
            nDaughters = 0;
            nDauGrdau = 0;
            for(int dauit=mother->fDau1; dauit<=mother->fDau2; dauit++)
            {	// <= is ok because fDau2 is id of last daughter in list
                TGenCand *gcDau = fpEvt->getGenCand(dauit);
                const unsigned int id = abs(gcDau->fID);
                // check if these particles come from the right mother
                if(gcDau->fMom1==mother->fNumber&&gcDau->fMom2==mother->fNumber)
                {
                    nDaughters++;
                    strDaughters += toString(id) + " ";
                    findDaughters grdau(gcDau,fpEvt);
                    strDauGrdau += toString(id) + " (" + grdau.getDaughters() + ") " ;
                    strDauGrdaus += toString(id) + " (" + grdau.getDauGrdau() + ") " ;
                    nDauGrdau += grdau.getNDaughters();
                }
                // else: particle comes from hadronic interactions or decays in flight
                // discard this for the moment
            }
        }
    };
    std::string getDaughters()
    {
        return strDaughters;
    };
    std::string getDauGrdau()
    {
        return strDauGrdau;
    };
    std::string getDauGrdaus()
    {
        return strDauGrdaus;
    };
    unsigned int getNDaughters()
    {
        return nDaughters;
    };
    unsigned int getNDauGrdau()
    {
        return nDauGrdau;
    };

private:
    TGenCand* mother;
    std::string strDaughters;
    std::string strDauGrdau;
    std::string strDauGrdaus;
    unsigned int nDaughters; // no. of direct daughters
    unsigned int nDauGrdau; // no. of daughters of the direct daughters
};

//==========================================================================

lambdaReader::lambdaReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName)
{
    cout << "==> lambdaReader: constructor..." << endl;
}

// ----------------------------------------------------------------------
lambdaReader::~lambdaReader()
{
    cout << "==> lambdaReader: destructor..." << endl;
}

// ----------------------------------------------------------------------
void lambdaReader::startAnalysis()
{
    cout << "==> lambdaReader: Starting analysis..." << endl;
    // fill the map with some anchor values
    if(CUTReadDecayMaps)
    {
	cout << "Reading " << CUTDecayMap1 << " - " << myDecayMap1Gen.readFile(CUTDecayMap1.c_str()) << " entries read." << endl;
	cout << "Reading " << CUTDecayMap2 << " - " << myDecayMap2Gen.readFile(CUTDecayMap2.c_str()) << " entries read." << endl;
    }
    // book reduced trees
    bookReducedTree();
}

// ----------------------------------------------------------------------
void lambdaReader::endAnalysis()
{
    if (CUTPrintDecayMaps)
    {
	cout << "Decay map 1Gen: " << endl;
        myDecayMap1Gen.printMap();
	cout << "Decay map 2Gen: " << endl;
        myDecayMap2Gen.printMap();
    }
}

// ----------------------------------------------------------------------
void lambdaReader::eventProcessing()
{

    ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());
    ((TH1D*)fpHistFile->Get("h2"))->Fill(fpEvt->nCands());

    TAnaCand *pCand;

    // initialize a flag to mark in the genTree if there exists a corresponding candidate from signal
    fghasCand = 0;

    // MC only
    if (fIsMC) doGenLevelStuff();

    // main loop
    for (int iC = 0; iC < fpEvt->nCands(); ++iC)
    {
        pCand = fpEvt->getCand(iC);
        ((TH1D*)fpHistFile->Get("h3"))->Fill(pCand->fType);

        if (CUTLbCandidate == pCand->fType)
        {
	    // ---------------------------------------------------------
	    // First we get some handier objects for later calculations
            TAnaCand *tacJpsi, *tacLambda0;
	    if (pCand->fDau1 < 0 || pCand->fDau2 < 0)
	    {
		cout << "Problem: pCand->fDau1: " << pCand->fDau1 << " pCand->fDau2: " << pCand->fDau2
		    << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
		continue;
	    }
            tacJpsi = fpEvt->getCand(pCand->fDau1);
            tacLambda0 = fpEvt->getCand(pCand->fDau2);

            // get the muons
	    if (tacJpsi->fSig1 < 0 || tacJpsi->fSig2 < 0)
	    {
		cout << "Problem: tacJpsi->fSig1: " << tacJpsi->fSig1 << " pCand->fDaui->fSig2: " << tacJpsi->fSig2
		    << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
		continue;
	    }
	    const TAnaTrack *tatSigMu1 = fpEvt->getSigTrack(tacJpsi->fSig1);
	    const TAnaTrack *tatSigMu2 = fpEvt->getSigTrack(tacJpsi->fSig2);
	    if (tatSigMu1->fIndex < 0 || tatSigMu2->fIndex < 0)
	    {
		cout << "Problem: tatSigMu1->fIndex: " << tatSigMu1->fIndex << " tatSigMu2->fIndex: " << tatSigMu2->fIndex 
		    << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
		continue;
	    }
            const TAnaTrack *tatRecMu1 = fpEvt->getRecTrack(tatSigMu1->fIndex);
            const TAnaTrack *tatRecMu2 = fpEvt->getRecTrack(tatSigMu2->fIndex);

	    // require the muons to be of at least some quality
	    if(tatRecMu1->fMuID==-1||(tatRecMu1->fMuID & CUTMuId1)!=CUTMuId1) continue;
	    if(tatRecMu2->fMuID==-1||(tatRecMu2->fMuID & CUTMuId2)!=CUTMuId2) continue;

            // get the lambda daughters
	    if (tacLambda0->fSig1 < 0 || tacLambda0->fSig2 < 0)
	    {
		cout << "Problem: tacLambda0->fSig1: " << tacLambda0->fSig1 << " tacLambda0->fSig2: " << tacLambda0->fSig2
		    << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
		continue;
	    }
	    const TAnaTrack *tatSigPi = fpEvt->getSigTrack(tacLambda0->fSig1);
	    const TAnaTrack *tatSigPr = fpEvt->getSigTrack(tacLambda0->fSig2);
	    if (tatSigPi->fIndex < 0 || tatSigPr->fIndex < 0)
	    {
		cout << "Problem: tatSigPr->fIndex: " << tatSigPr->fIndex << " tatSigPi->fIndex: " << tatSigPi->fIndex
		    << " Run: " << fRun << " Event: " << fEvent << " LS: " << fLS << " -- skipping " << endl;
		continue;
	    }
            const TAnaTrack *tatRecPi = fpEvt->getRecTrack(tatSigPi->fIndex);
            const TAnaTrack *tatRecPr = fpEvt->getRecTrack(tatSigPr->fIndex);

            // get some TLorentzVector for eta and phi (and reuse it later)
            TLorentzVector tlvRecMu1, tlvRecMu2, tlvRecPr, tlvRecPi;
            tlvRecMu1.SetVectM(tatRecMu1->fPlab, MMUON);
            tlvRecMu2.SetVectM(tatRecMu2->fPlab, MMUON);
            tlvRecPr.SetVectM(tatRecPr->fPlab, MPROTON);
            tlvRecPi.SetVectM(tatRecPi->fPlab, MPION);

            TLorentzVector tlvLambdaB, tlvLambda0, tlvJpsi;
            tlvLambdaB.SetVectM(pCand->fPlab, MLAMBDA_B);
            tlvLambda0.SetVectM(tacLambda0->fPlab, MLAMBDA_0);
            tlvJpsi.SetVectM(tacJpsi->fPlab, MJPSI);

	    // sigtracks
	    TLorentzVector tlvSigMu1, tlvSigMu2, tlvSigPr, tlvSigPi;
	    tlvSigMu1.SetVectM(tatSigMu1->fPlab, MMUON);
	    tlvSigMu2.SetVectM(tatSigMu2->fPlab, MMUON);
	    tlvSigPr.SetVectM(tatSigPr->fPlab,MPROTON);
	    tlvSigPi.SetVectM(tatSigPi->fPlab,MPION);

	    // ---------------------------------------------------------
	    // now we do some calculations

            // calculate alpha (pointing angles in rads)
            const TVector3 vecJpsiL0 = tacLambda0->fVtx.fPoint - tacJpsi->fVtx.fPoint;
	    const TVector3 vecPV = fpEvt->getPV(pCand->fPvIdx)->fPoint;
	    const TVector3 vecPVLb = tacJpsi->fVtx.fPoint - vecPV;
            const double alphal0 = vecJpsiL0.Angle(tacLambda0->fPlab);
	    const double alphalb = vecPVLb.Angle(pCand->fPlab);

	    // calculate beta factor as 3d vector
	    const TVector3 vecBetaL0 = vecJpsiL0 * (1. / tlvLambda0.E());
	    const TVector3 vecBetaLb = vecPVLb   * (1. / tlvLambdaB.E());

            // K_s hypothesis
            TLorentzVector tlvPrAsPiHypo;
            tlvPrAsPiHypo.SetVectM(tatRecPr->fPlab, MPION);
            const TLorentzVector tlvKsHypoth = tlvRecPi+tlvPrAsPiHypo;

	    fghasCand = 1; // obviously we have a candidate, so we want to mark this

	    // ---------------------------------------------------------
            // now comes the boring part: filling of the reduced tree

	    // candidate data from fit
            fmlb=pCand->fMass;
            fptlb=pCand->fPlab.Perp();
            fplb=pCand->fPlab.Mag();
            fetalb=tlvLambdaB.Eta();
            fphilb=tlvLambdaB.Phi();

            fml0=tacLambda0->fMass;
            fptl0=tacLambda0->fPlab.Perp();
            fpl0=tacLambda0->fPlab.Mag();
            fetal0=tlvLambda0.Eta();
            fphil0=tlvLambda0.Phi();

            fmjp=tacJpsi->fMass;
            fptjp=tacJpsi->fPlab.Perp();
            fpjp=tacJpsi->fPlab.Mag();
            fetajp=tlvJpsi.Eta();
            fphijp=tlvJpsi.Phi();

	    // reco tracks
            frpt1m=tatRecMu1->fPlab.Perp();
            freta1m=tlvRecMu1.Eta();
            frphi1m=tlvRecMu1.Phi();
            frq1m=tatRecMu1->fQ;
            if(tatRecMu1->fMuID==-1) frid1m=0;
            else frid1m=tatRecMu1->fMuID;
	    fchi21m=tatRecMu1->fChi2;
	    fprob1m=TMath::Prob(tatRecMu1->fChi2,tatRecMu1->fDof);
	    fndof1m=tatRecMu1->fDof;
	    fqual1m=tatRecMu1->fTrackQuality;

            frpt2m=tatRecMu2->fPlab.Perp();
            freta2m=tlvRecMu2.Eta();
            frphi2m=tlvRecMu2.Phi();
            frq2m=tatRecMu2->fQ;
            if(tatRecMu2->fMuID==-1) frid2m=0;
            else frid2m=tatRecMu2->fMuID;
	    fchi22m=tatRecMu2->fChi2;
	    fprob2m=TMath::Prob(tatRecMu2->fChi2,tatRecMu2->fDof);
	    fndof2m=tatRecMu2->fDof;
	    fqual2m=tatRecMu2->fTrackQuality;

            frptpr=tatRecPr->fPlab.Perp();
            fretapr=tlvRecPr.Eta();
            frphipr=tlvRecPr.Phi();
            frqpr=tatRecPr->fQ;
	    fchi2pr=tatRecPr->fChi2;
	    fprobpr=TMath::Prob(tatRecPr->fChi2,tatRecPr->fDof);
	    fndofpr=tatRecPr->fDof;
	    fqualpr=tatRecPr->fTrackQuality;

            frptpi=tatRecPi->fPlab.Perp();
            fretapi=tlvRecPi.Eta();
            frphipi=tlvRecPi.Phi();
            frqpi=tatRecPi->fQ;
	    fchi2pi=tatRecPi->fChi2;
	    fprobpi=TMath::Prob(tatRecPi->fChi2,tatRecPi->fDof);
	    fndofpi=tatRecPi->fDof;
	    fqualpi=tatRecPi->fTrackQuality;

	    // refitted tracks, "signal" tracks
	    fSpt1m=tlvSigMu1.Perp();
	    fSeta1m=tlvSigMu1.Eta();
	    fSphi1m=tlvSigMu1.Phi();

	    fSpt2m=tlvSigMu2.Perp();
	    fSeta2m=tlvSigMu2.Eta();
	    fSphi2m=tlvSigMu2.Phi();

	    fSptpr=tlvSigPr.Perp();
	    fSetapr=tlvSigPr.Eta();
	    fSphipr=tlvSigPr.Phi();

	    fSptpi=tlvSigPi.Perp();
	    fSetapi=tlvSigPi.Eta();
	    fSphipi=tlvSigPi.Phi();

	    // flight lengths and their errors
            fd3lb=pCand->fVtx.fD3d;
            fd3l0=tacLambda0->fVtx.fD3d;
            fd3jp=tacJpsi->fVtx.fD3d;
            fd3Elb=pCand->fVtx.fD3dE;
            fd3El0=tacLambda0->fVtx.fD3dE;
            fd3Ejp=tacJpsi->fVtx.fD3dE;

            fdxylb=pCand->fVtx.fDxy;
            fdxyl0=tacLambda0->fVtx.fDxy;
            fdxyjp=tacJpsi->fVtx.fDxy;
            fdxyElb=pCand->fVtx.fDxyE;
            fdxyEl0=tacLambda0->fVtx.fDxyE;
            fdxyEjp=tacJpsi->fVtx.fDxyE;

	    // qualities of vertex fits
            fchi2lb=pCand->fVtx.fChi2;
            fchi2l0=tacLambda0->fVtx.fChi2;
            fchi2jp=tacJpsi->fVtx.fChi2;
            fndoflb=pCand->fVtx.fNdof;
            fndofl0=tacLambda0->fVtx.fNdof;
            fndofjp=tacJpsi->fVtx.fNdof;
	    fproblb=TMath::Prob(fchi2lb,fndoflb);
	    fprobl0=TMath::Prob(fchi2l0,fndofl0);
	    fprobjp=TMath::Prob(fchi2jp,fndofjp);

	    // ctau - written as dist/p*m as this is numerically more stable than using
	    // the TLorentzVector::Gamma() functions to write beta*gamma
	    fctlb = pCand->fVtx.fD3d / pCand->fPlab.Mag() * MLAMBDA_B;
	    fctl0 = tacLambda0->fVtx.fD3d / tacLambda0->fPlab.Mag() * MLAMBDA_0;

	    // beta vector components
	    fbtlbx=vecBetaL0.x();
	    fbtlby=vecBetaL0.y();
	    fbtlbz=vecBetaL0.z();
	    fbtl0x=vecBetaLb.x();
	    fbtl0y=vecBetaLb.y();
	    fbtl0z=vecBetaLb.z();

	    // other cut values
	    // pointing angles
            falphalb=alphalb;
            falphal0=alphal0;

	    // doca
            fmaxdocalb=pCand->fMaxDoca;
            fmaxdocal0=tacLambda0->fMaxDoca;
            fmaxdocajp=tacJpsi->fMaxDoca;

	    // dR and angles
            fdRprpi=tlvRecPr.DeltaR(tlvRecPi);
            fdRmumu=tlvRecMu1.DeltaR(tlvRecMu2);
            fdRl0jp=tlvLambda0.DeltaR(tlvLambdaB);
            fanglbl0=tlvLambda0.Angle(tlvLambdaB.Vect());

            // K_s hypothesis
            fKshypo=tlvKsHypoth.M();

	    // Do truth matching
	    if (fIsMC)
	    {
		const int giMu1(fpEvt->getGenIndexWithDeltaR(tlvRecMu1,tatRecMu1->fQ));
		const int giMu2(fpEvt->getGenIndexWithDeltaR(tlvRecMu2,tatRecMu2->fQ));
		const int giPi(fpEvt->getGenIndexWithDeltaR(tlvRecPi,tatRecPi->fQ));
		const int giPr(fpEvt->getGenIndexWithDeltaR(tlvRecPr,tatRecPr->fQ));
		const int nGenCands(fpEvt->nGenCands());
		if( giMu1 < nGenCands && giMu2 < nGenCands && giPi < nGenCands && giPr < nGenCands &&
		    giMu1 > -1 && giMu2 > -1 && giPi > -1 && giPr > -1 )
		{
		    TGenCand *gcMu1 = fpEvt->getGenCand(giMu1);
		    TGenCand *gcMu2 = fpEvt->getGenCand(giMu2);
		    TGenCand *gcPion = fpEvt->getGenCand(giPi);
		    TGenCand *gcProton = fpEvt->getGenCand(giPr);

		    fIsMCmatch = (TMath::Abs(gcMu1->fID) == 13
				&& TMath::Abs(gcMu2->fID) == 13
				&& TMath::Abs(gcPion->fID) == 211
				&& TMath::Abs(gcProton->fID) == 2212);

		    if (fVerbose > 5)
		    {
			cout << "Truth matching: Mu1: " << gcMu1->fID
			     << " Mu2: " << gcMu2->fID
			     << " Pion: " << gcPion->fID
			     << " Proton: " << gcProton->fID
			     << " isMatch: " << fIsMCmatch << " isSig: " << fIsSig << endl;
		    }
		}
            }

            // Fill the trees
            fTree->Fill();
        }
    }
    if (fIsMC && fIsSig) fGenTree->Fill();
    return;
}

void lambdaReader::doGenLevelStuff()
{
    // Do some generator level stuff
    const int nGenCands(fpEvt->nGenCands());
    for(int gcit=0; gcit!=nGenCands; gcit++)
    {
	TGenCand *gc5122 = fpEvt->getGenCand(gcit);
	if (511==abs(gc5122->fID) || 521==abs(gc5122->fID) || 513==abs(gc5122->fID) || 523==abs(gc5122->fID))
	{
	    findDaughters fd(gc5122,fpEvt);
	}
	if(5122==abs(gc5122->fID))
        {
	    // reset tree variables
	    fgmlb = fgmlbsw = fgml0 = fgml0sw = 9999;
	    fgptpr = fgptpi = fgptmu1 = fgptmu2 = 9999;
	    fgetamu1 = fgetamu2 = 9999;
	    fgppr = fgppi = 9999;
	    fgptl0 = fgpl0 = 9999;
	    fgdRprpi = fgdRmumu = fgdRl0lb = 9999;
	    fganprpi = fganmumu = fganl0jp = 9999;
	    fgnDaughters = fgnGrandDaughters = fgmapRef1Gen = -1;
	    // analyse the decay and extract some data to store in the tree
	    findDaughters fd(gc5122,fpEvt);
	    fIsSig = ("5122 -> 443 (13 13 ) 3122 (211 2212 ) "==fd.getDauGrdau()) ? 1 : 0;
	    fnDaughters = fgnDaughters = fd.getNDaughters();
	    fnGrandDaughters = fgnGrandDaughters = fd.getNDauGrdau();
	    fmapRef1Gen = fgmapRef1Gen = myDecayMap1Gen.getPos(fd.getDaughters());
	    fmapRef2Gen = fgmapRef2Gen = myDecayMap2Gen.getPos(fd.getDauGrdau());
	    // now extract some kinematic data
	    TLorentzVector tlvGenJp, tlvGenL0, tlvGenPr, tlvGenPi, tlvGenMu1, tlvGenMu2, tlvGenLambda0, tlvGenLambdaB;
	    if (fIsSig)
	    {
		for(int gcDauit=gc5122->fDau1; gcDauit<=gc5122->fDau2; gcDauit++)
		{
		    TGenCand* gcDau=fpEvt->getGenCand(gcDauit);
                    if(443==abs(gcDau->fID))
                    {
			unsigned int mucounter(0);
                        tlvGenJp=gcDau->fP;
                        for(int gcGrDauit=gcDau->fDau1; gcGrDauit<=gcDau->fDau2; gcGrDauit++)
			{
                            TGenCand* gcGrDau=fpEvt->getGenCand(gcGrDauit);
                            if(13==abs(gcGrDau->fID))
                            {
				mucounter++;
				if (mucounter==1) tlvGenMu1=gcGrDau->fP;
				if (mucounter==2) tlvGenMu2=gcGrDau->fP;
			    }
			}
			if(mucounter>=2)
			{
			    fgptmu1 = tlvGenMu1.Pt();
			    fgptmu2 = tlvGenMu2.Pt();
			    fgetamu1 = tlvGenMu1.Eta();
			    fgetamu2 = tlvGenMu2.Eta();
			    fgdRmumu = tlvGenMu1.DeltaR(tlvGenMu2);
			    fganmumu = tlvGenMu1.Angle(tlvGenMu2.Vect());
			}
		    }
		    if(3122==abs(gcDau->fID))
		    {
			tlvGenL0=gcDau->fP;
                        bool prFound(false), piFound(false);
                        for(int gcGrDauit=gcDau->fDau1; gcGrDauit<=gcDau->fDau2; gcGrDauit++)
                        {
                            TGenCand* gcGrDau=fpEvt->getGenCand(gcGrDauit);
                            if(!prFound&&2212==abs(gcGrDau->fID)) // in case we have more than one pr we get the first
                            {
				tlvGenPr=gcGrDau->fP;
				prFound=true;
                            }
                            if(!piFound&&211==abs(gcGrDau->fID)) // in case we have more than one pi we get the first
                            {
                                tlvGenPi=gcGrDau->fP;
                                piFound=true;
                            }
			}
			if(prFound&&piFound) // we have a pr and a pi
			{
			    fgptpr = tlvGenPr.Pt();
                            fgppr = tlvGenPr.P();
                            fgptpi = tlvGenPi.Pt();
                            fgppi = tlvGenPi.P();
                            fgdRprpi = tlvGenPr.DeltaR(tlvGenPi);
                            fganprpi = tlvGenPr.Angle(tlvGenPi.Vect());
                            // make a L0
                            tlvGenLambda0 = tlvGenPr + tlvGenPi;
                            fgml0 = tlvGenLambda0.M();
                            fgptl0 = tlvGenLambda0.Pt();
                            fgpl0 = tlvGenLambda0.P();
                            // make a L0 with swapped mass assumptions for pr and pi
                            TLorentzVector tlvGenPrAsPi, tlvGenPiAsPr;
                            tlvGenPrAsPi.SetVectM(tlvGenPr.Vect(),MPION);
                            tlvGenPiAsPr.SetVectM(tlvGenPi.Vect(),MPROTON);
                            const TLorentzVector tlvGenLambdaSwapped=tlvGenPrAsPi+tlvGenPiAsPr;
                            fgml0sw = tlvGenLambdaSwapped.M();
                            // make a Lb
                            tlvGenLambdaB = tlvGenLambda0+tlvGenJp;
                            fgmlb = tlvGenLambdaB.M();
                            const TLorentzVector tlvGenLambdaBSwapped = tlvGenLambdaSwapped+tlvGenJp;
                            fgmlbsw = tlvGenLambdaBSwapped.M();

                            //
                            //if (fgml0>2) fpEvt->dumpGenBlock();
			}
		    }
		    fgdRl0lb = tlvGenJp.DeltaR(tlvGenLambda0);
		    fganl0jp = tlvGenJp.Angle(tlvGenLambda0.Vect());
		    fganl0lb = tlvGenLambdaB.Angle(tlvGenLambda0.Vect());
		    const double anl0mu1 = tlvGenLambda0.Angle(tlvGenMu1.Vect());
		    const double anl0mu2 = tlvGenLambda0.Angle(tlvGenMu2.Vect());
		    fganl0mumin = anl0mu1<anl0mu2 ? anl0mu1 : anl0mu2;
		    fganl0muPt = tlvGenMu1.Pt() > tlvGenMu2.Pt() ? anl0mu1 : anl0mu2;
		}
	    }
	}
    }
}

// ----------------------------------------------------------------------
void lambdaReader::initVariables()
{

}

// ----------------------------------------------------------------------
void lambdaReader::fillHist()
{

    cout << "fillHist()" << endl;
    ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());

    if (0 != fpCand1)
    {
        ((TH1D*)fpHistFile->Get("h10"))->Fill(fpCand1->fPlab.Perp());
        ((TH1D*)fpHistFile->Get("h11"))->Fill(fpCand1->fMass);
        ((TH1D*)fpHistFile->Get("h12"))->Fill(fpCand1->fVtx.fChi2);

        fTree->Fill();
    }

}

// ----------------------------------------------------------------------
void lambdaReader::bookHist()
{
    cout << "==> lambdaReader: bookHist " << endl;

    TH1D *h;
    h = new TH1D("h1", "Ntrk", 500, 0., 1000.);
    h = new TH1D("h2", "NCand", 20, 0., 20.);
    h = new TH1D("h3", "cand ID", 1000100, -100., 1000000.);
    h = new TH1D("h10", "pT", 40, 0., 20.);
    h = new TH1D("h11", "mass", 50, 1.6, 2.1);
    h = new TH1D("h12", "chi2", 50, 0., 10.);

    return;
}


// ----------------------------------------------------------------------
void lambdaReader::bookReducedTree()
{
    cout << "==> lambdaReader: bookReducedTree" << endl;

    // create the events tree ======================================
    fTree = new TTree("events", "events");

    // run info
    fTree->Branch("run",     &fRun,     "run/I");
    fTree->Branch("event",   &fEvent,   "event/I");
    fTree->Branch("LS",      &fLS,      "LS/I");

    // candidates
    fTree->Branch("mlb",     &fmlb,     "mlb/D");
    fTree->Branch("ml0",     &fml0,     "ml0/D");
    fTree->Branch("mjp",     &fmjp,     "mjp/D");
    fTree->Branch("ptlb",    &fptlb,    "ptlb/D");
    fTree->Branch("ptl0",    &fptl0,    "ptl0/D");
    fTree->Branch("pjp",     &fpjp,     "pjp/D");
    fTree->Branch("plb",     &fplb,     "plb/D");
    fTree->Branch("pl0",     &fpl0,     "pl0/D");
    fTree->Branch("ptjp",    &fptjp,    "ptjp/D");
    fTree->Branch("etalb",   &fetalb,   "etalb/D");
    fTree->Branch("etal0",   &fetal0,   "etal0/D");
    fTree->Branch("etajp",   &fetajp,   "etajp/D");
    fTree->Branch("philb",   &fphilb,   "philb/D");
    fTree->Branch("phil0",   &fphil0,   "phil0/D");
    fTree->Branch("phijp",   &fphijp,   "phijp/D");

    // signal tracks
    fTree->Branch("rpt1m",   &frpt1m,   "rpt1m/D");
    fTree->Branch("rpt2m",   &frpt2m,   "rpt2m/D");
    fTree->Branch("rptpr",   &frptpr,   "rptpr/D");
    fTree->Branch("rptpi",   &frptpi,   "rptpi/D");
    fTree->Branch("reta1m",  &freta1m,  "reta1m/D");
    fTree->Branch("reta2m",  &freta2m,  "reta2m/D");
    fTree->Branch("retapr",  &fretapr,  "retapr/D");
    fTree->Branch("retapi",  &fretapi,  "retapi/D");
    fTree->Branch("rphi1m",  &frphi1m,  "rphi1m/D");
    fTree->Branch("rphi2m",  &frphi2m,  "rphi2m/D");
    fTree->Branch("rphipr",  &frphipr,  "rphipr/D");
    fTree->Branch("rphipi",  &frphipi,  "rphipi/D");
    fTree->Branch("rid1m",   &frid1m,   "rid1m/I");
    fTree->Branch("rid2m",   &frid2m,   "rid2m/I");
    fTree->Branch("rq1m",    &frq1m,    "rq1m/I");
    fTree->Branch("rq2m",    &frq2m,    "rq2m/I");
    fTree->Branch("rqpr",    &frqpr,    "rqpr/I");
    fTree->Branch("rqpi",    &frqpi,    "rqpi/I");
    fTree->Branch("chi21m",  &fchi21m,  "chi21m/D");
    fTree->Branch("chi22m",  &fchi22m,  "chi22m/D");
    fTree->Branch("chi2pr",  &fchi2pr,  "chi2pr/D");
    fTree->Branch("chi2pi",  &fchi2pi,  "chi2pi/D");
    fTree->Branch("prob1m",  &fprob1m,  "prob1m/D");
    fTree->Branch("prob2m",  &fprob2m,  "prob2m/D");
    fTree->Branch("probpr",  &fprobpr,  "probpr/D");
    fTree->Branch("probpi",  &fprobpi,  "probpi/D");
    fTree->Branch("ndof1m",  &fndof1m,  "ndof1m/I");
    fTree->Branch("ndof2m",  &fndof2m,  "ndof2m/I");
    fTree->Branch("ndofpr",  &fndofpr,  "ndofpr/I");
    fTree->Branch("ndofpi",  &fndofpi,  "ndofpi/I");
    fTree->Branch("qual1m",  &fqual1m,  "qual1m/I");
    fTree->Branch("qual2m",  &fqual2m,  "qual2m/I");
    fTree->Branch("qualpr",  &fqualpr,  "qualpr/I");
    fTree->Branch("qualpi",  &fqualpi,  "qualpi/I");

    fTree->Branch("Kshypo",  &fKshypo,  "Kshypo/D");
    fTree->Branch("alphalb", &falphalb, "alphalb/D");
    fTree->Branch("alphal0", &falphal0, "alphal0/D");
    fTree->Branch("maxdocalb",&fmaxdocalb,    "maxdocalb/D");
    fTree->Branch("maxdocal0",&fmaxdocal0,    "maxdocal0/D");
    fTree->Branch("maxdocajp",&fmaxdocajp,    "maxdocajp/D");
    fTree->Branch("d3lb",    &fd3lb,    "d3lb/D");
    fTree->Branch("d3l0",    &fd3l0,    "d3l0/D");
    fTree->Branch("d3jp",    &fd3jp,    "d3jp/D");
    fTree->Branch("d3Elb",   &fd3Elb,   "d3Elb/D");
    fTree->Branch("d3El0",   &fd3El0,   "d3El0/D");
    fTree->Branch("d3Ejp",   &fd3Ejp,   "d3Ejp/D");
    fTree->Branch("dxylb",   &fdxylb,   "dxylb/D");
    fTree->Branch("dxyl0",   &fdxyl0,   "dxyl0/D");
    fTree->Branch("dxyjp",   &fdxyjp,   "dxyjp/D");
    fTree->Branch("dxyElb",  &fdxyElb,  "dxyElb/D");
    fTree->Branch("dxyEl0",  &fdxyEl0,  "dxyEl0/D");
    fTree->Branch("dxyEjp",  &fdxyEjp,  "dxyEjp/D");

    fTree->Branch("ctlb",    &fctlb,    "ctlb/D");
    fTree->Branch("ctl0",    &fctl0,    "ctl0/D");

    fTree->Branch("btlbx",   &fbtlbx,   "btlbx/D");
    fTree->Branch("btlby",   &fbtlby,   "btlby/D");
    fTree->Branch("btlbz",   &fbtlbz,   "btlbz/D");
    fTree->Branch("btl0x",   &fbtl0x,   "btl0x/D");
    fTree->Branch("btl0y",   &fbtl0y,   "btl0y/D");
    fTree->Branch("btl0z",   &fbtl0z,   "btl0z/D");

    fTree->Branch("chi2lb",  &fchi2lb,  "chi2lb/D");
    fTree->Branch("chi2l0",  &fchi2l0,  "chi2l0/D");
    fTree->Branch("chi2jp",  &fchi2jp,  "chi2jp/D");
    fTree->Branch("ndoflb",  &fndoflb,  "ndoflb/D");
    fTree->Branch("ndofl0",  &fndofl0,  "ndofl0/D");
    fTree->Branch("ndofjp",  &fndofjp,  "ndofjp/D");
    fTree->Branch("problb",  &fproblb,  "problb/D");
    fTree->Branch("probl0",  &fprobl0,  "probl0/D");
    fTree->Branch("probjp",  &fprobjp,  "probjp/D");

    fTree->Branch("dRprpi",  &fdRprpi,  "dRprpi/D");
    fTree->Branch("dRmumu",  &fdRmumu,  "dRmumu/D");
    fTree->Branch("dRl0jp",  &fdRl0jp,  "dRl0jp/D");

    fTree->Branch("anglbl0",&fanglbl0,"anglbl0/D");
    // gen info to main tree
    fTree->Branch("nDau",    &fnDaughters,"nDau/I");
    fTree->Branch("nGDau",   &fnGrandDaughters,"nGDau/I");
    fTree->Branch("nRef1G",  &fmapRef1Gen,"nRef1G/I");
    fTree->Branch("nRef2G",  &fmapRef2Gen,"nRef2G/I");
    fTree->Branch("isSig",   &fIsSig,     "isSig/I");
    fTree->Branch("isMCmatch",&fIsMCmatch, "isMCmatch/I");

    // sigtrack info
    fTree->Branch("Spt1m",   &fSpt1m,   "Spt1m/D");
    fTree->Branch("Seta1m",  &fSeta1m,  "Seta1m/D");
    fTree->Branch("Sphi1m",  &fSphi1m,  "Sphi1m/D");
    fTree->Branch("Spt2m",   &fSpt2m,   "Spt2m/D");
    fTree->Branch("Seta2m",  &fSeta2m,  "Seta2m/D");
    fTree->Branch("Sphi2m",  &fSphi2m,  "Sphi2m/D");
    fTree->Branch("Sptpr",   &fSptpr,   "Sptpr/D");
    fTree->Branch("Setapr",  &fSetapr,  "Setapr/D");
    fTree->Branch("Sphipr",  &fSphipr,  "Sphipr/D");
    fTree->Branch("Sptpi",   &fSptpi,   "Sptpi/D");
    fTree->Branch("Setapi",  &fSetapi,  "Setapi/D");
    fTree->Branch("Sphipi",  &fSphipi,  "Sphipi/D");

    // Generator tree ==========================================
    fGenTree = new TTree("genevents", "genevents");
    // run info
    fGenTree->Branch("run",     &fRun,     "run/I");
    fGenTree->Branch("event",   &fEvent,   "event/I");
    fGenTree->Branch("LS",      &fLS,      "LS/I");

    fGenTree->Branch("nDau",    &fgnDaughters,"nDau/I");
    fGenTree->Branch("nGDau",   &fgnGrandDaughters,"nGDau/I");
    fGenTree->Branch("nRef1G",  &fgmapRef1Gen,"nRef1G/I");
    fGenTree->Branch("nRef2G",  &fgmapRef2Gen,"nRef2G/I");
    fGenTree->Branch("hasCand", &fghasCand,"hasCand/I");

    fGenTree->Branch("mlb",     &fgmlb,    "gmlb/D");
    fGenTree->Branch("mlbSwap", &fgmlbsw,  "gmlbsw/D");
    fGenTree->Branch("ml0",     &fgml0,    "gml0/D");
    fGenTree->Branch("ml0Swap", &fgml0sw,  "gml0sw/D");

    fGenTree->Branch("ptmu1",   &fgptmu1,  "ptmu1/D");
    fGenTree->Branch("ptmu2",   &fgptmu2,  "ptmu2/D");
    fGenTree->Branch("etamu1",  &fgetamu1, "etamu1/D");
    fGenTree->Branch("etamu2",  &fgetamu2, "etamu2/D");

    fGenTree->Branch("ptpr",    &fgptpr,   "ptpr/D");
    fGenTree->Branch("ptpi",    &fgptpi,   "ptpi/D");
    fGenTree->Branch("ppr",     &fgppr,    "ppr/D");
    fGenTree->Branch("ppi",     &fgppi,    "ppi/D");

    fGenTree->Branch("ptl0",    &fgptl0,   "ptl0/D");
    fGenTree->Branch("pl0",     &fgpl0,    "pl0/D");

    fGenTree->Branch("dRprpi",  &fgdRprpi, "dRprpi/D");
    fGenTree->Branch("dRmumu",  &fgdRmumu, "dRmumu/D");
    fGenTree->Branch("dRl0lb",  &fgdRl0lb, "dRl0lb/D");

    fGenTree->Branch("anprpi",  &fganprpi, "anprpi/D");
    fGenTree->Branch("anmumu",  &fganmumu, "anmumu/D");
    fGenTree->Branch("anl0jp",  &fganl0jp, "anl0jp/D");
    fGenTree->Branch("anl0lb",  &fganl0lb, "anl0lb/D");
    fGenTree->Branch("anl0mumin",&fganl0mumin, "anl0mumin/D");
    fGenTree->Branch("anl0muPt",&fganl0muPt,"anl0muPt/D");
}

// ----------------------------------------------------------------------
void lambdaReader::readCuts(TString filename, int dump)
{
    fCutFile = filename;
    if (dump) cout << "==> lambdaReader: Reading " << fCutFile.Data() << " for cut settings" << endl;
    ifstream infile(fCutFile.Data());
    //std::string strCutName;
    //float fltCutValue;
    //int ok(0);

    TString fn(fCutFile.Data());

    if (dump)
    {
        cout << "====================================" << endl;
        cout << "==> lambdaReader: Cut file  " << fCutFile.Data() << endl;
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
	    if("CUTLbCandidate" == key)
	    {
		setCut(CUTLbCandidate, value);
		//setCut(CUTLbCandidate, value, hcuts, 1, "Candidate");
		continue;
	    }
	    if("CUTMuId1" == key)
	    {
		setCut(CUTMuId1, value);
		continue;
	    }
	    if("CUTMuId2" == key)
	    {
		setCut(CUTMuId2, value);
		continue;
	    }
	    if("CUTMjpMin" == key)
	    {
		setCut(CUTMjpMin, value, hcuts, 1, "m^{min}(J/#psi)");
		continue;
	    }
	    if("CUTMjpMax" == key)
	    {
		setCut(CUTMjpMax, value, hcuts, 2, "m^{max}(J/#psi)");
		continue;
	    }
	    if("CUTAlphal0Max" == key)
	    {
		setCut(CUTAlphal0Max, value, hcuts, 3, "#alpha^{max}(#Lambda^{b})");
		continue;
	    }
	    if("CUTMl0Max" == key)
	    {
		setCut(CUTMl0Max, value, hcuts, 4, "m^{max}(#Lambda^{0})");
		continue;
	    }
	    if("CUTD3dl0Min" == key)
	    {
		setCut(CUTD3dl0Min, value, hcuts, 5, "d3^{min}(#Lambda^{0})");
		continue;
	    }
	    if("CUTPtl0Min" == key)
	    {
		setCut(CUTPtl0Min, value, hcuts, 6, "p_{#perp}^{min}(#Lambda^{0})");
		continue;
	    }
	    if("CUTReadDecayMaps" == key)
	    {
		setCut(CUTReadDecayMaps, value);
		continue;
	    }
	    if("CUTPrintDecayMaps" == key)
	    {
		setCut(CUTPrintDecayMaps, value);
		continue;
	    }
	    if("CUTDecayMap1" == key)
	    {
		setCut(CUTDecayMap1, value);
		continue;
	    }
	    if("CUTDecayMap2" == key)
	    {
		setCut(CUTDecayMap2, value);
		continue;
	    }

	    cout << "==> lambdaReader: ERROR in cutfile: Don't know what to do with " << key << "=" << value << endl;
	}
    }

    if (dump)  cout << "------------------------------------" << endl;

}
