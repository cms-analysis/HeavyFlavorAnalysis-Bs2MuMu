#include "lambdaReader.hh"
#include <cstdlib>
#include <cmath>

#define require_true(COND,LABEL) if( !(COND) ) goto LABEL
#include "../interface/HFMasses.hh"

// test a <= x < b
static inline int in_interval(int x, int a, int b)
{
    return (a <= x) && (x < b);
} // in_interval()

lambdaReader::lambdaReader(TChain *tree, TString evtClassName) : treeReader01(tree, evtClassName),reduced_tree(NULL)
{
    fMomentumPtr = &fMomentum;
} // lambdaReader()

lambdaReader::~lambdaReader()
{ } // ~lambdaReader()

void lambdaReader::eventProcessing()
{
    int j;
    TAnaCand *pCand;
    TAnaTrack *mu1Trk, *mu2Trk;

    // Fill a reduced tree
    int nCands = fpEvt->nCands();
    runNr = fpEvt->fRunNumber;
    goodrun = CheckGoodRun(runNr);

    for (j=0; j<nCands; j++) {
        pCand = fpEvt->getCand(j);

        // Save in the tree
        fCandidate = pCand->fType;
        fMass = pCand->fMass;
	// write J/Psi candidate
        if(fCandidate==10443 && fabs(fMass-MJPSI)<0.3 ) {
	    // initialise the variables to unphysical values
	    jpsi_m= jpsi_pt= jpsi_dxy= jpsi_dxyE= jpsi_d3d= jpsi_d3dE = -9.;
	    jpsi_chi2vtx= jpsi_deltaR = -9.;
	    jpsi_dca= jpsi_mu1_pt= jpsi_mu2_pt= jpsi_mu1_eta= jpsi_mu2_eta = -9.;
	    fdistPV= fdistDp= jpsi_mu1_dR= jpsi_mu2_dR= fdeltaRMu = -9.;

            nPV=0;

            fMomentum = pCand->fPlab;

            jpsi_pt   = fMomentum.Perp();
            jpsi_m    = pCand->fMass;

            jpsi_dxy  = pCand->fVtx.fDxy;
            jpsi_dxyE = pCand->fVtx.fDxyE;
            jpsi_d3d  = pCand->fVtx.fD3d;
            jpsi_d3dE = pCand->fVtx.fD3dE;

            if((pCand->fSig2-pCand->fSig1) != 1) {
                std::cout <<  "Error detected for 10443: nTracks=" << pCand->fSig2-pCand->fSig1 << std::endl;
                continue;
            }
            mu1Trk = fpEvt->getSigTrack(pCand->fSig1);
            mu2Trk = fpEvt->getSigTrack(pCand->fSig2);

            jpsi_mu1_pt = mu1Trk->fPlab.Perp();
            jpsi_mu1_eta= fabs(mu1Trk->fPlab.Eta());
            jpsi_mu1_dR = fMomentum.DeltaR(mu1Trk->fPlab);

            jpsi_mu2_pt = mu2Trk->fPlab.Perp();
            jpsi_mu2_eta= fabs(mu2Trk->fPlab.Eta());
            jpsi_mu2_dR = fMomentum.DeltaR(mu2Trk->fPlab);

            //dcamin = pCand->fVar1;
            //dcamax = pCand->fVar2;

            jpsi_chi2vtx = pCand->fVtx.fChi2;
            nPV=fpEvt->nPV();
            /* if(nPV>0) {
                TVector3 pVtx=fpEvt->getPV(0)->fPoint;
                TVector3 svtx = pCand->fVtx.fPoint-pVtx;
                dalpha = fMomentum.Angle(svtx);
                //if (dalpha>PI/2) dalpha-=PI;
            } */

            jpsi_truth = checkTruth(pCand, 443);
	    reduced_tree->Fill();
            //if(dcamax<DCAMAX && dpdxy/dpdxyE>SIGVTX) reduced_tree->Fill();
        } 
    }
} // dplusReader::eventProcessing()

void lambdaReader::bookHist()
{
    // create the tree
    reduced_tree = new TTree("Tjpsi","J/Psi Candidate");

    // and add the branches
    reduced_tree->Branch("candidate",   &fCandidate,   "candidate/I");
    reduced_tree->Branch("run",         &runNr ,       "runNr/I");
    reduced_tree->Branch("goodrun",     &goodrun,      "goodrun/I");
    reduced_tree->Branch("nPV",         &nPV ,         "nPV/I");

    reduced_tree->Branch("jpsi_m",      &jpsi_m,       "jpsi_m/D");
    reduced_tree->Branch("jpsi_truth",  &jpsi_truth,   "jpsi_truth/I");
    reduced_tree->Branch("jpsi_p","TVector3",&fMomentumPtr);
    reduced_tree->Branch("jpsi_pt",     &jpsi_pt,      "jpsi_pt/D");
    reduced_tree->Branch("jpsi_dxy",    &jpsi_dxy,     "jpsi_dxy/D");
    reduced_tree->Branch("jpsi_dxyE",   &jpsi_dxyE,    "jpsi_dxyE/D");
    reduced_tree->Branch("jpsi_d3d",    &jpsi_d3d,     "jpsi_d3d/D");
    reduced_tree->Branch("jpsi_d3dE",   &jpsi_d3dE,    "jpsi_d3dE/D");

    reduced_tree->Branch("jpsi_chi2vtx",&jpsi_chi2vtx, "jpsi_chi2vtx/D");

    reduced_tree->Branch("jpsi_dca",    &jpsi_dca,     "jpsi_dca/D");
    reduced_tree->Branch("jpsi_mu1_pt", &jpsi_mu1_pt,  "jpsi_mu1_pt/D");
    reduced_tree->Branch("jpsi_mu2_pt", &jpsi_mu2_pt,  "jpsi_mu2_pt/D");
    reduced_tree->Branch("jpsi_mu1_eta", &jpsi_mu1_eta,"jpsi_mu1_eta/D");
    reduced_tree->Branch("jpsi_mu2_eta", &jpsi_mu2_eta,"jpsi_etamu2/D");
    reduced_tree->Branch("jpsi_mu1_dR", &jpsi_mu1_dR,  "jpsi_mu1_dR/D");
    reduced_tree->Branch("jpsi_mu2_dR", &jpsi_mu2_dR,  "jpsi_mu2_dR/D");

} // lambdaReader::bookHist()



int lambdaReader::checkTruth(TAnaCand *cand, int truth_type)
{
    TAnaTrack *sgTrack;
    TAnaTrack *recTrack;
    TGenCand *truthParticle;
    TGenCand *trackParticle;
    int j,succes = 0;
    int nSigs,nRecs,nGens;


    nSigs = fpEvt->nSigTracks();
    nRecs = fpEvt->nRecTracks();
    nGens = fpEvt->nGenCands();

    // get the first track and get the originating particle of type 'truth_type' => truthParticle
    require_true(in_interval(cand->fSig1, 0, nSigs),bail);
    sgTrack = fpEvt->getSigTrack(cand->fSig1);

    require_true(in_interval(sgTrack->fIndex, 0, nRecs),bail);
    recTrack = fpEvt->getRecTrack(sgTrack->fIndex);

    require_true(in_interval(recTrack->fGenIndex, 0, nGens),bail);
    truthParticle = fpEvt->getGenCand(recTrack->fGenIndex);

    while (abs(truthParticle->fID) != truth_type && in_interval(truthParticle->fMom1, 0, nGens)) {
        truthParticle = fpEvt->getGenCand(truthParticle->fMom1);
        //	  if(abs(truthParticle->fID)==411) std::cout << "Truth type = "<<truthParticle->fID<<"   "<<truth_type<<std::endl;
    }

    // check if our original is a valid.
    require_true(abs(truthParticle->fID) == truth_type,bail);

    // reconstruct the other tracks
    for (j = cand->fSig1+1; j <=cand->fSig2; j++) {

        require_true(in_interval(j, 0, nSigs),bail);
        sgTrack = fpEvt->getSigTrack(j);

        require_true(in_interval(sgTrack->fIndex, 0, nRecs),bail);
        recTrack = fpEvt->getRecTrack(sgTrack->fIndex);

        require_true(in_interval(recTrack->fGenIndex, 0, nGens),bail);
        trackParticle = fpEvt->getGenCand(recTrack->fGenIndex);

        while (abs(trackParticle->fID)!=truth_type && in_interval(trackParticle->fMom1, 0, nGens))
            trackParticle = fpEvt->getGenCand(trackParticle->fMom1);

        // check the particle
        require_true(trackParticle->fNumber==truthParticle->fNumber,bail);
    }

    // still here? then every track originated from the right particle
    succes = 1;

bail:
    return succes;
} // lambdaReader::checkTruth()


int lambdaReader::CheckGoodRun(int run) {
    for (std::vector<int>::iterator it=goodRuns.begin(); it!=goodRuns.end(); it++) {
        if(run==*it) return 1;
    }
    return 0;
}


void lambdaReader::startAnalysis() {
    ifstream infile("goodruns.txt",ifstream::in);
    int tmp_run;
    while(!infile.eof()) {
        infile >> tmp_run;
        goodRuns.push_back(tmp_run);
        std::cout <<"Good tracker run: "<<tmp_run<<std::endl;
    }
    infile.close();
}


void lambdaReader::readCuts(TString filename, int dump) {
    dump=1;
    char  buffer[200];
    fCutFile = filename;
    if (dump) std::cout << "==> lambdaReader: Reading " << fCutFile.Data() << " for cut settings" << std::endl;
    sprintf(buffer, "%s", fCutFile.Data());
    ifstream is(buffer);
    char CutName[100];
    float CutValue;

    TString fn(fCutFile.Data());

    if (dump) {
        std::cout << "====================================" << std::endl;
        std::cout << "==> lambdaReader: Cut file  " << fCutFile.Data() << std::endl;
        std::cout << "------------------------------------" << std::endl;
    }

    while (is.getline(buffer, 200, '\n')) {
        if (buffer[0] == '#') {
            continue;
        }
        if (buffer[0] == '/') {
            continue;
        }
        sscanf(buffer, "%s %f", CutName, &CutValue);

        if (!strcmp(CutName, "TYPE")) {
            TYPE = int(CutValue);
            if (dump) std::cout << "TYPE:           " << TYPE << std::endl;
        }

        if (!strcmp(CutName, "SIGVTX")) {
            SIGVTX = CutValue;
            if (dump) std::cout << "SIGVTX:           " << SIGVTX << std::endl;
        }

        if (!strcmp(CutName, "DCAMAX")) {
            DCAMAX = CutValue;
            if (dump) std::cout << "DCAMAX:           " << DCAMAX << std::endl;
        }
    }
}
