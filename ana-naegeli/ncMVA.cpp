/*
 *  ncMVA.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.02.12.
 *
 */

#include "ncMVA.h"
#include "ncAna.h"
#include "ncEvaluate.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <TLeaf.h>
#include <TFile.h>
#include <TEventList.h>
#include <TMVA/Factory.h>

using namespace std;

ncMVA::ncMVA() :
	fMCPath("/Users/cn/CMSData/Reduced/production-mix-general.root"),
	fDataPath("/Users/cn/CMSData/Reduced/data-2011.root"),
	fSigWeight(0.000156),
	fBkgWeight(0.2),
	fTrainFilename(NULL)
{
	fMVAOpts[0] = std::string("!H:V:TestRate=10:VarTransform=Norm:NeuronType=radial:NCycles=1291:HiddenLayers=N-5:LearningRate=0.318:DecayRate=0.022");
	fMVAOpts[1] = std::string("!H:V:TestRate=10:VarTransform=Norm:NeuronType=radial:NCycles=1107:HiddenLayers=N+7,N-1:LearningRate=0.086:DecayRate=0.0335");
} // ncMVA()

ncMVA::~ncMVA()
{
	if (fTrainFilename) {
		unlink(fTrainFilename->c_str()); // delete the temporary file
		delete fTrainFilename;
	}
} // ~ncMVA()

map<string,string>* ncMVA::getMVAVariables()
{
	static map<string,string> *mvaVars = NULL;
	
	if (!mvaVars) {
		mvaVars = new map<string,string>;
		
		mvaVars->insert(pair<string,string>("pt","pt"));
		mvaVars->insert(pair<string,string>("pt_mu1","pt_mu1"));
		mvaVars->insert(pair<string,string>("pt_mu2","pt_mu2"));
		mvaVars->insert(pair<string,string>("eta","eta"));
		
		mvaVars->insert(pair<string,string>("sig3d","d3 / d3e"));
		mvaVars->insert(pair<string,string>("d3","d3"));
		mvaVars->insert(pair<string,string>("alpha","alpha"));
		mvaVars->insert(pair<string,string>("normChi2","chi2/Ndof"));
		mvaVars->insert(pair<string,string>("ip","ip"));
		mvaVars->insert(pair<string,string>("sigip","ip / ipe"));
		
		mvaVars->insert(pair<string,string>("iso_mor12","iso_mor12"));
		mvaVars->insert(pair<string,string>("doca0","doca0"));
		mvaVars->insert(pair<string,string>("ntrk","ntrk"));
	}
	
	return mvaVars;
} // getMVAVariables()

void ncMVA::splitTree(TTree *tree, bool save)
{
	map<ncPair, ncPair> entries;
	map<ncPair, ncPair>::const_iterator it;
	ncPair key(0,0);
	Long64_t j;
	TLeaf *runLeaf = tree->FindLeaf("run");
	TLeaf *evtLeaf = tree->FindLeaf("event");
	unsigned e = 0; // even / odd marker
	TTree *theTrees[2];
	
	cout << "splitting the tree()" << endl;
	
	theTrees[0] = tree->CloneTree(0);
	theTrees[0]->SetName("Even");	
	theTrees[1] = tree->CloneTree(0);
	theTrees[1]->SetName("Odd");
	
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("run",1);
	tree->SetBranchStatus("event",1);
	
	for (j = 0; j < tree->GetEntries(); j++) {
		
		tree->GetEntry(j); // load the entry
		key.first = (Long64_t)runLeaf->GetValue();
		key.second = (Long64_t)evtLeaf->GetValue();
		if (entries.count(key) > 0)
			entries[key].second++;
		else
			entries[key] = ncPair(j,1);
	}
		
	// enable all branches again
	tree->SetBranchStatus("*",1);
	
	// iterate through all the events and sort them
	for (it = entries.begin(); it != entries.end(); ++it) {
		
		for (j = 0; j < it->second.second; j++) {
			
			// load and save the corresponding entry
			tree->GetEntry(it->second.first + j);
			theTrees[e++ % 2]->Fill();
		}
	}
	
	// save the split trees to the (temp) file
	if (save) {
		theTrees[0]->Write();
		theTrees[1]->Write();
	}
} // splitTree()

void ncMVA::prepTraining(unsigned channelIx)
{
	char *tmp = tempnam(".", "training-");
	TFile *tmpFile;
	TFile sigFile(fMCPath.c_str());
	TFile bkgFile(fDataPath.c_str());
	TCut preselCut = ncAna::cutMVAPresel() && ncAna::cutMuon() && ncAna::cutTrigger(true) && ncAna::cutChannel(channelIx);
	TCut cut;
	ncAna a; // default ncAna object
	TTree *sigTree;
	TTree *bkgTree;
	
	if (fTrainFilename) {
		unlink(fTrainFilename->c_str());
		delete fTrainFilename;
	}
	
	fTrainFilename = new string(Form("%s.root",tmp));
	tmpFile = new TFile(fTrainFilename->c_str(),"recreate");
	
	// change to the temp file
	tmpFile->cd();
	
	sigTree = (TTree*)sigFile.Get("T");
	bkgTree = (TTree*)bkgFile.Get("T");
	
	// extract the signal
	cout << "Truth matching signal..." << flush;
	cut = a.cutSigCand() && a.cutSanity() && a.cutSigTruth() && preselCut;
	sigTree = sigTree->CopyTree(cut.GetTitle());
	sigTree->Write("Sig");
	cout << "\tdone" << endl;
	
	// extract the sidebands
	cout << "Extracting sidebands..." << flush;
	cut = a.cutSigCand() && a.cutSanity() && a.cutHisto() && !a.cutBlindRegion() && preselCut;
	bkgTree = bkgTree->CopyTree(cut.GetTitle());
	bkgTree->Write("Bkg");
	cout << "\tdone" << endl;
	
	// split in even and odd for background
	splitTree(bkgTree);
	
	free(tmp);
	delete tmpFile;
} // prepTraining()

void ncMVA::runTraining(unsigned split, unsigned channelIx, bool prep)
{
	TFile *trainFile;
	TFile *tmvaFile;
	TString factOptions("Transformations=I");
	TString prepOptions("");
	TString methodTitle(Form("MLP_s%u_c%u",split,channelIx));
	map<string,string> *variables = getMVAVariables();
	map<string,string>::const_iterator it;
	TMVA::Factory *factory;
	const char *c;
	TTree *sigTree;
	TTree *bkgTree[2];
	
	cout << "==> ncMVA: Running training for split=" << split << " and channel=" << channelIx << endl;
	if(prep) prepTraining(channelIx);
	
	// reload the trees from the train file
	trainFile = new TFile(fTrainFilename->c_str());
	tmvaFile = new TFile(Form("TMVA_s%u_c%u.root",split,channelIx),"recreate");
	
	sigTree = (TTree*)trainFile->Get("Sig");
	bkgTree[0] = (TTree*)trainFile->Get("Even");
	bkgTree[1] = (TTree*)trainFile->Get("Odd");
	
	factory = new TMVA::Factory("ncMVA",tmvaFile,factOptions);
	factory->AddSignalTree(sigTree,fSigWeight);
	factory->AddBackgroundTree(bkgTree[(split+0) % 2],fBkgWeight);
	factory->AddBackgroundTree(bkgTree[(split+1) % 2],fBkgWeight);
	
	for (it = variables->begin(); it != variables->end(); ++it) {
		
		if (it->first.compare(it->second) == 0)	c = it->first.c_str();
		else									c = Form("%s := %s", it->first.c_str(), it->second.c_str());
		
		factory->AddVariable(c, 'F');
	}
	
	// prepare training and test
	prepOptions = Form("nTrain_Background=%lld:SplitMode=Block",bkgTree[split%2]->GetEntries());
	factory->PrepareTrainingAndTestTree("",prepOptions);
	
	// book neural network
	cout << "Booking method" << endl;
	factory->BookMethod(TMVA::Types::kMLP, methodTitle, fMVAOpts[channelIx]);
	
	// Training & Testing
	cout << "Training all methods" << endl;
	factory->TrainAllMethods();
	cout << "Testing all methods" << endl;
	factory->TestAllMethods();
	cout << "Evaluating all methods" << endl;
	factory->EvaluateAllMethods();
	
	delete trainFile;
	delete tmvaFile;
} // runTraining()

void ncMVA::runAllTrainings()
{
	runTraining(0, 0, true);
	runTraining(1, 0, false);
	runTraining(0, 1, true); // reload the trees for other channel
	runTraining(1, 1, false);
} // runTraining()

void ncMVA::evalFile(const char *filename, bool progressReport)
{
	set<pair<unsigned,unsigned> > s;
	
	s.insert(make_pair(0,0));
	s.insert(make_pair(0,1));
	s.insert(make_pair(1,0));
	s.insert(make_pair(1,1));
	
	evalFile(filename, &s, progressReport);
} // evalFile()

void ncMVA::evalFile(const char *filename, unsigned split, unsigned channelIx, bool progressReport)
{
	set<pair<unsigned,unsigned> > s;
	
	s.insert(make_pair(split,channelIx));
	
	evalFile(filename,&s,progressReport);
} // evalFile()

void ncMVA::evalFile(const char *filename, set<pair<unsigned,unsigned> >*mvas, bool progressReport)
{
	char *tname = tempnam(".", "eval");
	TFile tmpFile(tname,"create");
	TFile file(filename);
	TTree *tree = (TTree*)file.Get("T");
	int err;
	
	tmpFile.cd();
	tree = evalTree(tree, mvas, progressReport);
	tree->Write(tree->GetName(),TObject::kOverwrite);
	
	file.Close();
	tmpFile.Close();
	
	// swap the files...
	err = rename(tname, filename);
	if (err) fprintf(stderr, "Error renaming file, errno = %d (%s)\n",errno,strerror(errno));
	
	free(tname);
} // evalFile()

TTree *ncMVA::evalTree(TTree *tree, set<pair<unsigned,unsigned> > *mvas, bool progressReport)
{
	TEventList elist("elist");
	TTree *outTree = NULL;
	map<pair<unsigned,unsigned>, pair<Float_t,ncEvaluate*> > mlps;
	map<pair<unsigned,unsigned>, pair<Float_t,ncEvaluate*> >::iterator mit;
	set<pair<unsigned,unsigned> >::const_iterator sit;
	TBranch *branch;
	int32_t completed = 0, old = 0;
	int written = 0;
	Long64_t j,k;
	
	tree->Draw(">>elist",ncAna::cutMVAPresel() && ncAna::cutSanity());
	elist.Sort();
	
	outTree = tree->CloneTree(0);
	
	for (sit = mvas->begin(); sit != mvas->end(); ++sit) {
		string methodTitle(Form("MLP_s%u_c%u",sit->first,sit->second));
		mlps.insert(make_pair(*sit, make_pair(0.0, new ncEvaluate(Form("weights/ncMVA_%s",methodTitle.c_str()),tree,methodTitle.c_str()))));
		
		if ( (branch = outTree->FindBranch(Form("mlp_s%u_c%u",sit->first,sit->second))) != NULL )
			branch->SetAddress( &mlps[*sit] );
		else
			branch = outTree->Branch(Form("mlp_s%u_c%u",sit->first,sit->second),&mlps[*sit],Form("mlp_s%u_c%u/F",sit->first,sit->second));

	}
	
	if (progressReport) {
		written = printf("Evaluating: %d %%", completed);
		fflush(stdout);
	}
	
	for (j = 0; j < tree->GetEntries(); j++) {
		
		tree->GetEntry(j); // load the entry
		
		for (mit = mlps.begin(); mit != mlps.end(); ++mit) {
			
			if (tree->GetEntry(j) < 0)
				mit->second.first = -99.0; // failed the preselection cut
			else
				mit->second.first = mit->second.second->eval(j);
		}
		
		outTree->Fill();
		
		// update status
		completed = (int32_t)(100*(j+1) / tree->GetEntries());
		if (progressReport && completed > old) {
			old = completed;
			// erase old progress
			for (k = 0; k < written; k++)	fputc('\b', stdout);
			for (k = 0; k < written; k++)	fputc(' ', stdout);
			while (written-- > 0)			fputc('\b',stdout);
			
			written = printf("Evaluating: %d %%", completed);
			fflush(stdout);
		}
	}
	
	if (progressReport) {
		fputc('\n', stdout); // newline
	}
	
	// delete all the ncEvaluates
	for (mit = mlps.begin(); mit != mlps.end(); ++mit)
		delete mit->second.second;
	
	return outTree;
} // evalTree()
