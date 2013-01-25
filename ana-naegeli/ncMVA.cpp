/*
 *  ncMVA.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.02.12.
 *
 */

#include "ncMVA.h"
#include "ncAna.h"
#include "ncEvaluate.h"
#include "ncVarReader.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <TLeaf.h>
#include <TFile.h>
#include <TEventList.h>
#include <TMVA/Factory.h>

#define NBR_SPLITS 3

using namespace std;

/* FAVOURITE MVA CONFIGURATION STRINGS
 *	BARREL:
 *		!H:V:TestRate=10:VarTransform=Norm:NeuronType=radial:NCycles=667:HiddenLayers=N-1:LearningRate=0.497:DecayRate=0.0229
 *		!H:V:TestRate=10:VarTransform=Norm,Gauss,Deco:NeuronType=sigmoid:NCycles=1432:HiddenLayers=N-1:LearningRate=0.378:DecayRate=0.0034
 *	ENDCAP:
 *		!H:V:TestRate=10:VarTransform=Norm,Deco:NeuronType=sigmoid:NCycles=1018:HiddenLayers=N+7,N+4:LearningRate=0.065:DecayRate=0.0214
 */

ncMVA::ncMVA() :
	fMCPath("/Users/cn/CMSData/Reduced/production-mix-general.root"),
	fDataPath("/Users/cn/CMSData/Reduced/data-2011.root"),
	fVarPath("cuts/mvavars.def"),
	fSigWeight(0.000156),
	fBkgWeight(0.2),
	fTrainFilename(NULL)
{
	fMVAOpts[0] = std::string("!H:V:TestRate=10:VarTransform=Norm,Gauss,Deco:NeuronType=sigmoid:NCycles=1432:HiddenLayers=N-1:LearningRate=0.378:DecayRate=0.0034");
	fMVAOpts[1] = std::string("!H:V:TestRate=10:VarTransform=Norm,Deco:NeuronType=sigmoid:NCycles=1018:HiddenLayers=N+7,N+4:LearningRate=0.065:DecayRate=0.0214");
} // ncMVA()

ncMVA::~ncMVA()
{
	if (fTrainFilename) {
		unlink(fTrainFilename->c_str()); // delete the temporary file
		delete fTrainFilename;
	}
} // ~ncMVA()

std::set<ncCut> ncMVA::getMVAVariables()
{
	using std::cout; using std::endl;
	std::set<ncCut> result;
	ncVarReader reader;
	
	try {
		reader.loadFile(fVarPath.c_str());
	} catch (std::string err) {
		cout << "==> ncMVA::getMVAVariables(): ERROR '" << err << "'" << endl;
		goto bail;
	}
	
	if (reader.getNbr() > 0)
		result = *reader.getVars(0);
	
bail:
	return result;
} // getMVAVariables()

void ncMVA::splitTree(TTree *tree, bool save)
{
	map<ncPair, ncPair> entries;
	map<ncPair, ncPair>::const_iterator it;
	ncPair key(0,0);
	Long64_t j;
	TLeaf *runLeaf = tree->FindLeaf("run");
	TLeaf *evtLeaf = tree->FindLeaf("event");
	TTree *theTrees[NBR_SPLITS];
	uint32_t e;
	
	cout << "splitting the tree()" << endl;
	
	for (j = 0; j < NBR_SPLITS; j++) {
		theTrees[j] = tree->CloneTree(0);
		theTrees[j]->SetName(Form("Split_%d",(int)j));
	}
	
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
			e = (uint32_t)TMath::Abs(evtLeaf->GetValue()) % NBR_SPLITS;
			theTrees[e]->Fill();
		}
	}
	
	// save the split trees to the (temp) file
	for (j = 0; j < NBR_SPLITS && save; j++)
		theTrees[j]->Write();
} // splitTree()

void ncMVA::prepTraining(unsigned channelIx, std::set<ncCut> *vars)
{
	std::set<ncCut>::const_iterator it;
	char *tmp = tempnam(".", "training-");
	TFile *tmpFile;
	TFile *sigFile = TFile::Open(fMCPath.c_str());
	TFile *bkgFile = TFile::Open(fDataPath.c_str());
	TCut cut;
	ncAna a; // default ncAna object to get default cuts
	TTree *sigTree;
	TTree *bkgTree;
	TCut preselCut = a.cutSigCand() && a.cutMuon() && a.cutTrigger(true, channelIx == 0) && a.cutChannel(channelIx) && a.cutPreselection();
	
	for (it = vars->begin(); it != vars->end(); ++it) {
		preselCut = preselCut && TCut(Form("%e < %s && %s < %e",it->getCut().first,it->getFormula(),it->getFormula(),it->getCut().second));
	}
	
	if (fTrainFilename) {
		unlink(fTrainFilename->c_str());
		delete fTrainFilename;
	}
	
	fTrainFilename = new string(Form("%s.root",tmp));
	tmpFile = new TFile(fTrainFilename->c_str(),"recreate");
	
	// change to the temp file
	tmpFile->cd();
	
	sigTree = (TTree*)sigFile->Get("T");
	bkgTree = (TTree*)bkgFile->Get("T");
	
	// extract the signal
	cout << "Truth matching signal..." << flush;
	cut = preselCut && a.cutSigTruth();
	sigTree = sigTree->CopyTree(cut.GetTitle());
	sigTree->Write("Sig");
	cout << "\tdone" << endl;
	
	// extract the sidebands
	cout << "Extracting sidebands..." << flush;
	cut = preselCut && !a.cutBlindRegion();
	bkgTree = bkgTree->CopyTree(cut.GetTitle());
	bkgTree->Write("Bkg");
	cout << "\tdone" << endl;
	
	// split in even and odd for background
	splitTree(bkgTree);
	
	free(tmp);
	delete tmpFile;
	delete sigFile;
	delete bkgFile;
} // prepTraining()

void ncMVA::runTraining(unsigned split, unsigned channelIx, bool prep)
{
	TFile *trainFile;
	TFile *tmvaFile;
	TString factOptions("Transformations=I;G;D;G,D;D,G");
	TString prepOptions("");
	TString methodTitle(Form("MLP_s%u_c%u",split,channelIx));
	set<ncCut> variables = getMVAVariables();
	set<ncCut>::const_iterator it;
	TMVA::Factory *factory;
	const char *c;
	TTree *sigTree;
	TTree *bkgTree[NBR_SPLITS];
	int j;
	
	cout << "==> ncMVA: Running training for split=" << split << " and channel=" << channelIx << endl;
	if(prep) prepTraining(channelIx,&variables);
	
	// reload the trees from the train file
	trainFile = new TFile(fTrainFilename->c_str());
	tmvaFile = new TFile(Form("TMVA_s%u_c%u.root",split,channelIx),"recreate");
	
	sigTree = (TTree*)trainFile->Get("Sig");
	for (j = 0; j < NBR_SPLITS; j++)
		bkgTree[j] = (TTree*)trainFile->Get(Form("Split_%d",j));
	
	factory = new TMVA::Factory("ncMVA",tmvaFile,factOptions);
	factory->AddSignalTree(sigTree,fSigWeight);
	factory->AddBackgroundTree(bkgTree[(split+0) % NBR_SPLITS],fBkgWeight);
	factory->AddBackgroundTree(bkgTree[(split+1) % NBR_SPLITS],fBkgWeight);
	
	for (it = variables.begin(); it != variables.end(); ++it) {
		
		if (strcmp(it->getName(), it->getFormula()) == 0)	c = it->getName();
		else												c = Form("%s := %s", it->getName(), it->getFormula());
		
		factory->AddVariable(c, 'F');
	}
	
	// prepare training and test
	prepOptions = Form("nTrain_Background=%lld:SplitMode=Block",bkgTree[split%NBR_SPLITS]->GetEntries());
	factory->PrepareTrainingAndTestTree("",prepOptions);
	
	// book neural network
	cout << "Booking method" << endl;
	factory->BookMethod(TMVA::Types::kMLP, methodTitle, TString(fMVAOpts[channelIx].c_str()));
	
	// Training & Testing
	cout << "Training all methods" << endl;
	factory->TrainAllMethods();
	cout << "Testing all methods" << endl;
	factory->TestAllMethods();
	cout << "Evaluating all methods" << endl;
	factory->EvaluateAllMethods();
	
	delete factory;
	delete trainFile;
	delete tmvaFile;
} // runTraining()

void ncMVA::runAllTrainings()
{
	runTraining(0, 0, true);
	runTraining(1, 0, false);
	runTraining(2, 0, false);
	runTraining(0, 1, true); // reload the trees for other channel
	runTraining(1, 1, false);
	runTraining(2, 1, false);
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
	ncAna a;
	
	tree->Draw(">>elist",a.cutPreselection(0));
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
