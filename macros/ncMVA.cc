/*
 *  ncMVA.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.02.12.
 *
 */

#include "ncMVA.hh"
#include "ncEvaluate.hh"

#include <TFile.h>
#include <TEventList.h>
#include <TMVA/Factory.h>

using namespace std;

TTree* ncEvalAll(TTree *tree)
{
	Long64_t j;
	Float_t mlp;
	ncEvaluate v("weights/nnFact_MLP.weights.xml",tree);
	TEventList elist("elist");
	TTree *outTree = NULL;
	TBranch *branch;
	set<string>::const_iterator it;
	
	tree->Draw(">>elist", ncEvaluate::getPreselCut());
	elist.Sort();
	
	outTree = tree->CloneTree(0);
	if( (branch = outTree->FindBranch("mlp")) != NULL)
		branch->SetAddress(&mlp);
	else
		branch = outTree->Branch("mlp",&mlp,"mlp/F");
	
	for (j = 0; j < tree->GetEntries(); j++) {
		tree->GetEntry(j); // load the entry
		
		if (elist.GetIndex(j) < 0)
			mlp = -99.0;		// failed the preselection cut
		else
			mlp = v.eval(j);	// evaluate
		
		outTree->Fill();
	}
	
	return outTree;
} // ncEvalAll()

void ncEvalAll(const char *treeFileName)
{
	TFile file(treeFileName,"update");
	TTree *tree = (TTree*)file.Get("T");
	
	tree = ncEvalAll(tree);
	tree->Write(tree->GetName(),TObject::kOverwrite);
	
	file.Close();
} // ncEvalAll()

void ncRunTraining(TTree *signalTree, double signalWeight, TTree *bkgTree, double bkgWeight)
{
	TFile trainFile("TMVA.root","recreate");
	TString factOptions("");
	TString prepOptions("");
	TString methodTitle("MLP");
	TCut presel = ncEvaluate::getPreselCut();
	map<string,string> *variables = ncEvaluate::getDefaultVariables();
	map<string,string>::const_iterator it;
	const char *c;
	
	// FIXME: Understand booking options
	TString methodOptions("!H:V:VarTransform=Norm:NeuronType=sigmoid:NCycles=500:HiddenLayers=N,N-1");
	
	// create factory
	TMVA::Factory *factory = new TMVA::Factory("nnFact",&trainFile,factOptions);
	
	// create training data
	factory->AddSignalTree(signalTree,signalWeight);
	factory->AddBackgroundTree(bkgTree,bkgWeight);
	
	for (it = variables->begin(); it != variables->end(); ++it) {
		
		if (it->first.compare(it->second) == 0)	c = it->first.c_str();
		else									c = Form("%s := %s", it->first.c_str(), it->second.c_str());
		
		factory->AddVariable(c, 'F');
	}
	
	// Prepare the test tree
	factory->PrepareTrainingAndTestTree(presel, prepOptions);
	
	// book neural network
	factory->BookMethod(TMVA::Types::kMLP, methodTitle, methodOptions);
	
	// Training & Testing
	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();
	
	trainFile.Close();
} // ncRunTraining()

void ncRunDefaultTraining(const char *mcFile, const char *dataFile)
{
	char tmpfilename[] = "tmp_training.root";
	TFile signalFile(mcFile);
	TFile bkgFile(dataFile);
	TFile tmpFile(tmpfilename,"recreate");
	TCut candCut("candidate == 301313");
	TCut truthCut("true_decay == 1");
	TCut histoCut("4.9 < mass && mass < 5.9");
	TCut signalWindow("5.2 < mass && mass < 5.45");
	TCut cut;
	
	TTree *signalTree = (TTree*)signalFile.Get("T");
	TTree *bkgTree = (TTree*)bkgFile.Get("T");
	
	tmpFile.cd();
	cout << "Truth matching signal..." << flush;
	cut = candCut && truthCut;
	signalTree = signalTree->CopyTree(cut.GetTitle());
	cout << "\tdone" << endl;
	cout << "Extracting sidebands..." << flush;
	cut = candCut && histoCut && !signalWindow;
	bkgTree = bkgTree->CopyTree(cut.GetTitle());
	cout << "\tdone" << endl;
	
	// FIXME: peaking background?
	double bkgWeight = 0.200000;
	double signalWeight = 0.0016215;
	
	ncRunTraining(signalTree, signalWeight, bkgTree, bkgWeight);
	
	unlink(tmpfilename);
} // runTraining()
