/*
 *  ncMVA.cc
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 10.02.12.
 *
 */

#include "ncMVA.hh"
#include "ncEvaluate.hh"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <TFile.h>
#include <TEventList.h>
#include <TMVA/Factory.h>

using namespace std;

TTree* ncEvalAll(TTree *tree, bool verbose)
{
	Long64_t j;
	Float_t mlp[2];
	ncEvaluate v0("weights/ncMVA_MLP_0.weights.xml",tree,0);
	ncEvaluate v1("weights/ncMVA_MLP_1.weights.xml",tree,1);
	TEventList elist("elist");
	TTree *outTree = NULL;
	TBranch *branch;
	set<string>::const_iterator it;
	int32_t completed = 0, old = 0;
	int k,written = 0;
	
	tree->Draw(">>elist", ncEvaluate::getPreselCut());
	elist.Sort();
	
	outTree = tree->CloneTree(0);
	for (k = 0; k < 2; k++) {
		if( (branch = outTree->FindBranch(Form("mlp_%d",k))) != NULL)
		   branch->SetAddress(&mlp[k]);
		else
		   branch = outTree->Branch(Form("mlp_%d",k),&mlp[k],Form("mlp_%d/F",k));
	}
	
	if(verbose) {
		written = printf("Evaluating: %d %%", completed);
		fflush(stdout);
	}
	for (j = 0; j < tree->GetEntries(); j++) {
		
		tree->GetEntry(j); // load the entry
		
		if (elist.GetIndex(j) < 0) {
			mlp[0] = mlp[1] = -99.0;		// failed the preselection cut
		} else {
			mlp[0] = v0.eval(j);	// evaluate
			mlp[1] = v1.eval(j);	// evaluate
		}
		
		outTree->Fill();
		
		// update status
		completed = (int32_t)(100 * (j+1) / tree->GetEntries());
		if (completed > old) {
			old = completed;
			// erase old progress
			for (k = 0; k < written; k++)	fputc('\b', stdout);
			for (k = 0; k < written; k++)	fputc(' ', stdout);
			while (written-- > 0)			fputc('\b', stdout);
			
			// show progress
			if (verbose) {
				written = printf("Evaluating: %d %%", completed);
				fflush(stdout);
			}
		}		
	}
	if (verbose)
		fputc('\n', stdout); // newline...
	
	return outTree;
} // ncEvalAll()

void ncEvalAll(const char *treeFileName, bool verbose)
{
	char *tname = tempnam(".", "eval");
	TFile tmpFile(tname,"create");
	TFile file(treeFileName);
	TTree *tree = (TTree*)file.Get("T");
	int err;
	
	tmpFile.cd();
	tree = ncEvalAll(tree);
	tree->Write(tree->GetName(),TObject::kOverwrite);
	
	file.Close();
	tmpFile.Close();
	
	// swap the files...
	err = rename(tname, treeFileName);
	if(err) fprintf(stderr, "Error renaming file, errno = %d (%s)\n",errno,strerror(errno));
	
	free(tname);
} // ncEvalAll()

void ncRunTraining(TTree *signalTree, double signalWeight, TTree *bkgTree, double bkgWeight, int channelIx, TCut channelCut, string *mOptions)
{
	TFile trainFile(Form("TMVA_%d.root",channelIx),"recreate");
	TString factOptions("Transformations=I");
	TString prepOptions("");
	TString methodTitle(Form("MLP_%d",channelIx));
	TCut muonCut("tight_mu1 && tight_mu2");
	TCut presel = ncEvaluate::getPreselCut() && muonCut && TCut("triggered_bs") && channelCut;
	map<string,string> *variables = ncEvaluate::getDefaultVariables();
	map<string,string>::const_iterator it;
	const char *c;
	
	cout << "===> Running training for channel " << channelIx << endl;
	
	// FIXME: Optimize the booking options
	// TString allOptions("!H:V:VarTransform=Norm:NCycles=500:HiddenLayers=N,N-1:NeuronType=sigmoid:NeuronInputType=sum:LearningRate=0.02:DecayRate=0.01:TestRate=10:Tau=3:BPMode=sequential");
	TString methodOptions("!H:V:VarTransform=Norm:NCycles=500:HiddenLayers=N+2,N+2:NeuronType=sigmoid:LearningRate=0.08:DecayRate=0.03");
	// Also good
	//	"!H:V:VarTransform=Norm:NCycles=800:HiddenLayers=N+2,N+2:NeuronType=sigmoid:LearningRate=0.01:DecayRate=0.005"
	
	// create factory
	TMVA::Factory *factory = new TMVA::Factory("ncMVA",&trainFile,factOptions);
	
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
	if (mOptions) methodOptions = *mOptions;
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
	TCut candCut("candidate == 301313 && d3e > 0 && ipe > 0");
	TCut truthCut("true_decay == 1");
	TCut histoCut("4.9 < mass && mass < 5.9");
	TCut signalWindow("5.2 < mass && mass < 5.45");
	TCut barrelCut("TMath::Abs(eta_mu1) < 1.4 && TMath::Abs(eta_mu2) < 1.4");
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
	double signalWeight = 0.000156;
	
	ncRunTraining(signalTree, signalWeight, bkgTree, bkgWeight, 0, barrelCut, NULL);
	ncRunTraining(signalTree, signalWeight, bkgTree, bkgWeight, 1, !barrelCut, NULL);
	
	unlink(tmpfilename);
} // runTraining()
