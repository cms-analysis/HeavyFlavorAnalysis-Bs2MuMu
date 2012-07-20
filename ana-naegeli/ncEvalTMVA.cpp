#include <iostream>
#include <list>
#include <vector>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TRandom3.h>

using namespace std;

struct meas_t {
	Double_t median;
	std::pair<Double_t,Double_t> err; // first = loErr, second = highErr
};

static bool fillArrays(Int_t *nBkg, Double_t *bkg, Int_t *nSig, Double_t *sig, TTree *tree)
{
	TLeaf *idLeaf,*mlpLeaf = NULL;
	Long64_t j;
	bool result;
	
	idLeaf = tree->FindLeaf("classID");
	for (unsigned s = 0; !mlpLeaf && s < 2; s++)
		for (unsigned c = 0; !mlpLeaf && c < 2; c++)
			mlpLeaf = tree->FindLeaf(Form("MLP_s%u_c%u",s,c));
	
	result = idLeaf && mlpLeaf;
	for (j = 0; result && j < tree->GetEntries(); j++) {
		
		tree->GetEntry(j);
		switch ((int)idLeaf->GetValue()) {
			case 0: // signal
				*sig++ = mlpLeaf->GetValue();
				(*nSig)++;
				break;
			case 1: // bkg
				*bkg++ = mlpLeaf->GetValue();
				(*nBkg)++;
				break;
			default:
				result = false;
				break;
		}
	}
	
	return result;
} // fillArrays()

void evalTMVA(const char *path, Double_t S, Double_t B)
{
	TFile *file = TFile::Open(path);
	TTree *testTree,*trainTree;
	Double_t *sigTest = NULL,*sigTrain = NULL;
	Double_t *bkgTest = NULL,*bkgTrain = NULL;
	Double_t *mlpData;
	Double_t probSig = 0, probBkg = 0, sgn = 0;
	Double_t effSig,effBkg,cut;
	Double_t bestCut = -1.0;
	Int_t nSigTest = 0,nSigTrain = 0,nBkgTest = 0,nBkgTrain = 0;
	Int_t maximum,j;
	bool result;
	
	// get the trees
	if (!file) goto bail;
	testTree = (TTree*)file->Get("TestTree");
	trainTree = (TTree*)file->Get("TrainTree");
	
	if (!testTree || !trainTree) {
		cerr << "Not all necessary trees found in root file: '" << path << "'" << endl;
		goto bail;
	}
	
	maximum = std::max(testTree->GetEntries(),trainTree->GetEntries());
	
	// storage
	sigTest = new Double_t[maximum];
	sigTrain = new Double_t[maximum];
	bkgTest = new Double_t[maximum];
	bkgTrain = new Double_t[maximum];
	
	// fill the arrays
	result = fillArrays(&nBkgTest,bkgTest,&nSigTest,sigTest,testTree) && fillArrays(&nBkgTrain,bkgTrain,&nSigTrain,sigTrain,trainTree);
	if (!result) {
		cerr << "Couldn't fill the arrays in root file: '" << path << "'" << endl;
		goto bail;
	}
	
	// compute Kolmogorov-Smirnov
	sort(sigTest, sigTest + nSigTest);
	sort(sigTrain, sigTrain + nSigTrain);
	probSig = TMath::KolmogorovTest(nSigTest,sigTest,nSigTrain,sigTrain,"");
	
	sort(bkgTest, bkgTest + nBkgTest);
	sort(bkgTrain, bkgTrain + nBkgTrain);
	probBkg = TMath::KolmogorovTest(nBkgTest,bkgTest,nBkgTrain,bkgTrain,"");
	
	// compute significance (based on test sample)
	// recycle the sigTrain to store data;
	mlpData = sigTrain;
	merge(sigTest, sigTest + nSigTest, bkgTest, bkgTest + nBkgTest, mlpData);
	for (j = 1; j < nSigTest + nBkgTest ; j++) {
		
		cut = (mlpData[j] + mlpData[j-1])/2.;
		effSig = nSigTest - (Double_t)(lower_bound(sigTest, sigTest + nSigTest, cut) - sigTest);
		effBkg = nBkgTest - (Double_t)(lower_bound(bkgTest, bkgTest + nBkgTest, cut) - bkgTest);
		
		effSig /= (Double_t)nSigTest;
		effBkg /= (Double_t)nBkgTest;
		
		// normalize to real values
		effSig = effSig*S / TMath::Sqrt(effSig*S+effBkg*B);
		if (effSig > sgn) {
			sgn = effSig;
			bestCut = cut;
		}
	}
	
	cout << sgn << "\t" << bestCut << "\t" << probSig << "\t" << probBkg << "\t " << path << endl;
	
bail:
	delete [] sigTest; delete [] bkgTest;
	delete [] sigTrain; delete [] bkgTrain;
	delete file;
} // evalTMVA()

void evalTMVA2(const char *path, Double_t S, Double_t B)
{
	TFile *file = TFile::Open(path);
	TTree *testTree,*trainTree;
	Long64_t maximum;
	Double_t *sigTest = NULL,*sigTrain = NULL;
	Double_t *bkgTest = NULL,*bkgTrain = NULL;
	Double_t *mlpData;
	Double_t probSig,probBkg;
	Double_t cut,effSig,effBkg;
	Double_t obsSig,obsBkg;
	Double_t med,bestCut = -1;
	Int_t nSigTest = 0,nSigTrain = 0;
	Int_t nBkgTest = 0,nBkgTrain = 0;
	Int_t j,k;
	TRandom3 rand;
	vector<Double_t> v(100);
	meas_t bestMed;
	
	bestMed.median = -1;
	bestMed.err.first = bestMed.err.second = 0;
	
	if (!file) {
		cerr << "ERROR: Unable to open file '" << path << "'" << endl;
		goto bail;
	}
	
	testTree = (TTree*)file->Get("TestTree");
	trainTree = (TTree*)file->Get("TrainTree");
	
	if (!testTree || !trainTree) {
		cerr << "ERROR: Not all necessary trees found in root file '" << path << "'" << endl;
		goto bail;
	}
	
	maximum = std::max(testTree->GetEntries(),trainTree->GetEntries());
	sigTest = new Double_t[maximum];
	sigTrain = new Double_t[maximum];
	bkgTest = new Double_t[maximum];
	bkgTrain = new Double_t[maximum];
	
	// fill the arrays
	if(!(fillArrays(&nBkgTest,bkgTest,&nSigTest,sigTest,testTree) && fillArrays(&nBkgTrain,bkgTrain,&nSigTrain,sigTrain,trainTree))) {
		cerr << "ERROR: Could not fill the arrays from the root file." << endl;
		goto bail;
	}

	// Compute Kolmogorov-Smirnov
	sort(sigTest, sigTest + nSigTest);
	sort(sigTrain, sigTrain + nSigTrain);
	probSig = TMath::KolmogorovTest(nSigTest,sigTest,nSigTrain,sigTrain,"");
	
	sort(bkgTest, bkgTest + nBkgTest);
	sort(bkgTrain, bkgTrain + nBkgTrain);
	probBkg = TMath::KolmogorovTest(nBkgTest,bkgTest,nBkgTrain,bkgTrain,"");
	
	// compute the significance based on test sample plus error
	// recycle sigTrain
	mlpData = sigTrain;
	merge(sigTest, sigTest + nSigTest, bkgTest, bkgTest + nBkgTest, mlpData);
	for (j = 1; j < nSigTest + nBkgTest; j++) {
		
		cut = (mlpData[j] + mlpData[j-1])/2.;
		effSig = nSigTest - (Double_t)(lower_bound(sigTest, sigTest + nSigTest, cut) - sigTest);
		effBkg = nBkgTest - (Double_t)(lower_bound(bkgTest, bkgTest + nBkgTest, cut) - bkgTest);
		
		effSig /= (Double_t)nSigTest;
		effBkg /= (Double_t)nBkgTest;
		
		// normalize to observed values to get expected
		effSig *= S;
		effBkg *= B;
		
		for (k = 0; k < 100; k++) {
			obsSig = rand.Poisson(effSig);
			obsBkg = rand.Poisson(effBkg);
			
			v[k] = (obsSig + obsBkg > 0) ? obsSig / TMath::Sqrt(obsSig + obsBkg) : 0;
		}
		sort(v.begin(),v.end()); // sort the vector
		
		med = (v[v.size()/2-1]+v[v.size()/2]) / 2.;
		if (med > bestMed.median) {
			bestMed.median = med;
			bestMed.err.first = med - v[(int32_t)(0.16*v.size())];
			bestMed.err.second = v[(int32_t)(0.84*v.size())-1] - med;
			bestCut = cut;
		}
	}
	
	cout << bestMed.median << "+" << bestMed.err.second << "-" << bestMed.err.first << "\t" << bestCut << "\t" << probSig << "\t" << probBkg << "\t " << path << endl;
	
bail:
	delete [] sigTest; delete [] sigTrain;
	delete [] bkgTest; delete [] bkgTrain;
	delete file;
} // evalTMVA2()

static void usage(const char *name)
{
	cerr << "usage: " << name << " -s <signal_events> -b <bkg_events> <files>" << endl;
} // usage()

int main(int argc, const char *argv [])
{
	Double_t num_sig = 0,num_bkg = 0;
	string arg;
	list<string> files;
	list<string>::const_iterator it;
	string::size_type p;
	int j = 1;
	
	while (j < argc) {
		
		arg = string(argv[j++]);
		
		// search signal yield
		if( (p = arg.find("-s")) == 0) {
			
			arg = arg.substr(p + 2);
			
			if (arg.compare("") == 0) {
				if (j >= argc) {
					usage(argv[0]);
					goto bail;
				}
				arg = string(argv[j++]);
			}
			
			num_sig = atof(arg.c_str());
		} // search for background yield
		else if ( (p = arg.find("-b")) == 0) {
			arg = arg.substr(p + 2);
			
			if (arg.compare("") == 0) {
				if (j <= argc) {
					usage(argv[0]);
					goto bail;
				}
				arg = string(argv[j++]);
			}
			
			num_bkg = atof(arg.c_str());
		}
		else {
			files.push_back(string(arg));
		}
	}
	
	// dump stuff
	cerr << "S = " << num_sig << ", B = " << num_bkg << endl;
	if (num_sig == 0 || num_bkg == 0) {
		usage(argv[0]);
		goto bail;
	}
	
	for (it = files.begin(); it != files.end(); ++it) {
		cerr << "Processing '" << *it << "'" << endl;
//		evalTMVA(it->c_str(),num_sig,num_bkg);
		evalTMVA2(it->c_str(),num_sig,num_bkg);
	}
	
bail:
	return 0;
} // main()
