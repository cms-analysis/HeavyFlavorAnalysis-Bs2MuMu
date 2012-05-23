/*******************************************
 * helper routines to compute upper limits *
 *******************************************/

#include <TFile.h>
#include <TTree.h>
#include <RooDataSet.h>
#include <RooRealVar.h>

struct obs_t {
	int d0,d1;
	int s0,s1;
	void clear() {d0 = d1 = s0 = s1 = -1;}
	void print() {printf("d0=%d, d1=%d, s0=%d, s1=%d\n",d0,d1,s0,s1);}
};

bool operator<(const obs_t &o1, const obs_t& o2) {
	
	bool result = false;
	
	if		(o1.d0 < o2.d0)	result = true;
	else if	(o1.d0 > o2.d0)	result = false;
	else if (o1.d1 < o2.d1) result = true;
	else if (o1.d1 > o2.d1) result = false;
	else if (o1.s0 < o2.s0) result = true;
	else if (o1.s0 > o2.s0) result = false;
	else if (o1.s1 < o2.s1) result = true;
	else					result = false;
	
	return result;
}

void extract_into(RooDataSet *data, set<obs_t> &s)
{
	const RooArgSet *rooSet;
	obs_t o;
	
	for (int j = 0; j < data->numEntries(); j++) {
		rooSet = data->get(j);
		o.s0 = (int)round(((RooRealVar&)((*rooSet)["NsObs_0"])).getVal());
		o.s1 = (int)round(((RooRealVar&)((*rooSet)["NsObs_1"])).getVal());
		o.d0 = (int)round(((RooRealVar&)((*rooSet)["NdObs_0"])).getVal());
		o.d1 = (int)round(((RooRealVar&)((*rooSet)["NdObs_1"])).getVal());
		
		s.insert(o);
	}
	
	cout << s.size() << " independent entries found" << endl;
} // extract_into()

void write_set(set<obs_t> &s, const char *outfile)
{
	FILE *file = fopen(outfile,"w");
	
	for (set<obs_t>::const_iterator it = s.begin(); it != s.end(); ++it)
		fprintf(file,"%d %d %d %d\n",it->s0,it->d0,it->s1,it->d1);
	
	fclose(file);
} // write_set()

void extract(RooDataSet *data1, RooDataSet *data2, const char *outfile)
{	
	set<obs_t> s;
	if(data1) extract_into(data1,s);
	if(data2) extract_into(data2,s);
	write_set(s,outfile);
} // extract()

void extract(const char *data1, const char *data2, const char *outfile)
{
	TFile f1(data1);
	TFile f2(data2);
	
	RooDataSet *d1 = (RooDataSet*)f1.Get("signalsData");
	RooDataSet *d2 = (RooDataSet*)f2.Get("signalsData");
	
	extract(d1,d2,outfile);
} // extract()

map<obs_t,double> limit_table(const char *filename)
{
	map<obs_t,double> result;
	FILE *f = fopen(filename,"r");
	int dummy;
	char buffer[1024];
	obs_t o;
	double ul;
	
	while (true) {
		o.clear();
		if(!fgets(buffer, sizeof(buffer), f)) break;
		if(sscanf(buffer, "NbObs_0=%d, NsObs_0=%d, NdObs_0=%d", &dummy, &o.s0, &o.d0) < 3) {
			cerr << "Unable to parse line '" << buffer << "'" << endl;
			break;
		}
		
		if(!fgets(buffer, sizeof(buffer), f)) break;
		if(sscanf(buffer, "NbObs_1=%d, NsObs_1=%d, NdObs_1=%d", &dummy, &o.s1, &o.d1) < 3) {
			cerr << "Unable to parse line '" << buffer << "'" << endl;
			break;
		}
		
		if(!fgets(buffer, sizeof(buffer), f)) break;
		if(sscanf(buffer, "Bs->mumu upper limit with algorithm CLs: %lf",&ul) < 1) {
			cerr << "Unable to parse line '" << buffer << "'" << endl;
			break;
		}
		
		result[o] = ul;
	}
	
	fclose(f);
	
	cout << "map.size() = " << result.size() << endl;
	
	return result;
} // limit_table()

TTree *compute_expected(RooDataSet *data, map<obs_t,double> *table)
{
	TTree *tree = new TTree("L","");
	const RooArgSet *s;
	obs_t o;
	float ul;
	int j;
	int errs = 0;
	
	tree->Branch("NsObs_0",&o.s0,"NsObs_0/I");
	tree->Branch("NdObs_0",&o.d0,"NdObs_0/I");
	tree->Branch("NsObs_1",&o.s1,"NsObs_1/I");
	tree->Branch("NdObs_1",&o.d1,"NdObs_1/I");
	tree->Branch("ul",&ul,"ul/F");
	
	// iterate through data
	for (j = 0; j < data->numEntries(); j++) {
		
		s = data->get(j);
		o.d0 = ((RooRealVar&)(*s)["NdObs_0"]).getVal();
		o.d1 = ((RooRealVar&)(*s)["NdObs_1"]).getVal();
		o.s0 = ((RooRealVar&)(*s)["NsObs_0"]).getVal();
		o.s1 = ((RooRealVar&)(*s)["NsObs_1"]).getVal();
		
		if (table->count(o) > 0) {
			ul = (*table)[o];
			tree->Fill();
		} else
			errs++;
	}
	
	cout << "Lookup table is missing " << errs << " / " << data->numEntries() << " entries" << endl;
	
	return tree;
} // compute_expected()

TTree *compute_expected(RooDataSet *data, const char *filename)
{
	map<obs_t,double> table = limit_table(filename);
	return compute_expected(data,&table);
} // compute_expected()
