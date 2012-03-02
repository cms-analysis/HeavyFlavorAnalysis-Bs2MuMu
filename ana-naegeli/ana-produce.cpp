/*
 *  ana-produce.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 01.03.12.
 *
 */

#include <iostream>
#include <string>

#include <TCut.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

static string output_dir(".");
static string configfile("");

struct table_entry_t {
	string filename;
	TCut cut;
	string chainname;
	void clear() {filename = string(""); cut = TCut(""); chainname = string("");}
	void Print() {
		cout << "filename: '" << filename << "'" << endl;
		cout << "cut: '" << cut.GetTitle() << "'" << endl;
		cout << "chainname: '" << chainname << "'" << endl;
	}
};

static void usage() {cout << "ana-produce [-D <output_dir>] <configfile>" << endl;}

static bool parse_arguments(const char **first, const char **last)
{
	bool ok = false;
	const char *arg;
	
	while (first != last) {
		
		arg = *first++;
		if (arg[0] == '-') {
			if (strcmp(arg, "-D") == 0) {
				if (first == last) {
					usage();
					goto bail;
				}
				output_dir = string(*first++);
			} else {
				cout << "unknown option: " << arg << endl;
				usage();
				goto bail;
			}
		} else
			configfile = string(arg);
	}
	
	cout << "Reading configuration file: '" << configfile << "'" << endl;
	cout << "Output directory: " << output_dir << endl;
	
	ok = true;
bail:
	return ok;
} // parse_arguments()

static bool read_config(const char *configfile, vector<table_entry_t> *table)
{
	FILE *file = fopen(configfile,"r");
	table_entry_t entry;
	const char *input_line;
	string line;
	string::iterator it,it2;
	bool ok = false;
	size_t len;
	
	if (!file) {
		cout << "Unable to open file '" << configfile << "'" << endl;
		goto bail;
	}
	
	// look for an opening curly brace
	while ( (input_line = fgetln(file, &len)) ) {
		
		entry.clear();
		
		if (input_line[0] == '#') // comment line
			continue;
		
		if (input_line[len-1] == '\n' || input_line[len-1] == '\r')
			len--;
		
		if (len == 0) // empty line
			continue;
		
		// append this line to the buffer
		line = string(input_line,len);
		
		it = find(line.begin(), line.end(), '{');
		if (it == line.end()) {
			cout << "No opening '{' found. Discarding the following line:" << endl;
			cout << line << endl;
			continue;
		}
		it++; // go past '{'
		
		// get the first colon for rootfile
		it2 = find(it, line.end(), ',');
		if (it2 == line.end()) {
			cout << "Invalid format of line:" << endl << line << endl;
			continue;
		};
		
		// save the name of the the outputfile
		entry.filename = string(it,it2);
		it = it2+1;
		
		// get the second semicolon for cut
		it2 = find(it, line.end(), ',');
		if (it2 == line.end()) {
			cout << "Invalid format of lines:" << endl << line << endl;
			continue;
		}
		
		// save the cut
		entry.cut = TCut( string(it, it2).c_str() );
		it = it2+1;
		
		// get the closeing brace
		it2 = find(it,line.end(),'}');
		if (it2 == line.end()) {
			cout << "no closing '}' found in line:" << endl << line << endl;
			continue;
		}
		
		entry.chainname = string(it,it2);
		
		// still here? save the entry
		table->push_back(entry);
	}
	
	ok = true;
bail:
	if (file)
		fclose(file);
	
	return ok;
} // read_config()

static void process_entry(table_entry_t e)
{
	FILE *chainfile = fopen(e.chainname.c_str(),"r");
	const char *input_line;
	char chain_entry[1024];
	size_t len;
	unsigned size;
	int parsed;
	string line;
	TChain chain("T");
	TFile *rFile = NULL;
	TTree *copyTree = NULL;
	
	cout << "====> Processing entry:" << endl;
	e.Print();
	
	// open the chain
	while ( (input_line = fgetln(chainfile, &len)) ) {
		
		if(input_line[0] == '#') continue;
		if(input_line[len-1] == '\n' || input_line[len-1] == '\r') len--;
		if(!len) continue;
		
		line = string(input_line,len);
		parsed = sscanf(line.c_str(), "%s %u", chain_entry, &size);
		if (parsed > 1)
			chain.Add(chain_entry,size);
		else if(parsed > 0)
			chain.Add(chain_entry);
	}
	
	// open the root file
	rFile = new TFile(Form("%s/%s",output_dir.c_str(),e.filename.c_str()),"recreate");
	copyTree = chain.CopyTree(e.cut,"",chain.GetEntries());
	
	copyTree->Write();
	if(chainfile)
		fclose(chainfile);
	
	delete rFile;
} // process_entry()

int main(int argc, const char *argv[])
{
	vector<table_entry_t> entries;
	if (argc < 2) {
		usage();
		goto bail;
	}
	
	// parse the arguments...
	if(!parse_arguments(&argv[1], &argv[argc])) goto bail;
	
	if(!read_config(configfile.c_str(),&entries)) goto bail;
	
	// process each entry
	for_each(entries.begin(), entries.end(), process_entry);
		
bail:
	return 0;
} // main()
