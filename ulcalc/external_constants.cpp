/*
 *  external_constants.cpp
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#include "external_constants.h"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

////////////////////////
// External constants //
////////////////////////
// NEW LHCb VALUE
static const measurement_t fs_by_fu(0.267,0.021,0.021);

// PDG Values
static const measurement_t f_s(0.113,0.013,0.013);
static const measurement_t f_u(0.401,0.013,0.013);
static const measurement_t f_bbaryon(0.085,0.022,0.022);

// branching fractions
static const measurement_t sm_bstomumu(3.2e-9,0,0);
static const measurement_t sm_bdtomumu(1.0e-10,0,0);
static const measurement_t bptojpsik(1.014e-3,0.034e-3,0.034e-3);
static const measurement_t jpsitomumu(0.0593,0.0006,0.0006);

// rare backgrounds
static const measurement_t BsToKK(2.7e-5,0,0);
static const measurement_t BsToKPi(5.0e-6,0,0);
static const measurement_t BsToPiPi(1.2e-6,0,0);
static const measurement_t BsToKMuNu(1.3e-4,0,0);
static const measurement_t BdToPiPi(5.2e-6,0,0);
static const measurement_t BdToKPi(1.9e-5,0,0);
static const measurement_t BdToKK(1.5e-7,0,0);
static const measurement_t BdToPiMuNu(1.3e-4,0,0);
static const measurement_t LambdaBToPPi(3.5e-6,0,0);
static const measurement_t LambdaBToPK(5.6e-6,0,0);
static const measurement_t LambdaBToPMuNu(1.3e-4,0,0);

using std::map;

std::string find_bmm_name(bmm_param_tag p)
{
	std::string result;
	switch (p) {
		case kAcc_bplus:
			result = "ACC_BPLUS";
			break;
		case kEff_mu_bplus:
			result = "EFF_MU_BPLUS";
			break;
		case kEff_trig_bplus:
			result = "EFF_TRIG_BPLUS";
			break;
		case kEff_cand_bplus:
			result = "EFF_CAND_BPLUS";
			break;
		case kEff_ana_bplus:
			result = "EFF_ANA_BPLUS";
			break;
		case kObs_bplus:
			result = "OBS_BPLUS";
			break;
		case kAcc_bmm:
			result = "ACC_BMM";
			break;
		case kEff_mu_bmm:
			result = "EFF_MU_BMM";
			break;
		case kEff_trig_bmm:
			result = "EFF_TRIG_BMM";
			break;
		case kEff_cand_bmm:
			result = "EFF_CAND_BMM";
			break;
		case kEff_ana_bmm:
			result = "EFF_ANA_BMM";
			break;
		case kExp_bmm:
			result = "EXP_BMM";
			break;
		case kProb_swind_bmm:
			result = "PROB_SWIND";
			break;
		case kProb_dwind_bmm:
			result = "PROB_DWIND";
			break;
		case kObsBkg_bmm:
			result = "OBS_BKG";
			break;
		case kLow_signal_window_bmm:
			result = "LOW_SIGNAL_WINDOW";
			break;
		case kHigh_signal_window_bmm:
			result = "HIGH_SIGNAL_WINDOW";
			break;
		case kTot_bplus:
			result = "TOT_BPLUS";
			break;
		case kEff_total_bmm:
			result = "EFF_TOT_BMM";
			break;
		case kObsB_bmm:
			result = "OBS_BMM";
			break;
		case kPeakBkgOn_bmm:
			result = "PEAK_BKG_ON";
			break;
		case kPeakBkgOff_bmm:
			result = "PEAK_BKG_OFF";
			break;
		case kTau_bmm:
			result = "TAU";
			break;
		default:
			std::cerr << "Unknown bmm_param: " << p << std::endl;
			abort();
			break;
	}
	
	return result;
} // find_bmm_name()

bmm_param_tag find_bmm_param_by_name(std::string name, bool *bsparam)
{
	*bsparam = true;
	bmm_param_tag result = kUnknownParam;
	
	if (name.compare("OBS_BKG") == 0) {
		result = kObsBkg_bmm;
	} else if (name.compare("LOW_BD") == 0) {
		result = kLow_signal_window_bmm;
		*bsparam = false;
	} else if (name.compare("HIGH_BD") == 0) {
		result = kHigh_signal_window_bmm;
		*bsparam = false;
	} else if (name.compare("LOW_BS") == 0) {
		result = kLow_signal_window_bmm;
	} else if (name.compare("HIGH_BS") == 0) {
		result = kHigh_signal_window_bmm;
	} else if (name.compare("PSS") == 0) {
		result = kProb_swind_bmm;
	} else if (name.compare("PSD") == 0) {
		result = kProb_swind_bmm;
		*bsparam = false;
	} else if (name.compare("PDS") == 0) {
		result = kProb_dwind_bmm;
	} else if (name.compare("PDD") == 0) {
		result = kProb_dwind_bmm;
		*bsparam = false;
	} else if (name.compare("TOT_BPLUS") == 0) {
		result = kTot_bplus;
	} else if (name.compare("EFF_TOT_BSMM") == 0) {
		result = kEff_total_bmm;
	} else if (name.compare("EFF_TOT_BDMM") == 0) {
		result = kEff_total_bmm;
		*bsparam = false;
	} else if (name.compare("OBS_BSMM") == 0) {
		result = kObsB_bmm;
	} else if (name.compare("OBS_BDMM") == 0) {
		result = kObsB_bmm;
		*bsparam = false;
	} else if (name.compare("ACC_BPLUS") == 0) {
		result = kAcc_bplus;
	} else if (name.compare("EFF_MU_BPLUS") == 0) {
		result = kEff_mu_bplus;
	} else if (name.compare("EFF_TRIG_BPLUS") == 0) {
		result = kEff_trig_bplus;
	} else if (name.compare("EFF_CAND_BPLUS") == 0) {
		result = kEff_cand_bplus;
	} else if (name.compare("EFF_ANA_BPLUS") == 0) {
		result = kEff_ana_bplus;
	} else if (name.compare("OBS_BPLUS") == 0) {
		result = kObs_bplus;
	} else if (name.compare("ACC_BSMM") == 0) {
		result = kAcc_bmm;
	} else if (name.compare("EFF_MU_BSMM") == 0) {
		result = kEff_mu_bmm;
	} else if (name.compare("EFF_TRIG_BSMM") == 0) {
		result = kEff_trig_bmm;
	} else if (name.compare("EFF_CAND_BSMM") == 0) {
		result = kEff_cand_bmm;
	} else if (name.compare("EFF_ANA_BSMM") == 0) {
		result = kEff_ana_bmm;
	} else if (name.compare("ACC_BDMM") == 0) {
		result = kAcc_bmm;
		*bsparam = false;
	} else if (name.compare("EFF_MU_BDMM") == 0) {
		result = kEff_mu_bmm;
		*bsparam = false;
	} else if (name.compare("EFF_TRIG_BDMM") == 0) {
		result = kEff_trig_bmm;
		*bsparam = false;
	} else if (name.compare("EFF_CAND_BDMM") == 0) {
		result = kEff_cand_bmm;
		*bsparam = false;
	} else if (name.compare("EFF_ANA_BDMM") == 0) {
		result = kEff_ana_bmm;
		*bsparam = false;
	} else if (name.compare("PEAK_BKG_OFF") == 0) {
		result = kPeakBkgOff_bmm;
	} else if (name.compare("PEAK_BKG_BS") == 0) {
		result = kPeakBkgOn_bmm;
	} else if (name.compare("PEAK_BKG_BD") == 0) {
		result = kPeakBkgOn_bmm;
		*bsparam = false;
	} else if (name.compare("TAU_BS") == 0) {
		result = kTau_bmm;
	} else if (name.compare("TAU_BD") == 0) {
		result = kTau_bmm;
		*bsparam = false;
	}
	
	return result;
} // find_bmm_param_by_name()

const measurement_t bf_ratio_bsmm()
{
	measurement_t c_s = sm_bstomumu / (bptojpsik * jpsitomumu);
	
	return c_s;
} // c_s_theory()

const measurement_t bf_ratio_bdmm()
{
	measurement_t c_d = sm_bdtomumu / (bptojpsik * jpsitomumu);
	return c_d;
} // c_d_theory()

const measurement_t f_ratio()
{
	return fs_by_fu;
} // f_ratio()

const measurement_t f_ratio_lb()
{
	return f_bbaryon / f_u;
} // f_ratio_lb()

const measurement_t bf_BsToMuMu()		{return sm_bstomumu;}
const measurement_t bf_BsToKK()			{return BsToKK;}
const measurement_t bf_BsToKPi()		{return BsToKPi;}
const measurement_t bf_BsToPiPi()		{return BsToPiPi;}
const measurement_t bf_BsToKMuNu()		{return BsToKMuNu;}
const measurement_t bf_BdToMuMu()		{return sm_bdtomumu;}
const measurement_t bf_BdToPiPi()		{return BdToPiPi;}
const measurement_t bf_BdToKPi()		{return BdToKPi;}
const measurement_t bf_BdToKK()			{return BdToKK;}
const measurement_t bf_BdToPiMuNu()		{return BdToPiMuNu;}
const measurement_t bf_LambdaBToPPi()	{return LambdaBToPPi;}
const measurement_t bf_LambdaBToPK()	{return LambdaBToPK;}
const measurement_t bf_LambdaBToPMuNu()	{return LambdaBToPMuNu;}
const measurement_t bf_Bu2JpsiKp()		{return bptojpsik;}
const measurement_t bf_PsiToMuMu()		{return jpsitomumu;}

double bstomumu()
{
	return sm_bstomumu.getVal();
} // bstomumu()

double bdtomumu()
{
	return sm_bdtomumu.getVal();
} // bdtomumu()

double std_dev_binomail(double lambda,double n)
{
	return sqrt(lambda*(1-lambda)/ n);
} // std_dev_binomail()

#ifdef __linux__
/* Linux does not have that nice function fgetln
 * Use a very simple approach in that case */
char *fgetln(FILE *f, size_t *len)
{
	static char buffer[2048];
	char *result;
	
	result = fgets(buffer,sizeof(buffer),f);
	if (result) *len = strlen(buffer);
	else *len = 0;
	
	return result;
} // fgetln()
#endif

void parse_cuts(const char *filename, std::map<double,TCut> *cuts_read)
{
	using namespace std;
	FILE *file = fopen(filename, "r");
	char *input_line;
	size_t line_len;
	string line;
	string delim(" \t");
	string::iterator it;
	string::iterator it_e;
	double eta;
	TCut cut;
	
	while ( (input_line = fgetln(file, &line_len)) ) {
		
		if (input_line[0] == '#') // comment line
			continue;
		
		// remove endline at the end of line (if any)
		if (input_line[line_len-1] == '\n' || input_line[line_len-1] == '\r')
			line_len--;
		
		if (line_len == 0)
			continue;
		
		// create a c++ string out of input line
		line = string(input_line,line_len);
		
		// find first newline or tab
		it = find_first_of(line.begin(), line.end(), delim.begin(), delim.end());
		if (it == line.end()) { // no space, invalid line
			cerr << "Invalid line in cuts file: " << line << endl;
			continue;
		}
		
		eta = atof(string(line,0,it-line.begin()).c_str());
		if (eta == 0) {
			cerr << "eta in cuts file has to be > 0. Ignoring line: " << line << endl;
			continue;
		}
		
		it = find(it, line.end(), '"');
		it += 1;
		it_e = find(it, line.end(), '"');
		
		// save the cuts
		(*cuts_read)[eta] = TCut(string(line, it-line.begin(), it_e-it).c_str());
	}
	
	fclose(file);
} // parse_cuts()

measurement_t compute_efftot_bplus(map<bmm_param,measurement_t> *bmm, int channel)
{
	using std::make_pair;
	measurement_t eff = (*bmm)[make_pair(kAcc_bplus,channel)] * (*bmm)[make_pair(kEff_mu_bplus,channel)] * (*bmm)[make_pair(kEff_trig_bplus,channel)] * (*bmm)[make_pair(kEff_cand_bplus,channel)] * (*bmm)[make_pair(kEff_ana_bplus,channel)];
	
	return eff;
} // compute_efftot_bplus()

measurement_t compute_efftot_bmm(map<bmm_param,measurement_t> *bmm, int channel)
{
	using std::make_pair;
	measurement_t eff(0,0,0);
	
	if(bmm->count(make_pair(kEff_total_bmm, channel)) > 0)
		eff = (*bmm)[make_pair(kEff_total_bmm, channel)];
	else
		eff = (*bmm)[make_pair(kAcc_bmm,channel)] * (*bmm)[make_pair(kEff_mu_bmm,channel)] * (*bmm)[make_pair(kEff_trig_bmm,channel)] * (*bmm)[make_pair(kEff_cand_bmm,channel)] * (*bmm)[make_pair(kEff_ana_bmm,channel)];
	
	return eff;
} // compute_efftot_bmm()
