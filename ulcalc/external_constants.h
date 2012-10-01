/*
 *  external_constants.h
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef EXTERNAL_CONSTANTS_H
#define EXTERNAL_CONSTANTS_H

// Standard headers
#include <cmath>
#include <map>
#include <string>
#include <utility>

// ROOT headers
#include <TCut.h>

// Systematic uncertainties
enum systematics_t {
	g_sys_acc_efftrack,		// tracking efficiency for each additional hadron (muon should cancel)
	g_sys_acc_ppro_barrel,	// acceptance: production process uncertainty in barrel (ratio)
	g_sys_acc_ppro_endcap,	// acceptance: production process uncertainty in endcap (ratio)
	g_sys_effana,			// systematic uncertainty on ration of selection efficiency
	g_sys_massscale,		// systematic uncertainty on P_{ij}
	g_sys_effmu_barrel,		// ratio on muon efficiency barrel
	g_sys_effmu_endcap,		// ratio on muon efficiency endcap
	g_sys_efftrig_barrel,	// ratio on trigger efficiency in barrel
	g_sys_efftrig_endcap,	// ratio on trigger efficiency in endcap
	g_sys_normfit,			// systematic on normalization fitting procedure
	g_sys_shapecombbkg		// systematic uncertainty on tau
};

// measurement class
class measurement_t {
	
	public:
		explicit measurement_t(double theVal = 0.0, double theErr = 0.0) { val = theVal; err = theErr; }
		
		double getVal() const { return val; }
		double getErr() const { return err; }
		void setVal(double newVal) { val = newVal; }
		void setErr(double newErr) { err = (newErr < 0.0) ? 0.0 : newErr; }
		
		measurement_t multiply(measurement_t m) const {
			double v = val * m.val;
			double e = sqrt(m.val*m.val*err*err + val*val*m.err*m.err);
			return measurement_t(v, e);
		}
		measurement_t divide(measurement_t m) const {
			double v = 0.0, e = 0.0;
			if (m.val != 0) {
				v = val / m.val;
				e = sqrt( err*err/(m.val*m.val) + (val*m.err / (m.val*m.val))*(val*m.err / (m.val*m.val)) );
			}
			return measurement_t(v, e);
		}
		measurement_t add(measurement_t m) const {
			double v = val + m.val;
			double e = sqrt(err*err + m.err*m.err);
			return measurement_t(v, e);
		}
		
		measurement_t operator*(measurement_t m) const { return multiply(m); }
		measurement_t operator/(measurement_t m) const { return divide(m); }
		measurement_t operator+(measurement_t m) const { return add(m); }
		
	private:
		double val;
		double err;
};

/* external constants */
const measurement_t c_s_theory();
const measurement_t c_d_theory();
const measurement_t f_ratio();
const measurement_t f_ratio_lb();

/* braching fractions */
const measurement_t bf_BsToMuMu();
const measurement_t bf_BsToKK();
const measurement_t bf_BsToKPi();
const measurement_t bf_BsToPiPi();
const measurement_t bf_BsToKMuNu();
const measurement_t bf_BdToMuMu();
const measurement_t bf_BdToPiPi();
const measurement_t bf_BdToKPi();
const measurement_t bf_BdToKK();
const measurement_t bf_BdToPiMuNu();
const measurement_t bf_LambdaBToPPi();
const measurement_t bf_LambdaBToPK();
const measurement_t bf_LambdaBToPMuNu();
const measurement_t bf_Bs2JpsiPhi();
const measurement_t bf_Bu2JpsiKp();
const measurement_t bf_Bd2JpsiKstar();
const measurement_t bf_Bd2JpsiKs();
const measurement_t bf_PsiToMuMu();
const measurement_t bf_Psi2SToMuMu();
const measurement_t bf_Ups1SToMuMu();
const measurement_t bf_Ups2SToMuMu();
const measurement_t bf_Ups3SToMuMu();

double bstomumu();
double bdtomumu();

/* cut to be applied for analysis */
extern const char *bmmGeneratorCuts;

/* Parameters estimated, saved in map */
enum bmm_param_tag {
	/* Bp -> J/psi Kp */
	kAcc_bplus = 1,
	kEff_mu_bplus,
	kEff_trig_bplus,
	kEff_cand_bplus,
	kEff_ana_bplus,
	kObs_bplus,
	kTot_bplus,
	/* B -> mumu */
	kAcc_bmm,
	kEff_mu_bmm,
	kEff_trig_bmm,
	kEff_cand_bmm,
	kEff_ana_bmm,
	kEff_total_bmm,	// Used for shortcut in Grid search
	kExp_bmm,
	kExpUncor_bmm,
	kProb_swind_bmm,
	kProb_dwind_bmm,
	kLow_signal_window_bmm,
	kHigh_signal_window_bmm,
	kObsBkg_bmm,
	kObsB_bmm,
	kPeakBkgOn_bmm,
	kPeakBkgOff_bmm,
	kTau_bmm,
	/* Unknown param */
	kUnknownParam
};
typedef std::pair<bmm_param_tag,int> bmm_param;

/* conversion routines */
std::string find_bmm_name(bmm_param_tag p);
bmm_param_tag find_bmm_param_by_name(std::string name, bool *bsparam);

/* Utility routines */
double std_dev_binomail(double lambda,double n);
void parse_cuts(const char *filename, std::map<double,TCut> *cuts_read);
measurement_t compute_efftot_bplus(std::map<bmm_param,measurement_t> *bmm, int channel);
measurement_t compute_efftot_bmm(std::map<bmm_param,measurement_t> *bmm, int channel);

#ifdef __linux__
char *fgetln(FILE *f, size_t *len);
#endif

#endif
