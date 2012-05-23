/*
 *  bmm_estimator.h
 *  final_calculator
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 17.03.11.
 *  Copyright 2011 Christoph Nägeli. All rights reserved.
 *
 */

#ifndef BMM_ESTIMATOR
#define BMM_ESTIMATOR

#include "external_constants.h"
#include "bplus_estimator.h"

// STL
#include <map>
#include <string>

// ROOT headers
#include <TCut.h>
#include <TTree.h>

template<class T>
class triplet {
	public:
		triplet() : a(0,0), b(0,0), c(0,0) {}
		T a,b,c;
		T getIx(int ix) {
			T r;
			switch (ix) {
				case 0:
					r = a;
					break;
				case 1:
					r = b;
					break;
				case 2:
					r = c;
					break;
				default:
					break;
			}
			return r;
		}
		void setIx(int ix,T t) {
			switch (ix) {
				case 0:
					a = t;
					break;
				case 1:
					b = t;
					break;
				case 2:
					c = t;
					break;
				default:
					break;
			}
		}
};

// routines for B->mumu
std::map<int,triplet<measurement_t> > *estimate_bmm_eff(std::map<bmm_param,measurement_t> *bmm, TTree *accTree, TTree *sigTree, TTree *rareTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, std::pair<double,double> bd_window, std::pair<double,double> bs_window, bool is_bstomumu, std::map<systematics_t,double> *systematics_table);
void estimate_bmm_obs(std::map<bmm_param,measurement_t> *bmm, TTree *dataTree, double minEta, double maxEta, uint32_t channelIx, TCut anaCut, std::pair<double,double> bd_window, std::pair<double,double> bs_window, bool is_bstomumu, std::map<int,triplet<measurement_t> > *rare_effs, measurement_t tot_bu);

#endif
