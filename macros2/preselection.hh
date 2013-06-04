#ifndef PRESELECTION
#define PRESELECTION

#include <string>
#include "TMath.h"
#include "TH1.h"
#include "RedTreeData.hh"

// ----------------------------------------------------------------------
std::string preselection();
TH1D* getPreselectionNumbers();

// ----------------------------------------------------------------------
bool preselection(RedTreeData &b, int channel);
void printRedTreeEvt(RedTreeData &b);

#endif

