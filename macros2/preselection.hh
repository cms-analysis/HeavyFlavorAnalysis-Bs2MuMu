#ifndef PRESELECTION
#define PRESELECTION

#include <string>
#include "TMath.h"
#include "RedTreeData.hh"

// ----------------------------------------------------------------------
std::string preselection();

// ----------------------------------------------------------------------
bool preselection(RedTreeData &b, int channel, bool rejectInvIso);

#endif

