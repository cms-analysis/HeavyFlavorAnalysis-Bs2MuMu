#include "AnalysisCuts.hh"

#include <iostream>
#include <iomanip>


ClassImp(AnalysisCuts)

// ----------------------------------------------------------------------
AnalysisCuts::AnalysisCuts() {
  for (int i = 0; i < MAXCUTS; ++i) {
    fCutName[i]     = TString("unused"); 
    fDescription[i] = TString("nothing said"); 
    fCutValue[i]    = 0;
    fAocValue[i]    = 0;
    fCumValue[i]    = 0;
  }
  fNcuts = 0; 
  fUpdated = 0; 
}


// ----------------------------------------------------------------------
AnalysisCuts::~AnalysisCuts() {
  if (fUpdated == 0) {
    cout << "xx> AnalysisCuts: You never called update(). All answers were wrong!" << endl;
  }
}

// ----------------------------------------------------------------------
void AnalysisCuts::addCut(const char *name, int &location) {
  addCut(name, "no description", location); 
}


// ----------------------------------------------------------------------
void AnalysisCuts::addCut(const char *name, const char *description, int &location) {
  fCutName[fNcuts]     = TString(name); 
  fDescription[fNcuts] = TString(description); 
  fCutLocation[fNcuts] = &location; 
  fAocValue[fNcuts]    = -1; 
  fCumValue[fNcuts]    = -1; 
  cout << "--> Added cut " << name << "  " << description << endl;
  fNcuts++; 
}

// ----------------------------------------------------------------------
int AnalysisCuts::getIndex(const char *name) {
  for (int i = 0; i < fNcuts; ++i) {
    if (!strcmp(fCutName[i].Data(), name)) return i; 
  }
  return -1; 
}

// ----------------------------------------------------------------------
const char* AnalysisCuts::getName(int icut) {
  if ((icut >= 0) && (icut < fNcuts)) {
    return fCutName[icut].Data(); 
  } else {
    return "xx> AnalysisCuts: cut index out of range"; 
  }
}

// ----------------------------------------------------------------------
const char* AnalysisCuts::getDescription(int icut) {
  if ((icut >= 0) && (icut < fNcuts)) {
    return fDescription[icut].Data(); 
  } else {
    return "xx> AnalysisCuts: cut index out of range"; 
  }
}

// ----------------------------------------------------------------------
void AnalysisCuts::update() {
  // -- Fill all cuts 
  for (int i = 0; i < fNcuts; ++i) {
    fCutValue[i] = *(fCutLocation[i]); 
  }

  // -- Cache derived information 
  for (int i = 0; i < fNcuts; ++i) {
    fCumValue[i] = 1; 
    for (int j = 0; j <= i; ++j) {
      if (singleCutTrue(j) == 0) {
	fCumValue[i] = 0; 
	break;
      }
    }

    fAocValue[i] = 1; 
    for (int j = 0; j < fNcuts; ++j) {
      if (j == i) continue; 
      if (singleCutTrue(j) == 0) {
	fAocValue[i] = 0; 
	break;
      }
    }
  }    

  fUpdated = 1; 
}

// ----------------------------------------------------------------------
int AnalysisCuts::singleCutTrue(int icut) {
  return fCutValue[icut];
}

// ----------------------------------------------------------------------
int AnalysisCuts::cumulativeCutTrue(int icut) {
  return fCumValue[icut];
}

// ----------------------------------------------------------------------
int AnalysisCuts::allOtherCutsTrue(int icut) {
  return fAocValue[icut];
}


// ----------------------------------------------------------------------
int AnalysisCuts::singleCutTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return singleCutTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
  }    
}

// ----------------------------------------------------------------------
int AnalysisCuts::cumulativeCutTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return cumulativeCutTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
  }    
}

// ----------------------------------------------------------------------
int AnalysisCuts::allOtherCutsTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return allOtherCutsTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
  }    
}
