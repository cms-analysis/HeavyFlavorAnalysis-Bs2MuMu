#include "AnalysisCuts.hh"

#include <iostream>
#include <iomanip>


ClassImp(AnalysisCuts)

// ----------------------------------------------------------------------
AnalysisCuts::AnalysisCuts() {
  for (int i = 0; i < MAXCUTS; ++i) {
    fCutName[i]     = TString("unused"); 
    fDescription[i] = TString("nothing said"); 
    fCutValue[i]    = false;
    fAocValue[i]    = false;
    fCumValue[i]    = false;
    fNm1Value[i]    = false;
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

// // ----------------------------------------------------------------------
// void AnalysisCuts::addCut(const char *name, int &location) {
//   addCut(name, "no description", location); 
// }


// // ----------------------------------------------------------------------
// void AnalysisCuts::addCut(const char *name, const char *description, int &location) {
//   fCutName[fNcuts]     = TString(name); 
//   fDescription[fNcuts] = TString(description); 
//   fCutLocation[fNcuts] = &location; 
//   fAocValue[fNcuts]    = -1; 
//   fCumValue[fNcuts]    = -1; 
//   cout << "--> Added cut " << name << "  " << description << endl;
//   fNcuts++; 
// }



// ----------------------------------------------------------------------
void AnalysisCuts::dumpAll() {
  for (int i = 0; i < fNcuts; ++i) {
    cout << Form("Cut %2d: %s (%s) at %d", i, fCutName[i].Data(), fDescription[i].Data(), (fCutLocation[i]?1:0)) 
	 << Form(" nm=%i si=%i cu=%i", (fNm1Value[i]?1:0),  (fCutValue[i]?1:0),  (fCumValue[i]?1:0))
	 << endl; 
  }
}

// ----------------------------------------------------------------------
void AnalysisCuts::addCut(const char *name, bool &location) {
  addCut(name, "no description", location); 
}


// ----------------------------------------------------------------------
void AnalysisCuts::addCut(const char *name, const char *description, bool &location) {
  fCutName[fNcuts]     = TString(name); 
  fDescription[fNcuts] = TString(description); 
  fCutLocation[fNcuts] = &location; 
  fAocValue[fNcuts]    = false; 
  fCumValue[fNcuts]    = false; 
  fNm1Value[fNcuts]    = false;
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
    fNm1Value[i] = true; 
    for (int j = 0; j < i; ++j) {
      if (singleCutTrue(j) == false) {
	fNm1Value[i] = false; 
	break;
      }
    }


    fCumValue[i] = true; 
    for (int j = 0; j <= i; ++j) {
      if (singleCutTrue(j) == false) {
	fCumValue[i] = false; 
	break;
      }
    }

    fAocValue[i] = true; 
    for (int j = 0; j < fNcuts; ++j) {
      if (j == i) continue; 
      if (singleCutTrue(j) == false) {
	fAocValue[i] = false; 
	break;
      }
    }
  }    

  if (0 == fUpdated)  cout << "Updating AnalysisCuts for the first time" << endl;
  fUpdated = 1; 

}

// ----------------------------------------------------------------------
bool AnalysisCuts::singleCutTrue(int icut) {
  return fCutValue[icut];
}

// ----------------------------------------------------------------------
bool AnalysisCuts::cumulativeCutTrue(int icut) {
  return fCumValue[icut];
}

// ----------------------------------------------------------------------
bool AnalysisCuts::nMinus1CutsTrue(int icut) {
  return fNm1Value[icut];
}

// ----------------------------------------------------------------------
bool AnalysisCuts::allOtherCutsTrue(int icut) {
  return fAocValue[icut];
}


// ----------------------------------------------------------------------
bool AnalysisCuts::singleCutTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return singleCutTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
    return false;
  }    
}

// ----------------------------------------------------------------------
bool AnalysisCuts::cumulativeCutTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return cumulativeCutTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
    return false;
  }    
}

// ----------------------------------------------------------------------
bool AnalysisCuts::nMinus1CutsTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return nMinus1CutsTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
    return false;
  }    
}

// ----------------------------------------------------------------------
bool AnalysisCuts::allOtherCutsTrue(const char *name) {
  int icut = getIndex(name);
  if ((icut >= 0) && (icut < fNcuts)) {
    return allOtherCutsTrue(icut);
  } else {
    cout << "xx> AnalysisCuts: Unknown cut " << name << endl;
    return false;
  }    
}
