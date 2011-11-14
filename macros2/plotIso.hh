#ifndef PLOTISO
#define PLOTISO

#include "plotClass.hh"

class plotIso: public plotClass {

public:

  plotIso(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotIso();
  
  
  void makeAll(int channels = 3);
  void dataVsMc(std::string file1, std::string dir1, std::string file2, std::string dir2, const char *selection);
  

  int fMode;

  ClassDef(plotIso,1) //Testing plotOverlays

};


#endif

