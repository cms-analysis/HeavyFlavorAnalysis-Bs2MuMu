#ifndef PLOTOVERLAYS
#define PLOTOVERLAYS

#include "plotClass.hh"

class plotOverlays: public plotClass {

public:

  plotOverlays(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotOverlays();


  void makeAll(int channels = 3);
  
  ClassDef(plotOverlays,1) //Testing plotOverlays

  void sbsDistributionOverlay(std::string file1, std::string dir1, 
			      std::string file2, std::string dir2, const char *selection="Ao", const char *region = ""); 
  void sbsDistributionOverlaySameFile(std::string file1, std::string dir1, 
				      std::string r1 = "B_iso5", std::string r2 = "E_iso5", std::string sel = "Ao", 
				      std::string L1="barrel", std::string L2="endcap"); 


};


#endif

