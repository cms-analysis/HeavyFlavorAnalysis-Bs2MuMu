#ifndef PLOTOVERLAYS
#define PLOTOVERLAYS

#include "plotClass.hh"

class plotOverlays: public plotClass {

public:

  plotOverlays(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotOverlays();


  void makeAll(int verbose = 0);
  
  void productionMechanism(std::string dsample, std::string msample, std::string region = "A");
  void fitDistribution(std::string var, std::string sample, std::string mcsample, std::string selection);

  void profileVsEta(const char *var, const char *chan, const char *dir, const char *selection);

  void invertedMuonID(std::string var, std::string cuts, double lo, double hi, int nbin);

  void sbsDistributionOverlay(std::string file1, std::string dir1, std::string region1,
			      std::string file2, std::string dir2, std::string region2 = "A", 
			      std::string selection="Ao"); 
  void sbsDistributionOverlaySameFile(std::string file1, std::string dir1, 
				      std::string r1 = "B_iso5", std::string r2 = "E_iso5", std::string sel = "Ao", 
				      std::string L1="barrel", std::string L2="endcap"); 
  std::vector<double> computeDelta(double lo1, double lo1E, double lo2, double lo2E); 

  ClassDef(plotOverlays,1) //Testing plotOverlays

};


#endif

