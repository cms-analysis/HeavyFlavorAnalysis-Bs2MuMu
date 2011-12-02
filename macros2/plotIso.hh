#ifndef PLOTISO
#define PLOTISO

#include "plotClass.hh"

struct data {
  std::string name; 
  double doca, r, pt; 
  double eff1, eff2, ratio; 
};


class plotIso: public plotClass {

public:

  plotIso(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotIso();
  
  
  void makeAll(double isoCut = 0.80);
  void dataVsMc(std::string file1, std::string dir1, std::string file2, std::string dir2, std::string selection);
  void anaTextFiles(int cutval = 80);
  void readFile(std::string fname, std::vector<data> &v);

  int ptBin(double pt);
  int rBin(double r);
  void relabel(TH2 *h);

  int fMode;
  double fIsoCut; 

  ClassDef(plotIso,1) //Testing plotOverlays

};


#endif

