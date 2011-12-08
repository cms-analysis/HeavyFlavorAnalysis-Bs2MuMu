#ifndef PLOTEFFICIENCIES
#define PLOTEFFICIENCIES

#include "plotClass.hh"

class plotEfficiencies: public plotClass {

public:

  plotEfficiencies(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotEfficiencies();

  void makeAll(int channels = 3);
  
  void tnpVsMC(double m1pt, double m2pt, int chan = 3, std::string what = "default");
  void triggerSignal(std::string cuts = "fls3d>8&&chi2/dof<2");
  void triggerNormalization(std::string cuts = "fls3d>8&&chi2/dof<2&&alpha<0.05&&iso>0.5");

  ClassDef(plotEfficiencies,1) //Testing plotEfficiencies

};


#endif

