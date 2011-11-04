#ifndef PLOTPU
#define PLOTPU

#include "plotClass.hh"

class plotPU: public plotClass {

public:

  plotPU(const char *files="anaBmm.default.files", const char *cuts = "default", const char *dir = "default", int mode = 11);
  ~plotPU();

  void effVsNpv(const char *var = "iso", double cut=0.75, const char *ylabel="#epsilon(I>0.75)", 
		const char *chan = "A", const char* dir="candAnaBu2JpsiK", const char *selection="Ao"); 

  void makeAll(int channels = 3);

  ClassDef(plotPU,1) //Testing plotPU

};


#endif

