#include "plotEfficiencies.hh"
#include "plotPU.hh"
#include "plotOverlays.hh"
#include "plotResults.hh"

void makeAll(std::string files = "anaBmm.v11.files", int channels = 31) {

  plotResults a4(files.c_str()); 
  a4.makeAll(4);
  
}
