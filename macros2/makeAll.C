#include "plotEfficiencies.hh"
#include "plotPU.hh"
#include "plotOverlays.hh"
#include "plotResults.hh"

void makeAll(std::string files = "anaBmm.v11.files", int channels = 31) {

  plotBDT a3(files.c_str()); 
  a3.makeAll(7);

  plotReducedOverlays a4(files.c_str()); 
  a4.makeAll("Presel");

  plotResults a5(files.c_str()); 
  a5.makeAll(4);
  
}
