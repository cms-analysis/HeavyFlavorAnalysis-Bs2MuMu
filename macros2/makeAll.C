void makeAll(std::string files = "anaBmm.v11.files", int channels = 31) {

  if (channels & 16) {
    plotEfficiencies a3(files.c_str()); 
    a3.makeAll(7); 
  }

  if (channels & 8) {
    plotPU a2(files.c_str()); 
    a2.makeAll(); 
  }

  if (channels & 4) {
    plotOverlays a1(files.c_str()); 
    a1.makeAll(); 
  }

  if (channels & 1) {
    plotResults a4(files.c_str()); 
    a4.makeAll(4);

    plotResults a0(files.c_str()); 
    a0.makeAll(1);

    if (channels & 2) {
      a0.makeAll(2);
    }
  }
  
}
